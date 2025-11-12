// __BEGIN_LICENSE__
//  Copyright (c) 2009-2025, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NGT platform is licensed under the Apache License, Version 2.0 (the
//  "License"); you may not use this file except in compliance with the
//  License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__

#include <vw/Cartography/PointImageManipulation.h>
#include <vw/Cartography/GeoTransform.h>
#include <vw/Cartography/Map2CamTrans.h>
#include <vw/Cartography/CameraBBox.h>
#include <vw/Image/MaskViews.h>
#include <vw/Camera/CameraModel.h>

namespace vw { namespace cartography {

  // This transform is not thread-safe, unless caching is turned off. Use
  // mapproj_trans_copy() to make a copy. See that function for more details.
  Map2CamTrans::Map2CamTrans(vw::camera::CameraModel const* cam,
                              GeoReference const& image_georef,
                              GeoReference const& dem_georef,
                              std::string const& dem_file,
                              vw::Vector2i const& image_size,
                              bool call_from_mapproject,
                              bool nearest_neighbor):
    m_cam(cam), m_image_georef(image_georef), m_dem_georef(dem_georef),
    m_dem_file(dem_file), m_dem(dem_file), m_image_size(image_size),
    m_call_from_mapproject(call_from_mapproject),
    m_nearest_neighbor(nearest_neighbor), m_has_nodata(false),
    m_nodata(std::numeric_limits<double>::quiet_NaN()),
    m_use_cache(true) {

    boost::shared_ptr<vw::DiskImageResource>
      dem_rsrc(vw::DiskImageResourcePtr(dem_file));

    m_has_nodata = dem_rsrc->has_nodata_read();
    if (m_has_nodata) m_nodata = dem_rsrc->nodata_read();

    m_invalid_pix = vw::camera::CameraModel::invalid_pixel();

    // This is the full masked DEM. There is also m_cropped_masked_dem,
    // which should be used per tile.
    m_masked_dem = create_mask(m_dem, m_nodata);

    // An estimate of the DEM height can help the reliability of intersecting
    // a ray with the DEM.
    m_height_guess = vw::cartography::demHeightGuess(m_masked_dem);

    // Set up interpolation interface to the data we loaded into memory. The
    // nearest neighbor logic is useful when mapprojecting a mask.
    if (m_nearest_neighbor)
      m_interp_dem = interpolate(m_masked_dem,
                                 NearestPixelInterpolation(), ZeroEdgeExtension());
    else
      m_interp_dem = interpolate(m_masked_dem,
                                 BicubicInterpolation(), ZeroEdgeExtension());
  }

  void Map2CamTrans::set_use_cache(bool use_cache) {
    m_use_cache = use_cache;
  }

  // This function is not thread-safe by default. See above.
  vw::Vector2 Map2CamTrans::reverse(const vw::Vector2 &p) const {

    // If we have data for the location already cached. This is only useful
    // when processing tiles. For individual samples need to turn off the caching
    // ahead of time as it only adds overhead.
    if (m_use_cache && m_img_cache_box.contains(p)) {
      // Interpolate the output value using the cached data
      PixelMask<Vector2> v = m_cache_interp_mask(p.x() - m_img_cache_box.min().x(),
                                                 p.y() - m_img_cache_box.min().y());
      // We can just return the value if it is valid
      if (is_valid(v)) return v.child();
      else             return m_invalid_pix;
    }

    // No cached data
    int b = BicubicInterpolation::pixel_buffer;
    if (m_nearest_neighbor)
      b = NearestPixelInterpolation::pixel_buffer;
    Vector2 lonlat  = m_image_georef.pixel_to_lonlat(p);
    Vector2 dem_pix = m_dem_georef.lonlat_to_pixel(lonlat);
    if ((dem_pix[0] < b - 1) || (dem_pix[0] >= m_dem.cols() - b) ||
        (dem_pix[1] < b - 1) || (dem_pix[1] >= m_dem.rows() - b)) {
      // No DEM data
      return m_invalid_pix;
    }

    PixelMask<float> h;
    bool no_cache = true;
    if (m_use_cache) {
      no_cache = false;
      Vector2 crop_pix = dem_pix - m_dem_cache_box.min(); // since we cropped the DEM
      if (m_dem_cache_box.empty() ||
          (crop_pix[0] < b - 1) || (crop_pix[0] >= m_cropped_dem.cols() - b) ||
          (crop_pix[1] < b - 1) || (crop_pix[1] >= m_cropped_dem.rows() - b)) {
        no_cache = true;
      } else {
        h = m_cropped_interp_dem(crop_pix[0], crop_pix[1]);
      }
    }
    
    // If there is no cached data, use the full DEM
    if (no_cache) 
      h = m_interp_dem(dem_pix[0], dem_pix[1]);

    if (!is_valid(h))
      return m_invalid_pix;
    
    vw::Vector3 llh = Vector3(lonlat[0], lonlat[1], h.child());
    Vector3 xyz = m_dem_georef.datum().geodetic_to_cartesian(llh);
    Vector2 pt;
    try {
      pt = m_cam->point_to_pixel(xyz);
      if (m_call_from_mapproject &&
           (pt[0] < b - 1 || pt[0] >= m_image_size[0] - b ||
            pt[1] < b - 1 || pt[1] >= m_image_size[1] - b)) {
        // Won't be able to interpolate into image in transform(...)
        return m_invalid_pix;
      }
    } catch(...) { // If a point failed to project
      return m_invalid_pix;
    }

    return pt;
  }

  // Given a raw pixel, intersect the ray with the DEM and return the pixel
  // on the mapprojected image. This throws an exception if the intersection
  // failed.
  vw::Vector2 Map2CamTrans::forward(const vw::Vector2 &p) const {

    // TODO(oalexan1): If this logic is used as an inner loop by a CERES solver,
    // the tolerance here may not be good enough, as CERES makes very small
    // steps, and the result from here may be jumpy.
    double height_error_tol = 1e-3; // 1 mm
    double max_abs_tol      = 1e-14; // abs cost function change
    double max_rel_tol      = 1e-14; // rel cost function change
    int    num_max_iter     = 50;
    bool   treat_nodata_as_zero = false;
    bool has_intersection = false;

    vw::Vector3 prev_xyz(0, 0, 0); // do not have an initial guess
    vw::Vector3 xyz = vw::cartography::camera_pixel_to_dem_xyz
       (m_cam->camera_center(p), m_cam->pixel_to_vector(p), m_masked_dem,
        m_dem_georef, treat_nodata_as_zero, has_intersection,
        height_error_tol, max_abs_tol, max_rel_tol, num_max_iter,
        prev_xyz, m_height_guess);

    // If the intersection failed, throw an exception
    if (!has_intersection || xyz == vw::Vector3())
      vw::vw_throw(vw::ArgumentErr() << "Failed to intersect the camera ray with the DEM.\n");

    vw::Vector3 llh = m_dem_georef.datum().cartesian_to_geodetic(xyz);
    return m_image_georef.lonlat_to_pixel(vw::Vector2(llh[0], llh[1]));
  }

  // This function is not thread-safe, by default. See above. This function is slow,
  // but later speeds things up. Better not use it for sparse pixel queries.
  void Map2CamTrans::cache_dem(vw::BBox2i const& bbox) const {

    // TODO: This may fail around poles. Need to do the standard X trick, traverse
    // the edges and diagonals of the box. Use here the function sample_float_bbox().
    BBox2 dbox;
    dbox.grow(m_dem_georef.lonlat_to_pixel(m_image_georef.pixel_to_lonlat(Vector2(bbox.min().x(),   bbox.min().y())))); // Top left
    dbox.grow(m_dem_georef.lonlat_to_pixel(m_image_georef.pixel_to_lonlat(Vector2(bbox.max().x()-1, bbox.min().y())))); // Top right
    dbox.grow(m_dem_georef.lonlat_to_pixel(m_image_georef.pixel_to_lonlat(Vector2(bbox.min().x(),   bbox.max().y()-1)))); // Bottom left
    dbox.grow(m_dem_georef.lonlat_to_pixel(m_image_georef.pixel_to_lonlat(Vector2(bbox.max().x()-1, bbox.max().y()-1)))); // Bottom right

    // A lot of care is needed here when going from real box to int
    // box, and if in doubt, better expand more rather than less.
    dbox.expand(1);
    m_dem_cache_box = grow_bbox_to_int(dbox);
    if (m_nearest_neighbor)
      m_dem_cache_box.expand(NearestPixelInterpolation::pixel_buffer); // for interp
    else
      m_dem_cache_box.expand(BicubicInterpolation::pixel_buffer); // for interp
    m_dem_cache_box.crop(bounding_box(m_dem));

    // Read the dem in memory for speed in the region of the expanded bounding box.
    m_cropped_dem = crop(m_dem, m_dem_cache_box);

    if (m_has_nodata) {
      m_cropped_masked_dem = create_mask(m_cropped_dem, m_nodata);
    } else { // Don't need to handle nodata
      m_cropped_masked_dem = pixel_cast<PixelMask<float>>(m_cropped_dem);
    }
    // Set up interpolation interface to the data we loaded into memory
    if (m_nearest_neighbor)
      m_cropped_interp_dem = interpolate(m_cropped_masked_dem,
                                         NearestPixelInterpolation(), ZeroEdgeExtension());
    else
      m_cropped_interp_dem = interpolate(m_cropped_masked_dem,
                                         BicubicInterpolation(), ZeroEdgeExtension());
  } // End function cache_dem

  // This function will be called whenever we start to apply the
  // transform in a tile. It computes and caches the point cloud at
  // each pixel in the tile, to be used later when we iterate over pixels.
  // This function is not thread-safe, see above. See cache_dem() for
  // a note on performance.
  vw::BBox2i Map2CamTrans::reverse_bbox(vw::BBox2i const& bbox) const {

    // Custom reverse_bbox() function which can handle invalid pixels.
    if (!m_cached_rv_box.empty()) return m_cached_rv_box;
    
    cache_dem(bbox);

    // Cache the reverse transform
    m_img_cache_box = BBox2i();
    BBox2i local_cache_box = bbox;
    if (m_nearest_neighbor)
      local_cache_box.expand(NearestPixelInterpolation::pixel_buffer); // for interpolation
    else
      local_cache_box.expand(BicubicInterpolation::pixel_buffer); // for interpolation

    m_cache.set_size(local_cache_box.width(), local_cache_box.height());
    vw::BBox2 out_box;
    for (int32 y=local_cache_box.min().y(); y<local_cache_box.max().y(); y++) {
      for (int32 x=local_cache_box.min().x(); x<local_cache_box.max().x(); x++) {
        Vector2 p = reverse(Vector2(x,y));
        m_cache(x - local_cache_box.min().x(), y - local_cache_box.min().y()) = p;
        if (p == m_invalid_pix) continue;
        if (bbox.contains(Vector2i(x, y))) out_box.grow(p);
      }
    }
    out_box = grow_bbox_to_int(out_box);

    // Must happen after all calls to reverse finished.
    m_img_cache_box = local_cache_box;

    if (m_nearest_neighbor)
      m_cache_interp_mask = interpolate(create_mask(m_cache, m_invalid_pix),
                                        NearestPixelInterpolation(), ZeroEdgeExtension());
    else
      m_cache_interp_mask = interpolate(create_mask(m_cache, m_invalid_pix),
                                        BicubicInterpolation(), ZeroEdgeExtension());

    // Need the check below as to not try to create images with
    // negative dimensions.
    if (out_box.empty())
      out_box = vw::BBox2i(0, 0, 0, 0);

    m_cached_rv_box = out_box;

    return m_cached_rv_box;
  }

  /// This applies the forward transformation to an entire bounding box of pixels.
  BBox2i Map2CamTrans::forward_bbox(BBox2i const& /*output_bbox*/) const {
    vw::vw_throw(vw::NoImplErr() << "forward_bbox() is not implemented for Map2CamTrans.");
    return BBox2i();
  }

  // Make a copy of Map2CamTrans. If later queuing a very dense number of pixels
  // in a tile, see stereo_tri.cc for how the transform values can be cached for
  // that tile. Caching can greatly slow things down, however, if only a sparse
  // number of pixels are queried.
  TransformPtr mapproj_trans_copy(TransformPtr trans) {
    Map2CamTrans* t_ptr = dynamic_cast<Map2CamTrans*>(trans.get());
    if (!t_ptr)
      vw_throw(vw::NoImplErr() << "Expecting a transform of type Map2CamTrans.");
    return TransformPtr(new Map2CamTrans(*t_ptr));
  }

  Datum2CamTrans::Datum2CamTrans(camera::CameraModel const* cam,
                                  GeoReference const& image_georef,
                                  GeoReference const& dem_georef,
                                  float dem_height,
                                  Vector2i const& image_size,
                                  bool call_from_mapproject,
                                  bool nearest_neighbor):
    m_cam(cam), m_image_georef(image_georef), m_dem_georef(dem_georef),
    m_dem_height(dem_height), m_image_size(image_size),
    m_call_from_mapproject(call_from_mapproject),
    m_nearest_neighbor(nearest_neighbor) {

    m_invalid_pix = camera::CameraModel::invalid_pixel();
  }

  Vector2 Datum2CamTrans::reverse(const Vector2 &p) const{

    Vector2 lonlat = m_image_georef.pixel_to_lonlat(p);
    Vector3 lonlatAlt(lonlat[0], lonlat[1], m_dem_height);
    Vector3 xyz = m_dem_georef.datum().geodetic_to_cartesian(lonlatAlt);

    int b = BicubicInterpolation::pixel_buffer;
    if (m_nearest_neighbor)
      b = NearestPixelInterpolation::pixel_buffer;
    Vector2 pt;
    try{
      pt = m_cam->point_to_pixel(xyz);
      if (m_call_from_mapproject &&
          (pt[0] < b - 1 || pt[0] >= m_image_size[0] - b ||
            pt[1] < b - 1 || pt[1] >= m_image_size[1] - b)) {
        // Won't be able to interpolate into image in transform(...)
        return m_invalid_pix;
      }
    } catch(...) { // If a point failed to project
      return m_invalid_pix;
    }

    return pt;
  }

  BBox2i Datum2CamTrans::reverse_bbox(BBox2i const& bbox) const {

    BBox2 out_box;
    for (int32 y=bbox.min().y(); y<bbox.max().y(); y++) {
      for (int32 x=bbox.min().x(); x<bbox.max().x(); x++) {

        Vector2 p = reverse(Vector2(x,y));
        if (p == m_invalid_pix)
          continue;
        out_box.grow(p);
      }
    }
    out_box = grow_bbox_to_int(out_box);

    // Need the check below as to not try to create images with negative dimensions.
    if (out_box.empty())
      out_box = BBox2i(0, 0, 0, 0);

    return out_box;
  }

}} // namespace vw::cartography
