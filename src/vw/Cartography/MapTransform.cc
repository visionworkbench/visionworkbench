// __BEGIN_LICENSE__
//  Copyright (c) 2009-2013, United States Government as represented by the
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
#include <vw/Cartography/MapTransform.h>
#include <vw/Image/MaskViews.h>
#include <vw/Camera/CameraModel.h>

namespace vw { namespace cartography {


  MapTransform::MapTransform( vw::camera::CameraModel const* cam,
                              GeoReference const& image_georef,
                              GeoReference const& dem_georef,
                              boost::shared_ptr<DiskImageResource> dem_rsrc,
                              vw::Vector2i image_size) :
    m_cam(cam), m_image_georef(image_georef), m_dem_georef(dem_georef),
    m_dem_rsrc(dem_rsrc), m_dem(dem_rsrc), m_image_size(image_size),
    m_has_nodata(false),
    m_nodata(std::numeric_limits<double>::quiet_NaN() ){
    m_has_nodata = dem_rsrc->has_nodata_read();
    if (m_has_nodata) m_nodata = dem_rsrc->nodata_read();

    if (m_image_size != Vector2(-1, -1))
      m_invalid_pix = Vector2(-1e6, -1e6); // for mapproject
    else
      m_invalid_pix = Vector2(-1, -1);     // for stereo
  }

  vw::Vector2
  MapTransform::reverse(const vw::Vector2 &p) const {

    if (m_img_cache_box.contains(p)){
      PixelMask<Vector2> v = m_cache_interp_mask(p.x() - m_img_cache_box.min().x(),
                                                 p.y() - m_img_cache_box.min().y());
      if (is_valid(v)) return v.child();
      else             return m_invalid_pix;
    }

    int b = BicubicInterpolation::pixel_buffer;
    Vector2 lonlat = m_image_georef.pixel_to_lonlat(p);
    Vector2 dem_pix = m_dem_georef.lonlat_to_pixel(lonlat);
    if (dem_pix[0] < b - 1 || dem_pix[0] >= m_dem.cols() - b ||
        dem_pix[1] < b - 1 || dem_pix[1] >= m_dem.rows() - b
        ){
      // No DEM data
      return m_invalid_pix;
    }

    Vector2 sdem_pix = dem_pix - m_dem_cache_box.min(); // since we cropped the DEM
    if (m_dem_cache_box.empty() ||
        sdem_pix[0] < b - 1 || sdem_pix[0] >= m_cropped_dem.cols() - b ||
        sdem_pix[1] < b - 1 || sdem_pix[1] >= m_cropped_dem.rows() - b
        ){
      // Cache miss. Will not happen often.
      BBox2i box;
      box.min() = floor(p) - Vector2(1, 1);
      box.max() = ceil(p)  + Vector2(1, 1);
      cache_dem(box);
      return reverse(p);
    }

    PixelMask<float> h = m_interp_dem(sdem_pix[0], sdem_pix[1]);
    if (!is_valid(h))
      return m_invalid_pix;

    Vector3 xyz = m_dem_georef.datum().geodetic_to_cartesian
      (Vector3(lonlat[0], lonlat[1], h.child()));
    Vector2 pt;
    try{
      pt = m_cam->point_to_pixel(xyz);
      if ( (m_image_size != Vector2i(-1, -1)) &&
           (pt[0] < b - 1 || pt[0] >= m_image_size[0] - b ||
            pt[1] < b - 1 || pt[1] >= m_image_size[1] - b)
           ){
        // Won't be able to interpolate into image in transform(...)
        return m_invalid_pix;
      }
    }catch(...){ // If a point failed to project
      return m_invalid_pix;
    }

    return pt;
  }

  void MapTransform::cache_dem(vw::BBox2i const& bbox) const{

    BBox2 dbox;
    dbox.grow( m_dem_georef.lonlat_to_pixel(m_image_georef.pixel_to_lonlat( Vector2(bbox.min().x(),bbox.min().y()) ) )); // Top left
    dbox.grow( m_dem_georef.lonlat_to_pixel(m_image_georef.pixel_to_lonlat( Vector2(bbox.max().x()-1,bbox.min().y()) ) )); // Top right
    dbox.grow( m_dem_georef.lonlat_to_pixel(m_image_georef.pixel_to_lonlat( Vector2(bbox.min().x(),bbox.max().y()-1) ) )); // Bottom left
    dbox.grow( m_dem_georef.lonlat_to_pixel(m_image_georef.pixel_to_lonlat( Vector2(bbox.max().x()-1,bbox.max().y()-1) ) )); // Bottom right

    // A lot of care is needed here when going from real box to int
    // box, and if in doubt, better expand more rather than less.
    dbox.expand(1);
    m_dem_cache_box = grow_bbox_to_int(dbox);
    m_dem_cache_box.expand(BicubicInterpolation::pixel_buffer); // for interp
    m_dem_cache_box.crop(bounding_box(m_dem));
    // Read the dem in memory for speed.
    m_cropped_dem = crop(m_dem, m_dem_cache_box);

    if (m_has_nodata){
      m_masked_dem = create_mask(m_cropped_dem, m_nodata);
    }else{
      m_masked_dem = pixel_cast< PixelMask<float> >(m_cropped_dem);
    }
    m_interp_dem =
      interpolate(m_masked_dem, BicubicInterpolation(), ZeroEdgeExtension());

  }

  // This function will be called whenever we start to apply the
  // transform in a tile. It computes and caches the point cloud at
  // each pixel in the tile, to be used later when we iterate over
  // pixels.
  vw::BBox2i
  MapTransform::reverse_bbox( vw::BBox2i const& bbox ) const {

    // Custom reverse_bbox() function which can handle invalid pixels.
    if (!m_cached_rv_box.empty()) return m_cached_rv_box;

    cache_dem(bbox);

    // Cache the reverse transform

    m_img_cache_box = BBox2i();
    BBox2i local_cache_box = bbox;
    local_cache_box.expand(BicubicInterpolation::pixel_buffer); // for interpolation
    m_cache.set_size(local_cache_box.width(), local_cache_box.height());
    vw::BBox2 out_box;
    for( int32 y=local_cache_box.min().y(); y<local_cache_box.max().y(); ++y ){
      for( int32 x=local_cache_box.min().x(); x<local_cache_box.max().x(); ++x ){
        Vector2 p = reverse( Vector2(x,y) );
        m_cache(x - local_cache_box.min().x(), y - local_cache_box.min().y()) = p;
        if (p == m_invalid_pix) continue;
        if (bbox.contains(Vector2i(x, y))) out_box.grow( p );
      }
    }
    out_box = grow_bbox_to_int( out_box );

    // Must happen after all calls to reverse finished.
    m_img_cache_box = local_cache_box;

    m_cache_interp_mask = interpolate(create_mask(m_cache, m_invalid_pix),
                                      BicubicInterpolation(), ZeroEdgeExtension());

    // Need the check below as to not try to create images with
    // negative dimensions.
    if (out_box.empty())
      out_box = vw::BBox2i(0, 0, 0, 0);

    m_cached_rv_box = out_box;
    return m_cached_rv_box;
  }

#if 0
  // Old version. Kept here for a while.
  MapTransform::MapTransform( vw::camera::CameraModel const* cam,
                              GeoReference const& image_georef,
                              GeoReference const& dem_georef,
                              boost::shared_ptr<DiskImageResource> dem_rsrc ) :
    m_cam(cam), m_image_georef(image_georef), m_dem_georef(dem_georef),
    m_dem(dem_rsrc) {

    if ( dem_rsrc->has_nodata_read() )
      m_point_cloud =
        geo_transform(
                      geodetic_to_cartesian(
                                            dem_to_geodetic( create_mask( m_dem, dem_rsrc->nodata_read()),
                                                             m_dem_georef ), m_dem_georef.datum() ),
                      m_dem_georef, m_image_georef,
                      ValueEdgeExtension<Vector3>( Vector3() ),
                      BicubicInterpolation());
    else
      m_point_cloud =
        geo_transform(
                      geodetic_to_cartesian(
                                            dem_to_geodetic( m_dem, m_dem_georef ),
                                            m_dem_georef.datum() ),
                      m_dem_georef, m_image_georef,
                      ValueEdgeExtension<Vector3>( Vector3() ),
                      BicubicInterpolation());
  }

  vw::Vector2
  MapTransform::reverse(const vw::Vector2 &p) const {

    double NaN = std::numeric_limits<double>::quiet_NaN();

    // If possible, we will interpolate into the cached point cloud
    ImageViewRef<Vector3> interp_point_cloud =
      interpolate(m_point_cloud_cache, BicubicInterpolation(),
                  ZeroEdgeExtension());

    // Avoid interpolating close to edges
    BBox2i shrank_bbox = m_cache_size;
    shrank_bbox.contract(BicubicInterpolation::pixel_buffer);

    Vector3 xyz =
      shrank_bbox.contains( p ) ?
      interp_point_cloud(p.x() - m_cache_size.min().x(),
                         p.y() - m_cache_size.min().y()):
      m_point_cloud(p.x(),p.y());

    // The case when xyz does not have data
    if (xyz == Vector3())
      return vw::Vector2(NaN, NaN);

    try{
      return m_cam->point_to_pixel(xyz);
    }catch(...){}
    return vw::Vector2(NaN, NaN);

  }

  // This function will be called whenever we start to apply the
  // transform in a tile. It computes and caches the point cloud at
  // each pixel in the tile, to be used later when we iterate over
  // pixels.
  vw::BBox2i
  MapTransform::reverse_bbox( vw::BBox2i const& bbox ) const {
    cache_dem( bbox );
    vw::BBox2i out_box = vw::TransformBase<MapTransform>::reverse_bbox( bbox );
    if (out_box.empty()) return vw::BBox2i(0, 0, 1, 1);
    return out_box;
  }

  void
  MapTransform::cache_dem( vw::BBox2i const& bbox ) const {
    m_cache_size = bbox;

    // For bicubic interpolation later
    m_cache_size.expand(BicubicInterpolation::pixel_buffer);

    m_point_cloud_cache = crop( m_point_cloud, m_cache_size );
  }
#endif

}} // namespace vw::cartography
