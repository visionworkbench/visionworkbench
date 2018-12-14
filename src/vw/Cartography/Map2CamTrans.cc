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
#include <vw/Cartography/Map2CamTrans.h>
#include <vw/Image/MaskViews.h>
#include <vw/Camera/CameraModel.h>

namespace vw { namespace cartography {

  Map2CamTrans::Map2CamTrans( vw::camera::CameraModel const* cam,
                              GeoReference const& image_georef,
                              GeoReference const& dem_georef,
                              std::string const& dem_file,
                              vw::Vector2i const& image_size,
                              bool call_from_mapproject,
                              bool nearest_neighbor):
    m_cam(cam), m_image_georef(image_georef), m_dem_georef(dem_georef),
    m_dem(dem_file), m_image_size(image_size),
    m_call_from_mapproject(call_from_mapproject), 
    m_nearest_neighbor(nearest_neighbor), m_has_nodata(false),
    m_nodata(std::numeric_limits<double>::quiet_NaN()){

    boost::shared_ptr<vw::DiskImageResource>
      dem_rsrc( vw::DiskImageResourcePtr(dem_file) );

    m_has_nodata = dem_rsrc->has_nodata_read();
    if (m_has_nodata) m_nodata = dem_rsrc->nodata_read();

    m_invalid_pix = vw::camera::CameraModel::invalid_pixel();
  }

  vw::Vector2
  Map2CamTrans::reverse(const vw::Vector2 &p) const {

    // If we have data for the location already cached
    if (m_img_cache_box.contains(p)){
      // Interpolate the output value using the cached data
      PixelMask<Vector2> v = m_cache_interp_mask(p.x() - m_img_cache_box.min().x(),
                                                 p.y() - m_img_cache_box.min().y());
      // We can just return the value if it is valid!
      if (is_valid(v)) return v.child();
      else             return m_invalid_pix;
    }

    int b = BicubicInterpolation::pixel_buffer;
    if (m_nearest_neighbor)
      b = NearestPixelInterpolation::pixel_buffer;
    Vector2 lonlat  = m_image_georef.pixel_to_lonlat(p);
    Vector2 dem_pix = m_dem_georef.lonlat_to_pixel(lonlat);
    if ((dem_pix[0] < b - 1) || (dem_pix[0] >= m_dem.cols() - b) ||
        (dem_pix[1] < b - 1) || (dem_pix[1] >= m_dem.rows() - b)
        ){
      // No DEM data
      return m_invalid_pix;
    }

    Vector2 sdem_pix = dem_pix - m_dem_cache_box.min(); // since we cropped the DEM
    if (m_dem_cache_box.empty() ||
        (sdem_pix[0] < b - 1) || (sdem_pix[0] >= m_cropped_dem.cols() - b) ||
        (sdem_pix[1] < b - 1) || (sdem_pix[1] >= m_cropped_dem.rows() - b)
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
      if ( m_call_from_mapproject &&
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

  void Map2CamTrans::cache_dem(vw::BBox2i const& bbox) const{

    // TODO: This may fail around poles. Need to do the standard X trick, traverse
    // the edges and diagonals of the box.
    BBox2 dbox;
    dbox.grow( m_dem_georef.lonlat_to_pixel(m_image_georef.pixel_to_lonlat( Vector2(bbox.min().x(),   bbox.min().y()  ) ) )); // Top left
    dbox.grow( m_dem_georef.lonlat_to_pixel(m_image_georef.pixel_to_lonlat( Vector2(bbox.max().x()-1, bbox.min().y()  ) ) )); // Top right
    dbox.grow( m_dem_georef.lonlat_to_pixel(m_image_georef.pixel_to_lonlat( Vector2(bbox.min().x(),   bbox.max().y()-1) ) )); // Bottom left
    dbox.grow( m_dem_georef.lonlat_to_pixel(m_image_georef.pixel_to_lonlat( Vector2(bbox.max().x()-1, bbox.max().y()-1) ) )); // Bottom right

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

    if (m_has_nodata){
      m_masked_dem = create_mask(m_cropped_dem, m_nodata);
    }else{ // Don't need to handle nodata
      m_masked_dem = pixel_cast< PixelMask<float> >(m_cropped_dem);
    }
    // Set up interpolation interface to the data we loaded into memory
    if (m_nearest_neighbor)
      m_interp_dem = interpolate(m_masked_dem, NearestPixelInterpolation(), ZeroEdgeExtension());
    else
      m_interp_dem = interpolate(m_masked_dem, BicubicInterpolation(), ZeroEdgeExtension());

  } // End function cache_dem

  // This function will be called whenever we start to apply the
  // transform in a tile. It computes and caches the point cloud at
  // each pixel in the tile, to be used later when we iterate over pixels.
  vw::BBox2i
  Map2CamTrans::reverse_bbox( vw::BBox2i const& bbox ) const {

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
/*
  std::ostream& operator<<( std::ostream& os, Map2CamTrans const& trans ) {
    std::ostringstream oss; // To use custom precision
    oss.precision(17);
    oss << "TODO: Fill in Map2CamTrans printer!\n";
    os << oss.str();
    return os;
  }
*/
//=======================================================================
  

  Datum2CamTrans::Datum2CamTrans( camera::CameraModel const* cam,
                                  GeoReference const& image_georef,
                                  GeoReference const& dem_georef,
                                  float dem_height,
                                  Vector2i const& image_size,
                                  bool call_from_mapproject,
                                  bool nearest_neighbor):
    m_cam(cam), m_image_georef(image_georef), m_dem_georef(dem_georef),
    m_dem_height(dem_height), m_image_size(image_size),
    m_call_from_mapproject(call_from_mapproject),
    m_nearest_neighbor(nearest_neighbor){

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
      if ( m_call_from_mapproject &&
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

  BBox2i Datum2CamTrans::reverse_bbox( BBox2i const& bbox ) const {

    BBox2 out_box;      
    for( int32 y=bbox.min().y(); y<bbox.max().y(); ++y ){
      for( int32 x=bbox.min().x(); x<bbox.max().x(); ++x ){
      
        Vector2 p = reverse( Vector2(x,y) );
        if (p == m_invalid_pix) 
          continue;
        out_box.grow( p );
      }
    }
    out_box = grow_bbox_to_int( out_box );

    // Need the check below as to not try to create images with negative dimensions.
    if (out_box.empty())
      out_box = BBox2i(0, 0, 0, 0);

    return out_box;
  }


  
  
}} // namespace vw::cartography
