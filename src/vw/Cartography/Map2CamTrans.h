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


#ifndef __VW_CARTOGRAPHY_MAP_TRANSFORM_H__
#define __VW_CARTOGRAPHY_MAP_TRANSFORM_H__

#include <vw/Image/ImageViewRef.h>
#include <vw/Image/Transform.h>
#include <vw/FileIO/DiskImageView.h>
#include <vw/Cartography/GeoReference.h>


/// \file Map2CamTrans.h
/// Given a pixel in a map-projected image, convert it to lonlat, then
/// convert to the DEM pixel, then to the DEM lonlat, then to the DEM
/// xyz, then project into the camera, and find the camera pixel.

/// This class is not thread-safe. For performance reasons it better not
/// be invoked for individual pixels. The best approach is to first
/// create a new instance of this class for each tile, do one wholesale
/// caching computation on the entire tile, then invoke it repeatedly
/// individual pixels in the tile.

/// The class can handle DEMs with holes.

namespace vw { namespace camera{
  class CameraModel; // forward declaration
}}

namespace vw { namespace cartography {

  class VW_API Map2CamTrans : public vw::TransformBase<Map2CamTrans> {
    vw::camera::CameraModel const* m_cam;
    GeoReference m_image_georef, m_dem_georef;
    vw::DiskImageView<float> m_dem;
    vw::Vector2i m_image_size;
    bool         m_call_from_mapproject;
    bool         m_has_nodata;
    double       m_nodata;
    Vector2      m_invalid_pix;


    // We will always be modifying these
    mutable vw::BBox2i                         m_dem_cache_box;
    mutable vw::ImageView<float>               m_cropped_dem;
    mutable ImageViewRef< PixelMask<float> >   m_masked_dem;
    mutable ImageViewRef< PixelMask<float> >   m_interp_dem;
    mutable ImageView<Vector2>                 m_cache;
    mutable ImageViewRef< PixelMask<Vector2> > m_cache_interp_mask;
    mutable vw::BBox2i                         m_img_cache_box;
    mutable vw::BBox2i                         m_cached_rv_box;

  public:
    Map2CamTrans( vw::camera::CameraModel const* cam,
                  GeoReference const& image_georef,
                  GeoReference const& dem_georef,
                  std::string  const& dem_file,
                  vw::Vector2i const& image_size,
                  bool call_from_mapproject
                  );

    /// Convert Map Projected Coordinate to camera coordinate
    vw::Vector2 reverse(const vw::Vector2 &p) const;

    // Not thread safe ... you must copy this object
    void       cache_dem   ( vw::BBox2i const& bbox ) const;
    vw::BBox2i reverse_bbox( vw::BBox2i const& bbox ) const;
  };
  
  VW_API std::ostream& operator<<(std::ostream& os, const Map2CamTrans& trans);


}} // namespace vw::cartography

#endif//__VW_CARTOGRAPHY_MAP_TRANSFORM_H__
