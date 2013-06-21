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

namespace vw { namespace camera{
  class CameraModel; // forward declaration
}}

namespace vw {
namespace cartography {

  // MapTransform. Used to test the validity of IP matching on map
  // projected images. However, this could be used for performing an RPC
  // map projection.
  class MapTransform : public vw::TransformBase<MapTransform> {
    vw::camera::CameraModel const* m_cam;
    GeoReference m_image_georef, m_dem_georef;
    vw::DiskImageView<float> m_dem;
    vw::ImageViewRef<vw::Vector3> m_point_cloud;

    // We will always be modifying these
    mutable vw::ImageView<vw::Vector3> m_point_cloud_cache;
    mutable vw::BBox2i m_cache_size;
  public:
    MapTransform( vw::camera::CameraModel const* cam,
                     GeoReference const& image_georef,
                     GeoReference const& dem_georef,
                     boost::shared_ptr<vw::DiskImageResource> dem_rsrc );

    // Convert Map Projected Coordinate to camera coordinate
    vw::Vector2 reverse(const vw::Vector2 &p) const;

    vw::BBox2i reverse_bbox( vw::BBox2i const& bbox ) const;

    // Not thread safe ... you must copy this object
    void cache_dem( vw::BBox2i const& bbox ) const;
  };

  // MapTransform2. Used to test the validity of IP matching on map
  // projected images. However, this could be used for performing an
  // RPC map projection.
  class MapTransform2 : public vw::TransformBase<MapTransform2> {
    vw::camera::CameraModel const* m_cam;
    GeoReference m_image_georef, m_dem_georef;
    vw::DiskImageView<float> m_dem;
    vw::Vector2i m_image_size;
    bool m_has_nodata;
    double m_nodata;
    Vector2 m_invalid_pix;

    // We will always be modifying these
    mutable vw::BBox2i m_dem_box;
    mutable vw::ImageView<float> m_cropped_dem;
    mutable ImageViewRef< PixelMask<float> > m_masked_dem;
    mutable ImageViewRef< PixelMask <float> > m_interp_dem;
    mutable ImageView<Vector2> m_cache;
    mutable ImageViewRef< PixelMask<Vector2> > m_cache_interp_mask;
    mutable vw::BBox2i m_cache_box;
    mutable vw::BBox2i m_cached_rv_box;

  public:
    MapTransform2( vw::camera::CameraModel const* cam,
                   GeoReference const& image_georef,
                   GeoReference const& dem_georef,
                   boost::shared_ptr<vw::DiskImageResource> dem_rsrc,
                   vw::Vector2i image_size = vw::Vector2i(-1, -1)
                   );

    // Convert Map Projected Coordinate to camera coordinate
    vw::Vector2 reverse(const vw::Vector2 &p) const;

    // Not thread safe ... you must copy this object
    vw::BBox2i reverse_bbox( vw::BBox2i const& bbox ) const;
  };

}} // namespace vw::cartography

#endif//__VW_CARTOGRAPHY_MAP_TRANSFORM_H__
