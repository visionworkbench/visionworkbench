// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_CARTOGRAPHY_CAMERABBOX_H__
#define __VW_CARTOGRAPHY_CAMERABBOX_H__

#include <vw/config.h>
#if defined(VW_HAVE_PKG_CAMERA) && (VW_HAVE_PKG_CAMERA==1)

#include <vw/Cartography/GeoReference.h>
#include <vw/Camera/CameraModel.h>

#include <boost/shared_ptr.hpp>

namespace vw {
namespace cartography {

  BBox2 camera_bbox( GeoReference const& georef, 
                     boost::shared_ptr<vw::camera::CameraModel> camera_model, 
                     int32 cols, int32 rows, float &scale );

  inline BBox2 camera_bbox( GeoReference const& georef, boost::shared_ptr<vw::camera::CameraModel> camera_model, int32 cols, int32 rows ) {
    float scale;
    return camera_bbox( georef, camera_model, cols, rows, scale );
  }

} // namespace cartography
} // namespace vw

#endif // VW_HAVE_PKG_CAMERA

#endif // __VW_CARTOGRAPHY_CAMERABBOX_H__
