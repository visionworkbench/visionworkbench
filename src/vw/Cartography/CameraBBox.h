// __BEGIN_LICENSE__
// 
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
// 
// Copyright 2006 Carnegie Mellon University. All rights reserved.
// 
// This software is distributed under the NASA Open Source Agreement
// (NOSA), version 1.3.  The NOSA has been approved by the Open Source
// Initiative.  See the file COPYING at the top of the distribution
// directory tree for the complete NOSA document.
// 
// THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
// KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
// LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
// SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
// A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
// THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
// DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
// 
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

  BBox2 camera_bbox( GeoReference const& georef, float min_alt, float max_alt,
                     boost::shared_ptr<vw::camera::CameraModel> camera_model, int32 cols, int32 rows, float &scale );

  inline BBox2 camera_bbox( GeoReference const& georef, boost::shared_ptr<vw::camera::CameraModel> camera_model, int32 cols, int32 rows ) {
    float scale;
    return camera_bbox( georef, 0.0, 0.0, camera_model, cols, rows, scale );
  }

} // namespace cartography
} // namespace vw

#endif // VW_HAVE_PKG_CAMERA

#endif // __VW_CARTOGRAPHY_CAMERABBOX_H__
