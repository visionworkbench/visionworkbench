// __BEGIN_LICENSE__
//
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
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

/// \file CameraModel.h
/// 
/// An abstract base class referring to a camera model class.
/// 
#ifndef __VW_CAMERA_CAMERAMODEL_H__
#define __VW_CAMERA_CAMERAMODEL_H__

#include <vw/Core/Exception.h>
#include <vw/Math/Vector.h>
#include <vw/Image/ImageView.h>

namespace vw { 
namespace camera {

  /// This is the abstract base class for a camera model object.  You
  /// should create a child that adheres to the interface set forth
  /// below.
  class CameraModel {

  public:

    virtual ~CameraModel() {}

    //------------------------------------------------------------------
    // Generic Camera Model Interface
    //------------------------------------------------------------------
    virtual Vector2 vector_to_pixel(Vector3 const& vec) const = 0;
    virtual Vector3 pixel_to_vector(Vector2 const& pix) const = 0;
    virtual Vector3 camera_center(Vector2 const& pix) const = 0;

  };  

}}	// namespace vw::camera

#endif // __VW_CAMERA_CAMERAMODEL_H__
