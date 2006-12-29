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

/// \file CameraTransform.h
/// 
/// A functor and convenience functions for warping images between
/// different camera models.
/// 
#ifndef __VW_CAMERA_TRANSFORM_H__
#define __VW_CAMERA_TRANSFORM_H__

#include <vw/Image/Transform.h>

namespace vw {
namespace camera{

  /// This transform functor can be used along with the machinery in
  /// vw/Transform.h to warp an image from one camera's perspective
  /// into anothers.  In particular, this can be used to remove lens
  /// distortion by transforming from a nonlinear source camera to a
  /// linearized destination camera.  
  /// 
  /// NOTE: The pixel_to_vector/vector_to_pixel calling sequence only 
  /// works if both cameras have the some camera center (focal point).
  /// If the camera centers do not match it will vw_throw an exception.
  ///
  template <class SrcCameraT, class DstCameraT> 
  struct CameraTransform : public TransformBase<CameraTransform<SrcCameraT, DstCameraT> > {
    CameraTransform(SrcCameraT const& src_camera, 
                    DstCameraT const& dst_camera) : 
      m_src_camera(src_camera), m_dst_camera(dst_camera) {
    }

    /// This defines the transformation from coordinates in our target
    /// image back to coordinatess in the original image.
    inline Vector2 reverse(const Vector2 &p) const {
      VW_ASSERT(m_src_camera.camera_center(p) == m_dst_camera.camera_center(p),
                LogicErr() << "CameraTransformFunctor: Camera transformation require that the camera center is always the same for both cameras.");
      
      // (1) Call src PixelToVector to find the vector emanating from
      //     the camera center.
      Vector3 vec = m_dst_camera.pixel_to_vector(p);
      
      // (2) take resulting vector and call dest camera's
      //     VectorToPixel on it
      return m_src_camera.point_to_pixel(vec+m_src_camera.camera_center(p));
    }
    
    /// This defines the transformation from coordinates in our target
    /// image back to coordinatess in the original image.
    inline Vector2 forward(const Vector2 &p) const {
      VW_ASSERT(m_src_camera.camera_center(p) == m_dst_camera.camera_center(p),
                LogicErr() << "CameraTransformFunctor: Camera transformation require that the camera center is always the same for both cameras.");
      
      // (1) Call src PixelToVector to find the vector emanating from
      //     the camera center.
      Vector3 vec = m_src_camera.pixel_to_vector(p);
      
      // (2) take resulting vector and call dest camera's
      //     VectorToPixel on it.
      return m_dst_camera.vector_to_pixel(vec+m_dst_camera.camera_center(p));
    }
       
  private:
    SrcCameraT m_src_camera;
    DstCameraT m_dst_camera;
  };
  
}} // namespace vw::camera

#endif // __VW_CAMERA_TRANSFORM_H__
