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

/// \file CameraTransform.h
/// 
/// A functor and convenience functions for warping images between
/// different camera models.
/// 
#ifndef __VW_CAMERA_TRANSFORM_H__
#define __VW_CAMERA_TRANSFORM_H__

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
  /// If the camera centers do not match, an exception will be thrown.
  ///
  template <class SrcCameraT, class DstCameraT> 
  struct CameraTransformFunctor {
    CameraTransformFunctor(SrcCameraT const& src_camera, 
                           DstCameraT const& dst_camera) : 
      m_src_camera(src_camera), m_dst_camera(dst_camera) {}
    
    /// This defines the transformation from coordinates in our target
    /// image back to coordinatess in the original image.
    inline Vector2 operator()(const Vector2 &p) const {
      VW_ASSERT(m_src_camera.camera_center(p) == m_dst_camera.camera_center(p),
                LogicErr() << "CameraTransformFunctor: Camera transformation require that the camera center is always the same for both cameras.");
      
      // (1) Call src PixelToVector to find the vector emanating from
      //     the camera center.
      Vector3 vec = m_dst_camera.pixel_to_vector(p);
      
      // (2) take resulting vector and call dest camera's
      //     VectorToPixel on it
      return m_src_camera.vector_to_pixel(vec);
    }
    
    /// This defines the transformation from coordinates in our target
    /// image back to coordinatess in the original image.
    inline Vector2 reverse(const Vector2 &p) const {
      VW_ASSERT(m_src_camera.camera_center(p) == m_dst_camera.camera_center(p),
                LogicErr() << "CameraTransformFunctor: Camera transformation require that the camera center is always the same for both cameras.");
      
      // (1) Call src PixelToVector to find the vector emanating from
      //     the camera center.
      Vector3 vec = m_src_camera.pixel_to_vector(p);
      
      // (2) take resulting vector and call dest camera's
      //     VectorToPixel on it
      return m_dst_camera.vector_to_pixel(vec);
    }
   
  private:
    SrcCameraT m_src_camera;
    DstCameraT m_dst_camera;
  };

  /// Apply a fixed linear camera transform to the image.  The
  /// transform is defined by two camera models.  The width and height
  /// of the output image can be specified.
  template <class ImageT, class SrcCameraT, class DstCameraT, class InterpT, class EdgeT>
  FixedTransformView<ImageT, CameraTransformFunctor<SrcCameraT, DstCameraT>, InterpT, EdgeT > 
  inline fixed_camera_transform( ImageViewBase<ImageT> const& v, SrcCameraT const& src_camera, DstCameraT const& dst_camera, 
                                 int width, int height, InterpT, EdgeT) {
    return FixedTransformView<ImageT, CameraTransformFunctor<SrcCameraT, DstCameraT>, InterpT, EdgeT >
      ( v.impl(), CameraTransformFunctor<SrcCameraT, DstCameraT>(src_camera, dst_camera), width, height );
  }

  // FixedTransform: User specifies source and destination camera models, output width/height, and Interpolation type.
  template <class ImageT, class SrcCameraT, class DstCameraT, class InterpT>
  FixedTransformView<ImageT, CameraTransformFunctor<SrcCameraT, DstCameraT>, InterpT > 
  inline fixed_camera_transform( ImageViewBase<ImageT> const& v, SrcCameraT const& src_camera, DstCameraT const& dst_camera, 
                                 int width, int height, InterpT) {
    return FixedTransformView<ImageT, CameraTransformFunctor<SrcCameraT, DstCameraT>, InterpT >
      ( v.impl(), CameraTransformFunctor<SrcCameraT, DstCameraT>(src_camera, dst_camera), width, height );
  }

  // FixedTransform: User specifies source and destination camera models and output width/height.
  template <class ImageT, class SrcCameraT, class DstCameraT>
  FixedTransformView<ImageT, CameraTransformFunctor<SrcCameraT, DstCameraT> > 
  inline fixed_camera_transform( ImageViewBase<ImageT> const& v, SrcCameraT const& src_camera, DstCameraT const& dst_camera, 
                                 int width, int height) {
    return FixedTransformView<ImageT, CameraTransformFunctor<SrcCameraT, DstCameraT> >
      ( v.impl(), CameraTransformFunctor<SrcCameraT, DstCameraT>(src_camera, dst_camera), width, height );
  }

  /// Apply a fixed camera transform to the image.  The transform is
  /// defined by a source and destination camera model.  The width and
  /// height of the output image is set equal to the dimensions of the
  /// input image.  When performining camera transformations, this is
  /// generally what you want, since the camera model itself makes
  /// some assumptions about the output image size.
  template <class ImageT, class SrcCameraT, class DstCameraT, class InterpT, class EdgeT>
  typename boost::disable_if<IsScalar<InterpT>, FixedTransformView<ImageT, CameraTransformFunctor<SrcCameraT, DstCameraT>, InterpT, EdgeT > >::type
  inline fixed_camera_transform( ImageViewBase<ImageT> const& v, SrcCameraT const& src_camera, DstCameraT const& dst_camera, InterpT, EdgeT) {
    return FixedTransformView<ImageT, CameraTransformFunctor<SrcCameraT, DstCameraT>, InterpT, EdgeT >
      ( v.impl(), CameraTransformFunctor<SrcCameraT, DstCameraT>(src_camera, dst_camera) );
  }

  template <class ImageT, class SrcCameraT, class DstCameraT, class InterpT>
  typename boost::disable_if<IsScalar<InterpT>, FixedTransformView<ImageT, CameraTransformFunctor<SrcCameraT, DstCameraT>, InterpT > >::type
  inline fixed_camera_transform( ImageViewBase<ImageT> const& v, SrcCameraT const& src_camera, DstCameraT const& dst_camera, InterpT) {
    return FixedTransformView<ImageT, CameraTransformFunctor<SrcCameraT, DstCameraT>, InterpT>
      ( v.impl(), CameraTransformFunctor<SrcCameraT, DstCameraT>(src_camera, dst_camera) );
  }

  template <class ImageT, class SrcCameraT, class DstCameraT>
  FixedTransformView<ImageT, CameraTransformFunctor<SrcCameraT, DstCameraT> >
  inline fixed_camera_transform( ImageViewBase<ImageT> const& v, SrcCameraT const& src_camera, DstCameraT const& dst_camera) {
    return FixedTransformView<ImageT, CameraTransformFunctor<SrcCameraT, DstCameraT> >
      ( v.impl(), CameraTransformFunctor<SrcCameraT, DstCameraT>(src_camera, dst_camera) );
  }

  /// Apply a free camera transform to an image.  The transform is
  /// defined by a source and destination camera model.
  template <class ImageT, class SrcCameraT, class DstCameraT, class InterpT, class EdgeT>
  typename boost::disable_if<IsScalar<InterpT>, FreeTransformView<ImageT, CameraTransformFunctor<SrcCameraT, DstCameraT>, InterpT, EdgeT > >::type
  inline free_camera_transform( ImageViewBase<ImageT> const& v, SrcCameraT const& src_camera, DstCameraT const& dst_camera, InterpT, EdgeT ) {
    return FreeTransformView<ImageT, CameraTransformFunctor<SrcCameraT, DstCameraT>, InterpT, EdgeT >
      ( v.impl(), CameraTransformFunctor<SrcCameraT, DstCameraT>(src_camera, dst_camera) );
  }

  template <class ImageT, class SrcCameraT, class DstCameraT, class InterpT >
  typename boost::disable_if<IsScalar<InterpT>, FreeTransformView<ImageT, CameraTransformFunctor<SrcCameraT, DstCameraT>, InterpT > >::type
  inline free_camera_transform( ImageViewBase<ImageT> const& v, SrcCameraT const& src_camera, DstCameraT const& dst_camera, InterpT) {
    return FreeTransformView<ImageT, CameraTransformFunctor<SrcCameraT, DstCameraT>, InterpT>
      ( v.impl(), CameraTransformFunctor<SrcCameraT, DstCameraT>(src_camera, dst_camera) );
  }

  template <class ImageT, class SrcCameraT, class DstCameraT>
  FreeTransformView<ImageT, CameraTransformFunctor<SrcCameraT, DstCameraT> > 
  inline free_camera_transform( ImageViewBase<ImageT> const& v, SrcCameraT const& src_camera, DstCameraT const& dst_camera) {
    return FreeTransformView<ImageT, CameraTransformFunctor<SrcCameraT, DstCameraT> >
      ( v.impl(), CameraTransformFunctor<SrcCameraT, DstCameraT>(src_camera, dst_camera) );
  }

  /// Apply a linear free transform to an image.  The transform is
  /// defined by a 3x3 homography.
  ///
  /// The shift of the origin due to the automatic alignment process is
  /// returned in (di, dj).
  template <class ImageT, class SrcCameraT, class DstCameraT, class InterpT, class EdgeT>
  FreeTransformView<ImageT, CameraTransformFunctor<SrcCameraT, DstCameraT>, InterpT, EdgeT > 
  inline free_camera_transform( ImageViewBase<ImageT> const& v, SrcCameraT const& src_camera, DstCameraT const& dst_camera, 
                                double &di, double &dj, InterpT, EdgeT) {
    return FreeTransformView<ImageT, CameraTransformFunctor<SrcCameraT, DstCameraT>, InterpT, EdgeT >
      ( v.impl(), CameraTransformFunctor<SrcCameraT, DstCameraT>(src_camera, dst_camera), di, dj );
  }

  template <class ImageT, class SrcCameraT, class DstCameraT, class InterpT>
  FreeTransformView<ImageT, CameraTransformFunctor<SrcCameraT, DstCameraT>, InterpT> 
  inline free_camera_transform( ImageViewBase<ImageT> const& v, SrcCameraT const& src_camera, DstCameraT const& dst_camera, 
                                double &di, double &dj, InterpT) {
    return FreeTransformView<ImageT, CameraTransformFunctor<SrcCameraT, DstCameraT>, InterpT >
      ( v.impl(), CameraTransformFunctor<SrcCameraT, DstCameraT>(src_camera, dst_camera), di, dj );
  }

  template <class ImageT, class SrcCameraT, class DstCameraT>
  FreeTransformView<ImageT, CameraTransformFunctor<SrcCameraT, DstCameraT> > 
  inline free_camera_transform( ImageViewBase<ImageT> const& v, SrcCameraT const& src_camera, DstCameraT const& dst_camera, 
                                double &di, double &dj) {
    return FreeTransformView<ImageT, CameraTransformFunctor<SrcCameraT, DstCameraT> >
      ( v.impl(), CameraTransformFunctor<SrcCameraT, DstCameraT>(src_camera, dst_camera), di, dj );
  }

}} // namespace vw::camera

#endif // __VW_CAMERA_TRANSFORM_H__
