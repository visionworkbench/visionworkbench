// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NASA Vision Workbench is licensed under the Apache License,
//  Version 2.0 (the "License"); you may not use this file except in
//  compliance with the License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__


/// \file CameraTransform.h
///
/// A functor and convenience functions for warping images between
/// different camera models.
///
#ifndef __VW_CAMERA_TRANSFORM_H__
#define __VW_CAMERA_TRANSFORM_H__

#include <vw/Camera/CameraModel.h>
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
  /// works if both cameras have the same camera center (focal point).
  /// If the camera centers do not match it will vw_throw an exception.
  ///
  template <class SrcCameraT, class DstCameraT>
  struct CameraTransform : public TransformBase<CameraTransform<SrcCameraT, DstCameraT> > {
    CameraTransform(SrcCameraT const& src_camera,
                    DstCameraT const& dst_camera) :
      m_src_camera(src_camera), m_dst_camera(dst_camera) {
    }

    /// This defines the transformation from coordinates in our target
    /// image back to coordinates in the original image.
    inline Vector2 reverse(const Vector2 &p) const {
      VW_ASSERT(m_src_camera.camera_center(p) == m_dst_camera.camera_center(p),
                LogicErr() << "CameraTransformFunctor: Camera transformation require that the camera center is always the same for both cameras.");

      // (1) Call src PixelToVector to find the vector emanating from the camera center.
      Vector3 vec = m_dst_camera.pixel_to_vector(p);

      // (2) take resulting vector and call dest camera's VectorToPixel on it
      return m_src_camera.point_to_pixel(vec+m_dst_camera.camera_center(p));
    }

    /// This defines the transformation from coordinates in our source
    /// image to coordinatess in the target image.
    inline Vector2 forward(const Vector2 &p) const {
      VW_ASSERT(m_src_camera.camera_center(p) == m_dst_camera.camera_center(p),
                LogicErr() << "CameraTransformFunctor: Camera transformation require that the camera center is always the same for both cameras.");

      // (1) Call src PixelToVector to find the vector emanating from the camera center.
      Vector3 vec = m_src_camera.pixel_to_vector(p);

      // (2) take resulting vector and call dest camera's VectorToPixel on it.
      return m_dst_camera.point_to_pixel(vec+m_src_camera.camera_center(p));
    }

  private:
    SrcCameraT m_src_camera;
    DstCameraT m_dst_camera;
  };

  struct SmartPtrCameraTransform : public TransformBase<SmartPtrCameraTransform> {
    SmartPtrCameraTransform(boost::shared_ptr<CameraModel> const& src_camera,
                            boost::shared_ptr<CameraModel> const& dst_camera) :
      m_src_camera(src_camera), m_dst_camera(dst_camera) {
    }

    /// This defines the transformation from coordinates in our target
    /// image back to coordinatess in the original image.
    inline Vector2 reverse(const Vector2 &p) const {
      VW_ASSERT(m_src_camera->camera_center(p) == m_dst_camera->camera_center(p),
                LogicErr() << "CameraTransformFunctor: Camera transformation require that the camera center is always the same for both cameras.");

      // (1) Call src PixelToVector to find the vector emanating from
      //     the camera center.
      Vector3 vec = m_dst_camera->pixel_to_vector(p);

      // (2) take resulting vector and call dest camera's
      //     VectorToPixel on it
      return m_src_camera->point_to_pixel(vec+m_dst_camera->camera_center(p));
    }

    /// This defines the transformation from coordinates in our source
    /// image to coordinatess in the target image.
    inline Vector2 forward(const Vector2 &p) const {
      VW_ASSERT(m_src_camera->camera_center(p) == m_dst_camera->camera_center(p),
                LogicErr() << "CameraTransformFunctor: Camera transformation require that the camera center is always the same for both cameras.");

      // (1) Call src PixelToVector to find the vector emanating from
      //     the camera center.
      Vector3 vec = m_src_camera->pixel_to_vector(p);

      // (2) take resulting vector and call dest camera's
      //     VectorToPixel on it.
      return m_dst_camera->point_to_pixel(vec+m_src_camera->camera_center(p));
    }

  private:
    boost::shared_ptr<CameraModel> m_src_camera;
    boost::shared_ptr<CameraModel> m_dst_camera;
  };


  /// Transform an image from one camera model to another, explicitly
  /// specifying the edge extension and interpolation modes. Specify the size of the output image.
  template <class ImageT, class SrcCameraT, class DstCameraT, class EdgeT, class InterpT>
  TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, InterpT>, CameraTransform<SrcCameraT, DstCameraT> >
  inline camera_transform( ImageViewBase<ImageT> const& image, SrcCameraT const& src_camera, DstCameraT const& dst_camera, Vector2i size, EdgeT const& edge_func, InterpT const& interp_func )
  {
    CameraTransform<SrcCameraT, DstCameraT> ctx( src_camera, dst_camera );
    return transform( image, ctx, size[0], size[1], edge_func, interp_func );
  }

  /// Transform an image from one camera model to another, explicitly
  /// specifying the edge extension and interpolation modes.
  template <class ImageT, class SrcCameraT, class DstCameraT, class EdgeT, class InterpT>
  TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, InterpT>, CameraTransform<SrcCameraT, DstCameraT> >
  inline camera_transform( ImageViewBase<ImageT> const& image, SrcCameraT const& src_camera, DstCameraT const& dst_camera, EdgeT const& edge_func, InterpT const& interp_func )
  {
    CameraTransform<SrcCameraT, DstCameraT> ctx( src_camera, dst_camera );
    return transform( image, ctx, edge_func, interp_func );
  }

  /// Transform an image from one camera model to another using
  /// bilinear interpolation, explicitly specifying the edge extension mode.
  template <class ImageT, class SrcCameraT, class DstCameraT, class EdgeT>
  TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, BilinearInterpolation>, CameraTransform<SrcCameraT, DstCameraT> >
  inline camera_transform( ImageViewBase<ImageT> const& image, SrcCameraT const& src_camera, DstCameraT const& dst_camera, EdgeT const& edge_func )
  {
    CameraTransform<SrcCameraT, DstCameraT> ctx( src_camera, dst_camera );
    return transform( image, ctx, edge_func, BilinearInterpolation() );
  }

  /// Transform an image from one camera model to another, using
  /// zero (black) edge-extension and bilinear interpolation.
  template <class ImageT, class SrcCameraT, class DstCameraT>
  TransformView<InterpolationView<EdgeExtensionView<ImageT, ZeroEdgeExtension>, BilinearInterpolation>, CameraTransform<SrcCameraT, DstCameraT> >
  inline camera_transform( ImageViewBase<ImageT> const& image, SrcCameraT const& src_camera, DstCameraT const& dst_camera )
  {
    CameraTransform<SrcCameraT, DstCameraT> ctx( src_camera, dst_camera );
    return transform( image, ctx, ZeroEdgeExtension(), BilinearInterpolation() );
  }

  /// Transform an image from one camera model to another, using zero (black) edge-extension 
  /// and bilinear interpolation.  Specify the size of the output image.
  template <class ImageT, class SrcCameraT, class DstCameraT>
  TransformView<InterpolationView<EdgeExtensionView<ImageT, ZeroEdgeExtension>, BilinearInterpolation>, CameraTransform<SrcCameraT, DstCameraT> >
  inline camera_transform( ImageViewBase<ImageT> const& image, SrcCameraT const& src_camera, DstCameraT const& dst_camera,
                           Vector2i size)
  {
    CameraTransform<SrcCameraT, DstCameraT> ctx( src_camera, dst_camera );
    return transform( image, ctx, size[0], size[1], ZeroEdgeExtension(), BilinearInterpolation() );
  }

  /// Transform an image from one camera model to another, using zero (black) edge-extension 
  /// and bilinear interpolation.  Specify the bbox of the output image.
  template <class ImageT, class SrcCameraT, class DstCameraT>
  CropView<TransformView<InterpolationView<EdgeExtensionView<ImageT, ZeroEdgeExtension>, BilinearInterpolation>, CameraTransform<SrcCameraT, DstCameraT> > >
  inline camera_transform( ImageViewBase<ImageT> const& image, SrcCameraT const& src_camera, DstCameraT const& dst_camera,
                           BBox2 const &bbox)
  {
    CameraTransform<SrcCameraT, DstCameraT> ctx( src_camera, dst_camera );
    return transform( image, ctx, bbox, ZeroEdgeExtension(), BilinearInterpolation() );
  }

  /// Transform an image from a camera model to a linearized
  /// (i.e. undistorted) version of itself, explicitly specifying the
  /// edge extension and interpolation modes.
  template <class ImageT, class SrcCameraT, class EdgeT, class InterpT>
  TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, InterpT>, CameraTransform<SrcCameraT, typename SrcCameraT::linearized_type> >
  inline linearize_camera_transform( ImageViewBase<ImageT> const& image, SrcCameraT const& src_camera, EdgeT const& edge_func, InterpT const& interp_func )
  {
    typename SrcCameraT::linearized_type dst_camera =
      linearize_camera( src_camera, Vector2i(image.impl().cols(), image.impl().rows()),
                        Vector2i(image.impl().cols(), image.impl().rows()) );
    CameraTransform<SrcCameraT, typename SrcCameraT::linearized_type> ctx( src_camera, dst_camera );
    return transform( image, ctx, edge_func, interp_func );
  }

  /// Transform an image from a camera model to a linearized
  /// (i.e. undistorted) version of itself using bilinear
  /// interpolation, explicitly specifying the edge extension mode.
  template <class ImageT, class SrcCameraT, class EdgeT>
  TransformView<InterpolationView<EdgeExtensionView<ImageT, EdgeT>, BilinearInterpolation>, CameraTransform<SrcCameraT, typename SrcCameraT::linearized_type> >
  inline linearize_camera_transform( ImageViewBase<ImageT> const& image, SrcCameraT const& src_camera, EdgeT const& edge_func )
  {
    typename SrcCameraT::linearized_type dst_camera =
      linearize_camera( src_camera, Vector2i(image.impl().cols(), image.impl().rows()),
                        Vector2i(image.impl().cols(), image.impl().rows()) );
    CameraTransform<SrcCameraT, typename SrcCameraT::linearized_type> ctx( src_camera, dst_camera );
    return transform( image, ctx, edge_func, BilinearInterpolation() );
  }

  /// Transform an image from a camera model to a linearized
  /// (i.e. undistorted) version of itself, using zero (black)
  /// edge-extension and bilinear interpolation.
  template <class ImageT, class SrcCameraT>
  TransformView<InterpolationView<EdgeExtensionView<ImageT, ZeroEdgeExtension>, BilinearInterpolation>, CameraTransform<SrcCameraT, typename SrcCameraT::linearized_type> >
  inline linearize_camera_transform( ImageViewBase<ImageT> const& image, SrcCameraT const& src_camera )
  {
    typename SrcCameraT::linearized_type dst_camera =
      linearize_camera( src_camera, Vector2i(image.impl().cols(), image.impl().rows()),
                        Vector2i(image.impl().cols(), image.impl().rows()) );
    CameraTransform<SrcCameraT, typename SrcCameraT::linearized_type> ctx( src_camera, dst_camera );
    return transform( image, ctx, ZeroEdgeExtension(), BilinearInterpolation() );
  }



}} // namespace vw::camera

#endif // __VW_CAMERA_TRANSFORM_H__
