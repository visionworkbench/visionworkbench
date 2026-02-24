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

/// \file Algorithms.h
///
/// Basic algorithms operating on images. This includes only lazy view
/// implementations. See Manipulation.h for channel_cast(), PixelMath.h for
/// round(). etc.
///
#ifndef __VW_IMAGE_ALGORITHMS_H__
#define __VW_IMAGE_ALGORITHMS_H__

#include <vw/Image/AlgorithmFunctions.h>
#include <vw/Image/PerPixelViews.h>
#include <vw/Image/Statistics.h>
#include <vw/Image/EdgeExtension.h>
#include <vw/Image/PerPixelAccessorViews.h>

namespace vw {

  // *******************************************************************
  // clamp()
  // *******************************************************************

  template <class PixelT>
  class ChannelClampFunctor;

  /// Clamp the values in an image to fall within the range [low,high].
  template <class ImageT, class LowT, class HighT>
  UnaryPerPixelView<ImageT,UnaryCompoundFunctor<ChannelClampFunctor<typename ImageT::pixel_type>, typename ImageT::pixel_type > >
  inline clamp( ImageViewBase<ImageT> const& image, LowT low, HighT high );

  /// Clamp the values in an image to fall within the range [0,high].
  /// The low end of the range is actually determined by the
  /// ChannelRange type trait but is generally zero.
  template <class ImageT, class HighT>
  UnaryPerPixelView<ImageT,UnaryCompoundFunctor<ChannelClampFunctor<typename ImageT::pixel_type>, typename ImageT::pixel_type> >
  inline clamp( ImageViewBase<ImageT> const& image, HighT high );

  /// Clamp the values in an image to fall within the range [min,max],
  /// where min and max are determined by the ChannelRange type trait
  /// and are generally equal to 0.0 and 1.0 for floating point types
  /// and 0 and the largest positive value for integral types.
  template <class ImageT>
  UnaryPerPixelView<ImageT,UnaryCompoundFunctor<ChannelClampFunctor<typename ImageT::pixel_type>, typename ImageT::pixel_type> >
  inline clamp( ImageViewBase<ImageT> const& image );

  // *******************************************************************
  // normalize()
  // *******************************************************************

  /// \cond INTERNAL
  template <class PixelT>
  class ChannelNormalizeFunctor;
  /// \endcond

  template <class PixelT>
  class ChannelNormalizeRetainAlphaFunctor;

  /// Renormalize the values in an image to fall within the range
  /// [low,high), but leave the values in the alpha channel untouched.
  template <class ImageT>
  UnaryPerPixelView<ImageT, ChannelNormalizeRetainAlphaFunctor<typename ImageT::pixel_type> >
  inline normalize_retain_alpha( ImageViewBase<ImageT> const& image,
                                 typename ImageChannelType<ImageT>::type old_low,
                                 typename ImageChannelType<ImageT>::type old_high,
                                 typename ImageChannelType<ImageT>::type new_low,
                                 typename ImageChannelType<ImageT>::type new_high  );

  /// Renormalize the values in an image to fall within the range [low,high).
  template <class ImageT>
  UnaryPerPixelView<ImageT, UnaryCompoundFunctor<ChannelNormalizeFunctor<typename ImageT::pixel_type>, typename ImageT::pixel_type> >
  inline normalize( ImageViewBase<ImageT> const& image,
                    typename ImageChannelType<ImageT>::type old_low,
                    typename ImageChannelType<ImageT>::type old_high,
                    typename ImageChannelType<ImageT>::type new_low,
                    typename ImageChannelType<ImageT>::type new_high  );

  /// Renormalize the values in an image to fall within the range [low,high).
  template <class ImageT>
  UnaryPerPixelView<ImageT, UnaryCompoundFunctor<ChannelNormalizeFunctor<typename ImageT::pixel_type>, typename ImageT::pixel_type> >
  inline normalize( ImageViewBase<ImageT> const& image,
                    typename ImageChannelType<ImageT>::type low, typename ImageChannelType<ImageT>::type high );

  /// Renormalize the values in an image to fall within the range
  /// [0,high).  The low end of the range is actually determined by
  /// the ChannelRange type trait but is generally zero.
  template <class ImageT>
  UnaryPerPixelView<ImageT, UnaryCompoundFunctor<ChannelNormalizeFunctor<typename ImageT::pixel_type>, typename ImageT::pixel_type> >
  inline normalize( ImageViewBase<ImageT> const& image, typename ImageChannelType<ImageT>::type high );

  /// Renormalize the values in an image to fall within the range
  /// [min,max), where min and max are determined by the ChannelRange
  /// type trait and are generally equal to 0.0 and 1.0 for floating
  /// point types and 0 and the largest positive value for integral types.
  template <class ImageT>
  UnaryPerPixelView<ImageT, UnaryCompoundFunctor<ChannelNormalizeFunctor<typename ImageT::pixel_type>, typename ImageT::pixel_type> >
  inline normalize( ImageViewBase<ImageT> const& image );


  // *******************************************************************
  // threshold()
  // *******************************************************************

  // A per-pixel thresholding filter with adjustable threshold and
  // high and low values.
  template <class PixelT>
  class ChannelThresholdFunctor;

  /// Threshold the values in an image, generating a two-valued output
  /// image with values low and high.
  template <class ImageT, class ThreshT, class LowT, class HighT>
  UnaryPerPixelView<ImageT, UnaryCompoundFunctor<ChannelThresholdFunctor<typename ImageT::pixel_type>, typename ImageT::pixel_type> >
  inline threshold( ImageViewBase<ImageT> const& image, ThreshT thresh, LowT low, HighT high );

  /// Threshold the values in an image, generating a two-valued output
  /// image with values 0 and high.  The low value is actually
  /// determined by the ChannelRange type trait but is generally zero.
  template <class ImageT, class ThreshT, class HighT>
  UnaryPerPixelView<ImageT,UnaryCompoundFunctor<ChannelThresholdFunctor<typename ImageT::pixel_type>, typename ImageT::pixel_type> >
  inline threshold( ImageViewBase<ImageT> const& image, ThreshT thresh, HighT high );

  /// Threshold the values in an image, generating a two-valued output
  /// where the values are determined by the ChannelRange type trait
  /// and are generally equal to 0.0 and 1.0 for floating point types
  /// and 0 and the largest positive value for integral types.
  template <class ImageT, class ThreshT>
  UnaryPerPixelView<ImageT,UnaryCompoundFunctor<ChannelThresholdFunctor<typename ImageT::pixel_type>, typename ImageT::pixel_type> >
  inline threshold( ImageViewBase<ImageT> const& image, ThreshT thresh );

  /// Threshold the values in an image against zero, generating a
  /// two-valued output where the values are determined by the
  /// ChannelRange type trait and are generally equal to 0.0 and 1.0
  /// for floating point types and 0 and the largest positive value for
  /// integral types.
  template <class ImageT>
  UnaryPerPixelView<ImageT,UnaryCompoundFunctor<ChannelThresholdFunctor<typename ImageT::pixel_type>, typename ImageT::pixel_type> >
  inline threshold( ImageViewBase<ImageT> const& image );

  // ******************************************************************
  // MeanFillTransparent
  // ******************************************************************

  // This is a preprocess step that set the value of transparent
  // pixels to the mean of the nearby opaque pixels. This will not
  // produce a visible difference to the image as it only modifies
  // completely transparent pixels. The reason for this is to remove a
  // "bath tub ring" that happens when interpolating/resampling an
  // image with transparent sections.

  template <class ImageT>
  class MeanFillTransparent;

  template <class SourceT>
  MeanFillTransparent<SourceT>
  inline mean_fill_transparent( ImageViewBase<SourceT> const& src ) {
    return MeanFillTransparent<SourceT>( src.impl() );
  }

  // ******************************************************************
  // ComputeNormals
  // ******************************************************************

  // Compute a vector normal to the surface of a DEM for each given
  // pixel.  The normal is computed by forming a plane with three points
  // in the vicinity of the requested pixel, and then finding the vector
  // normal to that plane.  The user must specify the scale in the [u,v]
  // directions so that the direction of the vector in physical space
  // can be properly ascertained.  This is often contained in the (0,0)
  // and (1,1) entry of the georeference transform.
  class ComputeNormalsFunc : public ReturnFixedType<PixelMask<Vector3f> >
  {
    float m_u_scale, m_v_scale;

  public:
    ComputeNormalsFunc(float u_scale, float v_scale) :
      m_u_scale(u_scale), m_v_scale(v_scale) {}

    BBox2i work_area() const { return BBox2i(Vector2i(0, 0), Vector2i(1, 1)); }

    template <class PixelAccessorT>
    PixelMask<Vector3f> operator() (PixelAccessorT const& accessor_loc) const {
      PixelAccessorT acc = accessor_loc;

      // Pick out the three altitude values.
      if (is_transparent(*acc))
        return PixelMask<Vector3f>();
      float alt1 = *acc;

      acc.advance(1,0);
      if (is_transparent(*acc))
        return PixelMask<Vector3f>();
      float alt2 = *acc;

      acc.advance(-1,1);
      if (is_transparent(*acc))
        return PixelMask<Vector3f>();
      float alt3 = *acc;

      // Form two orthogonal vectors in the plane containing the three
      // altitude points
      Vector3f n1(m_u_scale, 0, alt2-alt1);
      Vector3f n2(0, m_v_scale, alt3-alt1);

      // Return the vector normal to the local plane.
      return normalize(cross_prod(n1,n2));
    }
  }; // End class ComputeNormalsFunc


  template <class ViewT>
  UnaryPerPixelAccessorView<EdgeExtensionView<ViewT,ConstantEdgeExtension>, ComputeNormalsFunc> 
  compute_normals(ImageViewBase<ViewT> const& image, float u_scale, float v_scale) {
    return UnaryPerPixelAccessorView<EdgeExtensionView<ViewT,ConstantEdgeExtension>, ComputeNormalsFunc>
                (edge_extend(image.impl(), ConstantEdgeExtension()), ComputeNormalsFunc (u_scale, v_scale));
  }

  // ******************************************************************
  // DotProduct
  // ******************************************************************

  /// Perform the dot product between each pixel and a constant vector.
  class DotProdFunc : public ReturnFixedType<PixelMask<PixelGray<float> > > {
    Vector3f m_vec;
  public:
    DotProdFunc(Vector3f const& vec) : m_vec(normalize(vec)) {}
    PixelMask<PixelGray<float> > operator() (PixelMask<Vector3f> const& pix) const {
      if (is_transparent(pix))
        return PixelMask<PixelGray<float> >();
      else
        return dot_prod(pix.child(),m_vec);
    }
  };

  template <class ViewT>
  UnaryPerPixelView<ViewT, DotProdFunc> dot_prod(ImageViewBase<ViewT> const& view, Vector3f const& vec) {
    return UnaryPerPixelView<ViewT, DotProdFunc>(view.impl(), DotProdFunc(vec));
    }

  // ******************************************************************
  // TwoThresholdFill
  // ******************************************************************

  /// Apply a double threshold to an image.
  /// - Pixels are set if they are above the high threshold.  In addition, a flood-fill is performed
  ///   from the high threshold pixels using pixels above the low threshold.
  /// - No built-in masked pixel handling.
  template <class ImageT>
  class TwoThresholdFill;

  /// Applies a flood fill from pixels which are above the high threshold through pixels above the low threshold.
  template <class ImageT>
  TwoThresholdFill<ImageT>
  two_threshold_fill(ImageViewBase<ImageT> const& image, int expand_size, double low_threshold, double high_threshold,
                     uint8 output_false = 0, uint8 output_true = 1) {
    return TwoThresholdFill<ImageT>(image.impl(), expand_size, low_threshold, high_threshold, output_false, output_true);
  }


} // namespace vw

#include "Algorithms.tcc"

#endif // __VW_IMAGE_ALGORITHMS_H__
