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
/// implementations.
///
#ifndef __VW_IMAGE_ALGORITHMS_H__
#define __VW_IMAGE_ALGORITHMS_H__

#include <vw/Image/AlgorithmFunctions.h>
#include <vw/Image/PerPixelViews.h>
#include <vw/Image/Statistics.h>
#include <vw/Image/Manipulation.h>
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
  /// and 0 and the largest positve value for integral types.
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
  /// point types and 0 and the largest positve value for integral types.
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
  /// and 0 and the largest positve value for integral types.
  template <class ImageT, class ThreshT>
  UnaryPerPixelView<ImageT,UnaryCompoundFunctor<ChannelThresholdFunctor<typename ImageT::pixel_type>, typename ImageT::pixel_type> >
  inline threshold( ImageViewBase<ImageT> const& image, ThreshT thresh );

  /// Threshold the values in an image against zero, generating a
  /// two-valued output where the values are determined by the
  /// ChannelRange type trait and are generally equal to 0.0 and 1.0
  /// for floating point types and 0 and the largest positve value for
  /// integral types.
  template <class ImageT>
  UnaryPerPixelView<ImageT,UnaryCompoundFunctor<ChannelThresholdFunctor<typename ImageT::pixel_type>, typename ImageT::pixel_type> >
  inline threshold( ImageViewBase<ImageT> const& image );

  // *******************************************************************
  // clear_nonopaque_pixels()
  //
  // This filter is useful for eliminating fringe effects along the
  // edges of images with some transparent or nodata values that have
  // be transformed with bilinear or bicubic interpolation.
  // *******************************************************************
  template <class PixelT>
  class ClearNonOpaqueFunctor;

  /// Zero out any pixels that aren't completely opaque.
  template <class ImageT>
  UnaryPerPixelView<ImageT,ClearNonOpaqueFunctor<typename ImageT::pixel_type> >
  inline clear_nonopaque_pixels( ImageViewBase<ImageT> const& image );

  // *******************************************************************
  // remap_pixel_value()
  //
  // This filter can be used to map one pixel value to another.  This
  // can be useful in many situations, for example when you need to
  // remap the nodata value used in a DEM.
  // *******************************************************************
  template <class PixelT>
  class RemapPixelFunctor;

  /// Zero out any pixels that aren't completely opaque.
  template <class ImageT>
  UnaryPerPixelView<ImageT,RemapPixelFunctor<typename ImageT::pixel_type> >
  inline remap_pixel_value( ImageViewBase<ImageT> const& image,
                            typename PixelChannelType<typename ImageT::pixel_type>::type src_val,
                            typename PixelChannelType<typename ImageT::pixel_type>::type dst_val);

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
  // FillHoles
  // ******************************************************************
  
  // A class to fill holes in images. Given the number
  // hole_fill_len, which gives the size of holes, look left, right,
  // up, and down, as far as hole_fill_len/2. Must have valid pixels
  // both left and right, otherwise both up and down, else do
  // nothing. The motivation here is to fill in only pixels
  // "surrounded" by valid pixels.
  
  // Use one of the two approaches:
  //
  // 1. Interpolate using the left and right values. Interpolate using
  // the up and down values. Average the results. Fast.
  //
  // 2. Find the weighted average of the points in the window
  // of size hole_fill_len centered at the current point. Slow.

  // After holes are filled, do several passes to average the results.
  
  // The image being passed in should be a PixelMask.
  template <class ImageT>
  class FillHoles;

  // Fill holes. This algorithm is very naive and works not so well.
  // The grassfire-based in-painting algorithm from ASP works much better.
  template <class ImageT> FillHoles<ImageT>
  fill_holes( ImageViewBase<ImageT> const& img, int hole_fill_mode, int hole_fill_num_smooth_iter, int hole_fill_len) {
    return FillHoles<ImageT>( img.impl(), hole_fill_mode, hole_fill_num_smooth_iter, hole_fill_len );
  }

  // ----------------------------------------------------------------------------

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
  class ComputeNormalsFunc;

  template <class ViewT>
  UnaryPerPixelAccessorView<EdgeExtensionView<ViewT,ConstantEdgeExtension>, ComputeNormalsFunc> 
  compute_normals(ImageViewBase<ViewT> const& image, float u_scale, float v_scale) {
    return UnaryPerPixelAccessorView<EdgeExtensionView<ViewT,ConstantEdgeExtension>, ComputeNormalsFunc>
                (edge_extend(image.impl(), ConstantEdgeExtension()), ComputeNormalsFunc (u_scale, v_scale));
  }

  // ******************************************************************
  // DotProduct
  // ******************************************************************

  class DotProdFunc;

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


//===================================================================
// TODO: Move WindowFunction code to a new file!


/// Generic class for implementing a function that generates a pixel value based
///  on a window around a pixel.
/// - TODO: Could ConvolutionView use this class?
template <class ImageT, class FuncT, class EdgeT>
class WindowFunctionView : public ImageViewBase<WindowFunctionView<ImageT,FuncT,EdgeT> >
{
private:
  ImageT   m_image;
  EdgeT    m_edge;     ///< Edge extension type
  FuncT    m_functor;  ///< Functor that operates on each window.
  Vector2i m_window_size;
  int      m_half_width;
  int      m_half_height;

public:
  typedef typename ImageT::pixel_type pixel_type;  ///< The pixel type of the image view.
  typedef pixel_type                  result_type; ///< We compute the result, so we return by value.
  typedef ProceduralPixelAccessor<WindowFunctionView<ImageT, FuncT, EdgeT> > 
                                      pixel_accessor; ///< The view's pixel_accessor type.

  /// Constructs a ConvolutionView with the given image and kernel and 
  /// with the origin of the kernel located at the point (ci,cj).
  WindowFunctionView( ImageT const& image, Vector2i window_size,
                      FuncT  const& functor, 
                      EdgeT  const& edge = EdgeT() )
    : m_image(image), m_edge(edge), m_functor(functor), m_window_size(window_size) {
    m_half_width  = m_window_size[0]/2;
    m_half_height = m_window_size[1]/2;
  }

  inline int32 cols  () const { return m_image.cols  (); }
  inline int32 rows  () const { return m_image.rows  (); }
  inline int32 planes() const { return m_image.planes(); }

  /// Returns a pixel_accessor pointing to the origin.
  inline pixel_accessor origin() const { return pixel_accessor( *this ); }

  /// Returns the pixel at the given position in the given plane.
  inline result_type operator()( int32 x, int32 y, int32 p=0 ) const {
    BBox2i roi(x-m_half_width, y-m_half_height, m_window_size[0], m_window_size[1]);
    return m_functor(edge_extend(m_image, roi, m_edge));
  }

  // Edge extension is done in the prerasterize function so the returned type does not need edge extension
  // - Currently the class does the hard work in the () function but it would probably be more efficient to
  //   do the computation on a per-tile basis.
  typedef WindowFunctionView<CropView<ImageView<typename ImageT::pixel_type> >, 
                             FuncT, NoEdgeExtension> prerasterize_type;
  inline prerasterize_type prerasterize( BBox2i const& bbox ) const {
    // Compute the required base of support for the input bounding box
    BBox2i src_bbox( bbox.min().x() - m_half_width, 
                     bbox.min().y() - m_half_height,
                     bbox.width () + m_window_size[0]-1, 
                     bbox.height() + m_window_size[1]-1 );
    // Take an edge extended image view of the input support region
    ImageView<typename ImageT::pixel_type> src = edge_extend(m_image, src_bbox, m_edge);
    // Use the crop trick to fake that the support region is the same size as the entire image.
    return prerasterize_type( crop(src, -src_bbox.min().x(), -src_bbox.min().y(), m_image.cols(), m_image.rows()),
                              m_window_size, m_functor, NoEdgeExtension() );
  }

  template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const {
    vw::rasterize( prerasterize(bbox), dest, bbox );
  }
}; // End class WindowFunctionView

/// Functor to find the median of an input image.
/// - Usually used with WindowFunctionView.
template <typename ImageT>
struct WindowMedianFunctor {

  mutable std::vector<double> m_values; ///< Persistent storage location

  /// Constructor
  WindowMedianFunctor(Vector2i window_size) {
    m_values.resize(window_size[0]*window_size[1]);
  }

  /// Returns the median pixel of the provided input image.
  /// - Generally the input image will be a cropped view of a whole image.
  template <class T>
  typename ImageT::pixel_type operator()(ImageViewBase<T> const& image) const {

    // Loop through the kernel and collect the values
    int index = 0;
    for (int r=0; r<image.impl().rows(); ++r) {
      for (int c=0; c<image.impl().cols(); ++c) {
        if (is_valid(image.impl()(c,r))) {
          m_values[index] = image.impl()(c,r);
          ++index;
        }
      }
    }
    if (index == 0) { // All pixels invalid!
      typename ImageT::pixel_type result(0);
      invalidate(result);
      return result;
    }

    // Now that we have all the values, compute the median.
    double median;
    if (index == static_cast<int>(m_values.size())) // No invalid pixels
      median = math::destructive_median(m_values);
    else { // Invalid pixels
      // Resize the vector twice so we can call the median function
      size_t full_size = m_values.size();
      m_values.resize(index);
      median = math::destructive_median(m_values);
      m_values.resize(full_size);
    }
    return typename ImageT::pixel_type(median);
  }
}; // End class WindowMedianFunctor

/// Apply a median filter to an input image
template <class ImageT, class EdgeT>
WindowFunctionView<ImageT, WindowMedianFunctor<ImageT>, EdgeT> 
median_filter_view(ImageT const& image, Vector2i window_size, EdgeT edge) {
  typedef WindowFunctionView<ImageT, WindowMedianFunctor<ImageT>, EdgeT> return_type;
  WindowMedianFunctor<ImageT> functor(window_size);
  return return_type(image, window_size, functor, edge);
}
/// Overload to set default edge extension.
template <class ImageT>
WindowFunctionView<ImageT, WindowMedianFunctor<ImageT>, ConstantEdgeExtension> 
median_filter_view(ImageT const& image, Vector2i window_size) {
  typedef WindowFunctionView<ImageT, WindowMedianFunctor<ImageT>, ConstantEdgeExtension> return_type;
  WindowMedianFunctor<ImageT> functor(window_size);
  return return_type(image, window_size, functor, ConstantEdgeExtension());
}

} // namespace vw

#include "Algorithms.tcc"

#endif // __VW_IMAGE_ALGORITHMS_H__
