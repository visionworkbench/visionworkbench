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
/// Lazy view image algorithms that implement a sliding window of some sort
/// and are not implemented in their own file.
///
#ifndef __VW_IMAGE_WINDOWALGORITHMS_H__
#define __VW_IMAGE_WINDOWALGORITHMS_H__

#include <vw/Math/Functors.h>
#include <vw/Image/Manipulation.h>
#include <vw/Image/EdgeExtension.h>
#include <vw/Image/ImageViewRef.h>

namespace vw {



/// For each pixel compute the standard deviation of the pixels in a neighborhood.
/// - Could factor out a generic sliding buffer class, which could also implement the median filter.
template <class ImageT, class EdgeT>
class StdDevView : public ImageViewBase<StdDevView<ImageT,EdgeT> >
{
private:
  ImageT   m_image;
  EdgeT    m_edge;     ///< Edge extension type
  Vector2i m_window_size;
  int      m_half_width;
  int      m_half_height;

public:
  typedef typename ImageT::pixel_type pixel_type;  ///< The pixel type of the image view.
  typedef pixel_type                  result_type; ///< We compute the result, so we return by value.
  typedef ProceduralPixelAccessor<StdDevView<ImageT, EdgeT> > 
                                      pixel_accessor; ///< The view's pixel_accessor type.

  /// Constructor
  StdDevView( ImageT const& image, Vector2i window_size, EdgeT  const& edge = EdgeT() )
    : m_image(image), m_edge(edge), m_window_size(window_size) {
    m_half_width  = m_window_size[0]/2;
    m_half_height = m_window_size[1]/2;
  }

  inline int32 cols  () const { return m_image.cols  (); }
  inline int32 rows  () const { return m_image.rows  (); }
  inline int32 planes() const { return m_image.planes(); }

  /// Returns a pixel_accessor pointing to the origin.
  inline pixel_accessor origin() const { return pixel_accessor( *this ); }

  /// Not implemented, the class is designed for fast operation on tiles.
  inline result_type operator()( int32 x, int32 y, int32 p=0 ) const {
    vw_throw(NoImplErr() << "StdDevView::operator()(...) is not implemented");
    return pixel_type();
  }

  // Generate the output tile
  typedef CropView<ImageView<typename ImageT::pixel_type> > prerasterize_type;
  inline prerasterize_type prerasterize( BBox2i const& bbox ) const {

    // Rasterize the input support region
    BBox2i larger_bbox = bbox;
    larger_bbox.expand(m_half_width);
    ImageView<typename ImageT::pixel_type> src = crop(edge_extend(m_image, m_edge), larger_bbox);
    const int col_off = m_half_width;
    const int row_off = m_half_width;

    ImageView<typename ImageT::pixel_type> dst(bbox.width(), bbox.height());
    const size_t window_area = m_window_size[0]*m_window_size[1];
    for (int r=0; r<dst.rows(); ++r) {

      // At the start of each row repopulate the sliding buffer from scratch.
      StdDevSlidingFunctor stddev_functor(window_area);
      for (int x=-m_half_width; x<=m_half_width; ++x) { // Iterate over columns, then rows.
        for (int y=-m_half_height; y<=m_half_height; ++y) {
          if (is_valid(src(0+x+col_off, r+y+row_off)))
            stddev_functor.push(src(0+x+col_off, r+y+row_off));
          else
            stddev_functor.pop();
        }
      }
      dst(0,r) = stddev_functor.get_std_dev();

      // For successive pixels use sliding window optimizations
      for (int c=1; c<dst.cols(); ++c) {
        for (int y=-m_half_height; y<=m_half_height; ++y) {
          if (is_valid(src(c+m_half_width+col_off, r+y+row_off)))
            stddev_functor.push(src(c+m_half_width+col_off, r+y+row_off)); // Add pixels at leading edge of the window.
          else
            stddev_functor.pop();
        }
        dst(c,r) = stddev_functor.get_std_dev();
      }

    } // End loop through rows

    // Use the crop trick to fake that the support region is the same size as the entire image.
    return crop(dst, -bbox.min().x(), -bbox.min().y(), m_image.cols(), m_image.rows());
  }

  template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const {
    vw::rasterize( prerasterize(bbox), dest, bbox );
  }
}; // End class StdDevView

/// Apply a stddev filter to an input image
template <class ImageT, class EdgeT>
StdDevView<ImageT, EdgeT> 
stddev_filter_view(ImageT const& image, Vector2i window_size, EdgeT edge) {
  typedef StdDevView<ImageT, EdgeT> return_type;
  return return_type(image, window_size, edge);
}
/// Overload to set default edge extension.
template <class ImageT>
StdDevView<ImageT, ConstantEdgeExtension> 
stddev_filter_view(ImageT const& image, Vector2i window_size) {
  typedef StdDevView<ImageT, ConstantEdgeExtension> return_type;
  return return_type(image, window_size, ConstantEdgeExtension());
}



//============================================================================



/// Generic class for implementing a function that generates a pixel value based
///  on a window around a pixel.
/// - This allows for simple but unoptimized implementations.
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

  /// Constructor
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
    
    // Require a certain percentage of input pixels for a valid output result.
    const double VALID_PERCENTAGE = 0.9;
    int min_valid_values = static_cast<double>(image.impl().rows()*image.impl().cols()) * VALID_PERCENTAGE;
    if (index < min_valid_values) {
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


#endif // __VW_IMAGE_WINDOWALGORITHMS_H__
