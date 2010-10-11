// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file Filter.h
///
/// Image filtering functions and classes.
///
/// These are the commonly-used image filters.  For the most part
/// the actual work is done lazily by filtering views, of which
/// the most important are the convolution views declared in
/// Convolution.h.
#ifndef __VW_IMAGE_FILTER_H__
#define __VW_IMAGE_FILTER_H__

#include <vector>

#include <boost/type_traits.hpp>
#include <boost/mpl/logical.hpp>

#include <vw/Image/ImageView.h>
#include <vw/Image/Convolution.h>

namespace vw {

  /// \cond INTERNAL

  /// This template function computes the default numerical type for
  /// the kernels used by the kernel-based filters, such as
  /// GaussianFilter and DerivativeFilter.  The type defaults to
  /// single-precision float for all pixel types except those based on
  /// double, for which it defaults to double.
  template <class PixelT>
  struct DefaultKernelT {
    typedef typename boost::mpl::if_< boost::is_same<typename PixelChannelType<PixelT>::type, double>, double, float >::type type;
  };

  /// Computes a Gaussian kernel.
  /// Instantiated by default only for float and double kernels.
  template <class KernelT>
  void generate_gaussian_kernel( std::vector<KernelT>& kernel, double sigma, int32 size=0 );

  /// Computes a differentiation kernel.
  /// Instantiated by default only for float and double kernels.
  template <class KernelT>
  void generate_derivative_kernel( std::vector<KernelT>& kernel, int32 deriv, int32 size=0 );

  /// Computes an oriented Gaussian derivative kernel.
  /// Instantiated by default only for float and double kernels.
  template <class KernelT>
  void generate_gaussian_derivative_kernel( ImageView<KernelT>& kernel, double x_sigma, int32 x_deriv, double y_sigma, int32 y_deriv, double angle, int32 size );

  /// Computes an oriented Gaussian derivative kernel.
  inline ImageView<double> generate_gaussian_derivative_kernel( double x_sigma, int32 x_deriv, double y_sigma, int32 y_deriv, double angle, int32 size ) {
    ImageView<double> result;
    generate_gaussian_derivative_kernel( result, x_sigma, x_deriv, y_sigma, y_deriv, angle, size );
    return result;
  }

  /// Computes a Laplacian of Gaussian kernel.
  /// Instantiated by default only for float and double kernels.
  template <class KernelT>
  void generate_laplacian_of_gaussian_kernel( ImageView<KernelT>& kernel, double sigma, int32 size );

  /// Computes a Laplacian of Gaussian kernel.
  inline ImageView<double> generate_laplacian_of_gaussian_kernel( double sigma, int32 size ) {
    ImageView<double> result;
    generate_laplacian_of_gaussian_kernel( result, sigma, size );
    return result;
  }

  /// \endcond

  // General 2D convolution filter functions

  /// This function computes the convolution of an image with a kernel
  /// stored in another image. It assumes the origin of the kernel is
  /// at the point <CODE>(cx,cy)</CODE> and uses the given edge
  /// extension mode to extend the source image as needed.
  /// \see FFTConvolutionFilter
  template <class SrcT, class KernelT, class EdgeT>
  inline ConvolutionView<SrcT,KernelT,EdgeT>
  convolution_filter( ImageViewBase<SrcT> const& src, KernelT const& kernel, int32 cx, int32 cy, EdgeT edge ) {
    return ConvolutionView<SrcT,KernelT,EdgeT>( src.impl(), kernel, cx, cy, edge );
  }

  /// This is an overloaded function provided for convenience; see
  /// vw::convolution_filter. It assumes that the origin of the kernel
  /// is at the point <B>((kernel.cols()-1)/2,(kernel.rows()-1)/2)</B>.
  template <class SrcT, class KernelT, class EdgeT>
  inline ConvolutionView<SrcT,KernelT,EdgeT>
  convolution_filter( ImageViewBase<SrcT> const& src, KernelT const& kernel, EdgeT edge ) {
    return ConvolutionView<SrcT,KernelT,EdgeT>( src.impl(), kernel, edge );
  }

  /// This is an overloaded function provided for convenience; see
  /// vw::convolution_filter. It uses the default
  /// vw::ConstantEdgeExtension mode.
  template <class SrcT, class KernelT>
  inline ConvolutionView<SrcT,KernelT,ConstantEdgeExtension>
  convolution_filter( ImageViewBase<SrcT> const& src, KernelT const& kernel, int32 cx, int32 cy ) {
    return ConvolutionView<SrcT,KernelT,ConstantEdgeExtension>( src.impl(), kernel, cx, cy );
  }

  /// This is an overloaded function provided for convenience; see
  /// vw::convolution_filter. It assumes that the origin of the kernel
  /// is at the point <B>((kernel.cols()-1)/2,(kernel.rows()-1)/2)</B>
  /// and uses the default vw::ConstantEdgeExtension mode.
  template <class SrcT, class KernelT>
  inline ConvolutionView<SrcT,KernelT,ConstantEdgeExtension>
  convolution_filter( ImageViewBase<SrcT> const& src, KernelT const& kernel ) {
    return ConvolutionView<SrcT,KernelT,ConstantEdgeExtension>( src.impl(), kernel );
  }


  // Separable convolution filter functions

  /// This function computes the convolution of an image with a
  /// separable kernel. The two components of the kernel may be stored
  /// in any kind of range object, that is, any object adhering to the
  /// usual 1D container semantics with begin() and end() methods,
  /// including an std::vector or a vw::ImageView. It assumes the
  /// origin of the kernel is at the point <CODE>(cx,cy)</CODE> and
  /// uses the given edge extension mode to extend the source image as
  /// needed.
  template <class SrcT, class KRangeT, class EdgeT>
  inline SeparableConvolutionView<SrcT,typename KRangeT::value_type,EdgeT>
  separable_convolution_filter( ImageViewBase<SrcT> const& src, KRangeT const& x_kernel, KRangeT const& y_kernel, int32 cx, int32 cy, EdgeT edge ) {
    return SeparableConvolutionView<SrcT,typename KRangeT::value_type,EdgeT>( src.impl(), x_kernel, y_kernel, cx, cy, edge );
  }

  /// This is an overloaded function provided for convenience; see
  /// vw::separable_convolution_filter. It assumes that the origin of the
  /// kernel is at the point
  /// <B>(x_kernel.size()/2,y_kernel.size()/2)</B>. (Note that this
  /// requires that the kernel objects support the size() method.)
  template <class SrcT, class KRangeT, class EdgeT>
  inline SeparableConvolutionView<SrcT,typename KRangeT::value_type,EdgeT>
  separable_convolution_filter( ImageViewBase<SrcT> const& src, KRangeT const& x_kernel, KRangeT const& y_kernel, EdgeT edge ) {
    return SeparableConvolutionView<SrcT,typename KRangeT::value_type,EdgeT>( src.impl(), x_kernel, y_kernel, edge );
  }

  /// This is an overloaded function provided for convenience; see vw::separable_convolution_filter.
  /// It uses the default vw::ConstantEdgeExtension mode.
  template <class SrcT, class KRangeT>
  inline SeparableConvolutionView<SrcT,typename KRangeT::value_type,ConstantEdgeExtension>
  separable_convolution_filter( ImageViewBase<SrcT> const& src, KRangeT const& x_kernel, KRangeT const& y_kernel, int32 cx, int32 cy ) {
    return SeparableConvolutionView<SrcT,typename KRangeT::value_type,ConstantEdgeExtension>( src.impl(), x_kernel, y_kernel, cx, cy );
  }

  /// This is an overloaded function provided for convenience; see
  /// vw::separable_convolution_filter. It assumes that the origin of the
  /// kernel is at the point
  /// <B>((x_kernel.size()-1)/2,(y_kernel.size()-1)/2)</B> and uses
  /// the default vw::ConstantEdgeExtension mode. (Note that this
  /// requires that the kernel objects support the size() method.)
  template <class SrcT, class KRangeT>
  inline SeparableConvolutionView<SrcT,typename KRangeT::value_type,ConstantEdgeExtension>
  separable_convolution_filter( ImageViewBase<SrcT> const& src, KRangeT const& x_kernel, KRangeT const& y_kernel ) {
    return SeparableConvolutionView<SrcT,typename KRangeT::value_type,ConstantEdgeExtension>( src.impl(), x_kernel, y_kernel );
  }


  // Gaussian convolution functions

  /// This function applies a Gaussian smoothing filter to an image.
  /// It uses an axis-aligned Guassian kernel with standard deviations
  /// of x_sigma and y_sigma, kernel dimensions of x_dim and y_dim,
  /// and the given edge extension mode to extend the source image as
  /// needed.  Specifying a zero falue of x_dim or y_dim causes the
  /// corresponding dimension to be chose automatically as appropriate
  /// for the requested standard deviation.
  template <class SrcT, class EdgeT>
  SeparableConvolutionView<SrcT, typename DefaultKernelT<typename SrcT::pixel_type>::type, EdgeT>
  inline gaussian_filter( ImageViewBase<SrcT> const& src, double x_sigma, double y_sigma, int32 x_dim, int32 y_dim, EdgeT edge ) {
    std::vector<typename DefaultKernelT<typename SrcT::pixel_type>::type> x_kernel, y_kernel;
    generate_gaussian_kernel( x_kernel, x_sigma, x_dim );
    generate_gaussian_kernel( y_kernel, y_sigma, y_dim );
    return SeparableConvolutionView<SrcT, typename DefaultKernelT<typename SrcT::pixel_type>::type, EdgeT>( src.impl(), x_kernel, y_kernel, edge );
  }

  /// This is an overloaded function provided for convenience; see
  /// vw::gaussian_filter. It uses the default
  /// vw::ConstantEdgeExtension mode.
  template <class SrcT>
  SeparableConvolutionView<SrcT, typename DefaultKernelT<typename SrcT::pixel_type>::type, ConstantEdgeExtension>
  inline gaussian_filter( ImageViewBase<SrcT> const& src, double x_sigma, double y_sigma, int32 x_dim, int32 y_dim ) {
    return gaussian_filter( src, x_sigma, y_sigma, x_dim, y_dim, ConstantEdgeExtension() );
  }

  /// This is an overloaded function provided for convenience; see
  /// vw::gaussian_filter. It chooses default kernel dimensions
  /// appropriate for the requested standard deviations.
  template <class SrcT, class EdgeT>
  SeparableConvolutionView<SrcT, typename DefaultKernelT<typename SrcT::pixel_type>::type, EdgeT>
  inline gaussian_filter( ImageViewBase<SrcT> const& src, double x_sigma, double y_sigma, EdgeT edge ) {
    return gaussian_filter( src, x_sigma, y_sigma, 0, 0, edge );
  }

  /// This is an overloaded function provided for convenience; see
  /// vw::gaussian_filter. It chooses default kernel dimensions
  /// appropriate for the requested standard deviations and uses the
  /// default vw::ConstantEdgeExtension mode.
  template <class SrcT>
  SeparableConvolutionView<SrcT, typename DefaultKernelT<typename SrcT::pixel_type>::type, ConstantEdgeExtension>
  inline gaussian_filter( ImageViewBase<SrcT> const& src, double x_sigma, double y_sigma ) {
    return gaussian_filter( src, x_sigma, y_sigma, 0, 0, ConstantEdgeExtension() );
  }

  /// This is an overloaded function provided for convenience; see
  /// vw::gaussian_filter. It uses the same standard deviation in both
  /// directions and chooses default kernel dimensions appropriate for
  /// the requested standard deviations.
  template <class SrcT, class EdgeT>
  SeparableConvolutionView<SrcT, typename DefaultKernelT<typename SrcT::pixel_type>::type, EdgeT>
  inline gaussian_filter( ImageViewBase<SrcT> const& src, double sigma, EdgeT edge ) {
    return gaussian_filter( src, sigma, sigma, 0, 0, edge );
  }

  /// This is an overloaded function provided for convenience; see
  /// vw::gaussian_filter. It uses the same standard deviation in both
  /// directions, chooses default kernel dimensions appropriate for
  /// the requested standard deviations, and uses the default
  /// vw::ConstantEdgeExtension mode.
  template <class SrcT>
  SeparableConvolutionView<SrcT, typename DefaultKernelT<typename SrcT::pixel_type>::type, ConstantEdgeExtension>
  inline gaussian_filter( ImageViewBase<SrcT> const& src, double sigma ) {
    return gaussian_filter( src, sigma, sigma, 0, 0, ConstantEdgeExtension() );
  }


  // Image differentiation functions

  /// Applies a differentiation filter to an image.  This function
  /// performs differentiation of order x_deriv in the x direction and
  /// order y_deriv in the y direction.  The dimensions of the
  /// differentiation kernel are x_dim and y_dim, and must be odd
  /// numbers large enough to accomodate a kernel of the requested
  /// order.  You may also specify zero for either kernel dimension,
  /// in which case it is chosen automatically to be the smallest
  /// suitable dimension.  The source image is edge-exetnded using the
  /// given edge extension mode as needed.
  template <class SrcT, class EdgeT>
  SeparableConvolutionView<SrcT, typename DefaultKernelT<typename SrcT::pixel_type>::type, EdgeT>
  inline derivative_filter( ImageViewBase<SrcT> const& src, int32 x_deriv, int32 y_deriv, int32 x_dim, int32 y_dim, EdgeT edge ) {
    std::vector<typename DefaultKernelT<typename SrcT::pixel_type>::type> x_kernel, y_kernel;
    generate_derivative_kernel( x_kernel, x_deriv, x_dim );
    generate_derivative_kernel( y_kernel, y_deriv, y_dim );
    return SeparableConvolutionView<SrcT, typename DefaultKernelT<typename SrcT::pixel_type>::type, EdgeT>( src.impl(), x_kernel, y_kernel, edge );
  }

  /// This is an overloaded function provided for convenience; see
  /// vw::derivative_filter. It uses a kernel with the default
  /// dimensions for the requested differentiation operation.
  template <class SrcT, class EdgeT>
  SeparableConvolutionView<SrcT, typename DefaultKernelT<typename SrcT::pixel_type>::type, EdgeT>
  inline derivative_filter( ImageViewBase<SrcT> const& src, int32 x_deriv, int32 y_deriv, EdgeT edge ) {
    return derivative_filter( src, x_deriv, y_deriv, 0, 0, edge );
  }

  /// This is an overloaded function provided for convenience; see
  /// vw::derivative_filter. It uses the default
  /// vw::ConstantEdgeExtension mode.
  template <class SrcT>
  SeparableConvolutionView<SrcT, typename DefaultKernelT<typename SrcT::pixel_type>::type, ConstantEdgeExtension>
  inline derivative_filter( ImageViewBase<SrcT> const& src, int32 x_deriv, int32 y_deriv, int32 x_dim, int32 y_dim ) {
    return derivative_filter( src, x_deriv, y_deriv, x_dim, y_dim, ConstantEdgeExtension() );
  }

  /// This is an overloaded function provided for convenience; see
  /// vw::derivative_filter. It uses a kernel with the default
  /// dimensions for the requested differentiation operation and the
  /// default vw::ConstantEdgeExtension mode.
  template <class SrcT>
  SeparableConvolutionView<SrcT, typename DefaultKernelT<typename SrcT::pixel_type>::type, ConstantEdgeExtension>
  inline derivative_filter( ImageViewBase<SrcT> const& src, int32 x_deriv, int32 y_deriv ) {
    return derivative_filter( src, x_deriv, y_deriv, 0, 0, ConstantEdgeExtension() );
  }


  // Image Laplacian functions

  /// Applies a Laplacian filter to an image.  This function
  /// computes the Laplacian \f$ \nabla^2\equiv\frac{d^2}{dx^2}+\frac{d^2}{dy^2} \f$
  /// of an image using a \f$ 3\times3 \f$ discrete differentiation
  /// kernel.  The source image is edge-exetnded using the
  /// given edge extension mode as needed.
  template <class SrcT, class EdgeT>
  ConvolutionView<SrcT, ImageView<typename DefaultKernelT<typename SrcT::pixel_type>::type>, EdgeT>
  inline laplacian_filter( ImageViewBase<SrcT> const& src, EdgeT edge ) {
    ImageView<typename DefaultKernelT<typename SrcT::pixel_type>::type> kernel(3,3);
    kernel(0,0)=0; kernel(1,0)=1;  kernel(2,0)=0;
    kernel(0,1)=1; kernel(1,1)=-4; kernel(2,1)=1;
    kernel(0,2)=0; kernel(1,2)=1;  kernel(2,2)=0;
    return ConvolutionView<SrcT, ImageView<typename DefaultKernelT<typename SrcT::pixel_type>::type>, EdgeT>( src.impl(), kernel, 1, 1, edge );
  }

  /// This is an overloaded function provided for convenience; see
  /// vw::laplacian_filter. It uses the default
  /// vw::ConstantEdgeExtension mode.
  template <class SrcT>
  ConvolutionView<SrcT, ImageView<typename DefaultKernelT<typename SrcT::pixel_type>::type>, ConstantEdgeExtension>
  inline laplacian_filter( ImageViewBase<SrcT> const& src ) {
    return laplacian_filter( src, ConstantEdgeExtension() );
  }


  // Per-pixel filter functions

  /// Filters an image by applying a user-supplied function to each
  /// pixel.  The function argument can be either a regular function
  /// pointer or a functor (i.e. a function object).  Either way it
  /// should take a single pixel value as its only argument and return
  /// a new pixel value, optionally of a different type.  The pixel
  /// type of the resulting image will correspond to the return type
  /// of the function.  If a functor is used, it must obey the return
  /// type deduction rules required by boost::result_of.
  /// \see vw::per_pixel_channel_filter
  /// \see vw/Functors.h
  template <class SrcT, class FuncT>
  UnaryPerPixelView<SrcT,FuncT>
  inline per_pixel_filter( ImageViewBase<SrcT> const& src, FuncT const& func ) {
    return UnaryPerPixelView<SrcT,FuncT>( src.impl(), func );
  }

  template <class Src1T, class Src2T, class FuncT>
  BinaryPerPixelView<Src1T,Src2T,FuncT>
  inline per_pixel_filter( ImageViewBase<Src1T> const& src1, ImageViewBase<Src2T> const& src2, FuncT const& func ) {
    return BinaryPerPixelView<Src1T,Src2T,FuncT>( src1.impl(), src2.impl(), func );
  }

  // Per-pixel-channel filter functions

  /// Filters an image by applying a user-supplied function to each
  /// channel of each pixel.  The function argument can be either a
  /// regular function pointer or a functor (i.e. a function object).
  /// Either way it should take a single pixel channel value as its
  /// only argument and return a new channel value, optionally of a
  /// different type.  The pixel type of the resulting image will
  /// correspond to the return type of the function.  If a functor is
  /// used, it must obey the return type deduction rules required by
  /// boost::result_of.
  /// \see vw::per_pixel_filter
  /// \see vw/Functors.h
  template <class SrcT, class FuncT>
  UnaryPerPixelView<SrcT,UnaryCompoundFunctor<FuncT, typename SrcT::pixel_type> >
  inline per_pixel_channel_filter( ImageViewBase<SrcT> const& src, FuncT const& func ) {
    return UnaryPerPixelView<SrcT,UnaryCompoundFunctor<FuncT, typename SrcT::pixel_type> >( src.impl(), UnaryCompoundFunctor<FuncT,typename SrcT::pixel_type>(func) );
  }

} // namespace vw

#endif // __VW_IMAGE_FILTER_H__
