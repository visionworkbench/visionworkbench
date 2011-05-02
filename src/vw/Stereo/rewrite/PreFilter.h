// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#ifndef __VW_STEREO_REWRITE_PREFILTER_H__
#define __VW_STEREO_REWRITE_PREFILTER_H__

#include <vw/Image/Filter.h>

namespace vw {
namespace stereo {
namespace rewrite {

  // Enum of operations
  enum PreFilterType {
    NULLOP,
    LAPLACIAN_OF_GAUSSIAN,
    SUBTRACTED_MEAN
  };

  template <int type>
  struct preprocessing {
    template <class ImageT> ImageT
    filter( ImageViewBase<ImageT> const& image ) const { return image.impl(); }
  };

  template <>
  struct preprocessing<NULLOP> {
    template <class ImageT> ImageT
    filter( ImageViewBase<ImageT> const& image ) const { return image.impl(); }
  };

  template <>
  struct preprocessing<LAPLACIAN_OF_GAUSSIAN> {
    float kernel_width;
    preprocessing( float size ) : kernel_width(size) {}

    template <class ImageT>
    ConvolutionView<SeparableConvolutionView<ImageT, typename DefaultKernelT<typename ImageT::pixel_type>::type, ConstantEdgeExtension>, ImageView<typename DefaultKernelT<typename ImageT::pixel_type>::type>, ConstantEdgeExtension>
    filter( ImageViewBase<ImageT> const& image ) const {
      return laplacian_filter(gaussian_filter(image.impl(),kernel_width));
    }
  };

  template <>
  struct preprocessing<SUBTRACTED_MEAN> {
    float kernel_width;
    preprocessing( float size ) : kernel_width(size) {}

    template <class ImageT>
    BinaryPerPixelView<ImageT,SeparableConvolutionView<ImageT, typename DefaultKernelT<typename ImageT::pixel_type>::type, ConstantEdgeExtension>,vw::ArgArgDifferenceFunctor>
    filter( ImageViewBase<ImageT> const& image ) const {
      return image.impl() - gaussian_filter( image.impl(), kernel_width );
    }
  };

}}} // end namespace vw::stereo::rewrite

#endif//__VW_STEREO_REWRITE_PREFILTER_H__
