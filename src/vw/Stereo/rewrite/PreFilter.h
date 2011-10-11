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

  // This is a base class that is used in other code to make sure the
  // user is passing an actual pre-processing filter as opposed to say
  // an 'int'.
  template <class ImplT>
  struct PreFilterBase {
    inline ImplT& impl() { return static_cast<ImplT&>(*this); }
    inline ImplT const& impl() const { return static_cast<ImplT const&>(*this); }
  };

  struct NullOperation : public PreFilterBase<NullOperation> {
    template <class ImageT>
    EdgeExtensionView<ImageT,ConstantEdgeExtension>
    filter( ImageViewBase<ImageT> const& image ) const {
      return edge_extend(image.impl(),ConstantEdgeExtension());
    }
  };

  struct LaplacianOfGaussian : public PreFilterBase<LaplacianOfGaussian> {
    float kernel_width;
    LaplacianOfGaussian( float size ) : kernel_width(size) {}

    template <class ImageT>
    ConvolutionView<SeparableConvolutionView<ImageT, typename DefaultKernelT<typename ImageT::pixel_type>::type, ConstantEdgeExtension>, ImageView<typename DefaultKernelT<typename ImageT::pixel_type>::type>, ConstantEdgeExtension>
    filter( ImageViewBase<ImageT> const& image ) const {
      return laplacian_filter(gaussian_filter(image.impl(),kernel_width));
    }
  };

  struct SubtractedMean : public PreFilterBase<SubtractedMean> {
    float kernel_width;
    SubtractedMean( float size ) : kernel_width(size) {}

    template <class ImageT>
    BinaryPerPixelView<EdgeExtensionView<ImageT, ConstantEdgeExtension>,SeparableConvolutionView<ImageT, typename DefaultKernelT<typename ImageT::pixel_type>::type, ConstantEdgeExtension>,vw::ArgArgDifferenceFunctor>
    filter( ImageViewBase<ImageT> const& image ) const {
      return edge_extend(image.impl(),ConstantEdgeExtension()) - gaussian_filter( image.impl(), kernel_width );
    }
  };

}}} // end namespace vw::stereo::rewrite

#endif//__VW_STEREO_REWRITE_PREFILTER_H__
