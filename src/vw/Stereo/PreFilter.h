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


#ifndef __VW_STEREO_PREFILTER_H__
#define __VW_STEREO_PREFILTER_H__

#include <vw/Image/Filter.h>

namespace vw {
namespace stereo {

  enum PrefilterModeType {
    PREFILTER_NONE    = 0,
    PREFILTER_MEANSUB = 1,
    PREFILTER_LOG     = 2
  };

  // TODO: Should we un-CRTP these classes?

  // This is a base class that is used in other code to make sure the
  // user is passing an actual pre-processing filter as opposed to say an 'int'.
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

/// Create the selected prefilter class and apply the filter to an image.
/// - The output image is rasterized to an ImageView.
template <class ImageT>
ImageView<typename ImageT::pixel_type> 
prefilter_image(ImageViewBase<ImageT> const& image,
                PrefilterModeType prefilter_mode,
                float             prefilter_width) {

  if (prefilter_mode == PREFILTER_LOG){  // LOG
    stereo::LaplacianOfGaussian prefilter(prefilter_width);
    return prefilter.filter(image);
  }
  if (prefilter_mode == PREFILTER_MEANSUB){  // Subtracted mean
    stereo::SubtractedMean prefilter(prefilter_width);
    return prefilter.filter(image);
  }
  //Default: PREFILTER_NONE
  stereo::NullOperation prefilter;
  return prefilter.filter(image);
}


}} // end namespace vw::stereo

#endif//__VW_STEREO_PREFILTER_H__
