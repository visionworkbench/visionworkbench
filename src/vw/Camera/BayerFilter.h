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


/// \file BayerFilter.h
///
/// This function performs a bayer pattern decoding of a grayscale image.
///
#ifndef __VW_CAMERA_BAYER__
#define __VW_CAMERA_BAYER__

#include <vw/Core/CompoundTypes.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/EdgeExtension.h>

namespace vw {
namespace camera {

  template <class ViewT>
  ImageView<PixelRGB<typename CompoundChannelType<typename ViewT::pixel_type>::type > >
  inverse_bayer_filter(ImageViewBase<ViewT > const& view_) {
    typedef typename CompoundChannelType<typename ViewT::pixel_type>::type channel_type;
    const EdgeExtensionView<ViewT, ZeroEdgeExtension> view = edge_extend(view_.impl(), ZeroEdgeExtension());
    channel_type r,g,b;

    ImageView<PixelRGB<channel_type> > output_image(view.cols(), view.rows());

    for( int32 j=0; j<output_image.rows(); ++j ) {
      for( int32 i=0; i<output_image.cols(); ++i ) {
        if((i%2)!=(j%2)) { //green
          g=view(i+1,j+1).v();
          if(i%2) {
            b=(view(i,j+1).v()+view(i+2,j+1).v())/2;
            r=(view(i+1,j).v()+view(i+1,j+2).v())/2;
          } else {
            r=(view(i,j+1).v()+view(i+2,j+1).v())/2;
            b=(view(i+1,j).v()+view(i+1,j+2).v())/2;
          }
        } else if(i%2) { //red
          r=view(i+1,j+1).v();
          b=(view(i,j).v()+view(i+2,j).v()+view(i+2,j+2).v()+view(i,j+2).v())/4;
          g=(view(i+1,j).v()+view(i+1,j+2).v()+view(i+2,j+1).v()+view(i,j+1).v())/4;
        } else { //blue
          b=view(i+1,j+1).v();
          r=(view(i,j).v()+view(i+2,j).v()+view(i+2,j+2).v()+view(i,j+2).v())/4;
          g=(view(i+1,j).v()+view(i+1,j+2).v()+view(i+2,j+1).v()+view(i,j+1).v())/4;
        }
        output_image(i,j)=PixelRGB<channel_type>(r,g,b);
      }
    }
    return output_image;
  }

}} // namespace vw::camera

#endif // __VW_CAMERA_BAYER__
