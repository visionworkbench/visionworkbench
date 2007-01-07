// __BEGIN_LICENSE__
// 
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
// 
// Copyright 2006 Carnegie Mellon University. All rights reserved.
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

/// \file BayerFilter.h
/// 
/// This function performs a bayer pattern decoding of a grayscale image..
/// 
#ifndef __VW_CAMERA_BAYER__
#define __VW_CAMERA_BAYER__

#include <vw/Image/ImageViewBase.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/PixelTypes.h>
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

    for( unsigned j=0; j<output_image.rows(); ++j ) {
      for( unsigned i=0; i<output_image.cols(); ++i ) {
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
