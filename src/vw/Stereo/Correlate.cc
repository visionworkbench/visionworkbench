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


#include <vw/Image/ImageMath.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/ImageViewBase.h>
#include <vw/Image/PixelMask.h>
#include <vw/Image/PixelAccessors.h>

#include <vw/Stereo/Correlate.h>

using namespace vw;
using namespace vw::stereo;


int
stereo::adjust_weight_image(ImageView<float> &weight,
                            ImageView<PixelMask<Vector2f> > const& disparity_map_patch,
                            ImageView<float> const& weight_template) {
  float sum = 0;
  int32 num_good_pix = 0;
  typedef ImageView<float>::pixel_accessor IViewFAcc;
  typedef ImageView<PixelMask<Vector2f> >::pixel_accessor IViewDAcc;
  IViewFAcc weight_row_acc = weight.origin();
  IViewFAcc template_row_acc = weight_template.origin();
  IViewDAcc disp_row_acc = disparity_map_patch.origin();
  for (int32 j = 0; j < weight_template.rows(); ++j) {
    IViewFAcc weight_col_acc = weight_row_acc;
    IViewFAcc template_col_acc = template_row_acc;
    IViewDAcc disp_col_acc = disp_row_acc;
    for (int32 i = 0; i < weight_template.cols(); ++i ) {

      // Mask is zero if the disparity map's pixel is missing...
      if ( !is_valid(*disp_col_acc) )
        *weight_col_acc = 0;

      // ... otherwise we use the weight from the weight template
      else {
        *weight_col_acc = *template_col_acc;
        sum += *weight_col_acc;
        ++num_good_pix;
      }

      disp_col_acc.next_col();
      weight_col_acc.next_col();
      template_col_acc.next_col();
    }
    disp_row_acc.next_row();
    weight_row_acc.next_row();
    template_row_acc.next_row();
  }

  // Normalize the weight image
  if (sum == 0)
    vw_throw(LogicErr() << "subpixel_weight: Sum of weight image was zero.  This isn't supposed to happen!");
  else
    weight /= sum;
  return num_good_pix;
}

ImageView<float>
stereo::detail::compute_spatial_weight_image(int32 kern_width, int32 kern_height,
                                             float two_sigma_sqr) {
  int32 center_pix_x = kern_width/2;
  int32 center_pix_y = kern_height/2;
  float sum;

  sum = 0.0;
  ImageView<float> weight(kern_width, kern_height);
  for (int32 j = 0; j < kern_height; ++j) {
    for (int32 i = 0; i < kern_width; ++i ) {
      weight(i,j) = exp(-1*((i-center_pix_x)*(i-center_pix_x) +
                            (j-center_pix_y)*(j-center_pix_y)) / two_sigma_sqr);
      sum += weight(i,j);
    }
  }

  weight /= sum;

  return weight;
}
