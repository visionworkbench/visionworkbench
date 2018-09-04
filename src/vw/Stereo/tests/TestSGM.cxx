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


#include <test/Helpers.h>
#include <vw/Image/EdgeExtension.h>
#include <vw/FileIO/DiskImageView.h>
#include <vw/Stereo/SGM.h>

using namespace vw;
using namespace vw::stereo;

TEST( SGM, constant_offset ) {
 
  // For this perfect test case, the correct disparity is (2,1) for each pixel!
  int min_disp_x  = -4;
  int max_disp_x  =  4;
  int min_disp_y  =  -4;
  int max_disp_y  =   4;
  int kernel_size = 3;
  DiskImageView<PixelGray<uint8> > inputLeft ("left.tif");
  DiskImageView<PixelGray<uint8> > inputRight("left_const_offset.tif");
  BBox2i leftRoi (0,0,400, 400);
  
  int disp_x_range = max_disp_x - min_disp_x + 1;
  int disp_y_range = max_disp_y - min_disp_y + 1;
  BBox2i rightRoi = leftRoi + Vector2i(min_disp_x, min_disp_y);
  rightRoi.max() += Vector2i(disp_x_range, disp_y_range);

  ImageView<uint8> left  = crop(inputLeft, leftRoi);
  ImageView<uint8> right = crop(inputRight, rightRoi);
  
  bool use_mgm = false;
  size_t memory_limit_mb = 1024;
  boost::shared_ptr<SemiGlobalMatcher> matcher_ptr;
  CostFunctionType cost_type = CENSUS_TRANSFORM;
  Vector2i search_buffer(4,4);

  SemiGlobalMatcher::DisparityImage result = calc_disparity_sgm(cost_type, left, right, 
                                                BBox2i(0,0,left.cols(), left.rows()),
                                                Vector2i(disp_x_range, disp_y_range),
                                                Vector2i(kernel_size, kernel_size),
                                                use_mgm, SemiGlobalMatcher::SUBPIXEL_LC_BLEND,
                                                search_buffer, memory_limit_mb, matcher_ptr);
                                                
  result = result + PixelMask<Vector2i>(min_disp_x, min_disp_y);
  //write_image("SGM_output.tif", result);

  double num_pixels  = result.rows()*result.cols();
  size_t num_correct = 0;
  for (int row=0; row<result.rows(); ++row) {
    for (int col=0; col<result.cols(); ++col) {
      PixelMask<Vector2i> val = result(col,row);
      if((val[0]==2) && (val[1]==1))
        ++num_correct;
    }
  }
  double percent_correct = static_cast<double>(num_correct) / num_pixels;
  EXPECT_GT(percent_correct, 0.99);
}

