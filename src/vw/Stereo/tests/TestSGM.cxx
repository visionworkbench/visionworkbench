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
#include <vw/Stereo/SGM.h>

/*
// These are used to read and write tif images with vector pixels
namespace vw {
  template<> struct PixelFormatID<Vector<uint8,  2> >  { static const PixelFormatEnum value = VW_PIXEL_GENERIC_2_CHANNEL; };
}
*/
using namespace vw;
using namespace vw::stereo;

TEST( SGM, basic ) {
  
/*
  // For this perfect test case, the correct disparity is (2,1) for each pixel!
  int min_disp_x  = -4;
  int max_disp_x  =  4;
  int min_disp_y  =  -3;
  int max_disp_y  =   3;
  int kernel_size = 1;
  DiskImageView<PixelGray<uint8> > inputLeft ("left.tif");
  DiskImageView<PixelGray<uint8> > inputRight("left_shift.tif");
  BBox2i leftRoi (0,0,400, 400);
  BBox2i rightRoi(0,0,400, 400);
*/
/*
  DiskImageView<PixelGray<uint8> > inputLeft ("left4.tif");
  DiskImageView<PixelGray<uint8> > inputRight("right4.tif");
  BBox2i leftRoi (0,0,80, 80);
  BBox2i rightRoi(0,0,80+dx, 80+dy);
*/
/*
  int min_disp_x  = -2;
  int max_disp_x  = 10;
  int min_disp_y  =  -2;
  int max_disp_y  =  12;
  int kernel_size = 1;
  DiskImageView<PixelGray<uint8> > inputLeft ("left3.tif");
  DiskImageView<PixelGray<uint8> > inputRight("right3.tif");
  BBox2i leftRoi (0,0,200, 200);
  BBox2i rightRoi(0,0,200, 200);
*/

  // The results may actually be pretty good except for the border region.
  int min_disp_x  = -5;
  int max_disp_x  =  5;
  int min_disp_y  =  -5;
  int max_disp_y  =   5;
  int kernel_size = 1;
  DiskImageView<PixelGray<uint8> > inputLeft ("moc_left_u8.tif");
  DiskImageView<PixelGray<uint8> > inputRight("moc_right_u8.tif");
  BBox2i leftRoi (0,0,340, 340);
  BBox2i rightRoi(0,0,340, 340);

/*
  // Note that the matches have higher column values in the left image.
  // X disparities are in the 0-140 range, y is 0 disparity.
  // - Results are similar to the 
  int min_disp_x  = 30;
  int max_disp_x  = 140;
  int min_disp_y  = 0;
  int max_disp_y  = 0;
  int kernel_size = 1;
  DiskImageView<PixelGray<uint8> > inputLeft ("cones_right.pgm");
  DiskImageView<PixelGray<uint8> > inputRight("cones_left.pgm");
  BBox2i leftRoi (0,0,700, 200);
  BBox2i rightRoi(0,0,700+max_disp_x, 200+max_disp_y);
*/


  ImageView<uint8> left  = crop(inputLeft, leftRoi);
  ImageView<uint8> right = crop(inputRight, rightRoi);
  
  write_image("left_in.tif", left);
  write_image("right_in.tif", right);
  
  
  
  SemiGlobalMatcher matcher;
  matcher.setParameters(min_disp_x, min_disp_y, max_disp_x, max_disp_y, kernel_size);
  SemiGlobalMatcher::DisparityImage result = matcher.semi_global_matching_func(left, right);
  
  std::cout << "Writing output...\n";
  write_image("SGM_output.tif", result);
  
  ASSERT_TRUE(false);
}

