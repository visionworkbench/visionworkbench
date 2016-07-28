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

  int dx = 20;
  int dy = 3;
  /*
  // For this example, the disparity is (2,1) for each pixel!
  DiskImageView<PixelGray<uint8> > inputLeft ("left.tif");
  DiskImageView<PixelGray<uint8> > inputRight("left_shift.tif");
  BBox2i leftRoi (0,0,400, 400);
  BBox2i rightRoi(0,0,400+dx, 400+dy);
*/
/*
  DiskImageView<PixelGray<uint8> > inputLeft ("left4.tif");
  DiskImageView<PixelGray<uint8> > inputRight("right4.tif");
  BBox2i leftRoi (0,0,80, 80);
  BBox2i rightRoi(0,0,80+dx, 80+dy);
*/
/*
  DiskImageView<PixelGray<uint8> > inputLeft ("left3.tif");
  DiskImageView<PixelGray<uint8> > inputRight("right3.tif");
  BBox2i leftRoi (0,0,200, 200);
  BBox2i rightRoi(0,0,200+dx, 200+dy);
*/

  DiskImageView<PixelGray<uint8> > inputLeft ("moc_right.tif");
  DiskImageView<PixelGray<uint8> > inputRight("moc_left.tif");
  BBox2i leftRoi (0,0,300, 300);
  BBox2i rightRoi(0,0,300+dx, 300+dy);


/*// Note that the matches have higher column values in the left image.
  // X disparities are in the 0-140 range, y is 0 disparity.
  DiskImageView<PixelGray<uint8> > inputLeft ("cones_right.pgm");
  DiskImageView<PixelGray<uint8> > inputRight("cones_left.pgm");
  BBox2i leftRoi (0,0,700, 700);
  BBox2i rightRoi(0,0,700+dx, 700+dy);
*/


  ImageView<uint8> left  = crop(inputLeft, leftRoi);
  ImageView<uint8> right = crop(inputRight, rightRoi);
  
  //write_image("left_in.tif", left);
  //write_image("right_in.tif", right);
  
  ImageView<Vector<uint8, 2> > result = 
  semi_global_matching_func(left, right);
  std::cout << "Writing output...\n";
  write_image("SGM_output.tif", result);
  
  ASSERT_TRUE(false);
}

