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


// TestDisparity.h
#include <test/Helpers.h>

#include <vw/Stereo/Correlate.h>

using namespace vw;
using namespace vw::stereo;

typedef PixelMask<Vector2i> PixelDisp;

TEST( Correlate, CrossCorrConsistency ) {

  ImageView<PixelDisp> r2l(3,3), l2r(3,3);
  fill( r2l, PixelDisp(Vector2i(0,0)) );
  fill( l2r, PixelDisp(Vector2i(0,0)) );

  fill( crop(l2r,2,0,1,3), PixelDisp(Vector2i(2,2)) );
  l2r(0,0) = PixelDisp(Vector2i(1,1));
  r2l(1,1) = PixelDisp(Vector2i(-1,-1));
  l2r(1,0) = PixelDisp(Vector2i(1,1));

  ImageView<PixelDisp> l2r_copy = copy(l2r);
  cross_corr_consistency_check( l2r_copy, r2l, 0 );
  EXPECT_FALSE( is_valid(l2r_copy(2,0)) );
  EXPECT_FALSE( is_valid(l2r_copy(2,1)) );
  EXPECT_FALSE( is_valid(l2r_copy(2,2)) );
  EXPECT_TRUE( is_valid( l2r_copy(0,0)) );
  EXPECT_FALSE( is_valid(l2r_copy(1,0)) );

  l2r_copy = l2r;
  cross_corr_consistency_check( l2r_copy, r2l, 2 );
  EXPECT_FALSE( is_valid(l2r_copy(2,0)) );
  EXPECT_FALSE( is_valid(l2r_copy(2,1)) );
  EXPECT_FALSE( is_valid(l2r_copy(2,2)) );
  EXPECT_TRUE( is_valid( l2r_copy(0,0)) );
  EXPECT_TRUE( is_valid( l2r_copy(1,0)) );
}
