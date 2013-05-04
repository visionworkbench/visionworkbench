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

#include <vw/Stereo/PreFilter.h>

using namespace vw;
using namespace vw::stereo;

template <typename PixelT>
class PreProcess : public ::testing::Test {
protected:
  PreProcess() {}

  typedef ImageView<PixelT> image_type;
  image_type input;

  virtual void SetUp() {
    input.set_size(4,4);
    int32 count = 0;
    for ( PixelT* p = &input(0,0); p < &input(3,3)+1; p++ ) {
      count++;
      *p = PixelT(count);
    }
  }
};

typedef PreProcess<PixelGray<uint8> > PreProcessGRAYU8;
typedef PreProcess<PixelGray<int16> > PreProcessGRAYI16;
typedef PreProcess<PixelGray<float> > PreProcessGRAYF32;

TEST_F( PreProcessGRAYU8, NullFilter ) {
  image_type result =
    NullOperation().filter( input );
  EXPECT_VW_EQ( input, result );
}

TEST_F( PreProcessGRAYU8, LOGFilter ) {
  image_type result =
    LaplacianOfGaussian(25).filter( input );
  image_type truth =
    laplacian_filter(gaussian_filter( input, 25 ));
  EXPECT_VW_EQ( truth, result );
}

TEST_F( PreProcessGRAYU8, ZeroMeanFilter ) {
  image_type result =
    SubtractedMean(25).filter( input );
  image_type truth =
    input - gaussian_filter( input, 25 );
  EXPECT_VW_EQ( truth, result );
}

TEST_F( PreProcessGRAYI16, NullFilter ) {
  image_type result =
    NullOperation().filter( input );
  EXPECT_VW_EQ( input, result );
}

TEST_F( PreProcessGRAYI16, LOGFilter ) {
  image_type result =
    LaplacianOfGaussian(25).filter( input );
  image_type truth =
    laplacian_filter(gaussian_filter( input, 25 ));
  EXPECT_VW_EQ( truth, result );
}

TEST_F( PreProcessGRAYI16, ZeroMeanFilter ) {
  image_type result =
    SubtractedMean(25).filter( input );
  image_type truth =
    input - gaussian_filter( input, 25 );
  EXPECT_VW_EQ( truth, result );
}

TEST_F( PreProcessGRAYF32, NullFilter ) {
  image_type result =
    NullOperation().filter( input );
  EXPECT_VW_EQ( input, result );
}

TEST_F( PreProcessGRAYF32, LOGFilter ) {
  image_type result =
    LaplacianOfGaussian(25).filter( input );
  image_type truth =
    laplacian_filter(gaussian_filter( input, 25 ));
  EXPECT_VW_EQ( truth, result );
}

TEST_F( PreProcessGRAYF32, ZeroMeanFilter ) {
  image_type result =
    SubtractedMean(25).filter( input );
  image_type truth =
    input - gaussian_filter( input, 25 );
  EXPECT_VW_EQ( truth, result );
}
