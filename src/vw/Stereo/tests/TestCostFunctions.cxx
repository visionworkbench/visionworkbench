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

#include <vw/Image/ImageView.h>
#include <vw/Stereo/CostFunctions.h>

using namespace vw;
using namespace vw::stereo;

template <typename PixelT>
class CostFunction : public ::testing::Test {
protected:
  CostFunction() {}

  typedef PixelT pixel_type;
  typedef ImageView<PixelT> image_type;
  typedef typename PixelChannelType<PixelT>::type channel_type;
  image_type input1, input2;
  Vector2i kernel_size;

  virtual void SetUp() {
    kernel_size = Vector2i(1,1);
    input1.set_size(2,2); input2.set_size(2,2);
    input1(0,0) = input2(0,0) = PixelT(128);
    input1(1,0) = input2(0,1) = PixelT(10);
    input2(1,0) = input1(0,1) = PixelT(100);
    input1(1,1) = PixelT(0);
    input2(1,1) = PixelT(ChannelRange<channel_type>::max());
  }
};

typedef CostFunction<PixelGray<uint8> > CostFunctionGRAYU8;
typedef CostFunction<PixelGray<int16> > CostFunctionGRAYI16;
typedef CostFunction<PixelGray<float> > CostFunctionGRAYF32;

typedef CostFunction<uint8> CostFunctionU8;
typedef CostFunction<PixelRGB<uint8> > CostFunctionRGBU8;

TEST_F( CostFunctionGRAYU8, AbsDiff ) {
  typedef PixelChannelCast<pixel_type,AbsAccumulatorType<PixelChannelType<pixel_type>::type>::type>::type rpixel_type;
  typedef ImageView<rpixel_type> result_type;
  result_type result = AbsoluteCost<image_type,boost::is_integral<PixelChannelType<pixel_type>::type>::value>(input1, input2, kernel_size)( input1, input2 );
  EXPECT_VW_EQ( rpixel_type( 0 ), result(0,0) );
  EXPECT_VW_EQ( rpixel_type( 90 ), result(0,1) );
  EXPECT_VW_EQ( rpixel_type( 90 ), result(1,0) );
  EXPECT_VW_EQ( rpixel_type( ChannelRange<channel_type>::max() ), result(1,1) );
}

TEST_F( CostFunctionU8, AbsDiff ) {
  typedef PixelChannelCast<pixel_type,AbsAccumulatorType<PixelChannelType<pixel_type>::type>::type>::type rpixel_type;
  typedef ImageView<rpixel_type> result_type;
  result_type result = AbsoluteCost<image_type,boost::is_integral<PixelChannelType<pixel_type>::type>::value>(input1, input2, kernel_size)( input1, input2 );
  EXPECT_VW_EQ( rpixel_type( 0 ), result(0,0) );
  EXPECT_VW_EQ( rpixel_type( 90 ), result(0,1) );
  EXPECT_VW_EQ( rpixel_type( 90 ), result(1,0) );
  EXPECT_VW_EQ( rpixel_type( ChannelRange<channel_type>::max() ), result(1,1) );
}

TEST_F( CostFunctionRGBU8, AbsDiff ) {
  typedef PixelChannelCast<pixel_type,AbsAccumulatorType<PixelChannelType<pixel_type>::type>::type>::type rpixel_type;
  typedef ImageView<rpixel_type> result_type;
  result_type result = AbsoluteCost<image_type,boost::is_integral<PixelChannelType<pixel_type>::type>::value>(input1, input2, kernel_size)( input1, input2 );
  EXPECT_VW_EQ( rpixel_type( 0 ), result(0,0) );
  EXPECT_VW_EQ( rpixel_type( 90 ), result(0,1) );
  EXPECT_VW_EQ( rpixel_type( 90 ), result(1,0) );
  EXPECT_VW_EQ( rpixel_type( ChannelRange<channel_type>::max() ), result(1,1) );
}

TEST_F( CostFunctionGRAYI16, AbsDiff ) {
  typedef PixelChannelCast<pixel_type,AbsAccumulatorType<PixelChannelType<pixel_type>::type>::type>::type rpixel_type;
  typedef ImageView<rpixel_type> result_type;
  result_type result = AbsoluteCost<image_type,boost::is_integral<PixelChannelType<pixel_type>::type>::value>(input1, input2, kernel_size)( input1, input2 );
  EXPECT_VW_EQ( rpixel_type( 0 ), result(0,0) );
  EXPECT_VW_EQ( rpixel_type( 90 ), result(0,1) );
  EXPECT_VW_EQ( rpixel_type( 90 ), result(1,0) );
  EXPECT_VW_EQ( rpixel_type( ChannelRange<channel_type>::max() ), result(1,1) );
}

TEST_F( CostFunctionGRAYF32, AbsDiff ) {
  typedef PixelChannelCast<pixel_type,AbsAccumulatorType<PixelChannelType<pixel_type>::type>::type>::type rpixel_type;
  typedef ImageView<rpixel_type> result_type;
  result_type result = AbsoluteCost<image_type,boost::is_integral<PixelChannelType<pixel_type>::type>::value>(input1, input2, kernel_size)( input1, input2 );
  EXPECT_VW_EQ( rpixel_type( 0 ), result(0,0) );
  EXPECT_VW_EQ( rpixel_type( 90 ), result(0,1) );
  EXPECT_VW_EQ( rpixel_type( 90 ), result(1,0) );
  EXPECT_VW_EQ( rpixel_type( ChannelRange<channel_type>::max() ), result(1,1) );
}

TEST_F( CostFunctionGRAYU8, SquaredDiff ) {
  typedef PixelChannelCast<pixel_type,SqrDiffAccumulatorType<PixelChannelType<pixel_type>::type>::type>::type rpixel_type;
  typedef ImageView<rpixel_type> result_type;
  result_type result = SquaredCost<image_type,boost::is_integral<PixelChannelType<pixel_type>::type>::value>(input1, input2, kernel_size)( input1, input2 );
  EXPECT_VW_EQ( rpixel_type( 0 ), result(0,0) );
  EXPECT_VW_EQ( rpixel_type( 8100 ), result(0,1) );
  EXPECT_VW_EQ( rpixel_type( 8100 ), result(1,0) );
  EXPECT_VW_EQ( rpixel_type( pow(ChannelRange<channel_type>::max(),2) ), result(1,1) );
}

TEST_F( CostFunctionU8, SquaredDiff ) {
  typedef PixelChannelCast<pixel_type,SqrDiffAccumulatorType<PixelChannelType<pixel_type>::type>::type>::type rpixel_type;
  typedef ImageView<rpixel_type> result_type;
  result_type result = SquaredCost<image_type,boost::is_integral<PixelChannelType<pixel_type>::type>::value>(input1, input2, kernel_size)( input1, input2 );
  EXPECT_VW_EQ( rpixel_type( 0 ), result(0,0) );
  EXPECT_VW_EQ( rpixel_type( 8100 ), result(0,1) );
  EXPECT_VW_EQ( rpixel_type( 8100 ), result(1,0) );
  EXPECT_VW_EQ( rpixel_type( pow(ChannelRange<channel_type>::max(),2) ), result(1,1) );

  ImageView<uint8> input3(2,1), input4(2,1);
  input3(0,0) = input4(1,0) = 0;
  input3(1,0) = input4(0,0) = 255;
  result_type result2 =
    SquaredCost<image_type,true>(input3,input4,Vector2i(2,1))(input3,input4);
  EXPECT_EQ( 65025, result2(0,0) );
  EXPECT_EQ( 65025, result2(1,0) );
}

TEST_F( CostFunctionRGBU8, SquaredDiff ) {
  typedef PixelChannelCast<pixel_type,SqrDiffAccumulatorType<PixelChannelType<pixel_type>::type>::type>::type rpixel_type;
  typedef ImageView<rpixel_type> result_type;
  result_type result = SquaredCost<image_type,boost::is_integral<PixelChannelType<pixel_type>::type>::value>(input1, input2, kernel_size)( input1, input2 );
  EXPECT_VW_EQ( rpixel_type( 0 ), result(0,0) );
  EXPECT_VW_EQ( rpixel_type( 8100 ), result(0,1) );
  EXPECT_VW_EQ( rpixel_type( 8100 ), result(1,0) );
  EXPECT_VW_EQ( rpixel_type( pow(ChannelRange<channel_type>::max(),2) ), result(1,1) );
}

TEST_F( CostFunctionGRAYI16, SquaredDiff ) {
  typedef PixelChannelCast<pixel_type,SqrDiffAccumulatorType<PixelChannelType<pixel_type>::type>::type>::type rpixel_type;
  typedef ImageView<rpixel_type> result_type;
  result_type result = SquaredCost<image_type,boost::is_integral<PixelChannelType<pixel_type>::type>::value>(input1, input2, kernel_size)( input1, input2 );
  EXPECT_VW_EQ( rpixel_type( 0 ), result(0,0) );
  EXPECT_VW_EQ( rpixel_type( 8100 ), result(0,1) );
  EXPECT_VW_EQ( rpixel_type( 8100 ), result(1,0) );
  EXPECT_VW_EQ( rpixel_type( pow(ChannelRange<channel_type>::max(),2) ), result(1,1) );
}

TEST_F( CostFunctionGRAYF32, SquaredDiff ) {
  typedef PixelChannelCast<pixel_type,SqrDiffAccumulatorType<PixelChannelType<pixel_type>::type>::type>::type rpixel_type;
  typedef ImageView<rpixel_type> result_type;
  result_type result = SquaredCost<image_type,boost::is_integral<PixelChannelType<pixel_type>::type>::value>(input1, input2, kernel_size)( input1, input2 );
  EXPECT_VW_EQ( rpixel_type( 0 ), result(0,0) );
  EXPECT_VW_EQ( rpixel_type( 8100 ), result(0,1) );
  EXPECT_VW_EQ( rpixel_type( 8100 ), result(1,0) );
  EXPECT_VW_EQ( rpixel_type( pow(ChannelRange<channel_type>::max(),2) ), result(1,1) );
}

TEST_F( CostFunctionGRAYU8, CrossCorrelation ) {
  typedef PixelChannelCast<pixel_type,SqrDiffAccumulatorType<PixelChannelType<pixel_type>::type>::type>::type rpixel_type;
  typedef ImageView<rpixel_type> result_type;
  result_type result = NCCCost<image_type,boost::is_integral<PixelChannelType<pixel_type>::type>::value>(input1, input2, kernel_size)( input1, input2 );
  EXPECT_VW_EQ( rpixel_type( 16384/256 ), result(0,0) );
  EXPECT_VW_EQ( rpixel_type( 1000/256 ), result(0,1) );
  EXPECT_VW_EQ( rpixel_type( 1000/256 ), result(1,0) );
  EXPECT_VW_EQ( rpixel_type( 0 ), result(1,1) );
}

TEST_F( CostFunctionU8, CrossCorrelation ) {
  typedef PixelChannelCast<pixel_type,SqrDiffAccumulatorType<PixelChannelType<pixel_type>::type>::type>::type rpixel_type;
  typedef ImageView<rpixel_type> result_type;
  result_type result = NCCCost<image_type,boost::is_integral<PixelChannelType<pixel_type>::type>::value>(input1, input2, kernel_size)( input1, input2 );
  EXPECT_VW_EQ( rpixel_type( 16384/256 ), result(0,0) );
  EXPECT_VW_EQ( rpixel_type( 1000/256 ), result(0,1) );
  EXPECT_VW_EQ( rpixel_type( 1000/256 ), result(1,0) );
  EXPECT_VW_EQ( rpixel_type( 0 ), result(1,1) );
}

TEST_F( CostFunctionRGBU8, CrossCorrelation ) {
  typedef PixelChannelCast<pixel_type,SqrDiffAccumulatorType<PixelChannelType<pixel_type>::type>::type>::type rpixel_type;
  typedef ImageView<rpixel_type> result_type;
  result_type result = NCCCost<image_type,boost::is_integral<PixelChannelType<pixel_type>::type>::value>(input1, input2, kernel_size)( input1, input2 );
  EXPECT_VW_EQ( rpixel_type( 16384/256 ), result(0,0) );
  EXPECT_VW_EQ( rpixel_type( 1000/256 ), result(0,1) );
  EXPECT_VW_EQ( rpixel_type( 1000/256 ), result(1,0) );
  EXPECT_VW_EQ( rpixel_type( 0 ), result(1,1) );
}

TEST_F( CostFunctionGRAYI16, CrossCorrelation ) {
  typedef PixelChannelCast<pixel_type,SqrDiffAccumulatorType<PixelChannelType<pixel_type>::type>::type>::type rpixel_type;
  typedef ImageView<rpixel_type> result_type;
  result_type result = NCCCost<image_type,boost::is_integral<PixelChannelType<pixel_type>::type>::value>(input1, input2, kernel_size)( input1, input2 );
  EXPECT_VW_EQ( rpixel_type( 16384/256 ), result(0,0) );
  EXPECT_VW_EQ( rpixel_type( 1000/256 ), result(0,1) );
  EXPECT_VW_EQ( rpixel_type( 1000/256 ), result(1,0) );
  EXPECT_VW_EQ( rpixel_type( 0 ), result(1,1) );
}

TEST_F( CostFunctionGRAYF32, CrossCorrelation ) {
  typedef PixelChannelCast<pixel_type,SqrDiffAccumulatorType<PixelChannelType<pixel_type>::type>::type>::type rpixel_type;
  typedef ImageView<rpixel_type> result_type;
  result_type result = NCCCost<image_type,boost::is_integral<PixelChannelType<pixel_type>::type>::value>(input1, input2, kernel_size)( input1, input2 );
  EXPECT_VW_EQ( rpixel_type( 16384 ), result(0,0) );
  EXPECT_VW_EQ( rpixel_type( 1000 ), result(0,1) );
  EXPECT_VW_EQ( rpixel_type( 1000 ), result(1,0) );
  EXPECT_VW_EQ( rpixel_type( 0 ), result(1,1) );
}

TEST( CostFunction, AccumulatorType ) {
  // I'm mostly just checking to make sure they exist
  EXPECT_TRUE( (boost::is_same<float,AbsoluteCost<ImageView<float>, false >::accumulator_type>::value) );
  EXPECT_TRUE( (boost::is_same<double,SquaredCost<ImageView<float>, false >::accumulator_type>::value) );
  EXPECT_TRUE( (boost::is_same<double,NCCCost<ImageView<float>, false>::accumulator_type>::value) );
}
