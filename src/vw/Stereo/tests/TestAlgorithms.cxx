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

#include <vw/Stereo/Algorithms.h>
#include <vw/Stereo/CostFunctions.h>

#include <test/Helpers.h>

#include <algorithm>

using namespace vw;
using namespace vw::stereo;

template <class T1, class T2>
static bool has_pixel_type( T2 ) {
  return boost::is_same<T1,typename T2::pixel_type>::value;
}

template <class T1, template<class> class AccuT, class T2>
static bool has_accum_type( T2 ) {
  return boost::is_same<T1, typename AccuT<T2>::type>::value;
}

TEST( AlgorithmsTest, AbsAccumulatorType ) {
  EXPECT_TRUE( (has_accum_type<double,AbsAccumulatorType>(ImageView<double>())) );
  EXPECT_TRUE( (has_accum_type<float,AbsAccumulatorType>( ImageView<float>())) );
  EXPECT_TRUE( (has_accum_type<int16,AbsAccumulatorType>( ImageView<uint8>())) );
  EXPECT_TRUE( (has_accum_type<int32,AbsAccumulatorType>( ImageView<uint32>())) );
}

TEST( AlgorithmsTest, FastBoxFloat ) {
  ImageView<float> input(7,5);
  {
    size_t count = 1;
    for ( PixelIterator<ImageView<float> > i = input.begin();
          i != input.end(); i++ ) {
      *i = count;
      count++;
    }
  }

  typedef AbsAccumulatorType<float>::type accum_type;

  EXPECT_EQ( 3, fast_box_sum<accum_type>( input, Vector2i(5,3) ).cols() );
  EXPECT_EQ( 3, fast_box_sum<accum_type>( input, Vector2i(5,3) ).rows() );
  ASSERT_TRUE( has_pixel_type<float>( fast_box_sum<accum_type>( input, Vector2i(5,3) ) ) );
  EXPECT_EQ( 5, fast_box_sum<accum_type>( input, Vector2i(3,3) ).cols() );
  EXPECT_EQ( 3, fast_box_sum<accum_type>( input, Vector2i(3,3) ).rows() );

  ImageView<float> output = fast_box_sum<accum_type>(input, Vector2i(5,3));
  EXPECT_EQ( 150, output(0,0) );
  EXPECT_EQ( 165, output(1,0) );
  EXPECT_EQ( 180, output(2,0) );
  EXPECT_EQ( 360, output(0,2) );
  EXPECT_EQ( 375, output(1,2) );
  EXPECT_EQ( 390, output(2,2) );
}

TEST( AlgorithmsTest, FastBoxInt32 ) {
  ImageView<int32> input(7,5);
  {
    size_t count = 1;
    for ( PixelIterator<ImageView<int32> > i = input.begin();
          i != input.end(); i++ ) {
      *i = count;
      count++;
    }
  }

  typedef AbsAccumulatorType<int32>::type accum_type;

  EXPECT_EQ( 3, fast_box_sum<accum_type>( input, Vector2i(5,3) ).cols() );
  EXPECT_EQ( 3, fast_box_sum<accum_type>( input, Vector2i(5,3) ).rows() );
  ASSERT_TRUE( has_pixel_type<int32>( fast_box_sum<accum_type>( input, Vector2i(5,3) ) ) );
  EXPECT_EQ( 5, fast_box_sum<accum_type>( input, Vector2i(3,3) ).cols() );
  EXPECT_EQ( 3, fast_box_sum<accum_type>( input, Vector2i(3,3) ).rows() );

  ImageView<int32> output = fast_box_sum<accum_type>(input, Vector2i(5,3));
  EXPECT_EQ( 150, output(0,0) );
  EXPECT_EQ( 165, output(1,0) );
  EXPECT_EQ( 180, output(2,0) );
  EXPECT_EQ( 360, output(0,2) );
  EXPECT_EQ( 375, output(1,2) );
  EXPECT_EQ( 390, output(2,2) );
}

TEST( AlgorithmsTest, FastBoxDouble ) {
  ImageView<double> input(7,5);
  std::fill( input.begin(), input.end(), 1 );
  std::fill( &input(0,0), &input(6,0)+1, 2 );
  input(6,0) = 3;

  typedef AbsAccumulatorType<double>::type accum_type;

  EXPECT_EQ( 3, fast_box_sum<accum_type>( input, Vector2i(5,3) ).cols() );
  EXPECT_EQ( 3, fast_box_sum<accum_type>( input, Vector2i(5,3) ).rows() );
  ASSERT_TRUE( has_pixel_type<double>( fast_box_sum<accum_type>( input, Vector2i(5,3) ) ) );
  EXPECT_EQ( 5, fast_box_sum<accum_type>( input, Vector2i(3,3) ).cols() );
  EXPECT_EQ( 3, fast_box_sum<accum_type>( input, Vector2i(3,3) ).rows() );

  ImageView<double> output = fast_box_sum<accum_type>(input, Vector2i(3,3));
  EXPECT_EQ( 12, output(0,0) );
  EXPECT_EQ( 9,  output(0,1) );
  EXPECT_EQ( 9,  output(0,2) );
  EXPECT_EQ( 12, output(1,0) );
  EXPECT_EQ( 12, output(2,0) );
  EXPECT_EQ( 12, output(3,0) );
  EXPECT_EQ( 13, output(4,0) );
}

TEST( AlgorithmsTest, FastBoxChar ) {
  ImageView<uint8> input(7,5);
  std::fill( input.begin(), input.end(), 30 ); // Picking values that
                                               // would overflow input
  std::fill( &input(0,0), &input(6,0)+1, 40 );
  input(6,0) = 50;

  typedef AbsAccumulatorType<uint8>::type accum_type;

  EXPECT_EQ( 3, fast_box_sum<accum_type>( input, Vector2i(5,3) ).cols() );
  EXPECT_EQ( 3, fast_box_sum<accum_type>( input, Vector2i(5,3) ).rows() );
  ASSERT_TRUE( has_pixel_type<int16>( fast_box_sum<accum_type>( input, Vector2i(5,3) ) ) );
  EXPECT_EQ( 5, fast_box_sum<accum_type>( input, Vector2i(3,3) ).cols() );
  EXPECT_EQ( 3, fast_box_sum<accum_type>( input, Vector2i(3,3) ).rows() );

  ImageView<int16> output = fast_box_sum<accum_type>(input, Vector2i(5,3));
  EXPECT_EQ( 450, output(0,2) );
  EXPECT_EQ( 450, output(2,2) );
  EXPECT_EQ( 500, output(0,0) );
  EXPECT_EQ( 510, output(2,0) );
  EXPECT_EQ( 450, output(1,1) );

  // Really force the fact that we are not over flowing
  std::fill( input.begin(), input.end(), 255 );
  EXPECT_EQ( 6375, fast_box_sum<accum_type>( input, Vector2i(5,5) )(0,0) );
}

TEST( AlgorithmsTest, FastBoxPixelU8 ) {
  ImageView<PixelGray<uint8> > input(7,5);
  std::fill( input.begin(), input.end(), PixelGray<uint8>(27) );
  std::fill( &input(0,0), &input(6,0)+1, PixelGray<uint8>(40) );

  typedef AbsAccumulatorType<uint8>::type accum_type;

  EXPECT_EQ( 3, fast_box_sum<accum_type>( input, Vector2i(5,3) ).cols() );
  EXPECT_EQ( 3, fast_box_sum<accum_type>( input, Vector2i(5,3) ).rows() );
  ASSERT_TRUE( has_pixel_type<PixelGray<int16> >( fast_box_sum<accum_type>( input, Vector2i(5,3) ) ) );
  EXPECT_EQ( 5, fast_box_sum<accum_type>( input, Vector2i(3,3) ).cols() );
  EXPECT_EQ( 3, fast_box_sum<accum_type>( input, Vector2i(3,3) ).rows() );

  ImageView<PixelGray<int16> > output = fast_box_sum<accum_type>( input, Vector2i(5,3) );
  EXPECT_VW_EQ( PixelGray<int16>(470), output(0,0) );
  EXPECT_VW_EQ( PixelGray<int16>(470), output(1,0) );
  EXPECT_VW_EQ( PixelGray<int16>(470), output(2,0) );
  EXPECT_VW_EQ( PixelGray<int16>(405), output(0,1) );
  EXPECT_VW_EQ( PixelGray<int16>(405), output(1,1) );
  EXPECT_VW_EQ( PixelGray<int16>(405), output(2,1) );
  EXPECT_VW_EQ( PixelGray<int16>(405), output(0,2) );
  EXPECT_VW_EQ( PixelGray<int16>(405), output(2,2) );
}
