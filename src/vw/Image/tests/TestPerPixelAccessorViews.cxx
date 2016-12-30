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


// TestPerPixelAccessorView.h
#include <gtest/gtest_VW.h>
#include <test/Helpers.h>

#include <vw/Image/PerPixelAccessorViews.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/PixelIterator.h>
#include <vw/Core/Functors.h>

using namespace vw;

template <template<class> class TraitT, class T>
static bool bool_trait( T const& /*arg*/ ) {
  return TraitT<T>::value;
}

template <class T1, class T2>
static bool has_pixel_type( T2 ) {
  return boost::is_same<T1,typename T2::pixel_type>::value;
}

// Pixel accessor function
struct AccessorNegationFunctor : public ReturnFixedType<float> {
  BBox2i work_area() const { return BBox2i(0,0,1,1); }

  template <class PixelAccessorT>
  typename boost::remove_reference<typename PixelAccessorT::result_type>::type
  operator()(PixelAccessorT acc) const {
    return -(*acc);
  }
};

TEST( PerPixelAccessorViews, Unary ) {
  ImageView<float> im(2,2); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4;
  UnaryPerPixelAccessorView<ImageView<float>, AccessorNegationFunctor> ppv( im, AccessorNegationFunctor() );
  ASSERT_EQ( ppv.cols(), 2 );
  ASSERT_EQ( ppv.rows(), 2 );
  ASSERT_EQ( ppv.planes(), 1 );

  // Test individual pixel access
  EXPECT_EQ( ppv(0,0), -1 );
  EXPECT_EQ( ppv(1,0), -2 );
  EXPECT_EQ( ppv(0,1), -3 );
  EXPECT_EQ( ppv(1,1), -4 );

  // Test full rasterizaion
  ImageView<float> im2 = ppv;
  EXPECT_EQ( im2.cols(), 2 );
  EXPECT_EQ( im2.rows(), 2 );
  EXPECT_EQ( im2.planes(), 1 );
  EXPECT_EQ( im2(0,0), ppv(0,0) );
  EXPECT_EQ( im2(1,0), ppv(1,0) );
  EXPECT_EQ( im2(0,1), ppv(0,1) );
  EXPECT_EQ( im2(1,1), ppv(1,1) );

  // Test partial rasterization
  ImageView<float> im3(1,2);
  ASSERT_NO_THROW( ppv.rasterize( im3, BBox2i(1,0,1,2) ) );
  EXPECT_EQ( im3(0,0), ppv(1,0) );
  EXPECT_EQ( im3(0,1), ppv(1,1) );
  ImageView<float> im4(2,1);
  ASSERT_NO_THROW( ppv.rasterize( im4, BBox2i(0,1,2,1) ) );
  EXPECT_EQ( im4(0,0), ppv(0,1) );
  EXPECT_EQ( im4(1,0), ppv(1,1) );

  // Test the accessor / generic rasterization
  ImageView<float> im5(2,2);
  vw::rasterize( ppv, im5, BBox2i(0,0,2,2) );
  EXPECT_EQ( im5(0,0), ppv(0,0) );
  EXPECT_EQ( im5(1,0), ppv(1,0) );
  EXPECT_EQ( im5(0,1), ppv(0,1) );
  EXPECT_EQ( im5(1,1), ppv(1,1) );

  // Test the iterator
  ImageView<float>::iterator im2i = im2.begin();
  UnaryPerPixelAccessorView<ImageView<float>, AccessorNegationFunctor>::iterator ppvi = ppv.begin();
  for( int i=0; i<4; ++i ) {
    EXPECT_EQ( *im2i, *ppvi );
    ASSERT_NO_THROW( ++ppvi );
    ++im2i;
  }

  // Test the types
  ASSERT_TRUE( has_pixel_type<float>( ppv ) );
  ASSERT_FALSE( bool_trait<IsMultiplyAccessible>(ppv) );
  ASSERT_TRUE( bool_trait<IsImageView>(ppv) );
}

TEST(PixelIterator, PixelLineIterator) {

  // Create a line iterator from the middle left to the middle top
  PixelLineIterator line(Vector2i(0,10), PixelLineIterator::TR, Vector2i(20,20));
  
  EXPECT_VECTOR_EQ(Vector2i(0,10), *line); line++;
  EXPECT_VECTOR_EQ(Vector2i(1, 9), *line); line++;
  EXPECT_VECTOR_EQ(Vector2i(2, 8), *line); line++;
  EXPECT_VECTOR_EQ(Vector2i(3, 7), *line); line++;
  line++;
  line++;
  line++;
  line++;
  line++;
  line++;
  EXPECT_VECTOR_EQ(Vector2i(10, 0), *line); line++;
  EXPECT_FALSE(line.is_good());
}


