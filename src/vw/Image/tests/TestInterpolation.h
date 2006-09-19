// __BEGIN_LICENSE__
//
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
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

// TestEdgeExtend.h
#include <cxxtest/TestSuite.h>

#include <vw/Image/ImageView.h>
#include <vw/Image/Interpolation.h>

using namespace vw;

class TestInterpolation : public CxxTest::TestSuite
{
public:

  template <template<class> class TraitT, class T>
  static bool bool_trait( T const& arg ) {
    return TraitT<T>::value;
  }

  void testBilinearInterpolation()
  {
    ImageView<double> im(2,3); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4; im(0,2)=5; im(1,2)=6;
    InterpolationView<ImageView<double>, BilinearInterpolation> im2 = bilinear_interpolation(im);
    TS_ASSERT_EQUALS( im2.cols(), 2 );
    TS_ASSERT_EQUALS( im2.rows(), 3 );
    TS_ASSERT_EQUALS( im2(1,1), 4 );
    TS_ASSERT_EQUALS( im2(1,1.5), 5 );
    TS_ASSERT_EQUALS( im2(0.5,1), 3.5 );
    TS_ASSERT_EQUALS( im2(0.5,0.5), 2.5 );

    // Teste accessor
    TS_ASSERT_EQUALS( *(im2.origin().advance(1,1.5)), 5 );
    TS_ASSERT_EQUALS( *(im2.origin().advance(1,1)), 4 );

    // Test the traits
    TS_ASSERT( bool_trait<IsFloatingPointIndexable>( im2 ) ); 
    TS_ASSERT( bool_trait<IsFloatingPointIndexable>( bilinear_interpolation(im) ) );
    TS_ASSERT( !bool_trait<IsReferenceable>( bilinear_interpolation(im) ) );
    TS_ASSERT( !bool_trait<IsMultiplyAccessible>( bilinear_interpolation(im) ) );
    TS_ASSERT( !bool_trait<IsReferenceable>( bilinear_interpolation(im) ) );
  }


  void testBicubicInterpolation()
  {
    ImageView<double> im(2,3); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4; im(0,2)=5; im(1,2)=6;
    InterpolationView<ImageView<double>, BicubicInterpolation> im2 = bicubic_interpolation(im);
    TS_ASSERT_EQUALS( im2.cols(), 2 );
    TS_ASSERT_EQUALS( im2.rows(), 3 );
    TS_ASSERT_EQUALS( im2(1,1), 4 );
    TS_ASSERT_EQUALS( im2(1,1.5), 5.5 );
    TS_ASSERT_EQUALS( im2(0.5,1), 3.9375);
    TS_ASSERT_DELTA( im2(0.5,0.5), 2.7773, 0.001);

    // Teste accessor
    TS_ASSERT_EQUALS( *(im2.origin().advance(1,1.5)), 5.5 );
    TS_ASSERT_EQUALS( *(im2.origin().advance(1,1)), 4 );

    // Test the traits
    TS_ASSERT( bool_trait<IsFloatingPointIndexable>( im2 ) ); 
    TS_ASSERT( bool_trait<IsFloatingPointIndexable>( bicubic_interpolation(im) ) );
    TS_ASSERT( !bool_trait<IsReferenceable>( bicubic_interpolation(im) ) );
    TS_ASSERT( !bool_trait<IsMultiplyAccessible>( bicubic_interpolation(im) ) );
    TS_ASSERT( !bool_trait<IsReferenceable>( bicubic_interpolation(im) ) );
  }

  void testNearestPixelInterpolation()
  {
    ImageView<double> im(2,3); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4; im(0,2)=5; im(1,2)=6;
    InterpolationView<ImageView<double>, NearestPixelInterpolation> im2 = nearest_pixel_interpolation(im);
    TS_ASSERT_EQUALS( im2.cols(), 2 );
    TS_ASSERT_EQUALS( im2.rows(), 3 );
    TS_ASSERT_EQUALS( im2(1,1), 4 );
    TS_ASSERT_EQUALS( im2(1,1.5), 6 );
    TS_ASSERT_EQUALS( im2(1,1.2), 4 );
    TS_ASSERT_EQUALS( im2(0.7,1), 4 );
    TS_ASSERT_EQUALS( im2(0.4,1.8), 5 );

    // Teste accessor
    TS_ASSERT_EQUALS( *(im2.origin().advance(1,1.5)), 6 );
    TS_ASSERT_EQUALS( *(im2.origin().advance(1,1)), 4 );

    // Test the traits
    TS_ASSERT( bool_trait<IsFloatingPointIndexable>( im2 ) ); 
    TS_ASSERT( bool_trait<IsFloatingPointIndexable>( nearest_pixel_interpolation(im) ) );
    TS_ASSERT( !bool_trait<IsReferenceable>( nearest_pixel_interpolation(im) ) );
    TS_ASSERT( !bool_trait<IsMultiplyAccessible>( nearest_pixel_interpolation(im) ) );
    TS_ASSERT( !bool_trait<IsReferenceable>( nearest_pixel_interpolation(im) ) );
  }


};
