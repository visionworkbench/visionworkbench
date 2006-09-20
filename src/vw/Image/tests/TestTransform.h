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

// TestTransform.h
#include <cxxtest/TestSuite.h>

#include <vw/Image/ImageView.h>
#include <vw/Image/Transform.h>

using namespace vw;

class TestTransform : public CxxTest::TestSuite
{
public:

  template <template<class> class TraitT, class T>
  static bool bool_trait( T const& arg ) {
    return TraitT<T>::value;
  }

  void testBBoxComputation()
  {
    ImageView<double> im(2,3); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4; im(0,2)=5; im(1,2)=6;
    TransformView<InterpolationView<EdgeExtendView<ImageView<double>, ZeroEdgeExtend>, BilinearInterpolation>, TranslateTransform> im2 = fixed_transform(im, TranslateTransform(1,1));

//     std::cout << compute_transformed_bbox(im, TranslateTransform(1,1)) << "\n";
//     std::cout << compute_transformed_bbox_fast(im, TranslateTransform(1,1)) << "\n";

  }


  void testTranslateTransform()
  {
    ImageView<double> im(2,3); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4; im(0,2)=5; im(1,2)=6;
    TransformView<InterpolationView<EdgeExtendView<ImageView<double>, ZeroEdgeExtend>, BilinearInterpolation>, TranslateTransform> im2 = fixed_transform(im, TranslateTransform(1,1));
    TS_ASSERT_EQUALS( im2.cols(), 2 );
    TS_ASSERT_EQUALS( im2.rows(), 3 );
    TS_ASSERT_EQUALS( im2(1,1), 1 );
    TS_ASSERT_EQUALS( im2(0,0), 0 );
    TS_ASSERT_EQUALS( im2(1,2), 3 );
    TS_ASSERT_EQUALS( im2(-1,-1), 0 );

    // Test accessor
    TS_ASSERT_EQUALS( *(im2.origin().advance(1,1)), 1 );
    TS_ASSERT_EQUALS( *(im2.origin().advance(-1,-1)), 0 );

    // Test the traits
    TS_ASSERT( bool_trait<IsFloatingPointIndexable>( im2 ) ); 
    TS_ASSERT( bool_trait<IsFloatingPointIndexable>( fixed_transform(im, TranslateTransform(1,1)) ) );
    TS_ASSERT( !bool_trait<IsReferenceable>( fixed_transform(im, TranslateTransform(1,1)) ) );
    TS_ASSERT( !bool_trait<IsMultiplyAccessible>( fixed_transform(im, TranslateTransform(1,1)) ) );
  }

};
