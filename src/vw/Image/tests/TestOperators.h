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

// TestOperators.h
#include <cxxtest/TestSuite.h>

#include <vw/Image/Operators.h>
#include <vw/Image/ImageView.h>

using namespace vw;

class TestOperators : public CxxTest::TestSuite
{
public:

  template <template<class> class TraitT, class T>
  static bool bool_trait( T const& arg ) {
    return TraitT<T>::value;
  }

  void testNegation()
  {
    ImageView<double> im1(2,2); im1(0,0)=1; im1(1,0)=2; im1(0,1)=3; im1(1,1)=4;

    ImageView<double> im2 = -im1;
    TS_ASSERT_EQUALS( im2(0,0), -1 );
    TS_ASSERT_EQUALS( im2(1,0), -2 );
    TS_ASSERT_EQUALS( im2(0,1), -3 );
    TS_ASSERT_EQUALS( im2(1,1), -4 );

    // Test the traits
    TS_ASSERT( !bool_trait<IsReferenceable>( -im1 ) );
    TS_ASSERT( !bool_trait<IsMultiplyAccessible>( -im1 ) );
    TS_ASSERT( !bool_trait<IsReferenceable>( (-im1).origin() ) );
  }

  void testSum()
  {
    ImageView<double> im1(2,2); im1(0,0)=1; im1(1,0)=2; im1(0,1)=3; im1(1,1)=4;
    ImageView<double> im2(2,2); im2(0,0)=0; im2(1,0)=3; im2(0,1)=2; im2(1,1)=9;

    ImageView<double> im3 = im1 + im2;
    TS_ASSERT_EQUALS( im3(0,0), 1 );
    TS_ASSERT_EQUALS( im3(1,0), 5 );
    TS_ASSERT_EQUALS( im3(0,1), 5 );
    TS_ASSERT_EQUALS( im3(1,1), 13 );

    im3 = im1 + 2;
    TS_ASSERT_EQUALS( im3(0,0), 3 );
    TS_ASSERT_EQUALS( im3(1,0), 4 );
    TS_ASSERT_EQUALS( im3(0,1), 5 );
    TS_ASSERT_EQUALS( im3(1,1), 6 );

    im3 = 1 + im1;
    TS_ASSERT_EQUALS( im3(0,0), 2 );
    TS_ASSERT_EQUALS( im3(1,0), 3 );
    TS_ASSERT_EQUALS( im3(0,1), 4 );
    TS_ASSERT_EQUALS( im3(1,1), 5 );

    im3 += im1;
    TS_ASSERT_EQUALS( im3(0,0), 3 );
    TS_ASSERT_EQUALS( im3(1,0), 5 );
    TS_ASSERT_EQUALS( im3(0,1), 7 );
    TS_ASSERT_EQUALS( im3(1,1), 9 );

    im3 += 2;
    TS_ASSERT_EQUALS( im3(0,0), 5 );
    TS_ASSERT_EQUALS( im3(1,0), 7 );
    TS_ASSERT_EQUALS( im3(0,1), 9 );
    TS_ASSERT_EQUALS( im3(1,1), 11 );

    // Test the traits
    TS_ASSERT( !bool_trait<IsReferenceable>( im1 + im2 ) );
    TS_ASSERT( !bool_trait<IsMultiplyAccessible>( im1 + im2 ) );
    TS_ASSERT( !bool_trait<IsReferenceable>( (im1 + im2).origin() ) );
  }

  void testDifference()
  {
    ImageView<double> im1(2,2); im1(0,0)=1; im1(1,0)=2; im1(0,1)=3; im1(1,1)=4;
    ImageView<double> im2(2,2); im2(0,0)=0; im2(1,0)=3; im2(0,1)=2; im2(1,1)=9;

    ImageView<double> im3 = im1 - im2;
    TS_ASSERT_EQUALS( im3(0,0), 1 );
    TS_ASSERT_EQUALS( im3(1,0), -1 );
    TS_ASSERT_EQUALS( im3(0,1), 1 );
    TS_ASSERT_EQUALS( im3(1,1), -5 );

    im3 = im1 - 2;
    TS_ASSERT_EQUALS( im3(0,0), -1 );
    TS_ASSERT_EQUALS( im3(1,0), 0 );
    TS_ASSERT_EQUALS( im3(0,1), 1 );
    TS_ASSERT_EQUALS( im3(1,1), 2 );

    im3 = 1 - im1;
    TS_ASSERT_EQUALS( im3(0,0), 0 );
    TS_ASSERT_EQUALS( im3(1,0), -1 );
    TS_ASSERT_EQUALS( im3(0,1), -2 );
    TS_ASSERT_EQUALS( im3(1,1), -3 );

    im3 -= im1;
    TS_ASSERT_EQUALS( im3(0,0), -1 );
    TS_ASSERT_EQUALS( im3(1,0), -3 );
    TS_ASSERT_EQUALS( im3(0,1), -5 );
    TS_ASSERT_EQUALS( im3(1,1), -7 );

    im3 -= 2;
    TS_ASSERT_EQUALS( im3(0,0), -3 );
    TS_ASSERT_EQUALS( im3(1,0), -5 );
    TS_ASSERT_EQUALS( im3(0,1), -7 );
    TS_ASSERT_EQUALS( im3(1,1), -9 );

    // Test the traits
    TS_ASSERT( !bool_trait<IsReferenceable>( im1 - im2 ) );
    TS_ASSERT( !bool_trait<IsMultiplyAccessible>( im1 - im2 ) );
    TS_ASSERT( !bool_trait<IsReferenceable>( (im1 - im2).origin() ) );
  }

  void testProduct()
  {
    ImageView<double> im1(2,2); im1(0,0)=1; im1(1,0)=2; im1(0,1)=3; im1(1,1)=4;
    ImageView<double> im2(2,2); im2(0,0)=0; im2(1,0)=3; im2(0,1)=2; im2(1,1)=9;

    ImageView<double> im3 = im1 * im2;
    TS_ASSERT_EQUALS( im3(0,0), 0 );
    TS_ASSERT_EQUALS( im3(1,0), 6 );
    TS_ASSERT_EQUALS( im3(0,1), 6 );
    TS_ASSERT_EQUALS( im3(1,1), 36 );

    im3 = im1 * 2;
    TS_ASSERT_EQUALS( im3(0,0), 2 );
    TS_ASSERT_EQUALS( im3(1,0), 4 );
    TS_ASSERT_EQUALS( im3(0,1), 6 );
    TS_ASSERT_EQUALS( im3(1,1), 8 );

    im3 = 3 * im1;
    TS_ASSERT_EQUALS( im3(0,0), 3 );
    TS_ASSERT_EQUALS( im3(1,0), 6 );
    TS_ASSERT_EQUALS( im3(0,1), 9 );
    TS_ASSERT_EQUALS( im3(1,1), 12 );

    im3 *= im1;
    TS_ASSERT_EQUALS( im3(0,0), 3 );
    TS_ASSERT_EQUALS( im3(1,0), 12 );
    TS_ASSERT_EQUALS( im3(0,1), 27 );
    TS_ASSERT_EQUALS( im3(1,1), 48 );

    im3 *= 0.5;
    TS_ASSERT_EQUALS( im3(0,0), 1.5 );
    TS_ASSERT_EQUALS( im3(1,0), 6 );
    TS_ASSERT_EQUALS( im3(0,1), 13.5 );
    TS_ASSERT_EQUALS( im3(1,1), 24 );

    // Test the traits
    TS_ASSERT( !bool_trait<IsReferenceable>( im1 * im2 ) );
    TS_ASSERT( !bool_trait<IsMultiplyAccessible>( im1 * im2 ) );
    TS_ASSERT( !bool_trait<IsReferenceable>( (im1 * im2).origin() ) );
  }

  void testQuotient()
  {
    ImageView<double> im1(2,2); im1(0,0)=1; im1(1,0)=2; im1(0,1)=3; im1(1,1)=4;
    ImageView<double> im2(2,2); im2(0,0)=1; im2(1,0)=4; im2(0,1)=2; im2(1,1)=8;

    ImageView<double> im3 = im1 / im2;
    TS_ASSERT_EQUALS( im3(0,0), 1 );
    TS_ASSERT_EQUALS( im3(1,0), 0.5 );
    TS_ASSERT_EQUALS( im3(0,1), 1.5 );
    TS_ASSERT_EQUALS( im3(1,1), 0.5 );

    im3 = im1 / 2;
    TS_ASSERT_EQUALS( im3(0,0), 0.5 );
    TS_ASSERT_EQUALS( im3(1,0), 1 );
    TS_ASSERT_EQUALS( im3(0,1), 1.5 );
    TS_ASSERT_EQUALS( im3(1,1), 2 );

    im3 = 2 / im2;
    TS_ASSERT_EQUALS( im3(0,0), 2 );
    TS_ASSERT_EQUALS( im3(1,0), 0.5 );
    TS_ASSERT_EQUALS( im3(0,1), 1 );
    TS_ASSERT_EQUALS( im3(1,1), 0.25 );

    im3 /= im2;
    TS_ASSERT_EQUALS( im3(0,0), 2 );
    TS_ASSERT_EQUALS( im3(1,0), 0.125 );
    TS_ASSERT_EQUALS( im3(0,1), 0.5 );
    TS_ASSERT_EQUALS( im3(1,1), 0.03125 );

    im3 /= 0.5;
    TS_ASSERT_EQUALS( im3(0,0), 4 );
    TS_ASSERT_EQUALS( im3(1,0), 0.25 );
    TS_ASSERT_EQUALS( im3(0,1), 1 );
    TS_ASSERT_EQUALS( im3(1,1), 0.0625 );

    // Test the traits
    TS_ASSERT( !bool_trait<IsReferenceable>( im1 / im2 ) );
    TS_ASSERT( !bool_trait<IsMultiplyAccessible>( im1 / im2 ) );
    TS_ASSERT( !bool_trait<IsReferenceable>( (im1 / im2).origin() ) );
  }

};
