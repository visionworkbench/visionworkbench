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

// TestPerPixelView.h
#include <cxxtest/TestSuite.h>

#include <vw/Image/PerPixelViews.h>
#include <vw/Image/ImageView.h>
#include <vw/Core/Functors.h>

using namespace vw;

class TestPerPixelViews : public CxxTest::TestSuite
{
public:

  template <template<class> class TraitT, class T>
  static bool bool_trait( T const& arg ) {
    return TraitT<T>::value;
  }

  template <class T1, class T2>
  static bool has_pixel_type( T2 ) {
    return boost::is_same<T1,typename T2::pixel_type>::value;
  }

  template <class T>
  static T square(T arg) {
    return arg*arg;
  }

  void test_per_pixel_view_function() {
    ImageView<float> im(2,2); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4;
    UnaryPerPixelView<ImageView<float>, float(*)(float)> ppv( im, square );
    TS_ASSERT_EQUALS( ppv.cols(), 2 );
    TS_ASSERT_EQUALS( ppv.rows(), 2 );
    TS_ASSERT_EQUALS( ppv.planes(), 1 );
    TS_ASSERT_EQUALS( ppv(0,0), 1 );
    TS_ASSERT_EQUALS( ppv(1,0), 4 );
    TS_ASSERT_EQUALS( ppv(0,1), 9 );
    TS_ASSERT_EQUALS( ppv(1,1), 16 );
    TS_ASSERT( has_pixel_type<float>( ppv ) );
    TS_ASSERT( !bool_trait<IsMultiplyAccessible>(ppv) );
    TS_ASSERT( bool_trait<IsImageView>(ppv) );
  }

  void test_per_pixel_view_functor() {
    ImageView<float> im(2,2); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4;
    UnaryPerPixelView<ImageView<float>, ArgNegationFunctor> ppv( im );
    TS_ASSERT_EQUALS( ppv.cols(), 2 );
    TS_ASSERT_EQUALS( ppv.rows(), 2 );
    TS_ASSERT_EQUALS( ppv.planes(), 1 );
    TS_ASSERT_EQUALS( ppv(0,0), -1 );
    TS_ASSERT_EQUALS( ppv(1,0), -2 );
    TS_ASSERT_EQUALS( ppv(0,1), -3 );
    TS_ASSERT_EQUALS( ppv(1,1), -4 );
    TS_ASSERT( has_pixel_type<float>( ppv ) );
    TS_ASSERT( !bool_trait<IsMultiplyAccessible>(ppv) );
    TS_ASSERT( bool_trait<IsImageView>(ppv) );
  }

};
