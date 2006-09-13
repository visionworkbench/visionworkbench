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
#include <vw/Image/PixelTypes.h>

using namespace vw;

class TestPerPixelViews : public CxxTest::TestSuite
{
public:

  template <template<class> class TraitT, class T>
  static bool bool_trait( T const& arg ) {
    return TraitT<T>::value;
  }

  void testSelectChannel()
  {
    ImageView<PixelRGB<double> > im(1,2); im(0,0)=PixelRGB<double>(1,2,3); im(0,1)=PixelRGB<double>(4,5,6);
    ImageView<double> im2 = select_channel(im,1);
    TS_ASSERT_EQUALS( im2.cols(), 1 );
    TS_ASSERT_EQUALS( im2.rows(), 2 );
    TS_ASSERT_EQUALS( im2.planes(), 1 );
    TS_ASSERT_EQUALS( im2(0,0), 2 );
    TS_ASSERT_EQUALS( im2(0,1), 5 );

    // Make sure it's really shallow.
    TS_ASSERT_EQUALS( select_channel(im,1)(0,1), im(0,1)[1] );
    TS_ASSERT_EQUALS( &(select_channel(im,1)(0,1)), &(im(0,1)[1]) );

    // Test the traits
    TS_ASSERT( bool_trait<IsReferenceable>( select_channel(im,1) ) );
    TS_ASSERT( bool_trait<IsMultiplyAccessible>( select_channel(im,1) ) );
    TS_ASSERT( bool_trait<IsReferenceable>( select_channel(im,1).origin() ) );
  }

  void testPixelCast()
  {
    ImageView<PixelRGB<double> > im(1,2); im(0,0)=PixelRGB<double>(1,2,3); im(0,1)=PixelRGB<double>(4,5,6);
    ImageView<PixelRGBA<double> > im2 = pixel_cast<PixelRGBA<double> >(im);
    TS_ASSERT_EQUALS( im2.cols(), 1 );
    TS_ASSERT_EQUALS( im2.rows(), 2 );
    TS_ASSERT_EQUALS( im2.planes(), 1 );
    TS_ASSERT_EQUALS( im2(0,0)[0], im(0,0)[0] );
    TS_ASSERT_EQUALS( im2(0,0)[1], im(0,0)[1] );
    TS_ASSERT_EQUALS( im2(0,0)[2], im(0,0)[2] );
    TS_ASSERT_EQUALS( im2(0,1)[0], im(0,1)[0] );
    TS_ASSERT_EQUALS( im2(0,1)[1], im(0,1)[1] );
    TS_ASSERT_EQUALS( im2(0,1)[2], im(0,1)[2] );

    // Test the traits
    TS_ASSERT( !bool_trait<IsReferenceable>( pixel_cast<PixelRGBA<double> >(im) ) );
    TS_ASSERT( !bool_trait<IsMultiplyAccessible>( pixel_cast<PixelRGBA<double> >(im) ) );
    TS_ASSERT( !bool_trait<IsReferenceable>( pixel_cast<PixelRGBA<double> >(im).origin() ) );
  }

  void testChannelCast()
  {
    ImageView<PixelRGB<double> > im(1,2); im(0,0)=PixelRGB<double>(1,2,3); im(0,1)=PixelRGB<double>(4,5,6);
    ImageView<PixelRGB<float> > im2 = channel_cast<float>(im);
    TS_ASSERT_EQUALS( im2.cols(), 1 );
    TS_ASSERT_EQUALS( im2.rows(), 2 );
    TS_ASSERT_EQUALS( im2.planes(), 1 );
    TS_ASSERT_EQUALS( im2(0,0)[0], im(0,0)[0] );
    TS_ASSERT_EQUALS( im2(0,0)[1], im(0,0)[1] );
    TS_ASSERT_EQUALS( im2(0,0)[2], im(0,0)[2] );
    TS_ASSERT_EQUALS( im2(0,1)[0], im(0,1)[0] );
    TS_ASSERT_EQUALS( im2(0,1)[1], im(0,1)[1] );
    TS_ASSERT_EQUALS( im2(0,1)[2], im(0,1)[2] );

    // Test the traits
    TS_ASSERT( !bool_trait<IsReferenceable>( channel_cast<float>(im) ) );
    TS_ASSERT( !bool_trait<IsMultiplyAccessible>( channel_cast<float>(im) ) );
    TS_ASSERT( !bool_trait<IsReferenceable>( channel_cast<float>(im).origin() ) );
  }

};
