// __BEGIN_LICENSE__
// 
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
// 
// Copyright 2006 Carnegie Mellon University. All rights reserved.
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

// TestEdgeExtension.h
#include <cxxtest/TestSuite.h>

#include <vw/Image/ImageView.h>
#include <vw/Image/EdgeExtension.h>

#include <boost/utility/result_of.hpp>
#include <boost/type_traits.hpp>

#define TS_ASSERT_BBOX(B,X,Y,W,H)    \
  { BBox2i b = B;                    \
    TS_ASSERT_EQUALS(b.min().x(),X); \
    TS_ASSERT_EQUALS(b.min().y(),Y); \
    TS_ASSERT_EQUALS(b.width(),W);   \
    TS_ASSERT_EQUALS(b.height(),H);  \
  }

using namespace vw;

class TestEdgeExtension : public CxxTest::TestSuite
{
public:

  template <template<class> class TraitT, class T>
  static bool bool_trait( T const& arg ) {
    return TraitT<T>::value;
  }

  class SomeType {};

  void testNoEdgeExtension()
  {
    ImageView<double> im(2,3); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4; im(0,2)=5; im(1,2)=6;
    EdgeExtensionView<ImageView<double>, NoEdgeExtension> im2 = edge_extend(im, NoEdgeExtension() );
    TS_ASSERT_EQUALS( im2.cols(), 2 );
    TS_ASSERT_EQUALS( im2.rows(), 3 );
    TS_ASSERT_EQUALS( im2(1,1), 4 );

    // Test the accessor
    TS_ASSERT_EQUALS( *(im2.origin().advance(1,1)), 4 );

    // Test the traits
    TS_ASSERT( !bool_trait<IsMultiplyAccessible>( edge_extend(im, ZeroEdgeExtension() ) ) );
    TS_ASSERT( (boost::is_same<boost::result_of<NoEdgeExtension(ImageView<SomeType>,int,int,int)>::type,SomeType>::value) );

    // Test the prerasterization bbox
    NoEdgeExtension ee;
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(0,0,2,1)), 0,0,2,1 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(1,1,1,2)), 1,1,1,2 );
  }


  void testZeroEdgeExtension()
  {
    ImageView<double> im(2,3); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4; im(0,2)=5; im(1,2)=6;
    EdgeExtensionView<ImageView<double>, ZeroEdgeExtension> im2 = edge_extend(im, ZeroEdgeExtension() );
    TS_ASSERT_EQUALS( im2.cols(), 2 );
    TS_ASSERT_EQUALS( im2.rows(), 3 );

    TS_ASSERT_EQUALS( im2(-1,0), 0 );
    TS_ASSERT_EQUALS( im2(0,-1), 0 );
    TS_ASSERT_EQUALS( im2(2,0), 0 );
    TS_ASSERT_EQUALS( im2(0,3), 0 );
    TS_ASSERT_EQUALS( im2(-1,-1), 0 );
    TS_ASSERT_EQUALS( im2(2,3), 0 );

    TS_ASSERT_EQUALS( im2(5,0), 0 );
    TS_ASSERT_EQUALS( im2(0,6), 0 );
    TS_ASSERT_EQUALS( im2(-4,0), 0 );
    TS_ASSERT_EQUALS( im2(1,-4), 0 );

    // Test the accessor
    TS_ASSERT_EQUALS( *(im2.origin().advance(-1,-1)), 0 );
    TS_ASSERT_EQUALS( *(im2.origin().advance(1,1)), 4 );

    // Test the traits
    TS_ASSERT( !bool_trait<IsMultiplyAccessible>( edge_extend(im, ZeroEdgeExtension() ) ) );
    TS_ASSERT( (boost::is_same<boost::result_of<ZeroEdgeExtension(ImageView<SomeType>,int,int,int)>::type,SomeType>::value) );

    // Test the prerasterization bbox
    ZeroEdgeExtension ee;
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(0,0,2,1)), 0,0,2,1 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(1,1,1,2)), 1,1,1,2 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(-1,-1,2,2)), 0,0,1,1 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(1,-1,2,2)), 1,0,1,1 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(-1,2,2,2)), 0,2,1,1 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(1,2,2,2)), 1,2,1,1 );
    TS_ASSERT( ee.source_bbox(im,BBox2i(-2,-2,2,2)).empty() );
    TS_ASSERT( ee.source_bbox(im,BBox2i(2,3,2,2)).empty() );
  }

  void testConstantEdgeExtension()
  {
    ImageView<double> im(2,3); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4; im(0,2)=5; im(1,2)=6;
    EdgeExtensionView<ImageView<double>, ConstantEdgeExtension> im2 = edge_extend(im, ConstantEdgeExtension() );
    TS_ASSERT_EQUALS( im2.cols(), 2 );
    TS_ASSERT_EQUALS( im2.rows(), 3 );

    TS_ASSERT_EQUALS( im2(-1,0), 1 );
    TS_ASSERT_EQUALS( im2(1,-1), 2 );
    TS_ASSERT_EQUALS( im2(2,0), 2 );
    TS_ASSERT_EQUALS( im2(0,3), 5 );
    TS_ASSERT_EQUALS( im2(-1,-1), 1 );
    TS_ASSERT_EQUALS( im2(2,3), 6 );

    TS_ASSERT_EQUALS( im2(5,0), 2 );
    TS_ASSERT_EQUALS( im2(0,6), 5 );
    TS_ASSERT_EQUALS( im2(-4,0), 1 );
    TS_ASSERT_EQUALS( im2(1,-4), 2 );

    // Test the accessor
    TS_ASSERT_EQUALS( *(im2.origin().advance(-1,-1)), 1 );
    TS_ASSERT_EQUALS( *(im2.origin().advance(1,1)), 4 );

    // Test the traits
    TS_ASSERT( !bool_trait<IsMultiplyAccessible>( edge_extend(im, ZeroEdgeExtension() ) ) );
    TS_ASSERT( (boost::is_same<boost::result_of<ConstantEdgeExtension(ImageView<SomeType>,int,int,int)>::type,SomeType>::value) );

    // Test the prerasterization bbox
    ConstantEdgeExtension ee;
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(0,0,2,1)), 0,0,2,1 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(1,1,1,2)), 1,1,1,2 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(-1,-1,2,2)), 0,0,1,1 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(1,-1,2,2)), 1,0,1,1 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(-1,2,2,2)), 0,2,1,1 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(1,2,2,2)), 1,2,1,1 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(-2,-2,2,2)), 0,0,1,1 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(2,3,2,2)), 1,2,1,1 );
  }

  void testPeriodicEdgeExtension()
  {
    ImageView<double> im(2,3); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4; im(0,2)=5; im(1,2)=6;
    EdgeExtensionView<ImageView<double>, PeriodicEdgeExtension> im2 = edge_extend(im, PeriodicEdgeExtension() );
    TS_ASSERT_EQUALS( im2.cols(), 2 );
    TS_ASSERT_EQUALS( im2.rows(), 3 );

    TS_ASSERT_EQUALS( im2(-1,0), 2 );
    TS_ASSERT_EQUALS( im2(1,-1), 6 );
    TS_ASSERT_EQUALS( im2(2,0), 1 );
    TS_ASSERT_EQUALS( im2(0,3), 1 );
    TS_ASSERT_EQUALS( im2(-1,-1), 6 );
    TS_ASSERT_EQUALS( im2(2,3), 1 );

    TS_ASSERT_EQUALS( im2(3,0), 2 );
    TS_ASSERT_EQUALS( im2(4,0), 1 );
    TS_ASSERT_EQUALS( im2(5,0), 2 );

    TS_ASSERT_EQUALS( im2(0,4), 3 );
    TS_ASSERT_EQUALS( im2(0,5), 5 );
    TS_ASSERT_EQUALS( im2(0,6), 1 );

    TS_ASSERT_EQUALS( im2(-2,0), 1 );
    TS_ASSERT_EQUALS( im2(-3,0), 2 );
    TS_ASSERT_EQUALS( im2(-4,0), 1 );

    TS_ASSERT_EQUALS( im2(1,-2), 4 );
    TS_ASSERT_EQUALS( im2(1,-3), 2 );
    TS_ASSERT_EQUALS( im2(1,-4), 6 );

    // Test the accessor
    TS_ASSERT_EQUALS( *(im2.origin().advance(-1,-1)), 6 );
    TS_ASSERT_EQUALS( *(im2.origin().advance(1,1)), 4 );

    // Test the traits
    TS_ASSERT( !bool_trait<IsMultiplyAccessible>( edge_extend(im, PeriodicEdgeExtension() ) ) );
    TS_ASSERT( (boost::is_same<boost::result_of<PeriodicEdgeExtension(ImageView<SomeType>,int,int,int)>::type,SomeType>::value) );

    // Test the prerasterization bbox
    PeriodicEdgeExtension ee;
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(0,0,2,1)), 0,0,2,1 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(1,1,1,2)), 1,1,1,2 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(-1,-1,2,2)), 0,0,2,3 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(1,-1,2,2)), 0,0,2,3 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(-1,2,2,2)), 0,0,2,3 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(1,2,2,2)), 0,0,2,3 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(-2,-2,2,2)), 0,1,2,2 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(2,3,2,2)), 0,0,2,2 );
    im.set_size(4,4);
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(-5,-5,2,2)), 0,0,4,4 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(-4,-4,2,2)), 0,0,2,2 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(-3,-3,2,2)), 1,1,2,2 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(-2,-2,2,2)), 2,2,2,2 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(-1,-1,2,2)), 0,0,4,4 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(0,0,2,2)), 0,0,2,2 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(1,1,2,2)), 1,1,2,2 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(2,2,2,2)), 2,2,2,2 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(3,3,2,2)), 0,0,4,4 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(4,4,2,2)), 0,0,2,2 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(5,5,2,2)), 1,1,2,2 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(-5,-5,3,3)), 0,0,4,4 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(-4,-4,3,3)), 0,0,3,3 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(-3,-3,3,3)), 1,1,3,3 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(-2,-2,3,3)), 0,0,4,4 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(-1,-1,3,3)), 0,0,4,4 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(0,0,3,3)), 0,0,3,3 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(1,1,3,3)), 1,1,3,3 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(2,2,3,3)), 0,0,4,4 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(3,3,3,3)), 0,0,4,4 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(4,4,3,3)), 0,0,3,3 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(5,5,3,3)), 1,1,3,3 );
  }

  void testReflectEdgeExtension()
  {
    ImageView<double> im(2,3); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4; im(0,2)=5; im(1,2)=6;
    EdgeExtensionView<ImageView<double>, ReflectEdgeExtension> im2 = edge_extend(im, ReflectEdgeExtension() );
    TS_ASSERT_EQUALS( im2.cols(), 2 );
    TS_ASSERT_EQUALS( im2.rows(), 3 );

    TS_ASSERT_EQUALS( im2(-1,0), 2 );
    TS_ASSERT_EQUALS( im2(1,-1), 4 );
    TS_ASSERT_EQUALS( im2(2,0), 1 );
    TS_ASSERT_EQUALS( im2(0,3), 3 );
    TS_ASSERT_EQUALS( im2(-1,-1), 4 );
    TS_ASSERT_EQUALS( im2(2,3), 3 );

    TS_ASSERT_EQUALS( im2(3,0), 2 );
    TS_ASSERT_EQUALS( im2(4,0), 1 );
    TS_ASSERT_EQUALS( im2(5,0), 2 );

    TS_ASSERT_EQUALS( im2(0,4), 1 );
    TS_ASSERT_EQUALS( im2(0,5), 3 );
    TS_ASSERT_EQUALS( im2(0,6), 5 );

    TS_ASSERT_EQUALS( im2(-2,0), 1 );
    TS_ASSERT_EQUALS( im2(-3,0), 2 );
    TS_ASSERT_EQUALS( im2(-4,0), 1 );

    TS_ASSERT_EQUALS( im2(1,-2), 6 );
    TS_ASSERT_EQUALS( im2(1,-3), 4 );
    TS_ASSERT_EQUALS( im2(1,-4), 2 );

    // Test the accessor
    TS_ASSERT_EQUALS( *(im2.origin().advance(-1,-1)), 4 );
    TS_ASSERT_EQUALS( *(im2.origin().advance(1,1)), 4 );

    // Test the traits
    TS_ASSERT( !bool_trait<IsMultiplyAccessible>( edge_extend(im, ReflectEdgeExtension() ) ) );
    TS_ASSERT( (boost::is_same<boost::result_of<ReflectEdgeExtension(ImageView<SomeType>,int,int,int)>::type,SomeType>::value) );

    // Test the prerasterization bbox
    ReflectEdgeExtension ee;
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(0,0,2,1)), 0,0,2,1 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(1,1,1,2)), 1,1,1,2 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(-1,-1,2,2)), 0,0,2,2 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(1,-1,2,2)), 0,0,2,2 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(-1,2,2,2)), 0,1,2,2 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(1,2,2,2)), 0,1,2,2 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(-2,-2,2,2)), 0,1,2,2 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(2,3,2,2)), 0,0,2,2 );
    im.set_size(4,4);
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(-5,-5,2,2)), 1,1,2,2 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(-4,-4,2,2)), 2,2,2,2 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(-3,-3,2,2)), 2,2,2,2 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(-2,-2,2,2)), 1,1,2,2 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(-1,-1,2,2)), 0,0,2,2 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(0,0,2,2)), 0,0,2,2 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(1,1,2,2)), 1,1,2,2 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(2,2,2,2)), 2,2,2,2 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(3,3,2,2)), 2,2,2,2 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(4,4,2,2)), 1,1,2,2 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(5,5,2,2)), 0,0,2,2 );

    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(-5,-5,3,3)), 1,1,3,3 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(-4,-4,3,3)), 2,2,2,2 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(-3,-3,3,3)), 1,1,3,3 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(-2,-2,3,3)), 0,0,3,3 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(-1,-1,3,3)), 0,0,2,2 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(0,0,3,3)), 0,0,3,3 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(1,1,3,3)), 1,1,3,3 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(2,2,3,3)), 2,2,2,2 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(3,3,3,3)), 1,1,3,3 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(4,4,3,3)), 0,0,3,3 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(5,5,3,3)), 0,0,2,2 );

    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(-5,-5,4,4)), 1,1,3,3 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(-4,-4,4,4)), 1,1,3,3 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(-3,-3,4,4)), 0,0,4,4 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(-2,-2,4,4)), 0,0,3,3 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(-1,-1,4,4)), 0,0,3,3 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(0,0,4,4)), 0,0,4,4 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(1,1,4,4)), 1,1,3,3 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(2,2,4,4)), 1,1,3,3 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(3,3,4,4)), 0,0,4,4 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(4,4,4,4)), 0,0,3,3 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(5,5,4,4)), 0,0,3,3 );
  }

  void testLinearEdgeExtension()
  {
    ImageView<double> im(2,3); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4; im(0,2)=5; im(1,2)=6;
    EdgeExtensionView<ImageView<double>, LinearEdgeExtension> im2 = edge_extend(im, LinearEdgeExtension() );
    TS_ASSERT_EQUALS( im2.cols(), 2 );
    TS_ASSERT_EQUALS( im2.rows(), 3 );

    TS_ASSERT_EQUALS( im2(-1,0), 0 );
    TS_ASSERT_EQUALS( im2(1,-1), 0 );
    TS_ASSERT_EQUALS( im2(2,0), 3 );
    TS_ASSERT_EQUALS( im2(0,3), 7 );
    TS_ASSERT_EQUALS( im2(-1,-1), -2 );
    TS_ASSERT_EQUALS( im2(2,3), 9 );

    TS_ASSERT_EQUALS( im2(3,0), 4 );
    TS_ASSERT_EQUALS( im2(4,0), 5 );
    TS_ASSERT_EQUALS( im2(5,0), 6 );

    TS_ASSERT_EQUALS( im2(0,4), 9 );
    TS_ASSERT_EQUALS( im2(0,5), 11 );
    TS_ASSERT_EQUALS( im2(0,6), 13 );

    TS_ASSERT_EQUALS( im2(-2,0), -1 );
    TS_ASSERT_EQUALS( im2(-3,0), -2 );
    TS_ASSERT_EQUALS( im2(-4,0), -3 );

    TS_ASSERT_EQUALS( im2(1,-2), -2 );
    TS_ASSERT_EQUALS( im2(1,-3), -4 );
    TS_ASSERT_EQUALS( im2(1,-4), -6 );

    // Test the accessor
    TS_ASSERT_EQUALS( *(im2.origin().advance(-1,-1)), -2 );
    TS_ASSERT_EQUALS( *(im2.origin().advance(1,1)), 4 );

    // Test the traits
    TS_ASSERT( !bool_trait<IsMultiplyAccessible>( edge_extend(im, LinearEdgeExtension() ) ) );
    TS_ASSERT( (boost::is_same<boost::result_of<LinearEdgeExtension(ImageView<SomeType>,int,int,int)>::type,SomeType>::value) );

    // Test the prerasterization bbox
    LinearEdgeExtension ee;
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(0,0,2,1)), 0,0,2,1 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(1,1,1,2)), 1,1,1,2 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(-1,-1,2,2)), 0,0,2,2 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(1,-1,2,2)), 0,0,2,2 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(-1,2,2,2)), 0,1,2,2 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(1,2,2,2)), 0,1,2,2 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(-2,-2,2,2)), 0,0,2,2 );
    TS_ASSERT_BBOX( ee.source_bbox(im,BBox2i(2,3,2,2)), 0,1,2,2 );
  }

};
