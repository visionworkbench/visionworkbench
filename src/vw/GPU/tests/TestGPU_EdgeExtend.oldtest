// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// TestEdgeExtend.h
#include <cxxtest/TestSuite.h>

#include <vw/GPU.h>

using namespace vw;
using namespace GPU;

class TestEdgeExtend : public CxxTest::TestSuite
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
  }

};
