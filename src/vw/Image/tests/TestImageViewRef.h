// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// TestImageViewRef.h

#include <cxxtest/TestSuite.h>

#include <vw/Image/ImageView.h>
#include <vw/Image/ImageViewRef.h>

using namespace vw;

class TestImageViewRef : public CxxTest::TestSuite
{
public:

  void testImageViewRef()
  {
    const int cols=3, rows=2;
    ImageView<float> image(cols,rows);
    for( int r=0; r<rows; ++r )
      for( int c=0; c<cols; ++c )
        image(c,r) = (float)(r*cols+c);

    ImageViewRef<float> ref = image;
    TS_ASSERT_EQUALS( ref.cols(), image.cols() );
    TS_ASSERT_EQUALS( ref.rows(), image.rows() );
    TS_ASSERT_EQUALS( ref.planes(), image.planes() );

    // Test pixel indexing
    for( int r=0; r<rows; ++r ) {
      for( int c=0; c<cols; ++c ) {
        TS_ASSERT_EQUALS( ref(c,r), (float)(r*cols+c) );
      }
    }

    // Test full rasterization: optimized case
    ImageView<float> im2 = ref;
    for( int r=0; r<rows; ++r ) {
      for( int c=0; c<cols; ++c ) {
        TS_ASSERT_EQUALS( im2(c,r), ref(c,r) );
      }
    }

    // Test full rasterization: general case
    ImageView<double> im3 = ref;
    for( int r=0; r<rows; ++r ) {
      for( int c=0; c<cols; ++c ) {
        TS_ASSERT_DELTA( im3(c,r), ref(c,r), 1e-8 );
      }
    }

    // Test accessor / generic rasterization
    ImageView<float> im4(cols,rows);
    vw::rasterize( ref, im4, BBox2i(0,0,cols,rows) );
    for( int r=0; r<rows; ++r ) {
      for( int c=0; c<cols; ++c ) {
        TS_ASSERT_EQUALS( im4(c,r), ref(c,r) );
      }
    }

    // Test partial rasterization
    ImageView<float> im5(cols-1,rows-1);
    ref.rasterize( im5, BBox2i(1,1,cols-1,rows-1) );
    for( int r=0; r<rows-1; ++r ) {
      for( int c=0; c<cols-1; ++c ) {
        TS_ASSERT_EQUALS( im5(c,r), ref(c+1,r+1) );
      }
    }

    // Test iterator
    int val=0;
    for( ImageViewRef<float>::iterator i=ref.begin(), end=ref.end(); i!=end; ++i, ++val ) {
      TS_ASSERT_EQUALS( *i, (float)(val) );
    }
  }

};
