// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// TestManipulation.h
#include <gtest/gtest.h>

#include <vw/Image/ImageView.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Image/Manipulation.h>
#include <vw/Image/Filter.h>
#include <vw/Image/Algorithms.h>

using namespace vw;

class PrerasterizationTestView : public ImageViewBase<PrerasterizationTestView> {
  ImageView<uint8> image;
public:
  typedef ImageView<uint8>::pixel_type pixel_type;
  typedef ImageView<uint8>::result_type result_type;
  typedef ImageView<uint8>::pixel_accessor pixel_accessor;

  PrerasterizationTestView( int32 cols, int32 rows ) : image(cols,rows) {}
  int32 cols() const { return image.cols(); }
  int32 rows() const { return image.rows(); }
  int32 planes() const { return 1; }
  pixel_accessor origin() const { return image.origin(); }

  BBox2i bbox() const {
    BBox2i result;
    for( int32 y=0; y<image.rows(); ++y ) {
      for( int32 x=0; x<image.cols(); ++x ) {
        if( image(x,y) ) result.grow( Vector2(x,y) );
      }
    }
    if( result != BBox2i() ) {
      result.max() += Vector2i(1,1);
    }
    return result;
  }

  typedef PrerasterizationTestView prerasterize_type;
  inline prerasterize_type prerasterize( BBox2i const& bbox ) const {
    if( bbox.min().x() < 0 || bbox.min().y() < 0 || bbox.max().x() > image.cols() || bbox.max().y() > image.rows() ) {
      vw_throw(ArgumentErr() << "PrerasterizationTestView::prerasterize() called with bounding box that exceeds image dimensions");
    }
    for( int32 y=bbox.min().y(); y<bbox.max().y(); ++y ) {
      for( int32 x=bbox.min().x(); x<bbox.max().x(); ++x ) {
        image(x,y) = 1;
      }
    }
    return *this;
  }
  template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const {
    vw::rasterize( prerasterize(bbox), dest, bbox );
  }
};


template <template<class> class TraitT, class T>
static bool bool_trait( T const& /*arg*/ ) {
  return TraitT<T>::value;
}

template <class T1, class T2>
static bool has_pixel_type( T2 ) {
  return boost::is_same<T1,typename T2::pixel_type>::value;
}

template <class T1, class T2>
static bool has_result_type( T2 ) {
  return boost::is_same<T1,typename T2::result_type>::value;
}

template <class T>
static T const& make_const_ref( T const& arg ) {
  return arg;
}

TEST( Manipulation, Copy ) {
  ImageView<double> im(2,3); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4; im(0,2)=5; im(1,2)=6;
  ImageView<double> im2 = copy(im);
  ASSERT_EQ( im2.cols(), 2 );
  ASSERT_EQ( im2.rows(), 3 );
  EXPECT_EQ( im2(0,0), 1 );
  EXPECT_EQ( im2(1,0), 2 );
  EXPECT_EQ( im2(0,1), 3 );
  EXPECT_EQ( im2(1,1), 4 );
  EXPECT_EQ( im2(0,2), 5 );
  EXPECT_EQ( im2(1,2), 6 );

  // Make sure it's really deep.
  EXPECT_NE( im2.data(), im.data() );
  EXPECT_EQ( copy(im)(1,0), im(1,0) );
  EXPECT_NE( &(copy(im)(1,0)), &(im(1,0)) );

  // Test the traits
  ASSERT_TRUE( bool_trait<IsMultiplyAccessible>( copy(im) ) );
}

TEST( Manipulation, TransposeView ) {
  ImageView<double> im(2,3); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4; im(0,2)=5; im(1,2)=6;
  TransposeView<ImageView<double> > rmv(im);
  ASSERT_EQ( rmv.cols(), 3 );
  ASSERT_EQ( rmv.rows(), 2 );
  ASSERT_EQ( rmv.planes(), 1 );

  // Test individual pixel access
  EXPECT_EQ( rmv(0,0), 1 );
  EXPECT_EQ( rmv(1,0), 3 );
  EXPECT_EQ( rmv(2,0), 5 );
  EXPECT_EQ( rmv(0,1), 2 );
  EXPECT_EQ( rmv(1,1), 4 );
  EXPECT_EQ( rmv(2,1), 6 );

  // Test full rasterizaion
  ImageView<double> im2 = rmv;
  ASSERT_EQ( im2.cols(), rmv.cols() );
  ASSERT_EQ( im2.rows(), rmv.rows() );
  ASSERT_EQ( im2.planes(), rmv.planes() );
  for ( int r=0; r<im2.rows(); ++r )
    for ( int c=0; c<im2.cols(); ++c )
      EXPECT_EQ( im2(c,r), rmv(c,r) );

  // Test partial rasterization
  ImageView<double> im3(rmv.cols()-1,rmv.rows()-1);
  ASSERT_NO_THROW( rmv.rasterize( im3, BBox2i(1,1,rmv.cols()-1,rmv.rows()-1) ) );
  for ( int r=0; r<im3.rows(); ++r )
    for ( int c=0; c<im3.cols(); ++c )
      EXPECT_EQ( im3(c,r), rmv(c+1,r+1) );

  // Test prerasterization
  PrerasterizationTestView ptv(4,4);
  TransposeView<PrerasterizationTestView> rtv(ptv);
  ASSERT_NO_THROW( rtv.prerasterize(BBox2i(2,0,1,2)) );
  EXPECT_EQ( ptv.bbox(), BBox2i(0,2,2,1) );

  // Test the accessor / generic rasterization
  ImageView<double> im4(rmv.cols(),rmv.rows());
  vw::rasterize( rmv, im4, BBox2i(0,0,rmv.cols(),rmv.rows()) );
  for ( int r=0; r<im2.rows(); ++r )
    for ( int c=0; c<im3.cols(); ++c )
      EXPECT_EQ( im4(c,r), rmv(c,r) );

  // Test the iterator
  ImageView<double>::iterator im2i = im2.begin();
  TransposeView<ImageView<double> >::iterator rmvi = rmv.begin();
  for ( int i=0; i<im2.cols()*im2.rows(); ++i ) {
    EXPECT_GT( rmv.end() - rmvi, 0 );
    EXPECT_EQ( *im2i, *rmvi );
    ASSERT_NO_THROW( ++rmvi );
    ++im2i;
  }
  EXPECT_EQ( rmv.end() - rmvi, 0 );

  // Test the types
  ASSERT_TRUE( has_pixel_type<double>( rmv ) );
  ASSERT_TRUE( bool_trait<IsMultiplyAccessible>(rmv) );
  ASSERT_TRUE( bool_trait<IsImageView>(rmv) );
}

TEST( Manipulation, Transpose ) {
  ImageView<double> im(2,3); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4; im(0,2)=5; im(1,2)=6;
  ImageView<double> im2 = transpose(im);
  EXPECT_EQ( im2.cols(), 3 );
  EXPECT_EQ( im2.rows(), 2 );

  // Test pixel indexing
  EXPECT_EQ( im2(0,0), 1 );
  EXPECT_EQ( im2(1,0), 3 );
  EXPECT_EQ( im2(2,0), 5 );
  EXPECT_EQ( im2(0,1), 2 );
  EXPECT_EQ( im2(1,1), 4 );
  EXPECT_EQ( im2(2,1), 6 );

  // Make sure it's really shallow.
  EXPECT_EQ( transpose(im)(1,0), im(0,1) );
  EXPECT_EQ( &(transpose(im)(1,0)), &(im(0,1)) );

  // Test the traits
  ASSERT_TRUE( bool_trait<IsMultiplyAccessible>( transpose(im) ) );
}

TEST( Manipulation, Rotate180View ) {
  ImageView<double> im(2,3); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4; im(0,2)=5; im(1,2)=6;
  Rotate180View<ImageView<double> > rmv(im);
  ASSERT_EQ( rmv.cols(), 2 );
  ASSERT_EQ( rmv.rows(), 3 );
  ASSERT_EQ( rmv.planes(), 1 );

  // Test individual pixel access
  EXPECT_EQ( rmv(0,0), 6 );
  EXPECT_EQ( rmv(1,0), 5 );
  EXPECT_EQ( rmv(0,1), 4 );
  EXPECT_EQ( rmv(1,1), 3 );
  EXPECT_EQ( rmv(0,2), 2 );
  EXPECT_EQ( rmv(1,2), 1 );

  // Test full rasterizaion
  ImageView<double> im2 = rmv;
  ASSERT_EQ( im2.cols(), rmv.cols() );
  ASSERT_EQ( im2.rows(), rmv.rows() );
  ASSERT_EQ( im2.planes(), rmv.planes() );
  for ( int r=0; r<im2.rows(); ++r )
    for ( int c=0; c<im2.cols(); ++c )
      EXPECT_EQ( im2(c,r), rmv(c,r) );

  // Test partial rasterization
  ImageView<double> im3(rmv.cols()-1,rmv.rows()-1);
  ASSERT_NO_THROW( rmv.rasterize( im3, BBox2i(1,1,rmv.cols()-1,rmv.rows()-1) ) );
  for ( int r=0; r<im3.rows(); ++r )
    for ( int c=0; c<im3.cols(); ++c )
      EXPECT_EQ( im3(c,r), rmv(c+1,r+1) );

  // Test prerasterization
  PrerasterizationTestView ptv(4,4);
  Rotate180View<PrerasterizationTestView> rtv(ptv);
  ASSERT_NO_THROW( rtv.prerasterize(BBox2i(2,0,1,2)) );
  EXPECT_EQ( ptv.bbox(), BBox2i(1,2,1,2) );

  // Test the accessor / generic rasterization
  ImageView<double> im4(rmv.cols(),rmv.rows());
  vw::rasterize( rmv, im4, BBox2i(0,0,rmv.cols(),rmv.rows()) );
  for ( int r=0; r<im2.rows(); ++r )
    for ( int c=0; c<im3.cols(); ++c )
      EXPECT_EQ( im4(c,r), rmv(c,r) );

  // Test the iterator
  ImageView<double>::iterator im2i = im2.begin();
  Rotate180View<ImageView<double> >::iterator rmvi = rmv.begin();
  for ( int i=0; i<im2.cols()*im2.rows(); ++i ) {
    EXPECT_GT( rmv.end() - rmvi, 0 );
    EXPECT_EQ( *im2i, *rmvi );
    ASSERT_NO_THROW( ++rmvi );
    ++im2i;
  }
  EXPECT_EQ( rmv.end() - rmvi, 0 );

  // Test the types
  ASSERT_TRUE( has_pixel_type<double>( rmv ) );
  ASSERT_TRUE( bool_trait<IsMultiplyAccessible>(rmv) );
  ASSERT_TRUE( bool_trait<IsImageView>(rmv) );
}

TEST( Manipulation, Rotate180 ) {
  ImageView<double> im(2,3); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4; im(0,2)=5; im(1,2)=6;
  ImageView<double> im2 = rotate_180(im);
  ASSERT_EQ( im2.cols(), 2 );
  ASSERT_EQ( im2.rows(), 3 );
  EXPECT_EQ( im2(0,0), 6 );
  EXPECT_EQ( im2(1,0), 5 );
  EXPECT_EQ( im2(0,1), 4 );
  EXPECT_EQ( im2(1,1), 3 );
  EXPECT_EQ( im2(0,2), 2 );
  EXPECT_EQ( im2(1,2), 1 );

  // Make sure it's really shallow.
  EXPECT_EQ( rotate_180(im)(1,0), im(0,2) );
  EXPECT_EQ( &(rotate_180(im)(1,0)), &(im(0,2)) );

  // Test the traits
  ASSERT_TRUE( bool_trait<IsMultiplyAccessible>( rotate_180(im) ) );
}

TEST( Manipulation, 90CWView ) {
  ImageView<double> im(2,3); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4; im(0,2)=5; im(1,2)=6;
  Rotate90CWView<ImageView<double> > rmv(im);
  ASSERT_EQ( rmv.cols(), 3 );
  ASSERT_EQ( rmv.rows(), 2 );
  ASSERT_EQ( rmv.planes(), 1 );

  // Test individual pixel access
  EXPECT_EQ( rmv(0,0), 5 );
  EXPECT_EQ( rmv(1,0), 3 );
  EXPECT_EQ( rmv(2,0), 1 );
  EXPECT_EQ( rmv(0,1), 6 );
  EXPECT_EQ( rmv(1,1), 4 );
  EXPECT_EQ( rmv(2,1), 2 );

  // Test full rasterizaion
  ImageView<double> im2 = rmv;
  EXPECT_EQ( im2.cols(), rmv.cols() );
  EXPECT_EQ( im2.rows(), rmv.rows() );
  EXPECT_EQ( im2.planes(), rmv.planes() );
  for ( int r=0; r<im2.rows(); ++r )
    for ( int c=0; c<im2.cols(); ++c )
      EXPECT_EQ( im2(c,r), rmv(c,r) );

  // Test partial rasterization
  ImageView<double> im3(rmv.cols()-1,rmv.rows()-1);
  ASSERT_NO_THROW( rmv.rasterize( im3, BBox2i(1,1,rmv.cols()-1,rmv.rows()-1) ) );
  for ( int r=0; r<im3.rows(); ++r )
    for ( int c=0; c<im3.cols(); ++c )
      EXPECT_EQ( im3(c,r), rmv(c+1,r+1) );

  // Test prerasterization
  PrerasterizationTestView ptv(4,4);
  Rotate90CWView<PrerasterizationTestView> rtv(ptv);
  ASSERT_NO_THROW( rtv.prerasterize(BBox2i(2,0,1,2)) );
  EXPECT_EQ( ptv.bbox(), BBox2i(0,1,2,1) );

  // Test the accessor / generic rasterization
  ImageView<double> im4(rmv.cols(),rmv.rows());
  vw::rasterize( rmv, im4, BBox2i(0,0,rmv.cols(),rmv.rows()) );
  for ( int r=0; r<im2.rows(); ++r )
    for ( int c=0; c<im3.cols(); ++c )
      EXPECT_EQ( im4(c,r), rmv(c,r) );

  // Test the iterator
  ImageView<double>::iterator im2i = im2.begin();
  Rotate90CWView<ImageView<double> >::iterator rmvi = rmv.begin();
  for ( int i=0; i<im2.cols()*im2.rows(); ++i ) {
    EXPECT_GT( rmv.end() - rmvi, 0 );
    EXPECT_EQ( *im2i, *rmvi );
    ASSERT_NO_THROW( ++rmvi );
    ++im2i;
  }
  EXPECT_EQ( rmv.end() - rmvi, 0 );

  // Test the types
  ASSERT_TRUE( has_pixel_type<double>( rmv ) );
  ASSERT_TRUE( bool_trait<IsMultiplyAccessible>(rmv) );
  ASSERT_TRUE( bool_trait<IsImageView>(rmv) );
}

TEST( Manipulation, Rotate90CW ) {
  ImageView<double> im(2,3); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4; im(0,2)=5; im(1,2)=6;
  ImageView<double> im2 = rotate_90_cw(im);
  ASSERT_EQ( im2.cols(), 3 );
  ASSERT_EQ( im2.rows(), 2 );
  EXPECT_EQ( im2(0,0), 5 );
  EXPECT_EQ( im2(1,0), 3 );
  EXPECT_EQ( im2(2,0), 1 );
  EXPECT_EQ( im2(0,1), 6 );
  EXPECT_EQ( im2(1,1), 4 );
  EXPECT_EQ( im2(2,1), 2 );

  // Make sure it's really shallow.
  EXPECT_EQ( rotate_90_cw(im)(2,1), im(1,0) );
  EXPECT_EQ( &(rotate_90_cw(im)(2,1)), &(im(1,0)) );

  // Test the traits
  ASSERT_TRUE( bool_trait<IsMultiplyAccessible>( rotate_90_cw(im) ) );
}

TEST( Manipulation, Rotate90CCWView ) {
  ImageView<double> im(2,3); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4; im(0,2)=5; im(1,2)=6;
  Rotate90CCWView<ImageView<double> > rmv(im);
  ASSERT_EQ( rmv.cols(), 3 );
  ASSERT_EQ( rmv.rows(), 2 );
  ASSERT_EQ( rmv.planes(), 1 );

  // Test individual pixel access
  EXPECT_EQ( rmv(0,0), 2 );
  EXPECT_EQ( rmv(1,0), 4 );
  EXPECT_EQ( rmv(2,0), 6 );
  EXPECT_EQ( rmv(0,1), 1 );
  EXPECT_EQ( rmv(1,1), 3 );
  EXPECT_EQ( rmv(2,1), 5 );

  // Test full rasterizaion
  ImageView<double> im2 = rmv;
  ASSERT_EQ( im2.cols(), rmv.cols() );
  ASSERT_EQ( im2.rows(), rmv.rows() );
  ASSERT_EQ( im2.planes(), rmv.planes() );
  for ( int r=0; r<im2.rows(); ++r )
    for ( int c=0; c<im2.cols(); ++c )
      EXPECT_EQ( im2(c,r), rmv(c,r) );

  // Test partial rasterization
  ImageView<double> im3(rmv.cols()-1,rmv.rows()-1);
  ASSERT_NO_THROW( rmv.rasterize( im3, BBox2i(1,1,rmv.cols()-1,rmv.rows()-1) ) );
  for ( int r=0; r<im3.rows(); ++r )
    for ( int c=0; c<im3.cols(); ++c )
      EXPECT_EQ( im3(c,r), rmv(c+1,r+1) );

  // Test prerasterization
  PrerasterizationTestView ptv(4,4);
  Rotate90CCWView<PrerasterizationTestView> rtv(ptv);
  ASSERT_NO_THROW( rtv.prerasterize(BBox2i(2,0,1,2)) );
  EXPECT_EQ( ptv.bbox(), BBox2i(2,2,2,1) );

  // Test the accessor / generic rasterization
  ImageView<double> im4(rmv.cols(),rmv.rows());
  vw::rasterize( rmv, im4, BBox2i(0,0,rmv.cols(),rmv.rows()) );
  for ( int r=0; r<im2.rows(); ++r )
    for ( int c=0; c<im3.cols(); ++c )
      EXPECT_EQ( im4(c,r), rmv(c,r) );

  // Test the iterator
  ImageView<double>::iterator im2i = im2.begin();
  Rotate90CCWView<ImageView<double> >::iterator rmvi = rmv.begin();
  for ( int i=0; i<im2.cols()*im2.rows(); ++i ) {
    EXPECT_GT( rmv.end() - rmvi, 0 );
    EXPECT_EQ( *im2i, *rmvi );
    ASSERT_NO_THROW( ++rmvi );
    ++im2i;
  }
  EXPECT_EQ( rmv.end() - rmvi, 0 );

  // Test the types
  ASSERT_TRUE( has_pixel_type<double>( rmv ) );
  ASSERT_TRUE( bool_trait<IsMultiplyAccessible>(rmv) );
  ASSERT_TRUE( bool_trait<IsImageView>(rmv) );
}

TEST( Manipulation, Rotate90CCW ) {
  ImageView<double> im(2,3); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4; im(0,2)=5; im(1,2)=6;
  ImageView<double> im2 = rotate_90_ccw(im);
  ASSERT_EQ( im2.cols(), 3 );
  ASSERT_EQ( im2.rows(), 2 );
  EXPECT_EQ( im2(0,0), 2 );
  EXPECT_EQ( im2(1,0), 4 );
  EXPECT_EQ( im2(2,0), 6 );
  EXPECT_EQ( im2(0,1), 1 );
  EXPECT_EQ( im2(1,1), 3 );
  EXPECT_EQ( im2(2,1), 5 );

  // Make sure it's really shallow.
  EXPECT_EQ( rotate_90_ccw(im)(0,0), im(1,0) );
  EXPECT_EQ( &(rotate_90_ccw(im)(0,0)), &(im(1,0)) );

  // Test the traits
  ASSERT_TRUE( bool_trait<IsMultiplyAccessible>( rotate_90_ccw(im) ) );
}

TEST( Manipulation, FlipVertView ) {
  ImageView<double> im(2,3); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4; im(0,2)=5; im(1,2)=6;
  FlipVerticalView<ImageView<double> > rmv(im);
  ASSERT_EQ( rmv.cols(), 2 );
  ASSERT_EQ( rmv.rows(), 3 );
  ASSERT_EQ( rmv.planes(), 1 );

  // Test individual pixel access
  EXPECT_EQ( rmv(0,0), 5 );
  EXPECT_EQ( rmv(1,0), 6 );
  EXPECT_EQ( rmv(0,1), 3 );
  EXPECT_EQ( rmv(1,1), 4 );
  EXPECT_EQ( rmv(0,2), 1 );
  EXPECT_EQ( rmv(1,2), 2 );

  // Test full rasterizaion
  ImageView<double> im2 = rmv;
  ASSERT_EQ( im2.cols(), rmv.cols() );
  ASSERT_EQ( im2.rows(), rmv.rows() );
  ASSERT_EQ( im2.planes(), rmv.planes() );
  for ( int r=0; r<im2.rows(); ++r )
    for ( int c=0; c<im2.cols(); ++c )
      EXPECT_EQ( im2(c,r), rmv(c,r) );

  // Test partial rasterization
  ImageView<double> im3(rmv.cols()-1,rmv.rows()-1);
  ASSERT_NO_THROW( rmv.rasterize( im3, BBox2i(1,1,rmv.cols()-1,rmv.rows()-1) ) );
  for ( int r=0; r<im3.rows(); ++r )
    for ( int c=0; c<im3.cols(); ++c )
      EXPECT_EQ( im3(c,r), rmv(c+1,r+1) );

  // Test prerasterization
  PrerasterizationTestView ptv(4,4);
  FlipVerticalView<PrerasterizationTestView> rtv(ptv);
  ASSERT_NO_THROW( rtv.prerasterize(BBox2i(2,0,1,2)) );
  EXPECT_EQ( ptv.bbox(), BBox2i(2,2,1,2) );

  // Test the accessor / generic rasterization
  ImageView<double> im4(rmv.cols(),rmv.rows());
  vw::rasterize( rmv, im4, BBox2i(0,0,rmv.cols(),rmv.rows()) );
  for ( int r=0; r<im2.rows(); ++r )
    for ( int c=0; c<im3.cols(); ++c )
      EXPECT_EQ( im4(c,r), rmv(c,r) );

  // Test the iterator
  ImageView<double>::iterator im2i = im2.begin();
  FlipVerticalView<ImageView<double> >::iterator rmvi = rmv.begin();
  for ( int i=0; i<im2.cols()*im2.rows(); ++i ) {
    EXPECT_GT( rmv.end() - rmvi, 0 );
    EXPECT_EQ( *im2i, *rmvi );
    ASSERT_NO_THROW( ++rmvi );
    ++im2i;
  }
  EXPECT_EQ( rmv.end() - rmvi, 0 );

  // Test the types
  ASSERT_TRUE( has_pixel_type<double>( rmv ) );
  ASSERT_TRUE( bool_trait<IsMultiplyAccessible>(rmv) );
  ASSERT_TRUE( bool_trait<IsImageView>(rmv) );
}

TEST( Manipulation, FlipVert ) {
  ImageView<double> im(2,3); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4; im(0,2)=5; im(1,2)=6;
  ImageView<double> im2 = flip_vertical(im);
  ASSERT_EQ( im2.cols(), 2 );
  ASSERT_EQ( im2.rows(), 3 );
  EXPECT_EQ( im2(0,0), 5 );
  EXPECT_EQ( im2(1,0), 6 );
  EXPECT_EQ( im2(0,1), 3 );
  EXPECT_EQ( im2(1,1), 4 );
  EXPECT_EQ( im2(0,2), 1 );
  EXPECT_EQ( im2(1,2), 2 );

  // Make sure it's really shallow.
  EXPECT_EQ( flip_vertical(im)(1,0), im(1,2) );
  EXPECT_EQ( &(flip_vertical(im)(1,0)), &(im(1,2)) );

  // Test the traits
  ASSERT_TRUE( bool_trait<IsMultiplyAccessible>( flip_vertical(im) ) );
}

TEST( Manipulation, FlipHorzView ) {
  ImageView<double> im(2,3); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4; im(0,2)=5; im(1,2)=6;
  FlipHorizontalView<ImageView<double> > rmv(im);
  ASSERT_EQ( rmv.cols(), 2 );
  ASSERT_EQ( rmv.rows(), 3 );
  ASSERT_EQ( rmv.planes(), 1 );

  // Test individual pixel access
  EXPECT_EQ( rmv(0,0), 2 );
  EXPECT_EQ( rmv(1,0), 1 );
  EXPECT_EQ( rmv(0,1), 4 );
  EXPECT_EQ( rmv(1,1), 3 );
  EXPECT_EQ( rmv(0,2), 6 );
  EXPECT_EQ( rmv(1,2), 5 );

  // Test full rasterizaion
  ImageView<double> im2 = rmv;
  ASSERT_EQ( im2.cols(), rmv.cols() );
  ASSERT_EQ( im2.rows(), rmv.rows() );
  ASSERT_EQ( im2.planes(), rmv.planes() );
  for ( int r=0; r<im2.rows(); ++r )
    for ( int c=0; c<im2.cols(); ++c )
      EXPECT_EQ( im2(c,r), rmv(c,r) );

  // Test partial rasterization
  ImageView<double> im3(rmv.cols()-1,rmv.rows()-1);
  ASSERT_NO_THROW( rmv.rasterize( im3, BBox2i(1,1,rmv.cols()-1,rmv.rows()-1) ) );
  for ( int r=0; r<im3.rows(); ++r )
    for ( int c=0; c<im3.cols(); ++c )
      EXPECT_EQ( im3(c,r), rmv(c+1,r+1) );

  // Test prerasterization
  PrerasterizationTestView ptv(4,4);
  FlipHorizontalView<PrerasterizationTestView> rtv(ptv);
  ASSERT_NO_THROW( rtv.prerasterize(BBox2i(2,0,1,2)) );
  EXPECT_EQ( ptv.bbox(), BBox2i(1,0,1,2) );

  // Test the accessor / generic rasterization
  ImageView<double> im4(rmv.cols(),rmv.rows());
  vw::rasterize( rmv, im4, BBox2i(0,0,rmv.cols(),rmv.rows()) );
  for ( int r=0; r<im2.rows(); ++r )
    for ( int c=0; c<im3.cols(); ++c )
      EXPECT_EQ( im4(c,r), rmv(c,r) );

  // Test the iterator
  ImageView<double>::iterator im2i = im2.begin();
  FlipHorizontalView<ImageView<double> >::iterator rmvi = rmv.begin();
  for ( int i=0; i<im2.cols()*im2.rows(); ++i ) {
    EXPECT_GT( rmv.end() - rmvi, 0 );
    EXPECT_EQ( *im2i, *rmvi );
    ASSERT_NO_THROW( ++rmvi );
    ++im2i;
  }
  EXPECT_EQ( rmv.end() - rmvi, 0 );

  // Test the types
  ASSERT_TRUE( has_pixel_type<double>( rmv ) );
  ASSERT_TRUE( bool_trait<IsMultiplyAccessible>(rmv) );
  ASSERT_TRUE( bool_trait<IsImageView>(rmv) );
}

TEST( Manipulation, FlipHorz ) {
  ImageView<double> im(2,3); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4; im(0,2)=5; im(1,2)=6;
  ImageView<double> im2 = flip_horizontal(im);
  ASSERT_EQ( im2.cols(), 2 );
  ASSERT_EQ( im2.rows(), 3 );
  EXPECT_EQ( im2(0,0), 2 );
  EXPECT_EQ( im2(1,0), 1 );
  EXPECT_EQ( im2(0,1), 4 );
  EXPECT_EQ( im2(1,1), 3 );
  EXPECT_EQ( im2(0,2), 6 );
  EXPECT_EQ( im2(1,2), 5 );

  // Make sure it's really shallow.
  EXPECT_EQ( flip_horizontal(im)(1,0), im(0,0) );
  EXPECT_EQ( &(flip_horizontal(im)(1,0)), &(im(0,0)) );

  // Test the traits
  ASSERT_TRUE( bool_trait<IsMultiplyAccessible>( flip_horizontal(im) ) );
}

TEST( Manipulation, CropView ) {
  ImageView<double> im(4,5); im(1,1)=1; im(2,1)=2; im(1,2)=3; im(2,2)=4; im(1,3)=5; im(2,3)=6;
  CropView<ImageView<double> > cv(im,1,1,2,3);
  ASSERT_EQ( cv.cols(), 2 );
  ASSERT_EQ( cv.rows(), 3 );
  ASSERT_EQ( cv.planes(), 1 );

  // Test individual pixel access
  EXPECT_EQ( cv(0,0), 1 );
  EXPECT_EQ( cv(1,0), 2 );
  EXPECT_EQ( cv(0,1), 3 );
  EXPECT_EQ( cv(1,1), 4 );
  EXPECT_EQ( cv(0,2), 5 );
  EXPECT_EQ( cv(1,2), 6 );

  // Test full rasterization
  ImageView<double> im2 = cv;
  ASSERT_EQ( im2.cols(), cv.cols() );
  ASSERT_EQ( im2.rows(), cv.rows() );
  ASSERT_EQ( im2.planes(), cv.planes() );
  for ( int r=0; r<im2.rows(); ++r )
    for ( int c=0; c<im2.cols(); ++c )
      EXPECT_EQ( im2(c,r), cv(c,r) );

  // Test partial rasterization
  ImageView<double> im3(cv.cols()-1,cv.rows()-1);
  ASSERT_NO_THROW( cv.rasterize( im3, BBox2i(1,1,cv.cols()-1,cv.rows()-1) ) );
  for ( int r=0; r<im3.rows(); ++r )
    for ( int c=0; c<im3.cols(); ++c )
      EXPECT_EQ( im3(c,r), cv(c+1,r+1) );

  // Test prerasterization
  PrerasterizationTestView ptv(4,4);
  CropView<PrerasterizationTestView> rtv(ptv,0,1,3,3);
  ASSERT_NO_THROW( rtv.prerasterize(BBox2i(2,0,1,2)) );
  EXPECT_EQ( ptv.bbox(), BBox2i(2,1,1,2) );

  // Test the accessor / generic rasterization
  ImageView<double> im4(cv.cols(),cv.rows());
  vw::rasterize( cv, im4, BBox2i(0,0,cv.cols(),cv.rows()) );
  for ( int r=0; r<im2.rows(); ++r )
    for ( int c=0; c<im3.cols(); ++c )
      EXPECT_EQ( im4(c,r), cv(c,r) );

  // Test the iterator
  ImageView<double>::iterator im2i = im2.begin();
  CropView<ImageView<double> >::iterator cvi = cv.begin();
  for ( int i=0; i<im2.cols()*im2.rows(); ++i ) {
    EXPECT_GT( cv.end() - cvi, 0 );
    EXPECT_EQ( *im2i, *cvi );
    ASSERT_NO_THROW( ++cvi );
    ++im2i;
  }
  EXPECT_EQ( cv.end() - cvi, 0 );

  // Test the types
  ASSERT_TRUE( has_pixel_type<double>( cv ) );
  ASSERT_TRUE( bool_trait<IsMultiplyAccessible>(cv) );
  ASSERT_TRUE( bool_trait<IsImageView>(cv) );
}

TEST( Manipulation, Crop ) {
  ImageView<double> im(2,3); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4; im(0,2)=5; im(1,2)=6;
  ImageView<double> im2 = crop(im,1,1,1,2);
  ASSERT_EQ( im2.cols(), 1 );
  ASSERT_EQ( im2.rows(), 2 );
  EXPECT_EQ( im2(0,0), 4 );
  EXPECT_EQ( im2(0,1), 6 );

  ImageView<double> im3 = crop(im,BBox2i(1,1,1,2));
  ASSERT_EQ( im3.cols(), 1 );
  ASSERT_EQ( im3.rows(), 2 );
  EXPECT_EQ( im3(0,0), 4 );
  EXPECT_EQ( im3(0,1), 6 );

  // Make sure it's really shallow.
  EXPECT_EQ( crop(im,1,1,1,2)(0,0), im(1,1) );
  EXPECT_EQ( &(crop(im,1,1,1,2)(0,0)), &(im(1,1)) );

  // Test the traits
  ASSERT_TRUE( bool_trait<IsMultiplyAccessible>( crop(im,1,1,1,2) ) );
}

TEST( Manipulation, CropAssignment ) {
  ImageView<int> dest(3,3);
  fill(dest,1);
  ImageView<int> a(2,2);
  a(0,0) = 1; // A|13|
  a(0,1) = 2; //  |24|
  a(1,0) = 3;
  a(1,1) = 4;

  crop(dest,1,1,2,2) = a;
  EXPECT_EQ( 1, dest(0,0) );
  EXPECT_EQ( 1, dest(1,0) );
  EXPECT_EQ( 1, dest(2,0) );
  EXPECT_EQ( 1, dest(0,1) );
  EXPECT_EQ( 1, dest(1,1) );
  EXPECT_EQ( 3, dest(2,1) );
  EXPECT_EQ( 1, dest(0,2) );
  EXPECT_EQ( 2, dest(1,2) );
  EXPECT_EQ( 4, dest(2,2) );

  crop(dest,1,0,2,1) = crop(a,0,1,2,1);
  EXPECT_EQ( 1, dest(0,0) );
  EXPECT_EQ( 2, dest(1,0) );
  EXPECT_EQ( 4, dest(2,0) );
  EXPECT_EQ( 1, dest(0,1) );
  EXPECT_EQ( 1, dest(1,1) );
  EXPECT_EQ( 3, dest(2,1) );
  EXPECT_EQ( 1, dest(0,2) );
  EXPECT_EQ( 2, dest(1,2) );
  EXPECT_EQ( 4, dest(2,2) );
}

TEST( Manipulation, SubsampleView ) {
  ImageView<double> im(4,5); im(0,0)=1; im(2,0)=2; im(0,2)=3; im(2,2)=4; im(0,4)=5; im(2,4)=6;
  SubsampleView<ImageView<double> > ssv(im,2,2);
  ASSERT_EQ( ssv.cols(), 2 );
  ASSERT_EQ( ssv.rows(), 3 );
  ASSERT_EQ( ssv.planes(), 1 );

  // Test individual pixel access
  EXPECT_EQ( ssv(0,0), 1 );
  EXPECT_EQ( ssv(1,0), 2 );
  EXPECT_EQ( ssv(0,1), 3 );
  EXPECT_EQ( ssv(1,1), 4 );
  EXPECT_EQ( ssv(0,2), 5 );
  EXPECT_EQ( ssv(1,2), 6 );

  // Test full rasterizaion
  ImageView<double> im2 = ssv;
  ASSERT_EQ( im2.cols(), ssv.cols() );
  ASSERT_EQ( im2.rows(), ssv.rows() );
  ASSERT_EQ( im2.planes(), ssv.planes() );
  for ( int r=0; r<im2.rows(); ++r )
    for ( int c=0; c<im2.cols(); ++c )
      EXPECT_EQ( im2(c,r), ssv(c,r) );

  // Test partial rasterization
  ImageView<double> im3(ssv.cols()-1,ssv.rows()-1);
  ASSERT_NO_THROW( ssv.rasterize( im3, BBox2i(1,1,ssv.cols()-1,ssv.rows()-1) ) );
  for ( int r=0; r<im3.rows(); ++r )
    for ( int c=0; c<im3.cols(); ++c )
      EXPECT_EQ( im3(c,r), ssv(c+1,r+1) );

  // Test prerasterization
  PrerasterizationTestView ptv(4,4);
  SubsampleView<PrerasterizationTestView> rtv(ptv,2,2);
  ASSERT_NO_THROW( rtv.prerasterize(BBox2i(1,0,1,2)) );
  EXPECT_EQ( ptv.bbox(), BBox2i(2,0,2,4) );

  // Test the accessor / generic rasterization
  ImageView<double> im4(ssv.cols(),ssv.rows());
  vw::rasterize( ssv, im4, BBox2i(0,0,ssv.cols(),ssv.rows()) );
  for ( int r=0; r<im2.rows(); ++r )
    for ( int c=0; c<im3.cols(); ++c )
      EXPECT_EQ( im4(c,r), ssv(c,r) );

  // Test the iterator
  ImageView<double>::iterator im2i = im2.begin();
  SubsampleView<ImageView<double> >::iterator ssvi = ssv.begin();
  for ( int i=0; i<im2.cols()*im2.rows(); ++i ) {
    EXPECT_GT( ssv.end() - ssvi, 0 );
    EXPECT_EQ( *im2i, *ssvi );
    ASSERT_NO_THROW( ++ssvi );
    ++im2i;
  }
  EXPECT_EQ( ssv.end() - ssvi, 0 );

  // Test the types
  ASSERT_TRUE( has_pixel_type<double>( ssv ) );
  ASSERT_TRUE( bool_trait<IsMultiplyAccessible>(ssv) );
  ASSERT_TRUE( bool_trait<IsImageView>(ssv) );
}

TEST( Manipulation, Subsample ) {
  ImageView<double> im(2,3); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4; im(0,2)=5; im(1,2)=6;
  ImageView<double> im2 = subsample(im,1,2);
  ASSERT_EQ( im2.cols(), 2 );
  ASSERT_EQ( im2.rows(), 2 );
  EXPECT_EQ( im2(0,0), 1 );
  EXPECT_EQ( im2(1,0), 2 );
  EXPECT_EQ( im2(0,1), 5 );
  EXPECT_EQ( im2(1,1), 6 );

  // Make sure it's really shallow.
  EXPECT_EQ( subsample(im,2)(0,1), im(0,2) );
  EXPECT_EQ( &(subsample(im,2)(0,1)), &(im(0,2)) );

  // Test the traits
  ASSERT_TRUE( bool_trait<IsMultiplyAccessible>( subsample(im,1,2) ) );
}

TEST( Manipulation, SelectCol ) {
  ImageView<double> im(2,2); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4;
  ImageView<double> im2 = select_col(im,1);
  ASSERT_EQ( im2.cols(), 1 );
  ASSERT_EQ( im2.rows(), 2 );
  EXPECT_EQ( im2.planes(), 1 );
  EXPECT_EQ( im2(0,0), 2 );
  EXPECT_EQ( im2(0,1), 4 );

  // Make sure it's really shallow.
  EXPECT_EQ( select_col(im,1)(0,1), im(1,1) );
  EXPECT_EQ( &(select_col(im,1)(0,1)), &(im(1,1)) );

  // Test the traits
  ASSERT_TRUE( bool_trait<IsMultiplyAccessible>( select_col(im,1) ) );
}

TEST( Manipulation, ColAssignment ) {
  ImageView<int32> a(2,2), b(2,2);
  a(0,0) = 1; // A|13| B|57|
  a(0,1) = 2; //  |24|  |68|
  a(1,0) = 3;
  a(1,1) = 4;
  b(0,0) = 5;
  b(0,1) = 6;
  b(1,0) = 7;
  b(1,1) = 8;

  ImageView<int32> tmp(1,2);
  tmp(0,0) = 9;
  tmp(0,1) = 10;

  select_col(a,1) = tmp;
  EXPECT_EQ( 1,  a(0,0) );
  EXPECT_EQ( 9,  a(1,0) );
  EXPECT_EQ( 2,  a(0,1) );
  EXPECT_EQ( 10, a(1,1) );

  // Remember:
  // Matrices index (row, col)
  // however Images index (col, row)

  select_col(a, 0) = select_col(b, 0);
  EXPECT_EQ( 5,  a(0,0) );
  EXPECT_EQ( 9,  a(1,0) );
  EXPECT_EQ( 6,  a(0,1) );
  EXPECT_EQ( 10, a(1,1) );
}


TEST( Manipulation, SelectRow ) {
  ImageView<double> im(2,2); im(0,0)=1; im(1,0)=2; im(0,1)=3; im(1,1)=4;
  ImageView<double> im2 = select_row(im,1);
  ASSERT_EQ( im2.cols(), 2 );
  ASSERT_EQ( im2.rows(), 1 );
  EXPECT_EQ( im2.planes(), 1 );
  EXPECT_EQ( im2(0,0), 3 );
  EXPECT_EQ( im2(1,0), 4 );

  // Make sure it's really shallow.
  EXPECT_EQ( select_row(im,1)(1,0), im(1,1) );
  EXPECT_EQ( &(select_row(im,1)(1,0)), &(im(1,1)) );

  // Test the traits
  ASSERT_TRUE( bool_trait<IsMultiplyAccessible>( select_row(im,1) ) );
}

TEST( Manipulation, RowAssignment ) {
  ImageView<int32> a(2,2), b(2,2);
  a(0,0) = 1; // A|13| B|57|
  a(0,1) = 2; //  |24|  |68|
  a(1,0) = 3;
  a(1,1) = 4;
  b(0,0) = 5;
  b(0,1) = 6;
  b(1,0) = 7;
  b(1,1) = 8;

  ImageView<int32> tmp(2,1);
  tmp(0,0) = 9;
  tmp(1,0) = 10;

  select_row(a,1) = tmp;
  EXPECT_EQ( 1,  a(0,0) );
  EXPECT_EQ( 3,  a(1,0) );
  EXPECT_EQ( 9,  a(0,1) );
  EXPECT_EQ( 10, a(1,1) );

  select_row(a,0) = select_row(b, 0);
  EXPECT_EQ( 5,  a(0,0) );
  EXPECT_EQ( 7,  a(1,0) );
  EXPECT_EQ( 9,  a(0,1) );
  EXPECT_EQ( 10, a(1,1) );
}

TEST( Manipulation, SelectPlane ) {
  ImageView<double> im(1,2,2); im(0,0,0)=1; im(0,1,0)=2; im(0,0,1)=3; im(0,1,1)=4;
  ImageView<double> im2 = select_plane(im,1);
  ASSERT_EQ( im2.cols(), 1 );
  ASSERT_EQ( im2.rows(), 2 );
  ASSERT_EQ( im2.planes(), 1 );
  EXPECT_EQ( im2(0,0), 3 );
  EXPECT_EQ( im2(0,1), 4 );

  // Make sure it's really shallow.
  EXPECT_EQ( select_plane(im,0)(0,1), im(0,1) );
  EXPECT_EQ( &(select_plane(im,0)(0,1)), &(im(0,1)) );

  // Test the traits
  ASSERT_TRUE( bool_trait<IsMultiplyAccessible>( select_plane(im,1) ) );
}

TEST( Manipulation, SelectChannel ) {
  ImageView<PixelRGB<double> > im(1,2); im(0,0)=PixelRGB<double>(1,2,3); im(0,1)=PixelRGB<double>(4,5,6);
  ImageView<double> im2 = select_channel(im,1);
  ASSERT_EQ( im2.cols(), 1 );
  ASSERT_EQ( im2.rows(), 2 );
  ASSERT_EQ( im2.planes(), 1 );
  EXPECT_EQ( im2(0,0), 2 );
  EXPECT_EQ( im2(0,1), 5 );

  // Make sure it's really shallow.
  EXPECT_EQ( select_channel(im,1)(0,1), im(0,1)[1] );
  EXPECT_EQ( &(select_channel(im,1)(0,1)), &(im(0,1)[1]) );

  // Test the traits
  ASSERT_TRUE( bool_trait<IsMultiplyAccessible>( select_channel(im,1) ) );
  ASSERT_TRUE( has_pixel_type<double>( select_channel(im,1) ) );
  ASSERT_TRUE( has_result_type<double&>( select_channel(im,1) ) );
  ASSERT_TRUE( has_pixel_type<double>( select_channel(per_pixel_filter(im,&make_const_ref<PixelRGB<double> >),1) ) );
  ASSERT_TRUE( has_result_type<double const&>( select_channel(per_pixel_filter(im,&make_const_ref<PixelRGB<double> >),1) ) );

  // Test a non-writable view
  ImageViewRef<PixelRGB<double> > im3 = im;
  ASSERT_NO_THROW( ImageView<double> im4 = select_channel(im3,1) );
}

TEST( Manipulation, ChannelsToPlanes ) {
  ImageView<PixelRGB<double> > im(1,2); im(0,0)=PixelRGB<double>(1,2,3); im(0,1)=PixelRGB<double>(4,5,6);
  ImageView<double> im2 = channels_to_planes(im);
  ASSERT_EQ( im2.cols(), 1 );
  ASSERT_EQ( im2.rows(), 2 );
  ASSERT_EQ( im2.planes(), 3 );
  EXPECT_EQ( im2(0,0,0), 1 );
  EXPECT_EQ( im2(0,1,0), 4 );
  EXPECT_EQ( im2(0,0,1), 2 );
  EXPECT_EQ( im2(0,1,1), 5 );
  EXPECT_EQ( im2(0,0,2), 3 );
  EXPECT_EQ( im2(0,1,2), 6 );

  // Make sure it's really shallow.
  EXPECT_EQ( channels_to_planes(im)(0,1,1), im(0,1)[1] );
  EXPECT_EQ( &(channels_to_planes(im)(0,1,1)), &(im(0,1)[1]) );

  // Test the traits
  ASSERT_TRUE( bool_trait<IsMultiplyAccessible>( channels_to_planes(im) ) );
}

TEST( Manipulation, PlanesToChannels ) {
  ImageView<double> im(1,2,3); im(0,0,0)=1; im(0,0,1)=2; im(0,0,2)=3; im(0,1,0)=4; im(0,1,1)=5; im(0,1,2)=6;
  ImageView<PixelRGB<double> > im2 = planes_to_channels<PixelRGB<double> >(im);
  ASSERT_EQ( im2.cols(), 1 );
  ASSERT_EQ( im2.rows(), 2 );
  ASSERT_EQ( im2.planes(), 1 );
  EXPECT_EQ( im2(0,0)[0], 1 );
  EXPECT_EQ( im2(0,1)[0], 4 );
  EXPECT_EQ( im2(0,0)[1], 2 );
  EXPECT_EQ( im2(0,1)[1], 5 );
  EXPECT_EQ( im2(0,0)[2], 3 );
  EXPECT_EQ( im2(0,1)[2], 6 );

  // Test the traits
  ASSERT_TRUE( !bool_trait<IsMultiplyAccessible>( planes_to_channels<PixelRGB<double> >(im) ) );
}

TEST( Manipulation, ChannelCast ) {
  ImageView<PixelRGB<double> > im(1,2); im(0,0)=PixelRGB<double>(1,2,3); im(0,1)=PixelRGB<double>(4,5,6);
  ImageView<PixelRGB<float> > im2 = channel_cast<float>(im);
  ASSERT_EQ( im2.cols(), 1 );
  ASSERT_EQ( im2.rows(), 2 );
  ASSERT_EQ( im2.planes(), 1 );
  EXPECT_EQ( im2(0,0)[0], im(0,0)[0] );
  EXPECT_EQ( im2(0,0)[1], im(0,0)[1] );
  EXPECT_EQ( im2(0,0)[2], im(0,0)[2] );
  EXPECT_EQ( im2(0,1)[0], im(0,1)[0] );
  EXPECT_EQ( im2(0,1)[1], im(0,1)[1] );
  EXPECT_EQ( im2(0,1)[2], im(0,1)[2] );

  // Test the traits
  ASSERT_FALSE( bool_trait<IsMultiplyAccessible>( channel_cast<float>(im) ) );
}

TEST( Manipulation, ChannelCastRescale ) {
  ImageView<PixelRGB<uint8> > im(1,2); im(0,0)=PixelRGB<uint8>(91,48,227); im(0,1)=PixelRGB<uint8>(53,189,98);
  ImageView<PixelRGB<float> > im2 = channel_cast_rescale<float>(im);
  ASSERT_EQ( im2.cols(), 1 );
  ASSERT_EQ( im2.rows(), 2 );
  ASSERT_EQ( im2.planes(), 1 );
  EXPECT_NEAR( im2(0,0)[0], im(0,0)[0]/255.0, 1e-7 );
  EXPECT_NEAR( im2(0,0)[1], im(0,0)[1]/255.0, 1e-7 );
  EXPECT_NEAR( im2(0,0)[2], im(0,0)[2]/255.0, 1e-7 );
  EXPECT_NEAR( im2(0,1)[0], im(0,1)[0]/255.0, 1e-7 );
  EXPECT_NEAR( im2(0,1)[1], im(0,1)[1]/255.0, 1e-7 );
  EXPECT_NEAR( im2(0,1)[2], im(0,1)[2]/255.0, 1e-7 );

  // Test the traits
  ASSERT_FALSE( bool_trait<IsMultiplyAccessible>( channel_cast_rescale<float>(im) ) );
}

TEST( Manipulation, PixelCast ) {
  ImageView<PixelRGB<double> > im(1,2); im(0,0)=PixelRGB<double>(1,2,3); im(0,1)=PixelRGB<double>(4,5,6);
  ImageView<PixelRGBA<double> > im2 = pixel_cast<PixelRGBA<double> >(im);
  ASSERT_EQ( im2.cols(), 1 );
  ASSERT_EQ( im2.rows(), 2 );
  ASSERT_EQ( im2.planes(), 1 );
  EXPECT_EQ( im2(0,0)[0], im(0,0)[0] );
  EXPECT_EQ( im2(0,0)[1], im(0,0)[1] );
  EXPECT_EQ( im2(0,0)[2], im(0,0)[2] );
  EXPECT_EQ( im2(0,0)[3], 1.0 );
  EXPECT_EQ( im2(0,1)[0], im(0,1)[0] );
  EXPECT_EQ( im2(0,1)[1], im(0,1)[1] );
  EXPECT_EQ( im2(0,1)[2], im(0,1)[2] );
  EXPECT_EQ( im2(0,1)[3], 1.0 );

  // Test the traits
  ASSERT_FALSE( bool_trait<IsMultiplyAccessible>( pixel_cast<PixelRGBA<double> >(im) ) );
}

TEST( Manipulation, Rescale ) {
  ImageView<PixelRGB<uint8> > im(1,2); im(0,0)=PixelRGB<uint8>(91,48,227); im(0,1)=PixelRGB<uint8>(53,189,98);
  ImageView<PixelRGBA<float> > im2 = pixel_cast_rescale<PixelRGBA<float> >(im);
  ASSERT_EQ( im2.cols(), 1 );
  ASSERT_EQ( im2.rows(), 2 );
  ASSERT_EQ( im2.planes(), 1 );
  EXPECT_NEAR( im2(0,0)[0], im(0,0)[0]/255.0, 1e-7 );
  EXPECT_NEAR( im2(0,0)[1], im(0,0)[1]/255.0, 1e-7 );
  EXPECT_NEAR( im2(0,0)[2], im(0,0)[2]/255.0, 1e-7 );
  EXPECT_EQ( im2(0,0)[3], 1.0 );
  EXPECT_NEAR( im2(0,1)[0], im(0,1)[0]/255.0, 1e-7 );
  EXPECT_NEAR( im2(0,1)[1], im(0,1)[1]/255.0, 1e-7 );
  EXPECT_NEAR( im2(0,1)[2], im(0,1)[2]/255.0, 1e-7 );
  EXPECT_EQ( im2(0,1)[3], 1.0 );

  // Test the traits
  ASSERT_FALSE( bool_trait<IsMultiplyAccessible>( pixel_cast_rescale<PixelRGBA<double> >(im) ) );
}
