// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// TestCameraBBox.h
#include <gtest/gtest.h>
#include <test/Helpers.h>

#include <vw/Image/ImageViewRef.h>
#include <vw/Cartography/OrthoImageView.h>
#include <vw/Cartography/SimplePointImageManipulation.h>
#include <vw/Camera/PinholeModel.h>
#include <vw/Image/Transform.h>

using namespace vw;
using namespace vw::cartography;
using namespace vw::test;
using namespace vw::camera;

template <template<class> class TraitT, class T>
static bool bool_trait( T const& /*arg*/ ) {
  return TraitT<T>::value;
}

// Must have protobuf to be able to read camera
#if defined(VW_HAVE_PKG_PROTOBUF) && VW_HAVE_PKG_PROTOBUF==1

// Test Pattern Image
template <class PixelT>
class TestPatternView : public ImageViewBase<TestPatternView<PixelT> > {
  int32 m_cols, m_rows, m_planes;
public:
  typedef PixelT pixel_type;
  typedef PixelT result_type;
  typedef ProceduralPixelAccessor<TestPatternView> pixel_accessor;

  TestPatternView( int32 cols, int32 rows, int32 planes = 1 )
    : m_cols(cols), m_rows(rows), m_planes(planes) {}

  inline int32 cols() const { return m_cols; }
  inline int32 rows() const { return m_rows; }
  inline int32 planes() const { return m_planes; }

  inline pixel_accessor origin() const { return pixel_accessor( *this ); }

  inline result_type operator()( int32 col, int32 row, int32 /*plane*/=0 ) const {
    typedef typename PixelChannelType<result_type>::type ChannelT;
    int32 intermediate = (col%1024)*(row%1024);
    return result_type(ChannelT(ChannelRange<result_type>::max() *
                                float(intermediate)/(1024.0f*1024.0f)));
  }

  typedef TestPatternView prerasterize_type;
  inline prerasterize_type prerasterize( BBox2i /*bbox*/ ) const { return *this; }
  template <class DestT> inline void rasterize( DestT const& dest, BBox2i bbox ) const {
    vw::rasterize( prerasterize(bbox), dest, bbox );
  }
};

namespace vw {
  template <class PixelT>
  struct IsMultiplyAccessible<TestPatternView<PixelT> > : public true_type {};
}

template <class PixelT>
TestPatternView<PixelT> test_pattern_view( PixelT const& /*value*/, int32 cols, int32 rows, int32 planes=1 ) {
  return TestPatternView<PixelT>( cols, rows, planes );
}


// Test Framework
class OrthoImageTest :  public ::testing::Test {
protected:
  virtual void SetUp() {
    apollo = boost::shared_ptr<CameraModel>(new PinholeModel("apollo.pinhole"));
    moon.set_well_known_geogcs("D_MOON");

    DEM.set_size(20,20);
    for ( uint32 i = 0; i < 20; i++ )
      for ( uint32 j = 0; j <20; j++ )
        DEM(i,j) = 10000 - 100*(pow(10.-i,2.)+pow(10.-j,2));
    Matrix<double> geotrans = vw::math::identity_matrix<3>();
    geotrans(0,2) = 85;
    geotrans(0,0) = 0.5;
    geotrans(1,1) = -0.5;
    geotrans(1,2) = 5;
    moon.set_transform(geotrans);
  }

  boost::shared_ptr<CameraModel> apollo;
  GeoReference moon;
  ImageView<float> DEM;
};

TEST_F( OrthoImageTest, OrthoImageRun ) {
  // These tests need to be careful. ImageViewRefs can't be used since
  // they are not
  {
    ImageView<PixelGray<uint8> > ortho_right =
      resample(orthoproject( interpolate(DEM,BicubicInterpolation()), moon,
                             test_pattern_view(PixelGray<uint8>(),5725,5725),
                             apollo, BicubicInterpolation(),
                             ZeroEdgeExtension()),16);
    ImageView<PixelGray<uint8> > ortho_wrong =
      resample(orthoproject( DEM, moon,
                             test_pattern_view(PixelGray<uint8>(),5725,5725),
                             apollo, BicubicInterpolation(),
                             ZeroEdgeExtension()),16);
    // Verifying that interpolating the DEM removes other interpolation error
    EXPECT_LT( ortho_right(142,16), ortho_right(146,16) );
    EXPECT_LT( ortho_wrong(142,16), ortho_right(146,16) );
    EXPECT_LE( ortho_wrong(146,16), ortho_right(146,16) );
    EXPECT_EQ( ortho_wrong(144,16), ortho_right(144,16) ); // check same
  }

  ImageView<PixelGray<uint8> > ortho =
    orthoproject( interpolate(DEM,BicubicInterpolation()), moon,
                  test_pattern_view(PixelGray<uint8>(),5725,5725),
                  apollo, BicubicInterpolation(),
                  ZeroEdgeExtension() );


  // Verifying the color
  Vector2 lonlat = moon.pixel_to_lonlat(Vector2(12,4));
  Vector3 xyz = LonLatRadToXYZFunctor()(Vector3(lonlat[0], lonlat[1],
                                                1737400+DEM(12,4)));
  Vector2i camloc = apollo->point_to_pixel( xyz );
  EXPECT_GT( camloc[0], 0 );
  EXPECT_GT( camloc[1], 0 );
  EXPECT_LT( camloc[0], 5725 );
  EXPECT_LT( camloc[1], 5725 );
  EXPECT_EQ( ortho(12,4)[0], test_pattern_view(PixelGray<uint8>(), 5725, 5725)(camloc[0],camloc[1])[0] );
}


TEST_F( OrthoImageTest, OrthoTraits ) {
  OrthoImageView<ImageView<float>,TestPatternView<PixelGray<uint8> >,
    BicubicInterpolation, ZeroEdgeExtension>
  ortho_nonfloat( DEM, moon, TestPatternView<PixelGray<uint8> >(5725,5725),
                  apollo, BicubicInterpolation(), ZeroEdgeExtension() );

  ASSERT_FALSE( bool_trait<IsFloatingPointIndexable>( ortho_nonfloat ) );
  ASSERT_FALSE( bool_trait<IsMultiplyAccessible>( ortho_nonfloat ) );

  OrthoImageView<InterpolationView<ImageView<float>,BicubicInterpolation>,TestPatternView<PixelGray<uint8> >,
    BicubicInterpolation, ZeroEdgeExtension>
  ortho_float( InterpolationView<ImageView<float>,BicubicInterpolation>(DEM),
               moon, TestPatternView<PixelGray<uint8> >(5725,5725),
               apollo, BicubicInterpolation(), ZeroEdgeExtension() );

  ASSERT_TRUE( bool_trait<IsFloatingPointIndexable>( ortho_float ) );
  ASSERT_FALSE( bool_trait<IsMultiplyAccessible>( ortho_float ) );

  // Additional tests to verify that interpolate is not messing up orthoproject
  ASSERT_TRUE( bool_trait<IsFloatingPointIndexable>(orthoproject(
                  interpolate(DEM,BicubicInterpolation()), moon,
                  test_pattern_view(PixelGray<uint8>(),5725,5725),
                  apollo, BicubicInterpolation(),
                  ZeroEdgeExtension())) );
  ASSERT_FALSE( bool_trait<IsMultiplyAccessible>(orthoproject(
                  interpolate(DEM,BicubicInterpolation()), moon,
                  test_pattern_view(PixelGray<uint8>(),5725,5725),
                  apollo, BicubicInterpolation(),
                  ZeroEdgeExtension())) );
  // Edge extend is not floating point indexable
  ASSERT_FALSE( bool_trait<IsFloatingPointIndexable>(edge_extend(orthoproject(
                  interpolate(DEM,BicubicInterpolation()), moon,
                  test_pattern_view(PixelGray<uint8>(),5725,5725),
                  apollo, BicubicInterpolation(),
                  ZeroEdgeExtension()), ZeroEdgeExtension())) );
}

#endif
