// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NASA Vision Workbench is licensed under the Apache License,
//  Version 2.0 (the "License"); you may not use this file except in
//  compliance with the License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__


// TestCameraBBox.h
#include <gtest/gtest_VW.h>
#include <test/Helpers.h>

#include <vw/Image/ImageViewRef.h>
#include <vw/Cartography/OrthoImageView.h>
#include <vw/Cartography/PointImageManipulation.h>
#include <vw/Camera/PinholeModel.h>
#include <vw/Image/Transform.h>

// Must have protobuf to be able to read camera
#if defined(VW_HAVE_PKG_PROTOBUF) && defined(VW_HAVE_PKG_CAMERA) 

using namespace vw;
using namespace vw::cartography;
using namespace vw::test;
using namespace vw::camera;

template <template<class> class TraitT, class T>
static bool bool_trait( T const& /*arg*/ ) {
  return TraitT<T>::value;
}

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

  inline int32 cols  () const { return m_cols;   }
  inline int32 rows  () const { return m_rows;   }
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
                             apollo.get(), BicubicInterpolation(),
                             ZeroEdgeExtension()),16);
    ImageView<PixelGray<uint8> > ortho_wrong =
      resample(orthoproject( DEM, moon,
                             test_pattern_view(PixelGray<uint8>(),5725,5725),
                             apollo.get(), BicubicInterpolation(),
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
                  apollo.get(), BicubicInterpolation(),
                  ZeroEdgeExtension() );


  // Verifying the color - Maybe a problem since these are not accurate?
  Vector2 lonlat = moon.pixel_to_lonlat(Vector2(12,4));
  Vector3 xyz = LonLatRadToXYZEstimateFunctor()(Vector3(lonlat[0], lonlat[1],
                                                1737400+DEM(12,4)));
  Vector2i camloc = apollo->point_to_pixel( xyz );
  EXPECT_GT( camloc[0], 0 );
  EXPECT_GT( camloc[1], 0 );
  EXPECT_LT( camloc[0], 5725 );
  EXPECT_LT( camloc[1], 5725 );
  EXPECT_EQ( ortho(12,4)[0], test_pattern_view(PixelGray<uint8>(), 5725, 5725)(camloc[0],camloc[1])[0] );

  // Test the version of orthoproject which distinguishes no-data vs
  // no-processed-data.
  ImageView<PixelMask<float> > maskedDEM = DEM;
  for (int col =  0; col < maskedDEM.cols()/2; col++){ 
    for (int row =  0; row < maskedDEM.rows()/2; row++){
      // Make some pixels transparent in the input
      maskedDEM(col, row).invalidate();
    }
  }
  ImageView<PixelGrayA<uint8> > ortho_NP =
    orthoproject_markNoProcessedData(interpolate(maskedDEM, BicubicInterpolation()), moon,
                                     test_pattern_view(PixelGrayA<uint8>(),5725,5725),
                                     apollo.get(), BicubicInterpolation(),
                                     ZeroEdgeExtension() );
  EXPECT_EQ(ortho_NP(0, 0),  PixelGrayA<uint8>(0,0));    // transparent pixel, no-data
  EXPECT_EQ(ortho_NP(5, 6),  PixelGrayA<uint8>(0,255));  // black pixel, no-processed-data
  EXPECT_EQ(ortho_NP(19, 4), PixelGrayA<uint8>(27,255)); // non-transparent and non-black pixel
    
}


TEST_F( OrthoImageTest, OrthoTraits ) {
  OrthoImageView<ImageView<float>,TestPatternView<PixelGray<uint8> >,
    BicubicInterpolation, ZeroEdgeExtension, false>
  ortho_nonfloat( DEM, moon, TestPatternView<PixelGray<uint8> >(5725,5725),
                  apollo.get(), BicubicInterpolation(), ZeroEdgeExtension() );

  ASSERT_FALSE( bool_trait<IsFloatingPointIndexable>( ortho_nonfloat ) );
  ASSERT_FALSE( bool_trait<IsMultiplyAccessible>( ortho_nonfloat ) );

  OrthoImageView<InterpolationView<ImageView<float>,BicubicInterpolation>,TestPatternView<PixelGray<uint8> >,
    BicubicInterpolation, ZeroEdgeExtension, false>
  ortho_float( InterpolationView<ImageView<float>,BicubicInterpolation>(DEM),
               moon, TestPatternView<PixelGray<uint8> >(5725,5725),
               apollo.get(), BicubicInterpolation(), ZeroEdgeExtension() );

  ASSERT_TRUE( bool_trait<IsFloatingPointIndexable>( ortho_float ) );
  ASSERT_FALSE( bool_trait<IsMultiplyAccessible>( ortho_float ) );

  // Additional tests to verify that interpolate is not messing up orthoproject
  ASSERT_TRUE( bool_trait<IsFloatingPointIndexable>(orthoproject(
                  interpolate(DEM,BicubicInterpolation()), moon,
                  test_pattern_view(PixelGray<uint8>(),5725,5725),
                  apollo.get(), BicubicInterpolation(),
                  ZeroEdgeExtension())) );
  ASSERT_FALSE( bool_trait<IsMultiplyAccessible>(orthoproject(
                  interpolate(DEM,BicubicInterpolation()), moon,
                  test_pattern_view(PixelGray<uint8>(),5725,5725),
                  apollo.get(), BicubicInterpolation(),
                  ZeroEdgeExtension())) );
  // Edge extend is not floating point indexable
  ASSERT_FALSE( bool_trait<IsFloatingPointIndexable>(edge_extend(orthoproject(
                  interpolate(DEM,BicubicInterpolation()), moon,
                  test_pattern_view(PixelGray<uint8>(),5725,5725),
                  apollo.get(), BicubicInterpolation(),
                  ZeroEdgeExtension()), ZeroEdgeExtension())) );
}

#endif
