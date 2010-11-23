// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <gtest/gtest.h>
#include <test/Helpers.h>
#include <vw/Plate/PlateManager.h>
#include <vw/Plate/PlateCarreePlateManager.h>
#include <vw/Cartography/GeoReference.h>

using namespace std;
using namespace vw;
using namespace vw::platefile;
using namespace vw::test;
using namespace vw::cartography;

template <class PixelT>
class PlateCarreeExposed : public PlateCarreePlateManager<PixelT> {
public:
  PlateCarreeExposed(boost::shared_ptr<PlateFile> platefile) :
    PlateCarreePlateManager<PixelT>(platefile) {}

  void transform_image(cartography::GeoReference const& georef,
                       ImageViewRef<PixelT>& image,
                       TransformRef& txref, int& level ) const {
    PlateCarreePlateManager<PixelT>::transform_image(georef,image,txref,level);
  }

  void affected_tiles(BBox2i const& image_size,
                      TransformRef const& tx, int tile_size,
                      int level, std::list<TileInfo>& tiles ) const {
    PlateCarreePlateManager<PixelT>::affected_tiles(image_size,tx,tile_size,
                                                    level,tiles);
  }
};

class PlateCarreeExposedTest : public ::testing::Test {
protected:
  typedef PixelGrayA<uint8> PixelT;

  virtual void SetUp() {
    image = ImageView<PixelT>(20,20);
    std::fill(image.begin(),image.end(), PixelT(255,255) );
    image_ref = image;

    platename = UnlinkName("test.plate");
    platefile =
      boost::shared_ptr<PlateFile>( new PlateFile( Url(platename),
                                                   "", "", 256, "",
                                                   VW_PIXEL_GRAYA,
                                                   VW_CHANNEL_UINT8) );
    platemanager = boost::shared_ptr<PlateCarreeExposed<PixelT> >( new PlateCarreeExposed<PixelT>( platefile ) );
  }

  UnlinkName platename;
  ImageView<PixelT> image;
  ImageViewRef<PixelT> image_ref;
  boost::shared_ptr<PlateFile> platefile;
  boost::shared_ptr<PlateCarreeExposed<PixelT> > platemanager;
};

TEST_F( PlateCarreeExposedTest, WestWrap ) {
  Matrix3x3 affine = identity_matrix<3>();
  affine(1,1) = -1;
  affine(0,2) = -190;
  affine(1,2) = 10;
  GeoReference georef( Datum(), affine );

  TransformRef txref(ResampleTransform(1,1));
  int pyramid_level;
  platemanager->transform_image( georef, image_ref, txref, pyramid_level );

  EXPECT_EQ( 1, pyramid_level );
  EXPECT_EQ( 512, image_ref.cols() );
  EXPECT_EQ( 512, image_ref.rows() );

  std::list<TileInfo> tiles;
  platemanager->affected_tiles( bounding_box(image), txref, 256,
                                pyramid_level, tiles );

  // Checking the validity of the tiles requested for rasterization.
  EXPECT_EQ( 4, tiles.size() );
  BOOST_FOREACH( TileInfo const& tile, tiles ) {
    EXPECT_LT( tile.i, 2 );
    EXPECT_LT( tile.j, 2 );
    EXPECT_LE( tile.bbox.max().x(), 512 );
    EXPECT_LE( tile.bbox.max().y(), 512 );
    EXPECT_GE( tile.i, 0 );
    EXPECT_GE( tile.j, 0 );
    EXPECT_GE( tile.bbox.min().x(), 0 );
    EXPECT_GE( tile.bbox.min().y(), 0 );

    ImageView<PixelT > raster = crop(image_ref,tile.bbox);
    EXPECT_FALSE( is_transparent( raster ) );
  }
}

TEST_F( PlateCarreeExposedTest, EastWrap ) {
  Matrix3x3 affine = identity_matrix<3>();
  affine(1,1) = -1;
  affine(0,2) = 170;
  affine(1,2) = 10;
  GeoReference georef( Datum(), affine );

  TransformRef txref(ResampleTransform(1,1));
  int pyramid_level;
  platemanager->transform_image( georef, image_ref, txref, pyramid_level );

  EXPECT_EQ( 1, pyramid_level );
  EXPECT_EQ( 512, image_ref.cols() );
  EXPECT_EQ( 512, image_ref.rows() );

  std::list<TileInfo> tiles;
  platemanager->affected_tiles( bounding_box(image), txref, 256,
                                pyramid_level, tiles );

  // Checking the validity of the tiles requested for rasterization.
  EXPECT_EQ( 4, tiles.size() );
  BOOST_FOREACH( TileInfo const& tile, tiles ) {
    EXPECT_LT( tile.i, 2 );
    EXPECT_LT( tile.j, 2 );
    EXPECT_LE( tile.bbox.max().x(), 512 );
    EXPECT_LE( tile.bbox.max().y(), 512 );
    EXPECT_GE( tile.i, 0 );
    EXPECT_GE( tile.j, 0 );
    EXPECT_GE( tile.bbox.min().x(), 0 );
    EXPECT_GE( tile.bbox.min().y(), 0 );

    ImageView<PixelT > raster = crop(image_ref,tile.bbox);
    EXPECT_FALSE( is_transparent( raster ) );
  }
}

TEST_F( PlateCarreeExposedTest, SingleTile ) {
  Matrix3x3 affine = identity_matrix<3>();
  affine(1,1) = -1;
  affine(0,2) = 0;
  affine(1,2) = 20;
  GeoReference georef( Datum(), affine );

  TransformRef txref(ResampleTransform(1,1));
  int pyramid_level;
  platemanager->transform_image( georef, image_ref, txref, pyramid_level );

  EXPECT_EQ( 1, pyramid_level );

  std::list<TileInfo> tiles;
  platemanager->affected_tiles( bounding_box(image), txref, 256,
                                pyramid_level, tiles );

  // Checking the validity of the tiles requested for rasterization.
  ASSERT_EQ( 1, tiles.size() );
  EXPECT_EQ( 1, tiles.begin()->i );
  EXPECT_EQ( 0, tiles.begin()->j );
  EXPECT_VW_EQ( Vector2i(256,0), tiles.begin()->bbox.min() );
  EXPECT_VW_EQ( Vector2i(512,256), tiles.begin()->bbox.max() );

  ImageView<PixelT > raster = crop(image_ref,
                                   tiles.begin()->bbox);
  EXPECT_FALSE( is_transparent( raster ) );
}
