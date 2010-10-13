// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// TestGeoReference.h
#include <gtest/gtest.h>
#include <test/Helpers.h>

#include <vw/Cartography/ToastTransform.h>

using namespace vw;
using namespace vw::cartography;
using namespace vw::test;

static const int32 toast_resolution = 4*255+1;
static const int32 lonlat_resolution = 4*256;

class ToastTransformTest : public ::testing::Test {
protected:

  ToastTransformTest() : txform(GeoReference(),1) {}

  virtual void SetUp() {
    GeoReference georef;
    georef.set_pixel_interpretation(GeoReference::PixelAsPoint);
    Matrix3x3 M;
    M(0,0) = 360.0 / lonlat_resolution;
    M(0,2) = -180;
    M(1,1) = -M(0,0);
    M(1,2) = 90;
    M(2,2) = 1;
    georef.set_transform(M);
    basicref = georef;
    txform = ToastTransform(georef, toast_resolution);
    toast_res_vec = Vector2(toast_resolution,
                            toast_resolution);
    lonlat_res_vec = Vector2(lonlat_resolution,
                             lonlat_resolution);
  }

  ToastTransform txform;
  GeoReference basicref;
  Vector2 point;

  Vector2 toast_res_vec;
  Vector2 lonlat_res_vec;
};

TEST_F( ToastTransformTest, BasicForward ) {

  // North pole
  point = txform.forward(Vector2(0,0));
  EXPECT_VECTOR_NEAR( point, (toast_res_vec-Vector2(1,1))/2, 1e-5 );

  // North pole
  point = txform.forward(Vector2(lonlat_resolution/4,0));
  EXPECT_VECTOR_NEAR( point, (toast_res_vec-Vector2(1,1))/2, 1e-5 );

  // North pole
  point = txform.forward(Vector2(lonlat_resolution/2,0));
  EXPECT_VECTOR_NEAR( point, (toast_res_vec-Vector2(1,1))/2, 1e-5 );

  // North pole
  point = txform.forward(Vector2(3*lonlat_resolution/4,0));
  EXPECT_VECTOR_NEAR( point, (toast_res_vec-Vector2(1,1))/2, 1e-5 );

  // North pole
  point = txform.forward(Vector2(lonlat_resolution,0));
  EXPECT_VECTOR_NEAR( point, (toast_res_vec-Vector2(1,1))/2, 1e-5 );

  // (-180,45)
  point = txform.forward(Vector2(0,lonlat_resolution/8));
  EXPECT_VECTOR_NEAR( point, (toast_resolution-1)*Vector2(3./4,1./2), 1e-5 );

  // (-90,45)
  point = txform.forward(Vector2(lonlat_resolution/4,lonlat_resolution/8));
  EXPECT_VECTOR_NEAR( point, (toast_resolution-1)*Vector2(1./2,1./4), 1e-5 );

  // (0,45)
  point = txform.forward(Vector2(lonlat_resolution/2,lonlat_resolution/8));
  EXPECT_VECTOR_NEAR( point, (toast_resolution-1)*Vector2(1./4,1./2), 1e-5 );

  // (90,45)
  point = txform.forward(Vector2(3*lonlat_resolution/4,lonlat_resolution/8));
  EXPECT_VECTOR_NEAR( point, (toast_resolution-1)*Vector2(1./2,3./4), 1e-5 );

  // (180,45)
  point = txform.forward(Vector2(0,lonlat_resolution/8));
  EXPECT_VECTOR_NEAR( point, (toast_resolution-1)*Vector2(3./4,1./2), 1e-5 );

  // (-180,0)
  point = txform.forward(Vector2(0,lonlat_resolution/4));
  EXPECT_VECTOR_NEAR( point, (toast_resolution-1)*Vector2(1,1./2), 1e-5 );

  // (-135,0)
  point = txform.forward(Vector2(lonlat_resolution/8,lonlat_resolution/4));
  EXPECT_VECTOR_NEAR( point, (toast_resolution-1)*Vector2(3./4,1./4), 1e-5 );

  // (-90,0)
  point = txform.forward(Vector2(lonlat_resolution/4,lonlat_resolution/4));
  EXPECT_VECTOR_NEAR( point, (toast_resolution-1)*Vector2(1./2,0), 1e-5 );

  // (-45,0)
  point = txform.forward(Vector2(3*lonlat_resolution/8,lonlat_resolution/4));
  EXPECT_VECTOR_NEAR( point, (toast_resolution-1)*Vector2(1./4,1./4), 1e-5 );

  // (0,0)
  point = txform.forward(Vector2(lonlat_resolution/2,lonlat_resolution/4));
  EXPECT_VECTOR_NEAR( point, (toast_resolution-1)*Vector2(0,1./2), 1e-5 );

  // (45,0)
  point = txform.forward(Vector2(5*lonlat_resolution/8,lonlat_resolution/4));
  EXPECT_VECTOR_NEAR( point, (toast_resolution-1)*Vector2(1./4,3./4), 1e-5 );

  // (90,0)
  point = txform.forward(Vector2(3*lonlat_resolution/4,lonlat_resolution/4));
  EXPECT_VECTOR_NEAR( point, (toast_resolution-1)*Vector2(1./2,1), 1e-5 );

  // (135,0)
  point = txform.forward(Vector2(7*lonlat_resolution/8,lonlat_resolution/4));
  EXPECT_VECTOR_NEAR( point, (toast_resolution-1)*Vector2(3./4,3./4), 1e-5 );

  // (180,0)
  point = txform.forward(Vector2(0,lonlat_resolution/4));
  EXPECT_VECTOR_NEAR( point, (toast_resolution-1)*Vector2(1,1./2), 1e-5 );

  // (-180,-45)
  point = txform.forward(Vector2(0,3*lonlat_resolution/8));
  EXPECT_VECTOR_NEAR( point, (toast_resolution-1)*Vector2(1,3./4), 1e-5 );

  // (-90,-45)
  point = txform.forward(Vector2(lonlat_resolution/4,3*lonlat_resolution/8));
  EXPECT_VECTOR_NEAR( point, (toast_resolution-1)*Vector2(3./4,0), 1e-5 );

  // (0,-45)
  point = txform.forward(Vector2(lonlat_resolution/2,3*lonlat_resolution/8));
  EXPECT_VECTOR_NEAR( point, (toast_resolution-1)*Vector2(0,3./4), 1e-5 );

  // (90,-45)
  point = txform.forward(Vector2(3*lonlat_resolution/4,3*lonlat_resolution/8));
  EXPECT_VECTOR_NEAR( point, (toast_resolution-1)*Vector2(3./4,1), 1e-5 );

  // (180,-45)
  point = txform.forward(Vector2(0,3*lonlat_resolution/8));
  EXPECT_VECTOR_NEAR( point, (toast_resolution-1)*Vector2(1,3./4), 1e-5 );

  // South pole
  point = txform.forward(Vector2(0,lonlat_resolution/2));
  EXPECT_VECTOR_NEAR( point, toast_res_vec-Vector2(1,1), 1e-5 );

  // South pole
  point = txform.forward(Vector2(lonlat_resolution/4,lonlat_resolution/2));
  EXPECT_VECTOR_NEAR( point, (toast_resolution-1)*Vector2(1,0), 1e-5 );

  // South pole
  point = txform.forward(Vector2(lonlat_resolution/2,lonlat_resolution/2));
  EXPECT_VECTOR_NEAR( point, (toast_resolution-1)*Vector2(0,1), 1e-5 );

  // South pole
  point = txform.forward(Vector2(3*lonlat_resolution/4,lonlat_resolution/2));
  EXPECT_VECTOR_NEAR( point, toast_res_vec-Vector2(1,1), 1e-5 );

  // South pole
  point = txform.forward(Vector2(lonlat_resolution,lonlat_resolution/2));
  EXPECT_VECTOR_NEAR( point, toast_res_vec-Vector2(1,1), 1e-5 );
}

TEST_F( ToastTransformTest, BasicReverse ) {

  // Top left: (*,-90)
  point = txform.reverse(Vector2(0,0));
  EXPECT_NEAR( point.y(), lonlat_resolution/2, 1e-5 );

  // Top left edge: (-90,-45)
  point = txform.reverse(Vector2((toast_resolution-1)/4,0));
  EXPECT_VECTOR_NEAR( point, lonlat_resolution*Vector2(1./4,3./8), 1e-5 );

  // Top center: (-90,0)
  point = txform.reverse(Vector2((toast_resolution-1)/2,0));
  EXPECT_VECTOR_NEAR( point, lonlat_resolution*Vector2(1./4,1./4), 1e-5 );

  // Top right edge: (-90,-45)
  point = txform.reverse(Vector2(3*(toast_resolution-1)/4,0));
  EXPECT_VECTOR_NEAR( point, lonlat_resolution*Vector2(1./4,3./8), 1e-5 );

  // Top right: (*,-90)
  point = txform.reverse(Vector2(toast_resolution-1,0));
  EXPECT_NEAR( point.y(), lonlat_resolution/2, 1e-5 );

  // Upper left edge: (0,-45)
  point = txform.reverse(Vector2(0,(toast_resolution-1)/4));
  EXPECT_VECTOR_NEAR( point, lonlat_resolution*Vector2(1./2,3./8), 1e-5 );

  // Upper left middle: (-45,0)
  point = txform.reverse(Vector2((toast_resolution-1)/4,(toast_resolution-1)/4));
  EXPECT_VECTOR_NEAR( point, lonlat_resolution*Vector2(3./8,1./4), 1e-5 );

  // Upper middle: (-90,45)
  point = txform.reverse(Vector2((toast_resolution-1)/2,(toast_resolution-1)/4));
  EXPECT_VECTOR_NEAR( point, lonlat_resolution*Vector2(1./4,1./8), 1e-5 );

  // Upper right middle: (-135,0)
  point = txform.reverse(Vector2(3*(toast_resolution-1)/4,(toast_resolution-1)/4));
  EXPECT_VECTOR_NEAR( point, lonlat_resolution*Vector2(1./8,1./4), 1e-5 );

  // Upper right edge: (-180,-45)
  point = txform.reverse(Vector2(toast_resolution-1,(toast_resolution-1)/4));
  EXPECT_VECTOR_NEAR( point, lonlat_resolution*Vector2(0,3./8), 1e-5 );

  // Left center: (0,0)
  point = txform.reverse(Vector2(0,(toast_resolution-1)/2));
  EXPECT_VECTOR_NEAR( point, lonlat_resolution*Vector2(1./2,1./4), 1e-5 );

  // Left middle: (0,45)
  point = txform.reverse(Vector2((toast_resolution-1)/4,(toast_resolution-1)/2));
  EXPECT_VECTOR_NEAR( point, lonlat_resolution*Vector2(1./2,1./8), 1e-5 );

  // Center: (*, 90)
  point = txform.reverse(Vector2((toast_resolution-1)/2,(toast_resolution-1)/2));
  EXPECT_NEAR( point.y(), 0, 1e-5 );

  // Right middle: (-180,45)
  point = txform.reverse(Vector2(3*(toast_resolution-1)/4,(toast_resolution-1)/2));
  EXPECT_VECTOR_NEAR( point, lonlat_resolution*Vector2(0,1./8), 1e-5 );

  // Right center: (-180,0)
  point = txform.reverse(Vector2(toast_resolution-1,(toast_resolution-1)/2));
  EXPECT_VECTOR_NEAR( point, lonlat_resolution*Vector2(0,1./4), 1e-5 );

  // Lower left edge: (0,-45)
  point = txform.reverse(Vector2(0,3*(toast_resolution-1)/4));
  EXPECT_VECTOR_NEAR( point, lonlat_resolution*Vector2(1./2,3./8), 1e-5 );

  // Lower left middle: (45,0)
  point = txform.reverse(Vector2((toast_resolution-1)/4,3*(toast_resolution-1)/4));
  EXPECT_VECTOR_NEAR( point, lonlat_resolution*Vector2(5./8,1./4), 1e-5 );

  // Lower middle: (90,45)
  point = txform.reverse(Vector2((toast_resolution-1)/2,3*(toast_resolution-1)/4));
  EXPECT_VECTOR_NEAR( point, lonlat_resolution*Vector2(3./4,1./8), 1e-5 );

  // Lower right middle: (135,0)
  point = txform.reverse(Vector2(3*(toast_resolution-1)/4,3*(toast_resolution-1)/4));
  EXPECT_VECTOR_NEAR( point, lonlat_resolution*Vector2(7./8,1./4), 1e-5 );

  // Lower right edge: (-180,-45)
  point = txform.reverse(Vector2(toast_resolution-1,3*(toast_resolution-1)/4));
  EXPECT_VECTOR_NEAR( point, lonlat_resolution*Vector2(1,3./8), 1e-5 );

  // Bottom left: (*,-90)
  point = txform.reverse(Vector2(0,toast_resolution-1));
  EXPECT_NEAR( point.y(), lonlat_resolution/2, 1e-5 );

  // Bottom left edge: (90,-45)
  point = txform.reverse(Vector2((toast_resolution-1)/4,toast_resolution-1));
  EXPECT_VECTOR_NEAR( point, lonlat_resolution*Vector2(3./4,3./8), 1e-5 );

  // Bottom center: (90,0)
  point = txform.reverse(Vector2((toast_resolution-1)/2,toast_resolution-1));
  EXPECT_VECTOR_NEAR( point, lonlat_resolution*Vector2(3./4,1./4), 1e-5 );

  // Bottom right edge: (90,-45)
  point = txform.reverse(Vector2(3*(toast_resolution-1)/4,toast_resolution-1));
  EXPECT_VECTOR_NEAR( point, lonlat_resolution*Vector2(3./4,3./8), 1e-5 );

  // Bottom right: (*,-90)
  point = txform.reverse(Vector2(toast_resolution-1,toast_resolution-1));
  EXPECT_NEAR( point.y(), lonlat_resolution/2, 1e-5 );
}

TEST_F( ToastTransformTest, RandomForwardReverse ) {
  Vector2 lonlat;

  for( int i=0; i<1000; ++i ) {
    Vector2 obj( 1+(double)(lonlat_resolution-1)*rand()/RAND_MAX,
                 1+(double)(lonlat_resolution/2-1)*rand()/RAND_MAX );
    point = txform.forward(obj);
    lonlat = txform.reverse(point);
    EXPECT_VECTOR_NEAR( lonlat, obj, 1e-4 );
  }
}

TEST_F( ToastTransformTest, BasicBBoxCheck ) {
  // Just make sure that forward_bbox and reverse_bbox behave

  // form a bbox with the point in question at the center (PixelAsPoint, so
  // a 1x1 box at the point in question should represent the point).

  // (-90,45)
  BBox2i in_box(lonlat_resolution/4, lonlat_resolution/8,1,1);
  BBox2i out_box;
  EXPECT_NO_THROW( out_box = txform.forward_bbox(in_box); );
  EXPECT_VECTOR_NEAR(out_box.min(), (toast_resolution-1)*Vector2(1./2,1./4), 1e-5);

  // Lower right edge: (-180,-45)
  in_box = BBox2(toast_resolution-1,3*(toast_resolution-1)/4,1,1);
  EXPECT_NO_THROW( out_box = txform.reverse_bbox(in_box); );
  EXPECT_VECTOR_NEAR(out_box.min(), lonlat_resolution*Vector2(1,3./8), 1e-5 );
}

TEST_F( ToastTransformTest, StereographicBBox ) {
  // forward_bbox needs to override for the polar case. Let's make sure the
  // polar case works for stereographic, where the poles aren't defined.
  Datum datum("unknown", "unnamed", "Greenwich", 3376200, 3376200, 0);

  basicref.set_stereographic(90, 0, 1, 0, 0);
  basicref.set_datum(datum);
  ToastTransform tx(basicref, toast_resolution);

  // (-90,45)
  BBox2i in_box(lonlat_resolution/4, lonlat_resolution/8,1,1);
  BBox2i out_box;

  // We mostly only care that it doesn't throw an exception, but also do a
  // basic check to make sure it ended up inside the toast bbox.
  EXPECT_NO_THROW( out_box = tx.forward_bbox(in_box); );

  BBox2i global(0,0,toast_resolution,toast_resolution);
  EXPECT_TRUE( global.contains(out_box) );
}
