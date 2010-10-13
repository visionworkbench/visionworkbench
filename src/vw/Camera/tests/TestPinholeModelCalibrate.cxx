// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <gtest/gtest.h>

#include <vw/Camera/PinholeModelCalibrate.h>
#include <test/Helpers.h>

#include <cstdlib>
#include <ctime>
#include <functional>

using namespace vw;
using namespace vw::camera;

// Utility Functions
// -----------------------

void LessThanEqualDelta( double const& a,
                         double const& b, double delta=1e-4 ) {
  if ( a >= b )
    EXPECT_NEAR(a,b,delta);
}

double mean_error( PinholeModel const& camera,
                   std::vector<Vector2> const& pixels,
                   std::vector<Vector3> const& points ) {
  VW_ASSERT(pixels.size() > 0, LogicErr() << "mean_error: refusing to divide by zero");
  double mean = 0;
  for (uint32 i = 0; i < pixels.size(); i++)
    mean += norm_2(pixels[i] - camera.point_to_pixel(points[i]));
  return mean / pixels.size();
}

double mean_error( PinholeModel const& camera,
                   std::vector<Vector2> const& pixels,
                   std::vector<Vector3> const& points,
                   std::vector<int> const& indices ) {
  VW_ASSERT(indices.size() > 0, LogicErr() << "mean_error: refusing to divide by zero");
  double mean = 0;
  for (uint32 j = 0; j < indices.size(); j++)
    mean += norm_2(pixels[indices[j]] -
                   camera.point_to_pixel(points[indices[j]]));
  return mean / indices.size();
}

double mean_sqr_error(PinholeModel const& camera,
                      std::vector<Vector2> const& pixels,
                      std::vector<Vector3> const& points) {
  VW_ASSERT(pixels.size() > 0, LogicErr() << "mean_sqr_error: refusing to divide by zero");
  double mean = 0;
  for (uint32 i = 0; i < pixels.size(); i++)
    mean += norm_2_sqr(pixels[i] -
                       camera.point_to_pixel(points[i]));
  return mean / pixels.size();
}

double mean_sqr_error( PinholeModel const& camera,
                       std::vector<Vector2> const& pixels,
                       std::vector<Vector3> const& points,
                       std::vector<int> const& indices) {
  VW_ASSERT(indices.size() > 0, LogicErr() << "mean_sqr_error: refusing to divide by zero");
  double mean = 0;
  for (uint32 j = 0; j < indices.size(); j++)
    mean += norm_2_sqr(pixels[indices[j]] -
                       camera.point_to_pixel(points[indices[j]]));
  return mean / indices.size();
}

class RepeatableSequence {
  uint32 seed;
public:
  RepeatableSequence( uint32 init = 10 ) : seed(init) {}
  double operator()( double max = 100 ) {
    seed = seed*1103515245 + 12345;
    return max*static_cast<double>(seed) / std::numeric_limits<uint32>::max();
  }
};

// Actual Tests
// -----------------------

TEST( PinholeModelCalibrate, SerializeFour ) {
  for (int i = 0; i < 10; i++) {
    Vector3 cc(rand(), rand(), rand());
    Vector3 rv(rand(), rand(), rand());
    Matrix3x3 rm( math::axis_angle_to_matrix(rv));
    Vector2 focal( rand(), rand() );
    Vector2 offst( rand(), rand() );
    Vector3 u(1, 0, 0);
    Vector3 v(0, 1, 0);
    Vector3 w(0, 0, 1);
    Vector4 tsaiv(Vector4(rand(), rand(), rand(), rand()));
    TsaiLensDistortion tsai(tsaiv);
    PinholeModel m(cc, rm,
                   focal[0], focal[1],
                   offst[0], offst[1],
                   u, v, w, tsai);

    Vector<double> serial(serialize_pinholemodel<PinholeModelSerializeIntrinsic, PinholeModelSerializeRotation, PinholeModelSerializeTranslation, PinholeModelSerializeTSAI>(m));

    // check if correctly serialized
    EXPECT_VECTOR_NEAR( subvector(serial,0,2), focal, 1e-10 );
    EXPECT_VECTOR_NEAR( subvector(serial,2,2), offst, 1e-10 );
    EXPECT_VECTOR_NEAR( subvector(serial,4,3),
                        m.camera_pose().axis_angle(), 1e-10 );
    EXPECT_VECTOR_NEAR( subvector(serial,7,3),
                        m.camera_center(), 1e-10 );
    EXPECT_VECTOR_NEAR( subvector(serial,10,4),
                        tsaiv, 1e-10 );

    // test deserialization
    PinholeModel d;
    d.set_coordinate_frame(u, v, w);
    deserialize_pinholemodel<PinholeModelSerializeIntrinsic, PinholeModelSerializeRotation, PinholeModelSerializeTranslation, PinholeModelSerializeTSAI>(d, serial);
    EXPECT_VECTOR_NEAR( d.focal_length(), focal, 1e-10 );
    EXPECT_VECTOR_NEAR( d.point_offset(), offst, 1e-10 );
    EXPECT_VECTOR_NEAR( d.camera_pose().axis_angle(),
                        m.camera_pose().axis_angle(), 1e-14 );
    EXPECT_VECTOR_NEAR( d.camera_center(),
                        m.camera_center(), 1e-8 );
    Vector4 tsai_d( dynamic_cast<const TsaiLensDistortion*>(d.lens_distortion())->distortion_parameters() );
    EXPECT_VECTOR_NEAR( tsai_d, tsaiv, 1e-8 );
  }
}

TEST( PinholeModelCalibrate, SerializeThree ) {
  for (int i = 0; i < 10; i++) {
    Vector3 cc(rand(), rand(), rand());
    Vector3 rv(rand(), rand(), rand());
    Matrix3x3 rm(math::axis_angle_to_matrix(rv));
    Vector2 focal( rand(), rand() );
    Vector2 offst( rand(), rand() );
    Vector3 u(1, 0, 0);
    Vector3 v(0, 1, 0);
    Vector3 w(0, 0, 1);
    Vector4 tsaiv(Vector4(rand(), rand(), rand(), rand()));
    TsaiLensDistortion tsai(tsaiv);
    PinholeModel m(cc, rm,
                   focal[0], focal[1],
                   offst[0], offst[1],
                   u, v, w, tsai);

    Vector<double> serial(serialize_pinholemodel<PinholeModelSerializeIntrinsic, PinholeModelSerializeRotation, PinholeModelSerializeTranslation>(m));

    // check if correctly serialized
    EXPECT_VECTOR_NEAR( subvector(serial,0,2), focal, 1e-10 );
    EXPECT_VECTOR_NEAR( subvector(serial,2,2), offst, 1e-10 );
    EXPECT_VECTOR_NEAR( subvector(serial,4,3),
                        m.camera_pose().axis_angle(), 1e-10 );
    EXPECT_VECTOR_NEAR( subvector(serial,7,3),
                        m.camera_center(), 1e-10 );

    // test deserialization
    PinholeModel d;
    d.set_coordinate_frame(u, v, w);
    deserialize_pinholemodel<PinholeModelSerializeIntrinsic, PinholeModelSerializeRotation, PinholeModelSerializeTranslation>(d, serial);
    EXPECT_VECTOR_NEAR( d.focal_length(), focal, 1e-10 );
    EXPECT_VECTOR_NEAR( d.point_offset(), offst, 1e-10 );
    EXPECT_VECTOR_NEAR( d.camera_pose().axis_angle(),
                        m.camera_pose().axis_angle(), 1e-14 );
    EXPECT_VECTOR_NEAR( d.camera_center(),
                        m.camera_center(), 1e-8 );
  }
}

TEST( PinholeModelCalibrate, SerializeTwo ) {
  for (int i = 0; i < 10; i++) {
    Vector3 cc(rand(), rand(), rand());
    Vector3 rv(rand(), rand(), rand());
    Matrix3x3 rm(math::axis_angle_to_matrix(rv));
    Vector2 focal( rand(), rand() );
    Vector2 offst( rand(), rand() );
    Vector3 u(1, 0, 0);
    Vector3 v(0, 1, 0);
    Vector3 w(0, 0, 1);
    Vector4 tsaiv(Vector4(rand(), rand(), rand(), rand()));
    TsaiLensDistortion tsai(tsaiv);
    PinholeModel m(cc, rm,
                   focal[0], focal[1],
                   offst[0], offst[1],
                   u, v, w, tsai);

    Vector<double> serial(serialize_pinholemodel<PinholeModelSerializeIntrinsic, PinholeModelSerializeRotation>(m));

    // check if correctly serialized
    EXPECT_VECTOR_NEAR( subvector(serial,0,2), focal, 1e-10 );
    EXPECT_VECTOR_NEAR( subvector(serial,2,2), offst, 1e-10 );
    EXPECT_VECTOR_NEAR( subvector(serial,4,3),
                        m.camera_pose().axis_angle(), 1e-10 );

    // test deserialization
    PinholeModel d;
    d.set_coordinate_frame(u, v, w);
    deserialize_pinholemodel<PinholeModelSerializeIntrinsic, PinholeModelSerializeRotation>(d, serial);
    EXPECT_VECTOR_NEAR( d.focal_length(), focal, 1e-10 );
    EXPECT_VECTOR_NEAR( d.point_offset(), offst, 1e-10 );
    EXPECT_VECTOR_NEAR( d.camera_pose().axis_angle(),
                        m.camera_pose().axis_angle(), 1e-14 );
  }
}

TEST( PinholeModelCalibrate, SerializeOne ) {
  for (int i = 0; i < 10; i++) {
    Vector3 cc(rand(), rand(), rand());
    Vector3 rv(rand(), rand(), rand());
    Matrix3x3 rm(vw::math::axis_angle_to_matrix(rv));
    Vector2 focal( rand(), rand() );
    Vector2 offst( rand(), rand() );
    Vector3 u(1, 0, 0);
    Vector3 v(0, 1, 0);
    Vector3 w(0, 0, 1);
    Vector4 tsaiv(vw::Vector4(rand(), rand(), rand(), rand()));
    TsaiLensDistortion tsai(tsaiv);
    PinholeModel m(cc, rm,
                   focal[0], focal[1],
                   offst[0], offst[1],
                   u, v, w, tsai);

    vw::Vector<double> serial(serialize_pinholemodel<PinholeModelSerializeIntrinsic>(m));

    // check if correctly serialized
    EXPECT_VECTOR_NEAR( subvector(serial,0,2), focal, 1e-10 );
    EXPECT_VECTOR_NEAR( subvector(serial,2,2), offst, 1e-10 );

    // test deserialization
    PinholeModel d;
    d.set_coordinate_frame(u, v, w);
    deserialize_pinholemodel<PinholeModelSerializeIntrinsic>(d, serial);
    EXPECT_VECTOR_NEAR( d.focal_length(), focal, 1e-10 );
    EXPECT_VECTOR_NEAR( d.point_offset(), offst, 1e-10 );
  }
}

TEST( PinholeModelCalibrate, Calibrate ) {
  RepeatableSequence seq;

  Vector3 cc(1, 1, 1);
  Vector3 rv(1, 1, 1);
  Matrix3x3 rm(vw::math::axis_angle_to_matrix(rv));
  Vector2 focal( 1000, 500 );
  Vector2 offst( 200,  100 );
  Vector3 u(1, 0, 0);
  Vector3 v(0, 1, 0);
  Vector3 w(0, 0, 1);
  Vector4 tsaiv(1e-2, 1e-2, 1e-2, 1e-2);
  TsaiLensDistortion tsai(tsaiv);
  PinholeModel m(cc, rm,
                 focal[0], focal[1],
                 offst[0], offst[1],
                 u, v, w, tsai);

  for (int ni = 0; ni < 5; ni++) {
    std::vector<Vector3> points; // in 3D space
    std::vector<Vector2> pixels; // from projection through m
    int n = 30;
    for (int i = 0; i < n; i++) {
      Vector3 p( seq(), seq(), seq() );
      points.push_back(p);
      Vector2 noise( seq(10), seq(10) ); // add some noise!
      pixels.push_back(m.point_to_pixel(p) + noise);
    }

    double mean = mean_sqr_error(m, pixels, points);

    // see if the optimizer improves the results; the mean error
    // should decrease as the number of variables the optimizer gets
    // to play with increases. These tests can fail on real data
    // depending on your srand implementation without anything being
    // actually broken - which is why a "fake" srand/rand function
    // pair local to this class was implemented to generate pseudo
    // pseudo random values (i.e. sequence that is always the same)
    // for point positions and noise - thus the test should be
    // completely repeatable
    {
      PinholeModel c(m);
      pinholemodel_calibrate<PinholeModelSerializeIntrinsic>(c, pixels, points, 1000);
      double new_mean = mean_sqr_error(c, pixels, points);
      EXPECT_LT( new_mean, mean );
    }

    {
      PinholeModel c(m);
      pinholemodel_calibrate<PinholeModelSerializeTSAI>(c, pixels, points, 1000);
      double new_mean = mean_sqr_error(c, pixels, points);
      EXPECT_LT( new_mean, mean );
    }

    {
      PinholeModel c(m);
      pinholemodel_calibrate<PinholeModelSerializeRotation>(c, pixels, points, 1000);
      double new_mean = mean_sqr_error(c, pixels, points);
      EXPECT_LT( new_mean, mean );
    }

    PinholeModel c(m);
    {
      pinholemodel_calibrate<PinholeModelSerializeTranslation>(c, pixels, points, 1000);
      double new_mean = mean_sqr_error(c, pixels, points);
      EXPECT_LT( new_mean, mean );
      mean = new_mean;
    }

    {
      pinholemodel_calibrate<PinholeModelSerializeTranslation, PinholeModelSerializeRotation>(c, pixels, points, 1000);
      double new_mean = mean_sqr_error(c, pixels, points);
      LessThanEqualDelta(new_mean, mean);
      mean = new_mean;
    }

    {
      pinholemodel_calibrate<PinholeModelSerializeTranslation, PinholeModelSerializeRotation, PinholeModelSerializeIntrinsic>(c, pixels, points, 1000);
      double new_mean = mean_sqr_error(c, pixels, points);
      LessThanEqualDelta(new_mean, mean);
      mean = new_mean;
    }

    {
      pinholemodel_calibrate<PinholeModelSerializeTranslation, PinholeModelSerializeRotation, PinholeModelSerializeIntrinsic, PinholeModelSerializeTSAI>(c, pixels, points, 1000);
      double new_mean = mean_sqr_error(c, pixels, points);
      LessThanEqualDelta(new_mean, mean);
      mean = new_mean;
    }
  }
}

TEST( PinholeModelCalibrate, Ransac ) {
  RepeatableSequence seq;

  Vector3 cc(1, 1, 1);
  Vector3 rv(1, 1, 1);
  Matrix3x3 rm(vw::math::axis_angle_to_matrix(rv));
  Vector2 focal( 1000, 500 );
  Vector2 offst( 200,  100 );
  Vector3 u(1, 0, 0);
  Vector3 v(0, 1, 0);
  Vector3 w(0, 0, 1);
  Vector4 tsaiv(1e-2, 1e-2, 1e-2, 1e-2);
  TsaiLensDistortion tsai(tsaiv);
  PinholeModel m(cc, rm,
                 focal[0], focal[1],
                 offst[0], offst[1],
                 u, v, w, tsai);

  std::vector<Vector3> points; // in 3 space
  std::vector<Vector2> pixels; // from projection through m
  int n = 30;
  for (int i = 0; i < n; i++) {
    Vector3 p(seq(), seq(), seq());
    points.push_back(p);
    Vector2 noise(seq(20), seq(20)); // add some noise!
    pixels.push_back(m.point_to_pixel(p) + noise);
  }

  //double mean = mean_error(m, pixels, points);

  // these tests verify mainly that RANSAC actually does respect the
  // inlier_threshold value passed in (The maximum error in the image
  // plane for a resulting camera model can be equal to
  // inlier_threshold) and that the number of inliers is not "too
  // small" play with these parameters to make the test more or less
  // stringent and/or fast
  double inlier_threshold = 15;
  const unsigned ransac_inlier_threshold = 10; // how many inliers to we require
  const unsigned ransac_iter = 20; // number of ransac iterations
  const unsigned lm_iter = 5; // number of levenberg marquardt iterations at every ransac iteration
  {
    PinholeModel c(m);
    std::vector<int> inliers( pinholemodel_calibrate_ransac<PinholeModelSerializeIntrinsic>(c, pixels, points, inlier_threshold, ransac_iter, lm_iter) );
    EXPECT_LT(ransac_inlier_threshold, inliers.size()); // only critical if this fails repeatably, as ransac classification depends on random numbers
    EXPECT_LT(mean_error(c, pixels, points, inliers), inlier_threshold);
  }

  {
    PinholeModel c(m);
    std::vector<int> inliers( pinholemodel_calibrate_ransac<PinholeModelSerializeTSAI>(c, pixels, points, inlier_threshold, ransac_iter, lm_iter) );
    EXPECT_LT(ransac_inlier_threshold, inliers.size()); // only critical if this fails repeatably, as ransac classification depends on random numbers
    EXPECT_LT(mean_error(c, pixels, points, inliers), inlier_threshold);
  }

  {
    PinholeModel c(m);
    std::vector<int> inliers( pinholemodel_calibrate_ransac<PinholeModelSerializeRotation>(c, pixels, points, inlier_threshold, ransac_iter, lm_iter) );
    EXPECT_LT(ransac_inlier_threshold, inliers.size()); // only critical if this fails repeatably, as ransac classification depends on random numbers
    EXPECT_LT(mean_error(c, pixels, points, inliers), inlier_threshold);
  }

  {
    PinholeModel c(m);
    std::vector<int> inliers( pinholemodel_calibrate_ransac<PinholeModelSerializeTranslation>(c, pixels, points, inlier_threshold, ransac_iter, lm_iter) );
    EXPECT_LT(ransac_inlier_threshold, inliers.size()); // only critical if this fails repeatably, as ransac classification depends on random numbers
    EXPECT_LT(mean_error(c, pixels, points, inliers), inlier_threshold);
  }

  {
    PinholeModel c(m);
    std::vector<int> inliers( pinholemodel_calibrate_ransac<PinholeModelSerializeTranslation, PinholeModelSerializeRotation>(c, pixels, points, inlier_threshold, ransac_iter, lm_iter) );
    EXPECT_LT(ransac_inlier_threshold, inliers.size()); // only critical if this fails repeatably, as ransac classification depends on random numbers
    EXPECT_LT(mean_error(c, pixels, points, inliers), inlier_threshold);
  }

  {
    PinholeModel c(m);
    std::vector<int> inliers( pinholemodel_calibrate_ransac<PinholeModelSerializeTranslation, PinholeModelSerializeRotation, PinholeModelSerializeIntrinsic>(c, pixels, points, inlier_threshold, ransac_iter, lm_iter) );
    EXPECT_LT(ransac_inlier_threshold, inliers.size()); // only critical if this fails repeatably, as ransac classification depends on random numbers
    EXPECT_LT(mean_error(c, pixels, points, inliers), inlier_threshold);
  }

  {
    PinholeModel c(m);
    std::vector<int> inliers( pinholemodel_calibrate_ransac<PinholeModelSerializeTranslation, PinholeModelSerializeRotation, PinholeModelSerializeIntrinsic, PinholeModelSerializeTSAI>(c, pixels, points, inlier_threshold, ransac_iter, lm_iter) );
    EXPECT_LT(ransac_inlier_threshold, inliers.size()); // only critical if this fails repeatably, as ransac classification depends on random numbers
    EXPECT_LT(mean_error(c, pixels, points, inliers), inlier_threshold);
  }
}
