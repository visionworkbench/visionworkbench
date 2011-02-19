// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <gtest/gtest.h>

#include <boost/random.hpp>

#include <vw/Camera/PinholeModel.h>
#include <vw/Camera/CameraGeometry.h>
#include <vw/Math/EulerAngles.h>
#include <test/Helpers.h>
#include <vw/Math/RANSAC.h>

using namespace vw;
using namespace vw::camera;

class CameraGeometryTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    // Create Camera
    Matrix<double,3,3> pose = math::euler_to_rotation_matrix(1.3,2.0,-.7,"xyz");
    pinhole = PinholeModel( Vector3(-1,4,2),
                            pose, 600, 700,
                            500, 500, NullLensDistortion() );

    // Building measurements
    boost::minstd_rand random_gen(42u);
    boost::normal_distribution<double> normal(0,20);
    boost::variate_generator<boost::minstd_rand&,
      boost::normal_distribution<double> > generator( random_gen, normal );
    for ( uint8 i = 0; i < 50; i++ ) {
      world_m.push_back( Vector4( generator(), generator(),
                                  generator() + 60.0, 1.0 ) );
      subvector(world_m.back(),0,3) = inverse(pose)*subvector(world_m.back(),0,3);
      Vector2 pixel = pinhole.point_to_pixel( subvector(world_m.back(),0,3) );
      image_m.push_back( Vector3( pixel[0], pixel[1], 1.0 ) );
    }

    // Building even noiser data
    for ( uint8 i = 0; i < 50; i++ ) {
      Vector4 noise( generator()*0.025, generator()*0.025,
                     generator()*0.025, 0.0 );
      if ( norm_2( noise ) > 1 )
        noise /= norm_2(noise);
      Vector3 im_noise( generator()*0.01, generator()*0.01, 0.0 );
      noisy_world_m.push_back( world_m[i] + noise );
      noisy_image_m.push_back( image_m[i] + im_noise );
    }
  }

  PinholeModel pinhole;
  std::vector<Vector<double> > world_m, image_m;
  std::vector<Vector<double> > noisy_world_m, noisy_image_m;
};

TEST_F( CameraGeometryTest, LinearSolve ) {
  std::vector<Vector<double> > world_small, image_small;
  for ( uint32 i = 0; i < 6; i++ ) {
    world_small.push_back( world_m[i] );
    image_small.push_back( image_m[i] );
  }
  Matrix<double> P =
    CameraMatrixFittingFunctor()(world_small,image_small);
  ASSERT_EQ( P.rows(), 3u );
  ASSERT_EQ( P.cols(), 4u );

  for ( uint8 i = 0; i < 10; i++ ) {
    Vector3 p_result = P*noisy_world_m[i];
    p_result /= p_result[2];
    Vector2 cam_result =
      pinhole.point_to_pixel(subvector(noisy_world_m[i],0,3));
    EXPECT_VECTOR_NEAR( subvector(p_result,0,2),
                        cam_result, 1 );
  }
}

TEST_F( CameraGeometryTest, IteratorSolve ) {
  Matrix<double> P =
    CameraMatrixFittingFunctor()(noisy_world_m,
                                 noisy_image_m );
  ASSERT_EQ( P.rows(), 3u );
  ASSERT_EQ( P.cols(), 4u );

  for ( uint8 i = 0; i < 10; i++ ) {
    Vector3 p_result = P*world_m[i];
    p_result /= p_result[2];
    Vector2 cam_result =
      pinhole.point_to_pixel(subvector(world_m[i],0,3));
    EXPECT_VECTOR_NEAR( subvector(p_result,0,2),
                        cam_result, 20); //
  }
}

TEST_F( CameraGeometryTest, DISABLED_RansacSolve ) {
  // Add some extra wrong data points
  // 10% are just completely wrong
  for ( uint32 i=0,j=0; i < 5 && j <50; i++, j+=10 ) {
    world_m[j][0] = rand()%200-100;
    world_m[j][1] = rand()%200-100;
    world_m[j][2] = rand()%200;
  }

  math::RandomSampleConsensus<CameraMatrixFittingFunctor,CameraMatrixErrorMetric> ransac( CameraMatrixFittingFunctor(), CameraMatrixErrorMetric(), 2 );
  Matrix<double> P(ransac(world_m,image_m));
  ASSERT_EQ( P.rows(), 3u );
  ASSERT_EQ( P.cols(), 4u );

  for ( uint8 i = 0; i < 10; i++ ) {
    Vector3 p_result = P*world_m[i];
    p_result /= p_result[2];
    Vector2 cam_result =
      pinhole.point_to_pixel(subvector(world_m[i],0,3));
    EXPECT_VECTOR_NEAR( subvector(p_result,0,2),
                        cam_result, 500);
    // I don't understand why this are so far off.
  }
}

// Test Fundamental Matrix ----------------------------------------------
class FundamentalMatrixStaticTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    // Push back standard points
    world_points.clear();
    world_points.push_back( Vector3(-1.197167786710194e-01, 1.535837754235135, 6.160134902124728e-01) );
    world_points.push_back( Vector3(-2.398868915842913e-01, 2.592493439223233, 2.106859900249888e-01) );
    world_points.push_back( Vector3(5.700164332624225e-02, 2.680805101933800, 8.513699491560502e-01) );
    world_points.push_back( Vector3(-5.638564814956753e-03, 2.159721813657774, 1.300915066875919) );
    world_points.push_back( Vector3(4.668814127918872e-01, 1.557199966822651, 5.615713956684834e-01) );
    world_points.push_back( Vector3(3.549516562975479e-01, 1.488383485452233, 5.962236919954063e-01) );
    world_points.push_back( Vector3(-2.272124188441233e-02, 2.096124383794463, 2.316413346082943e-01) );
    world_points.push_back( Vector3(-6.089950299762996e-01, 1.642522467534197, 1.124964063561710) );
    world_points.push_back( Vector3(-4.324481312564093e-02, 1.682823080340085, 9.892865764503735e-01) );
    world_points.push_back( Vector3(3.063284697835341e-01, 1.712468548845789, 8.583171749675742e-01) );
    world_points.push_back( Vector3(3.701320821868042e-01, 1.742080228159405, 6.799617290423526e-01) );
    world_points.push_back( Vector3(-1.235800856367202e-01, 2.416933695634718e+00, 7.931139899916800e-02) );
    world_points.push_back( Vector3(-1.361779654316855e-01, 2.273596665248890e+00, 4.734887351320431e-01) );
    world_points.push_back( Vector3(1.921310249433029e-01, 2.579181736347006e+00, 3.578803440576578e-01) );
    world_points.push_back( Vector3(-3.652748558723435e-01,2.520761700565743e+00,3.542000261284033e-01) );
    world_points.push_back( Vector3(1.060947027908127e-01,1.806012449723466e+00,-2.293187222584325e-02) );
    world_points.push_back( Vector3(-3.157784173745294e-01,2.437412553760792e+00,3.616192474219446e-01) );
    world_points.push_back( Vector3(-5.006093526417452e-02,2.100012579179123e+00,8.922727263631631e-01) );
    world_points.push_back( Vector3(5.079908576059367e-02,1.913836937602315e+00,2.158429714772892e-01) );
    world_points.push_back( Vector3(1.942716829862292e-01,1.303897751202779e+00,1.083393790683977e-01) );

    // Construct standard cameras
    pinhole1 =
      PinholeModel( Vector3(), Matrix3x3(1,0,0,0,0,1,0,-1,0),
                    700, 700, 640, 480, Vector3(1,0,0),
                    Vector3(0,1,0), Vector3(0,0,1),
                    NullLensDistortion() );
    pinhole2 =
      PinholeModel( Vector3(-.7,-.7,0),
                    Matrix3x3(1,0,0,0,0,1,0,-1,0)*math::rotation_y_axis(M_PI/6),
                    700, 700, 640, 480, Vector3(1,0,0),
                    Vector3(0,1,0), Vector3(0,0,1),
                    NullLensDistortion() );

    // Getting measurements
    measure1.clear(); measure2.clear();
    BOOST_FOREACH( Vector3 const& point, world_points ) {
      measure1.push_back( pinhole1.point_to_pixel( point ) );
      measure2.push_back( pinhole2.point_to_pixel( point ) );
    }

    expected_F = Matrix3x3(-8.111434829751515e-12,1.093101634752133e-05,-5.246849944201554e-03,-1.493213410928517e-05,-1.685262742074159e-12,1.235728252446591e-02,7.167400675360256e-03,-1.464760650157697e-02,1.099333713562344e+00);
  }

  std::vector<Vector3> world_points;
  PinholeModel pinhole1, pinhole2;
  std::vector<Vector<double> > measure1, measure2;
  Matrix<double> expected_F; // Result from Epipolar Geometry Toolbox
};

// Test the test
TEST_F( FundamentalMatrixStaticTest, SanityCheck ) {
  EXPECT_EQ( 2, rank(expected_F) );
  EXPECT_NEAR( 1, norm_frobenius(expected_F), 0.1 );

  // EGT guess is worse because they perform denormalization out of
  // step against Hartley's recommendation.
  for ( unsigned i = 0; i < measure1.size(); i++ )  {
    EXPECT_LT( FundamentalMatrixSampsonErrorMetric()(expected_F, Vector3( measure1[i][0], measure1[i][1], 1),  Vector3( measure2[i][0], measure2[i][1], 1) ), 7 );
  }
}

TEST_F( FundamentalMatrixStaticTest, EightPointAlgorithm ) {
  Matrix<double> F = FundamentalMatrix8PFittingFunctor()( measure1, measure2 );

  EXPECT_EQ( 2, rank(F) );
  EXPECT_NEAR( 1, norm_frobenius(F), 0.1 );

  for ( unsigned i = 0; i < measure1.size(); i++ )  {
    EXPECT_LT( FundamentalMatrixSampsonErrorMetric()(F, Vector3( measure1[i][0], measure1[i][1], 1),  Vector3( measure2[i][0], measure2[i][1], 1) ), 1e-6 );
  }
}

TEST_F( FundamentalMatrixStaticTest, MLAlgorithm ) {
  boost::minstd_rand random_gen(42u);
  boost::normal_distribution<double> normal(0,4);
  boost::variate_generator<boost::minstd_rand&,
    boost::normal_distribution<double> > generator( random_gen, normal );

  // Adding Noise to measurements
  for ( unsigned i = 0; i < measure1.size(); i++ ) {
    measure1[i] += Vector2( generator(), generator() );
    measure2[i] += Vector2( generator(), generator() );
  }

  // Creating Seed
  Matrix<double> seed =
    FundamentalMatrix8PFittingFunctor()( measure1, measure2 );

  // Actual measurement
  Matrix<double> F =
    FundamentalMatrixMLFittingFunctor()( measure1, measure2, seed );

  EXPECT_EQ( 2, rank(F) );
  EXPECT_NEAR( 1, norm_frobenius(F), 0.4 );

  for ( unsigned i = 0; i < measure1.size(); i++ )  {
    EXPECT_LT( FundamentalMatrixSampsonErrorMetric()(seed, Vector3( measure1[i][0], measure1[i][1], 1), Vector3( measure2[i][0], measure2[i][1], 1) ), 0.2 );
    EXPECT_LT( FundamentalMatrixSampsonErrorMetric()(F, Vector3( measure1[i][0], measure1[i][1], 1),  Vector3( measure2[i][0], measure2[i][1], 1) ), 0.15 );
  }
}
