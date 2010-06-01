// __BEGIN_LICENSE__
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
  CameraGeometryTest() {}

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
  for ( uint i = 0; i < 6; i++ ) {
    world_small.push_back( world_m[i] );
    image_small.push_back( image_m[i] );
  }
  Matrix<double> P =
    CameraMatrixFittingFunctor()(world_small,image_small);
  ASSERT_EQ( P.rows(), 3 );
  ASSERT_EQ( P.cols(), 4 );

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
  ASSERT_EQ( P.rows(), 3 );
  ASSERT_EQ( P.cols(), 4 );

  for ( uint8 i = 0; i < 10; i++ ) {
    Vector3 p_result = P*world_m[i];
    p_result /= p_result[2];
    Vector2 cam_result =
      pinhole.point_to_pixel(subvector(world_m[i],0,3));
    EXPECT_VECTOR_NEAR( subvector(p_result,0,2),
                        cam_result, 20); //
  }
}

TEST_F( CameraGeometryTest, RansacSolve ) {
  // Add some extra wrong data points
  // 10% are just completely wrong
  for ( uint i=0,j=0; i < 5 && j <50; i++, j+=10 ) {
    world_m[j][0] = rand()%200-100;
    world_m[j][1] = rand()%200-100;
    world_m[j][2] = rand()%200;
  }

  math::RandomSampleConsensus<CameraMatrixFittingFunctor,CameraMatrixErrorMetric> ransac( CameraMatrixFittingFunctor(), CameraMatrixErrorMetric(), 2 );
  Matrix<double> P(ransac(world_m,image_m));
  ASSERT_EQ( P.rows(), 3 );
  ASSERT_EQ( P.cols(), 4 );

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


