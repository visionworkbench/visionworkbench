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


#include <gtest/gtest_VW.h>

#include <vw/Math/Vector.h>
#include <vw/Math/EulerAngles.h>
#include <vw/Math/LinearAlgebra.h>
#include <vw/Math/MatrixSparseSkyline.h>
#include <test/Helpers.h>

#include <vw/Camera/PinholeModel.h>
#include <vw/BundleAdjustment/AdjustRef.h>
#include <vw/BundleAdjustment/AdjustRobustRef.h>
#include <vw/BundleAdjustment/AdjustSparse.h>
#include <vw/BundleAdjustment/AdjustRobustSparse.h>

// This test is for non robust bundle adjustment.

using namespace vw;
using namespace vw::camera;
using namespace vw::ba;

// Building a Model
// ------------------------
class TestBAModel : public ba::ModelBase< TestBAModel, 6, 3 > {

  typedef Vector<double, 6> camera_vector_t;
  typedef Vector<double, 3> point_vector_t;

  std::vector< boost::shared_ptr<PinholeModel> > m_cameras;
  boost::shared_ptr<ControlNetwork> m_cnet;
  std::vector<camera_vector_t> m_cam_vec, m_cam_target_vec;
  std::vector<point_vector_t> m_point_vec, m_point_target_vec;
  size_t m_num_pixel_observations;

public:
  // Constructor
  TestBAModel( std::vector< boost::shared_ptr<PinholeModel> > const& cameras,
               boost::shared_ptr<ControlNetwork> network ) : m_cameras(cameras), m_cnet(network) {

    // Compute the number of observations from the bundle.
    m_num_pixel_observations = 0;
    for (size_t i = 0; i < network->size(); ++i)
      m_num_pixel_observations += (*network)[i].size();

    // Setting up the camera vectors
    m_cam_vec.resize( m_cameras.size() );
    m_cam_target_vec.resize( m_cam_vec.size() );
    for ( size_t j = 0; j < m_cameras.size(); j++ ) {
      m_cam_vec[j] = camera_vector_t();
      m_cam_target_vec[j] = m_cam_vec[j];
    }

    // Setting up the point vectors
    m_point_vec.resize( m_cnet->size() );
    m_point_target_vec.resize( m_point_vec.size() );
    for ( size_t i = 0; i < m_cnet->size(); i++ ) {
      m_point_vec[i] = (*m_cnet)[i].position();
      m_point_target_vec[i] = m_point_vec[i];
    }

  }

  // -- REQUIRED STUFF ---------------------------------------

  // Access to the cameras
  Vector2 cam_pixel(size_t /*i*/, size_t j,
		    camera_vector_t const& cam_j,
		    point_vector_t const& point_i ) const {
    // Quaternions are the last half of this equation
    AdjustedCameraModel cam( m_cameras[j],
                             subvector(cam_j,0,3),
                             math::euler_to_quaternion(cam_j[3],cam_j[4],cam_j[5],"xyz") );

    return cam.point_to_pixel( point_i );
  }

  inline Matrix<double,6,6> cam_inverse_covariance( size_t /*j*/ ) {
    Matrix<double,6,6> result;
    result.set_identity();
    result *= 2;
    return result;
  }
  inline Matrix<double,3,3> point_inverse_covariance( size_t /*i*/ ) {
    Matrix<double,3,3> result;
    result.set_identity();
    result *= 2;
    return result;
  }

  size_t num_cameras() const { return m_cam_vec.size(); }
  size_t num_points() const { return m_point_vec.size(); }
  camera_vector_t cam_params( size_t j ) const { return m_cam_vec[j]; }
  point_vector_t point_params( size_t i ) const { return m_point_vec[i]; }
  camera_vector_t cam_target( size_t j ) const { return m_cam_target_vec[j]; }
  point_vector_t point_target( size_t i ) const { return m_point_target_vec[i]; }
  size_t num_pixel_observations() const { return m_num_pixel_observations; }
  void set_cam_params(size_t j, camera_vector_t const& cam_j) { m_cam_vec[j] = cam_j; }
  void set_point_params(size_t i, point_vector_t const& point_i) { m_point_vec[i] = point_i; }

  boost::shared_ptr<ControlNetwork> control_network(void) {
    return m_cnet; }
};

// Generating Data
// ----------------------

// Function for surface
inline float surface_func( float const& x, float const& y ) {
  return 4*sin( 10*x*y );
}

// Generate Camera data
inline void generate_camera_data( std::vector<boost::shared_ptr<PinholeModel> > & cameras,
                                  boost::shared_ptr<ControlNetwork> & cnet ) {
  Matrix<double,3,3> pose1 = math::euler_to_rotation_matrix(M_PI-0.01,0,
                                                            0.212,"xyz");
  Matrix<double,3,3> pose2 = math::euler_to_rotation_matrix(0.001,M_PI+0.003,
                                                            0.7,"xyz");
  Matrix<double,3,3> pose3 = math::euler_to_rotation_matrix(M_PI+0.002,0.01,
                                                            1.4,"xyz");

  cameras.push_back( boost::shared_ptr<PinholeModel>( new PinholeModel(Vector3(-10,-12,20), pose1, 600, 600, 500, 500) ) );
  cameras.push_back( boost::shared_ptr<PinholeModel>( new PinholeModel(Vector3(-17, 19,22), pose2, 600, 600, 500, 500) ) );
  cameras.push_back( boost::shared_ptr<PinholeModel>( new PinholeModel(Vector3(  0,  0,19), pose3, 500, 450, 500, 500) ) );
  cameras.push_back( boost::shared_ptr<PinholeModel>( new PinholeModel(Vector3( 21,-17,21), pose1, 600, 600, 500, 500) ) );
  cameras.push_back( boost::shared_ptr<PinholeModel>( new PinholeModel(Vector3( 12, 10,18), pose2, 600, 600, 500, 500) ) );

  BBox2i image(0,0,1000,1000);
  cnet = boost::shared_ptr<ControlNetwork>( new ControlNetwork("Test case:") );

  for ( int i = -10; i <= 10; i += 2 ) {
    for ( int j = -10; j <= 10; j += 3 ) {
      Vector3 position(i,j,surface_func(i,j));

      ControlPoint cpoint( ControlPoint::TiePoint );
      cpoint.set_position( position );
      cpoint.set_sigma( Vector3(2,2,2) );

      for ( uint32 ci = 0; ci < cameras.size(); ci++ ) {
        Vector2i pixel = cameras[ci]->point_to_pixel( position );
        if ( image.contains(pixel) ) {
          ControlMeasure cmeasure( pixel[0], pixel[1],
                                   1, 1, ci );
          cpoint.add_measure( cmeasure );
        }
      }
      cnet->add_control_point( cpoint );
    }
  }
}

// Standard test that provides perfect data
class NullTest : public ::testing::Test {
protected:
  NullTest() {}

  virtual void SetUp() {
    generate_camera_data( cameras, cnet );
  }

  std::vector<boost::shared_ptr<PinholeModel> > cameras;
  boost::shared_ptr<ControlNetwork> cnet;
};

// Test that provides corrupted data
class ComparisonTest : public ::testing::Test {
protected:
  ComparisonTest() {}

  virtual void SetUp() {
    generate_camera_data( cameras, cnet );

    // Corrupting cameras a bit
    cameras[0]->set_camera_center( cameras[0]->camera_center() +
                                   Vector3(1.2,0,-1) );
    cameras[1]->set_camera_center( cameras[1]->camera_center() +
                                   Vector3(0.5,2,1) );
    cameras[2]->set_camera_center( cameras[2]->camera_center() +
                                   Vector3(-0.2,-1,3.0) );
    cameras[3]->set_camera_pose( cameras[3]->camera_pose() +
                                 math::euler_to_quaternion(0.1,-0.1,0.0,"xyz") );
    for ( uint32 i = 0; i < cnet->size(); i++ ) {
      if ( i % 2 ) {
        (*cnet)[i].set_position( (*cnet)[i].position()+Vector3(-1,0.5,-0.7) );
      } else {
        (*cnet)[i].set_position( (*cnet)[i].position()+Vector3(0.4,-1.5,0.3) );
      }
    }
  }

  std::vector<boost::shared_ptr<PinholeModel> > cameras;
  boost::shared_ptr<ControlNetwork> cnet;
};

// Null Tests
// -----------------------
TEST_F( NullTest, AdjustRef ) {
  TestBAModel model( cameras, cnet );
  AdjustRef< TestBAModel, L2Error > adjuster( model, L2Error(), false, false);

  // Running BA
  double abs_tol = 1e10, rel_tol = 1e10;
  for ( uint32 i = 0; i < 5; i++ )
    adjuster.update(abs_tol,rel_tol);

  // Checking solutions
  Vector<double,6> zero_vector;
  for ( uint32 i = 0; i < 5; i++ ) {
    Vector<double> solution = model.cam_params(i);
    EXPECT_VECTOR_NEAR( solution, zero_vector, 1e-1 );
  }
}

TEST_F( NullTest, AdjustSparse ) {
  TestBAModel model( cameras, cnet );
  AdjustSparse< TestBAModel, L2Error > adjuster( model, L2Error(), false, false);

  // Running BA
  double abs_tol = 1e10, rel_tol = 1e10;
  for ( uint32 i = 0; i < 5; i++ )
    adjuster.update(abs_tol,rel_tol);

  // Checking solutions
  Vector<double,6> zero_vector;
  for ( uint32 i = 0; i < 5; i++ ) {
    Vector<double> solution = model.cam_params(i);
    EXPECT_VECTOR_NEAR( solution, zero_vector, 1e-1 );
  }
}

TEST_F( NullTest, AdjustRobustRef ) {
  TestBAModel model( cameras, cnet );
  AdjustRobustRef< TestBAModel, L2Error > adjuster( model, L2Error(), false, false);

  // Running BA
  double abs_tol = 1e10, rel_tol = 1e10;
  for ( uint32 i = 0; i < 5; i++ )
    adjuster.update(abs_tol,rel_tol);

  // Checking solutions
  Vector<double,6> zero_vector;
  for ( uint32 i = 0; i < 5; i++ ) {
    Vector<double> solution = model.cam_params(i);
    EXPECT_VECTOR_NEAR( solution, zero_vector, 1e-1 );
  }
}

TEST_F( NullTest, AdjustRobustSparse ) {
  TestBAModel model( cameras, cnet );
  AdjustRobustSparse< TestBAModel, L2Error > adjuster( model, L2Error(), false, false);

  // Running BA
  double abs_tol = 1e10, rel_tol = 1e10;
  for ( uint32 i = 0; i < 5; i++ )
    adjuster.update(abs_tol,rel_tol);

  // Checking solutions
  Vector<double,6> zero_vector;
  for ( uint32 i = 0; i < 5; i++ ) {
    Vector<double> solution = model.cam_params(i);
    EXPECT_VECTOR_NEAR( solution, zero_vector, 1e-1 );
  }
}

// Comparison Tests
// -----------------------
TEST_F( ComparisonTest, Ref_VS_Sparse ) {
  std::vector<Vector<double> > ref_solution;
  std::vector<Vector<double> > spr_solution;

  { // Performing Ref BA
    TestBAModel model( cameras, cnet );
    AdjustRef< TestBAModel, L2Error > adjuster( model, L2Error(), false, false);

    // Running BA
    double abs_tol = 1e10, rel_tol = 1e10;
    for ( unsigned i = 0; i < 10; i++ )
      adjuster.update(abs_tol,rel_tol);

    // Storing result
    for ( uint32 i = 0; i < 5; i++ )
      ref_solution.push_back( model.cam_params(i) );
  }

  { // Performing Sparse BA
    TestBAModel model( cameras, cnet );
    AdjustSparse< TestBAModel, L2Error > adjuster( model, L2Error(), false, false);

    // Running BA
    double abs_tol = 1e10, rel_tol = 1e10;
    for ( unsigned i = 0; i < 10; i++ )
      adjuster.update(abs_tol,rel_tol);

    // Storing result
    for ( uint32 i = 0; i < 5; i++ )
      spr_solution.push_back( model.cam_params(i) );
  }

  // Comparison
  for ( uint32 i = 0; i < 5; i++ )
    ASSERT_VECTOR_NEAR( ref_solution[i],
                        spr_solution[i],
                        1e-3 );
}

// For whatever reason .. RobustRef and RobustSparse diverge
// quickly. This is probably do to unwise application of floats or
// arithmetic ordering.
TEST_F( ComparisonTest, DISABLED_RobustRef_VS_RobustSparse ) {
  std::vector<Vector<double> > ref_solution;
  std::vector<Vector<double> > spr_solution;

  { // Performing Ref BA
    TestBAModel model( cameras, cnet );
    AdjustRobustRef< TestBAModel, L2Error > adjuster( model, L2Error(),
                                                      false, false);

    // Running BA
    double abs_tol = 1e10, rel_tol = 1e10;
    for ( unsigned i = 0; i < 2; i++ )
      adjuster.update(abs_tol,rel_tol);

    // Storing result
    for ( uint32 i = 0; i < 5; i++ )
      ref_solution.push_back( model.cam_params(i) );
  }

  { // Performing Sparse BA
    TestBAModel model( cameras, cnet );
    AdjustRobustSparse< TestBAModel, L2Error > adjuster( model, L2Error(),
                                                         false, false);

    // Running BA
    double abs_tol = 1e10, rel_tol = 1e10;
    for ( unsigned i = 0; i < 2; i++ )
      adjuster.update(abs_tol,rel_tol);

    // Storing result
    for ( uint32 i = 0; i < 5; i++ )
      spr_solution.push_back( model.cam_params(i) );
  }

  // Comparison
  for ( uint32 i = 0; i < 5; i++ )
    ASSERT_VECTOR_NEAR( ref_solution[i],
                        spr_solution[i],
                        1e-2 );
}
