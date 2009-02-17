// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <math.h>

#include <ctime>
#include <cstdlib>

// TestMatrix.h
#include <cxxtest/TestSuite.h>
#include <vw/Math/Quaternion.h>

using namespace vw;

class TestQuaternion : public CxxTest::TestSuite
{
public:
  void test_axis_angle() {
    Vector<double,3> x_axis(1,0,0);
    Vector<double,3> y_axis(0,1,0);
    Vector<double,3> z_axis(0,0,1);

    // No rotation;  Quaternion should be (1,0,0,0)
    Quaternion<double> q1(x_axis, 0);
    TS_ASSERT_DELTA(math::norm_2(q1 - Quaternion<double>(1,0,0,0)), 0, 1e-10);

    // 90-degree rotation about X axis;  Quaternion should be (.707,.707,0,0)
    Quaternion<double> q2(x_axis, .5*M_PI);
    TS_ASSERT_DELTA(math::norm_2(q2 - Quaternion<double>(1./sqrt(2.),1./sqrt(2.),0,0)), 0, 1e-10);

    // 180-degree rotation about X axis;  Quaternion should be (0,1,0,0)
    Quaternion<double> q3(x_axis, M_PI);
    TS_ASSERT_DELTA(math::norm_2(q3 - Quaternion<double>(0,1,0,0)), 0, 1e-10);

    // 180-degree rotation about Y axis;  Quaternion should be (0,0,1,0)
    Quaternion<double> q4(y_axis, M_PI);
    TS_ASSERT_DELTA(math::norm_2(q4 - Quaternion<double>(0,0,1,0)), 0, 1e-10);

    // 180-degree rotation about Z axis;  Quaternion should be (0,0,0,1)
    Quaternion<double> q5(z_axis, M_PI);
    TS_ASSERT_DELTA(math::norm_2(q5 - Quaternion<double>(0,0,0,1)), 0, 1e-10);

    // Arbitrary axis
    Vector<double,3> axis(1.2, -2.3, 3.4);
    axis /= math::norm_2(axis);

    for (double angle = -4*M_PI; angle < 4*M_PI; angle += .1) {
      Quaternion<double> q(axis, angle);
      Vector<double,3> test_axis;
      double test_angle;
      q.axis_angle(test_axis, test_angle);
      //TS_TRACE(stringify("in: "        ) + stringify(axis) + "," + stringify(angle));
      //TS_TRACE(stringify("quaternion: ") + stringify(q));
      //TS_TRACE(stringify("out: "       ) + stringify(test_axis) + "," + stringify(test_angle))

      if (math::dot_prod(test_axis, axis) < 0) {
        test_axis = -test_axis;
        test_angle = -test_angle;
      }

      TS_ASSERT_DELTA(fmod(angle - test_angle, M_PI*2), 0, 1e-10);
      if (fabs(fmod(angle, M_PI*2)) > 1e-10) {
        TS_ASSERT_DELTA(math::norm_2(test_axis - axis), 0, 1e-10);
      }
    }
  }

  void test_rotation()
  {
    // Create a test vector to rotate
    Vector3 test_vec(0.1, 0.84, 0.23);

    double theta = 32*M_PI/180;  // 32 degree rotation
    Matrix3x3 rotation = math::identity_matrix<3>();
    rotation(0,0) = cos(theta);
    rotation(0,1) = -sin(theta);
    rotation(1,0) = sin(theta);
    rotation(1,1) = cos(theta);

    Quaternion<double> quat(rotation);
    Matrix3x3 rotation2 = quat.rotation_matrix();

    Vector3 sol1 = rotation * test_vec;
    Vector3 sol2 = rotation2 * test_vec;
    Vector3 sol3 = quat.rotate(test_vec);
    Vector3 sol4 = -quat.rotate(-test_vec);

    TS_ASSERT_DELTA(sol1(0),sol2(0),1e-9);
    TS_ASSERT_DELTA(sol1(1),sol2(1),1e-9);
    TS_ASSERT_DELTA(sol1(2),sol2(2),1e-9);
    TS_ASSERT_DELTA(sol1(0),sol3(0),1e-9);
    TS_ASSERT_DELTA(sol1(1),sol3(1),1e-9);
    TS_ASSERT_DELTA(sol1(2),sol3(2),1e-9);
    TS_ASSERT_DELTA(sol1(0),sol4(0),1e-9);
    TS_ASSERT_DELTA(sol1(1),sol4(1),1e-9);
    TS_ASSERT_DELTA(sol1(2),sol4(2),1e-9);
  }


  void test_rotationmatrix() {
    vw::Matrix3x3 I;
    I.set_identity();

    vw::Matrix3x3 T, dQ;


    // test ww
    double R_data_ww[9] = { 0.36, 0.48, -0.8, -0.8, 0.6, 0, 0.48, 0.64, 0.6 };
    vw::Matrix3x3 R_ww(R_data_ww);

    // check if R*transpose(R) is identity
    T = R_ww * transpose(R_ww) - I;
    for (Matrix3x3::const_iterator it = T.begin(); it != T.end(); ++it)
      TS_ASSERT_DELTA(*it, 0, 1e-14);

    // converting the matrix to quaternion and back should not change it
    dQ = R_ww - vw::Quaternion<double>(R_ww).rotation_matrix();
    for (Matrix3x3::const_iterator it = dQ.begin(); it != dQ.end(); ++it)
      TS_ASSERT_DELTA(*it, 0, 1e-14);


    // test xx
    double R_data_xx[9] = { 0.11276615807984,   -0.0161365698014631,  0.993490515660293,
			    -0.633815522454304, -0.771198883306014,   0.0594151990953122,
			    0.765220018744864, -0.636389733970151,  -0.097192743605217,
    };
    vw::Matrix3x3 R_xx(R_data_xx);

    // check if is identity
    T = R_xx * transpose(R_xx) - I;
    for (Matrix3x3::const_iterator it = T.begin(); it != T.end(); ++it)
      TS_ASSERT_DELTA(*it, 0, 1e-14);

    // converting the matrix to quaternion and back should not change it
    dQ = R_xx - vw::Quaternion<double>(R_xx).rotation_matrix();
    for (Matrix3x3::const_iterator it = dQ.begin(); it != dQ.end(); ++it)
      TS_ASSERT_DELTA(*it, 0, 1e-14);


    // test yy
    double R_data_yy[9] = { -0.771198883306014, -0.633815522454304,  0.0594151990953122,
			    -0.0161365698014631, 0.11276615807984,   0.993490515660293,
			    -0.636389733970151,  0.765220018744864, -0.097192743605217,
    };
    vw::Matrix3x3 R_yy(R_data_yy);

    // check if is identity
    T = R_yy * transpose(R_yy) - I;
    for (Matrix3x3::const_iterator it = T.begin(); it != T.end(); ++it)
      TS_ASSERT_DELTA(*it, 0, 1e-14);

    // converting the matrix to quaternion and back should not change it
    dQ = R_yy - vw::Quaternion<double>(R_yy).rotation_matrix();
    for (Matrix3x3::const_iterator it = dQ.begin(); it != dQ.end(); ++it)
      TS_ASSERT_DELTA(*it, 0, 1e-14);


    // test zz
    double R_data_zz[9] = { -0.771198883306014,  0.0594151990953122, -0.633815522454304,
			    -0.636389733970151, -0.097192743605217,   0.765220018744864,
			    -0.0161365698014631, 0.993490515660293,   0.11276615807984 
    };
    vw::Matrix3x3 R_zz(R_data_zz);

    T = R_zz * transpose(R_zz) - I;
    for (Matrix3x3::const_iterator it = T.begin(); it != T.end(); ++it)
      TS_ASSERT_DELTA(*it, 0, 1e-14);

    // converting the matrix to quaternion and back should not change it
    dQ = R_zz - vw::Quaternion<double>(R_zz).rotation_matrix();
    for (Matrix3x3::const_iterator it = dQ.begin(); it != dQ.end(); ++it)
      TS_ASSERT_DELTA(*it, 0, 1e-14);
  }

  void test_axis_angle_quaternion() {
    vw::Vector3 x(M_PI/2., 0, 0);
    vw::Vector3 y(0, M_PI/2., 0);
    vw::Vector3 z(0, 0, M_PI/2.);
    
    // Test constructor taking normalized axis and angle
    Quaternion<double> q_x(x / vw::math::norm_2(x), vw::math::norm_2(x));
    vw::Vector3 test(2.13, -812.183, 18.31);
    // rotate vector test by 90° around x axis and compare
    vw::Vector3 diff_1( vw::Vector3(2.13, -18.31, -812.183) - q_x.rotate(test) );
    TS_ASSERT_DELTA(diff_1.x(), 0, 1e-12);
    TS_ASSERT_DELTA(diff_1.y(), 0, 1e-12);
    TS_ASSERT_DELTA(diff_1.z(), 0, 1e-12);

    // Test standalone method taking nonnormalized axis with angle = |axis|
    vw::Vector3 diff_2( vw::Vector3(2.13, -18.31, -812.183) -  vw::math::axis_angle_to_quaternion(x).rotate(test));
    TS_ASSERT_DELTA(diff_2.x(), 0, 1e-12);
    TS_ASSERT_DELTA(diff_2.y(), 0, 1e-12);
    TS_ASSERT_DELTA(diff_2.z(), 0, 1e-12);
    
    // Test conversion back to vector using two-step method
    vw::Vector3 axis_x;
    double angle_x;
    q_x.axis_angle(axis_x, angle_x);
    TS_ASSERT_DELTA(x.x(), angle_x*axis_x.x(), 1e-10);
    TS_ASSERT_DELTA(x.y(), angle_x*axis_x.y(), 1e-10);
    TS_ASSERT_DELTA(x.z(), angle_x*axis_x.z(), 1e-10);

    // Test conversion back to vector using "direct" method
    vw::Vector3 back_x( q_x.axis_angle() );
    TS_ASSERT_DELTA(x.x(), back_x.x(), 1e-10);
    TS_ASSERT_DELTA(x.y(), back_x.y(), 1e-10);
    TS_ASSERT_DELTA(x.z(), back_x.z(), 1e-10);


    
    // Test constructor taking normalized axis and angle
    Quaternion<double> q_y(y / vw::math::norm_2(y), vw::math::norm_2(y));
    vw::Vector3 test_y(-58.12, -81.19, -48);
    // rotate vector test by 90° around y axis and compare
    diff_1 = vw::Vector3( vw::Vector3(-48, -81.19, 58.12) - q_y.rotate(test_y) );
    TS_ASSERT_DELTA(diff_1.x(), 0, 1e-12);
    TS_ASSERT_DELTA(diff_1.y(), 0, 1e-12);
    TS_ASSERT_DELTA(diff_1.z(), 0, 1e-12);

    // Test standalone method taking nonnormalized axis with angle = |axis|
    diff_2 = vw::Vector3( vw::Vector3(-48, -81.19, 58.12) -  vw::math::axis_angle_to_quaternion(y).rotate(test_y));
    TS_ASSERT_DELTA(diff_2.x(), 0, 1e-12);
    TS_ASSERT_DELTA(diff_2.y(), 0, 1e-12);
    TS_ASSERT_DELTA(diff_2.z(), 0, 1e-12);
    
    // Test conversion back to vector using two-step method
    vw::Vector3 axis_y;
    double angle_y;
    q_y.axis_angle(axis_y, angle_y);
    TS_ASSERT_DELTA(y.x(), angle_y*axis_y.x(), 1e-10);
    TS_ASSERT_DELTA(y.y(), angle_y*axis_y.y(), 1e-10);
    TS_ASSERT_DELTA(y.z(), angle_y*axis_y.z(), 1e-10);

    // Test conversion back to vector using "direct" method
    vw::Vector3 back_y( q_y.axis_angle() );
    TS_ASSERT_DELTA(y.x(), back_y.x(), 1e-10);
    TS_ASSERT_DELTA(y.y(), back_y.y(), 1e-10);
    TS_ASSERT_DELTA(y.z(), back_y.z(), 1e-10);



    // Test constructor taking normalized axis and angle
    Quaternion<double> q_z(z / vw::math::norm_2(z), vw::math::norm_2(z));
    vw::Vector3 test_z(183.12, 845.43, 73.2);
    // rotate vector test by 90° around z axis and compare
    diff_1 = vw::Vector3( vw::Vector3(-845.43, 183.12, 73.2) - q_z.rotate(test_z) );
    TS_ASSERT_DELTA(diff_1.x(), 0, 1e-12);
    TS_ASSERT_DELTA(diff_1.y(), 0, 1e-12);
    TS_ASSERT_DELTA(diff_1.z(), 0, 1e-12);

    // Test standalone method taking nonnormalized axis with angle = |axis|
    diff_2 = vw::Vector3( vw::Vector3(-845.43, 183.12, 73.2) -  vw::math::axis_angle_to_quaternion(z).rotate(test_z));
    TS_ASSERT_DELTA(diff_2.x(), 0, 1e-12);
    TS_ASSERT_DELTA(diff_2.y(), 0, 1e-12);
    TS_ASSERT_DELTA(diff_2.z(), 0, 1e-12);
    
    // Test conversion back to vector using two-step method
    vw::Vector3 axis_z;
    double angle_z;
    q_z.axis_angle(axis_z, angle_z);
    TS_ASSERT_DELTA(z.x(), angle_z*axis_z.x(), 1e-10);
    TS_ASSERT_DELTA(z.y(), angle_z*axis_z.y(), 1e-10);
    TS_ASSERT_DELTA(z.z(), angle_z*axis_z.z(), 1e-10);

    // Test conversion back to vector using "direct" method
    vw::Vector3 back_z( q_z.axis_angle() );
    TS_ASSERT_DELTA(z.x(), back_z.x(), 1e-10);
    TS_ASSERT_DELTA(z.y(), back_z.y(), 1e-10);
    TS_ASSERT_DELTA(z.z(), back_z.z(), 1e-10);
  }
  
  void test_axis_angle_to_from_matrix() {
    double matrix_raw[][9] = { 
      { 1, 0, 0, 0, -1, 0, 0, 0, -1 },
      { -1, 0, 0, 0, 1, 0, 0, 0, -1 },
      { -1, 0, 0, 0, -1, 0, 0, 0, 1 },
      { 0.766460869268653, -0.6332673325641, -0.107285699825701, 0.54750932113484, 0.556843741813312, 0.624626760933132, -0.335814352210275, -0.537491890832317, 0.773518705746082 }, 
      { 0.959857331298381, 0.224736332046459,  0.167831715152572, -0.213344180548401, 0.973418345385044, -0.0833125770455766, -0.182093833445511, 0.0441622681273, 0.982288923838079 }, 
      { 0.986507883482759, 0.131321875125383, -0.0977586872861997, -0.162531630378502, 0.857235846002895, -0.488600218434465, 0.0196383541284319, 0.497896846185286, 0.867013878554358 },
      { 0.721137182192153, -0.639876267870011, -0.265555128507638, 0.692124442131289, 0.648587293563264, 0.316699035730664, -0.0304125149307545, -0.412180645403235, 0.910594418218429 },
      { 0.764713536132754, 0.19163445642993, 0.615214956551078, -0.378788528301041, 0.906061660893168, 0.188604129029681, -0.521279635591945, -0.377264498420297, 0.765466550378794 }
    };

    double vector_raw[][3] = {
      { M_PI, 0, 0},
      { 0, M_PI, 0},
      { 0, 0, M_PI},
      { -0.688156207637604, 0.135324745474648, 0.699204666272306 },
      { 0.0646497978603147, 0.177467295232607, -0.222175728453655 },
      { 0.518493243020897, -0.0617027408465663, -0.154446538673176 },
      { -0.415601155986352, -0.134076370289505, 0.759495770262097 },
      { -0.312930331688783, 0.628491512746504, -0.315448931352088 },
    };
    for (int i = 0; i < sizeof(matrix_raw)/sizeof(matrix_raw[1]); i++) {
      vw::Vector3 v(vector_raw[i]);
      vw::Matrix3x3 m(matrix_raw[i]);

      vw::Vector3 diff_v(matrix_to_axis_angle(m) - v);
      TS_ASSERT_DELTA(diff_v.x(), 0, 1e-14);
      TS_ASSERT_DELTA(diff_v.y(), 0, 1e-14);
      TS_ASSERT_DELTA(diff_v.z(), 0, 1e-14);
      
      vw::Matrix3x3 mat_v(axis_angle_to_matrix(v) - m);
      TS_ASSERT_DELTA(mat_v(0, 0), 0, 1e-14);
      TS_ASSERT_DELTA(mat_v(0, 1), 0, 1e-14);
      TS_ASSERT_DELTA(mat_v(0, 2), 0, 1e-14);
      TS_ASSERT_DELTA(mat_v(1, 0), 0, 1e-14);
      TS_ASSERT_DELTA(mat_v(1, 1), 0, 1e-14);
      TS_ASSERT_DELTA(mat_v(1, 2), 0, 1e-14);
      TS_ASSERT_DELTA(mat_v(2, 0), 0, 1e-14);
      TS_ASSERT_DELTA(mat_v(2, 1), 0, 1e-14);
      TS_ASSERT_DELTA(mat_v(2, 2), 0, 1e-14);
    }
  }

  void test_axis_angle_quaternion_and_matrix_random() {
    std::srand(std::time(0));

    for (int i = 0; i < 1000; i++) {
      vw::Vector3 v(std::rand(), std::rand(), std::rand());
      v /= vw::math::norm_2(v); // normalize
      v *= 2.0 * M_PI * std::rand() / RAND_MAX; // rotation angle between 0 and 2 pi

      Quaternion<double> q(axis_angle_to_quaternion(v));

      // conversion back to axis_angle should give original vector
      vw::Vector3 axis;
      double angle;
      q.axis_angle(axis, angle);
      axis *= angle;
      TS_ASSERT_DELTA(axis.x(), v.x(), 1e-14);
      TS_ASSERT_DELTA(axis.y(), v.y(), 1e-14);
      TS_ASSERT_DELTA(axis.z(), v.z(), 1e-14);

      // conversion back to axis_angle should give original vector
      TS_ASSERT_DELTA(q.axis_angle().x(), v.x(), 1e-14);
      TS_ASSERT_DELTA(q.axis_angle().y(), v.y(), 1e-14);
      TS_ASSERT_DELTA(q.axis_angle().z(), v.z(), 1e-14);    

      vw::Matrix3x3 m( axis_angle_to_matrix(v) );
      vw::Matrix3x3 m_q(q.rotation_matrix());
    
      // test if m created directly from axis_angle and m_q created via quaternion are equal
      TS_ASSERT_DELTA(m(0,0), m_q(0,0), 1e-14);    
      TS_ASSERT_DELTA(m(0,1), m_q(0,1), 1e-14);    
      TS_ASSERT_DELTA(m(0,2), m_q(0,2), 1e-14);    
      TS_ASSERT_DELTA(m(1,0), m_q(1,0), 1e-14);    
      TS_ASSERT_DELTA(m(1,1), m_q(1,1), 1e-14);    
      TS_ASSERT_DELTA(m(1,2), m_q(1,2), 1e-14);    
      TS_ASSERT_DELTA(m(2,0), m_q(2,0), 1e-14);    
      TS_ASSERT_DELTA(m(2,1), m_q(2,1), 1e-14);    
      TS_ASSERT_DELTA(m(2,2), m_q(2,2), 1e-14);    
    
      // test if m is unitary
      vw::Matrix3x3 prod(m*vw::math::transpose(m));
      TS_ASSERT_DELTA(prod(0,0), 1, 1e-14);    
      TS_ASSERT_DELTA(prod(0,1), 0, 1e-14);    
      TS_ASSERT_DELTA(prod(0,2), 0, 1e-14);    
      TS_ASSERT_DELTA(prod(1,0), 0, 1e-14);    
      TS_ASSERT_DELTA(prod(1,1), 1, 1e-14);    
      TS_ASSERT_DELTA(prod(1,2), 0, 1e-14);    
      TS_ASSERT_DELTA(prod(2,0), 0, 1e-14);    
      TS_ASSERT_DELTA(prod(2,1), 0, 1e-14);    
      TS_ASSERT_DELTA(prod(2,2), 1, 1e-14);    

      // test if m is rotation matrix (and not reflection)
      TS_ASSERT_DELTA(det(m), 1, 1e-14);    
      
      // test if we can convert m back to a rotation vector
      vw::Vector3 matrix_back(matrix_to_axis_angle(m));
      TS_ASSERT_DELTA(matrix_back.x(), v.x(), 1e-14);
      TS_ASSERT_DELTA(matrix_back.y(), v.y(), 1e-14);
      TS_ASSERT_DELTA(matrix_back.z(), v.z(), 1e-14);    

   
      // Quaternion and matrix rotation should give equal results
      vw::Vector3 test(std::rand(), std::rand(), std::rand());   
      test /= RAND_MAX;
      vw::Vector3 diff (m*test - q.rotate(test));
      TS_ASSERT_DELTA(diff.x(), 0, 1e-14);
      TS_ASSERT_DELTA(diff.y(), 0, 1e-14);
      TS_ASSERT_DELTA(diff.z(), 0, 1e-14);

      // rotation in one direction by (transposed) matrix and other by standard quaternion should give original vector
      vw::Vector3 null_op(vw::math::transpose(m)*q.rotate(test));
      TS_ASSERT_DELTA(null_op.x(), test.x(), 1e-14);
      TS_ASSERT_DELTA(null_op.y(), test.y(), 1e-14);
      TS_ASSERT_DELTA(null_op.z(), test.z(), 1e-14);
    }    
  }
  
}; // class TestQuaternion
