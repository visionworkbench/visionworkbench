// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/Core/Log.h>
#include <vw/Camera/PinholeModel.h>
#include <vw/Math/EulerAngles.h>

#if defined(VW_HAVE_PKG_LAPACK) && VW_HAVE_PKG_LAPACK==1
#include <vw/Math/LinearAlgebra.h>
#endif

// For std::setprecision
#include <iomanip>

#if defined(VW_HAVE_PKG_PROTOBUF) && VW_HAVE_PKG_PROTOBUF==1
#include <vw/Camera/TsaiFile.pb.h>
#endif

// Reads in a file containing parameters of a pinhole model with
// a tsai lens distortion model.
void vw::camera::PinholeModel::read_file(std::string const& filename) {
#if defined(VW_HAVE_PKG_PROTOBUF) && VW_HAVE_PKG_PROTOBUF==1
  std::fstream input( filename.c_str(), std::ios::in | std::ios::binary );
  if ( !input )
    vw_throw( IOErr() << "Pinhole::read_file: Could not open " << filename << "\n" );
  TsaiFile file;
  file.ParseFromIstream( &input );
  input.close();

  // Making sure protobuf seems correct
  VW_ASSERT( file.focal_length_size() == 2,
             IOErr() << "Pinhole::read_file: Unexpected amount of focal lengths." );
  VW_ASSERT( file.center_point_size() == 2,
             IOErr() << "Pinhole::read_file: Unexpected amount of center points." );
  VW_ASSERT( file.u_direction_size() == 3,
             IOErr() << "Pinhole::read_file: Unexpected size of u vector." );
  VW_ASSERT( file.v_direction_size() == 3,
             IOErr() << "Pinhole::read_file: Unexpected size of v vector." );
  VW_ASSERT( file.w_direction_size() == 3,
             IOErr() << "Pinhole::read_file: Unexpected size of w vector." );
  VW_ASSERT( file.camera_center_size() == 3,
             IOErr() << "Pinhole::read_file: Unexpected size of camera vector." );
  VW_ASSERT( file.camera_rotation_size() == 9,
             IOErr() << "Pinhole::read_file: Unexpected size of rotation matrix." );

  for ( char i = 0; i < 3; i++ ) {
    m_u_direction[i] = file.u_direction(i);
    m_v_direction[i] = file.v_direction(i);
    m_w_direction[i] = file.w_direction(i);
    m_camera_center[i] = file.camera_center(i);
  }
  m_fu = file.focal_length(0);
  m_fv = file.focal_length(1);
  m_cu = file.center_point(0);
  m_cv = file.center_point(1);

  {
    Matrix<double>::iterator iter = m_rotation.begin();
    for ( char i = 0; i < 9; i++ ) {
      *iter = file.camera_rotation(i);
      iter++;
    }
  }

  this->rebuild_camera_matrix();

  Vector<double> distortion_params( file.distortion_vector_size() );
  for ( int i = 0; i < file.distortion_vector_size(); i++ )
    distortion_params[i] = file.distortion_vector(i);
  if ( file.distortion_name() == "NULL" ) {
    VW_ASSERT( distortion_params.size() == 0,
               IOErr() << "Pinhole::read_file: Unexpected distortion vector." );
    m_distortion.reset( new NullLensDistortion());
  } else if ( file.distortion_name() == "TSAI" ) {
    VW_ASSERT( distortion_params.size() == 4,
               IOErr() << "Pinhole::read_file: Unexpected distortion vector." );
    m_distortion.reset( new TsaiLensDistortion(distortion_params));
  } else if ( file.distortion_name() == "BROWNCONRADY" ) {
    VW_ASSERT( distortion_params.size() == 8,
               IOErr() << "Pinhole::read_file: Unexpected distortion vector." );
    m_distortion.reset( new BrownConradyDistortion(distortion_params));
  }
#else
  // If you hit this point, you need to install Google Protobuffers to
  // be in order to write.
  vw_throw( IOErr() << "Pinhole::read_file: Camera IO not supported without Google Protobuffers" );
#endif
}


// Write parameters of an exiting PinholeModel into a .tsai file for later use.
void vw::camera::PinholeModel::write_file(std::string const& filename) const {
#if defined(VW_HAVE_PKG_PROTOBUF) && VW_HAVE_PKG_PROTOBUF==1
  TsaiFile file;
  file.add_focal_length( m_fu );
  file.add_focal_length( m_fv );
  file.add_center_point( m_cu );
  file.add_center_point( m_cv );

  for ( char i = 0; i < 3; i++ ) {
    file.add_u_direction(m_u_direction[i]);
    file.add_v_direction(m_v_direction[i]);
    file.add_w_direction(m_w_direction[i]);
    file.add_camera_center(m_camera_center[i]);
  }

  for ( Matrix<double>::const_iterator iter = m_rotation.begin();
        iter != m_rotation.end(); iter++ ) {
    file.add_camera_rotation(*iter);
  }

  file.set_distortion_name( m_distortion->name() );
  Vector<double> distort_vec = m_distortion->distortion_parameters();
  for ( unsigned i = 0; i < distort_vec.size(); i++ )
    file.add_distortion_vector(distort_vec(i));

  std::ofstream output(filename.c_str());
  if( !output.is_open() )
    vw_throw( IOErr() << "PinholeModel::write_file: Could not open file\n" );
  file.SerializeToOstream( &output );
  output.close();
#else
  // If you hit this point, you need to install Google Protobuffers to
  // be in order to write.
  vw_throw( IOErr() << "Pinhole::write_file: Camera IO not supported without Google Protobuffers" );
#endif
}

vw::Vector2 vw::camera::PinholeModel::point_to_pixel(vw::Vector3 const& point) const {

  //  Multiply the pixel location by the camera matrix.
  double denominator = m_camera_matrix(2,0)*point(0) + m_camera_matrix(2,1)*point(1) + m_camera_matrix(2,2)*point(2) + m_camera_matrix(2,3);
  Vector2 pixel = Vector2( (m_camera_matrix(0,0)*point(0) + m_camera_matrix(0,1)*point(1) + m_camera_matrix(0,2)*point(2) + m_camera_matrix(0,3)) / denominator,
      (m_camera_matrix(1,0)*point(0) + m_camera_matrix(1,1)*point(1) + m_camera_matrix(1,2)*point(2) + m_camera_matrix(1,3)) / denominator);

  //  Apply the lens distortion model
  return m_distortion->distorted_coordinates(*this, pixel);
}

bool vw::camera::PinholeModel::projection_valid(Vector3 const& point) const {
  // z coordinate after extrinsic transformation
  double z = m_extrinsics(2, 0)*point(0) + m_extrinsics(2, 1)*point(1) + m_extrinsics(2, 2)*point(2) + m_extrinsics(2,3);
  return z > 0;
}

vw::Vector3 vw::camera::PinholeModel::pixel_to_vector (vw::Vector2 const& pix) const {

  // Apply the inverse lens distortion model
  vw::Vector2 undistorted_pix = m_distortion->undistorted_coordinates(*this, pix);

  // Compute the direction of the ray emanating from the camera center.
  vw::Vector3 p(0,0,1);
  subvector(p,0,2) = undistorted_pix;
  return normalize( m_inv_camera_transform * p);
}

void vw::camera::PinholeModel::set_camera_matrix( Matrix<double,3,4> const& p ) {
#if defined(VW_HAVE_PKG_LAPACK) && VW_HAVE_PKG_LAPACK==1
  // Solving for camera center
  Matrix<double> cam_nullsp = nullspace(p);
  Vector<double> cam_center = select_col(cam_nullsp,0);
  cam_center /= cam_center[3];
  m_camera_center = subvector(cam_center,0,3);

  // Solving for intrinsics with RQ decomposition
  Matrix<double> M = submatrix(p,0,0,3,3);
  Matrix<double> R,Q;
  rqd( M, R, Q );
  Matrix<double> sign_fix(3,3);
  sign_fix.set_identity();
  if ( R(0,0) < 0 )
    sign_fix(0,0) = -1;
  if ( R(1,1) < 0 )
    sign_fix(1,1) = -1;
  if ( R(2,2) < 0 )
    sign_fix(2,2) = -1;
  R = R*sign_fix;
  Q = sign_fix*Q;
  R /= R(2,2);

  // Pulling out intrinsic and last extrinsic
  Matrix<double,3,3> uvwRotation;
  select_row(uvwRotation,0) = m_u_direction;
  select_row(uvwRotation,1) = m_v_direction;
  select_row(uvwRotation,2) = m_w_direction;
  m_rotation = inverse(uvwRotation*Q);
  m_fu = R(0,0);
  m_fv = R(1,1);
  m_cu = R(0,2);
  m_cv = R(1,2);

  if ( fabs(R(0,1)) >= 1.2 )
    vw_out(vw::WarningMessage,"camera") << "Significant skew not modelled by pinhole camera\n";

  // Rebuild
  rebuild_camera_matrix();
#else
  vw_throw( NoImplErr() << "PinholeModel::set_Camera_Matrix is unavailable without LAPACK" );
#endif
}

void vw::camera::PinholeModel::rebuild_camera_matrix() {

  /// The intrinsic portion of the camera matrix is stored as
  ///
  ///    [  fx   0   cx  ]
  /// K= [  0    fy  cy  ]
  ///    [  0    0   1   ]
  ///
  /// with fx, fy the focal length of the system (in horizontal and
  /// vertical pixels), and (cx, cy) the pixel coordinates of the
  /// central pixel (the principal point on the image plane).

  m_intrinsics(0,0) = m_fu;
  m_intrinsics(0,1) = 0;
  m_intrinsics(0,2) = m_cu;
  m_intrinsics(1,0) = 0;
  m_intrinsics(1,1) = m_fv;
  m_intrinsics(1,2) = m_cv;
  m_intrinsics(2,0) = 0;
  m_intrinsics(2,1) = 0;
  m_intrinsics(2,2) = 1;

  // The extrinsics are normally built as the matrix:  [ R | -R*C ].
  // To allow for user-specified coordinate frames, the
  // extrinsics are now build to include the u,v,w rotation
  //
  //               | u_0  u_1  u_2  |
  //     Extr. =   | v_0  v_1  v_2  | * [ R | -R*C]
  //               | w_0  w_1  w_2  |
  //
  // The vectors u,v, and w must be orthonormal.

  /*   check for orthonormality of u,v,w              */
  VW_LINE_ASSERT( dot_prod(m_u_direction, m_v_direction) == 0 );
  VW_LINE_ASSERT( dot_prod(m_u_direction, m_w_direction) == 0 );
  VW_LINE_ASSERT( dot_prod(m_v_direction, m_w_direction) == 0 );
  VW_LINE_ASSERT( fabs( norm_2(m_u_direction) - 1 ) < 0.001 );
  VW_LINE_ASSERT( fabs( norm_2(m_v_direction) - 1 ) < 0.001 );
  VW_LINE_ASSERT( fabs( norm_2(m_w_direction) - 1 ) < 0.001 );

  Matrix<double,3,3> uvwRotation;

  select_row(uvwRotation,0) = m_u_direction;
  select_row(uvwRotation,1) = m_v_direction;
  select_row(uvwRotation,2) = m_w_direction;

  Matrix<double,3,3> m_rotation_inverse = transpose(m_rotation);
  submatrix(m_extrinsics,0,0,3,3) = uvwRotation * m_rotation_inverse;
  select_col(m_extrinsics,3) = uvwRotation * -m_rotation_inverse * m_camera_center;

  m_camera_matrix = m_intrinsics * m_extrinsics;
  m_inv_camera_transform = inverse(uvwRotation*m_rotation_inverse) * inverse(m_intrinsics);
}

// scale_camera
//  Used to modify camera in the event to user resizes the image
vw::camera::PinholeModel vw::camera::scale_camera(vw::camera::PinholeModel const& camera_model, float const& scale) {
  double fu, fv, cu, cv;
  camera_model.intrinsic_parameters(fu, fv, cu, cv);
  fu *= scale;
  fv *= scale;
  cu *= scale;
  cv *= scale;
  boost::shared_ptr<LensDistortion> lens = camera_model.lens_distortion()->copy();
  lens->scale( scale );
  return vw::camera::PinholeModel( camera_model.camera_center(),
                                   camera_model.camera_pose().rotation_matrix(),
                                   fu, fv, cu, cv,
                                   camera_model.coordinate_frame_u_direction(),
                                   camera_model.coordinate_frame_v_direction(),
                                   camera_model.coordinate_frame_w_direction(),
                                   *lens );
}

//   /// Given two pinhole camera models, this method returns two new camera
//   /// models that have been epipolar rectified.
//   template <>
//   void epipolar(PinholeModel<NoLensDistortion> const& src_camera0,
//                 PinholeModel<NoLensDistortion> const& src_camera1,
//                 PinholeModel<NoLensDistortion> &dst_camera0,
//                 PinholeModel<NoLensDistortion> &dst_camera1);

vw::camera::PinholeModel vw::camera::linearize_camera(vw::camera::PinholeModel const& camera_model) {
  double fu, fv, cu, cv;
  camera_model.intrinsic_parameters(fu, fv, cu, cv);
  vw::camera::NullLensDistortion distortion;
  return vw::camera::PinholeModel(camera_model.camera_center(),
      camera_model.camera_pose().rotation_matrix(),
      fu, fv, cu, cv,
      camera_model.coordinate_frame_u_direction(),
      camera_model.coordinate_frame_v_direction(),
      camera_model.coordinate_frame_w_direction(),
      distortion);
}

std::ostream& vw::camera::operator<<(std::ostream& str, vw::camera::PinholeModel const& model) {
  double fu, fv, cu, cv;
  model.intrinsic_parameters(fu, fv, cu, cv);

  str << "Pinhole camera: \n";
  str << "\tCamera Center: " << model.camera_center() << "\n";
  str << "\tRotation Matrix: " << model.camera_pose() << "\n";
  str << "\tIntrinsics:\n";
  str << "\t  f_u: " << fu << "    f_v: " << fv << "\n";
  str << "\t  c_u: " << cu << "    c_v: " << cv << "\n";
  str << model.lens_distortion() << "\n";

  return str;
}
