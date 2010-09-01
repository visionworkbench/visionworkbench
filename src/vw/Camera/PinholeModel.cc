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
using google::protobuf::RepeatedFieldBackInserter;
#endif

#include <boost/filesystem/convenience.hpp>
namespace fs = boost::filesystem;

// Reads in a file containing parameters of a pinhole model with
// a tsai lens distortion model.
void vw::camera::PinholeModel::read_file(std::string const& filename){
  this->read(filename);
}
void vw::camera::PinholeModel::read(std::string const& filename) {
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

  typedef VectorProxy<double,3> Vector3P;
  m_u_direction = Vector3P(file.mutable_u_direction()->mutable_data());
  m_v_direction = Vector3P(file.mutable_v_direction()->mutable_data());
  m_w_direction = Vector3P(file.mutable_w_direction()->mutable_data());
  m_camera_center = Vector3P(file.mutable_camera_center()->mutable_data());
  m_fu = file.focal_length(0);
  m_fv = file.focal_length(1);
  m_cu = file.center_point(0);
  m_cv = file.center_point(1);
  m_rotation = MatrixProxy<double,3,3>(file.mutable_camera_rotation()->mutable_data());
  m_pixel_pitch = file.pixel_pitch();

  this->rebuild_camera_matrix();

  if ( file.distortion_name() == "NULL" ) {
    VW_ASSERT( file.distortion_vector_size() == 0,
               IOErr() << "Pinhole::read_file: Unexpected distortion vector." );
    m_distortion.reset( new NullLensDistortion());
  } else if ( file.distortion_name() == "TSAI" ) {
    VW_ASSERT( file.distortion_vector_size() == 4,
               IOErr() << "Pinhole::read_file: Unexpected distortion vector." );
    m_distortion.reset( new TsaiLensDistortion(VectorProxy<double,4>(file.mutable_distortion_vector()->mutable_data())));
  } else if ( file.distortion_name() == "BROWNCONRADY" ) {
    VW_ASSERT( file.distortion_vector_size() == 8,
               IOErr() << "Pinhole::read_file: Unexpected distortion vector." );
    m_distortion.reset( new BrownConradyDistortion(VectorProxy<double,8>(file.mutable_distortion_vector()->mutable_data())));
  } else if ( file.distortion_name() == "AdjustableTSAI" ) {
    VW_ASSERT( file.distortion_vector_size() > 3,
               IOErr() << "Pinhole::read_file: Unexpected distortion vector." );
    m_distortion.reset( new AdjustableTsaiLensDistortion(VectorProxy<double>(file.distortion_vector_size(),file.mutable_distortion_vector()->mutable_data())));
  }
#else
  // If you hit this point, you need to install Google Protobuffers to
  // be in order to write.
  vw_throw( IOErr() << "Pinhole::read_file: Camera IO not supported without Google Protobuffers" );
#endif
}


// Write parameters of an exiting PinholeModel into a .tsai file for later use.
void vw::camera::PinholeModel::write_file(std::string const& filename) const {
  write(filename);
}
void vw::camera::PinholeModel::write(std::string const& filename) const {
#if defined(VW_HAVE_PKG_PROTOBUF) && VW_HAVE_PKG_PROTOBUF==1
  std::string output_file =
    fs::path(filename).replace_extension(".pinhole").string();

  TsaiFile file;
  file.add_focal_length( m_fu );
  file.add_focal_length( m_fv );
  file.add_center_point( m_cu );
  file.add_center_point( m_cv );
  file.set_pixel_pitch( m_pixel_pitch );

  std::copy(m_u_direction.begin(), m_u_direction.end(),
            RepeatedFieldBackInserter(file.mutable_u_direction()));
  std::copy(m_v_direction.begin(), m_v_direction.end(),
            RepeatedFieldBackInserter(file.mutable_v_direction()));
  std::copy(m_w_direction.begin(), m_w_direction.end(),
            RepeatedFieldBackInserter(file.mutable_w_direction()));
  std::copy(m_camera_center.begin(), m_camera_center.end(),
            RepeatedFieldBackInserter(file.mutable_camera_center()));

  std::copy(m_rotation.begin(), m_rotation.end(),
            RepeatedFieldBackInserter(file.mutable_camera_rotation()));

  file.set_distortion_name( m_distortion->name() );
  Vector<double> distort_vec = m_distortion->distortion_parameters();
  std::copy(distort_vec.begin(),distort_vec.end(),
            RepeatedFieldBackInserter(file.mutable_distortion_vector()));

  std::ofstream output(output_file.c_str());
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
  double denominator = m_camera_matrix(2,0)*point(0) + m_camera_matrix(2,1)*point(1) +
    m_camera_matrix(2,2)*point(2) + m_camera_matrix(2,3);
  Vector2 pixel = Vector2( (m_camera_matrix(0,0)*point(0) + m_camera_matrix(0,1)*point(1) +
                            m_camera_matrix(0,2)*point(2) + m_camera_matrix(0,3)) / denominator,
                           (m_camera_matrix(1,0)*point(0) + m_camera_matrix(1,1)*point(1) +
                            m_camera_matrix(1,2)*point(2) + m_camera_matrix(1,3)) / denominator);

  //  Apply the lens distortion model
  return m_distortion->distorted_coordinates(*this, pixel)/m_pixel_pitch;
}

bool vw::camera::PinholeModel::projection_valid(Vector3 const& point) const {
  // z coordinate after extrinsic transformation
  double z = m_extrinsics(2, 0)*point(0) + m_extrinsics(2, 1)*point(1) +
    m_extrinsics(2, 2)*point(2) + m_extrinsics(2,3);
  return z > 0;
}

vw::Vector3 vw::camera::PinholeModel::pixel_to_vector (vw::Vector2 const& pix) const {
  // Apply the inverse lens distortion model
  vw::Vector2 undistorted_pix = m_distortion->undistorted_coordinates(*this, pix*m_pixel_pitch);

  // Compute the direction of the ray emanating from the camera center.
  vw::Vector3 p(0,0,1);
  subvector(p,0,2) = undistorted_pix;
  return normalize( m_inv_camera_transform * p);
}

void vw::camera::PinholeModel::intrinsic_parameters(double& f_u, double& f_v,
                                                    double& c_u, double& c_v) const {
  f_u = m_fu;  f_v = m_fv;  c_u = m_cu;  c_v = m_cv;
}

void vw::camera::PinholeModel::set_intrinsic_parameters(double f_u, double f_v,
                                                        double c_u, double c_v) {
  m_fu = f_u;  m_fv = f_v;  m_cu = c_u;  m_cv = c_v;
  rebuild_camera_matrix();
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
  Vector2 focal = camera_model.focal_length();
  Vector2 offset = camera_model.point_offset();
  focal *= scale;
  offset *= scale;
  boost::shared_ptr<LensDistortion> lens = camera_model.lens_distortion()->copy();
  lens->scale( scale );
  return vw::camera::PinholeModel( camera_model.camera_center(),
                                   camera_model.camera_pose().rotation_matrix(),
                                   focal[0], focal[1], offset[0], offset[1],
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
  Vector2 focal = camera_model.focal_length();
  Vector2 offset = camera_model.point_offset();
  vw::camera::NullLensDistortion distortion;
  return vw::camera::PinholeModel(camera_model.camera_center(),
      camera_model.camera_pose().rotation_matrix(),
      focal[0], focal[1], offset[0], offset[1],
      camera_model.coordinate_frame_u_direction(),
      camera_model.coordinate_frame_v_direction(),
      camera_model.coordinate_frame_w_direction(),
      distortion);
}

std::ostream& vw::camera::operator<<(std::ostream& str, vw::camera::PinholeModel const& model) {
  str << "Pinhole camera: \n";
  str << "\tCamera Center: " << model.camera_center() << "\n";
  str << "\tRotation Matrix: " << model.camera_pose() << "\n";
  str << "\tIntrinsics:\n";
  str << "\t  focal: " << model.focal_length() << "\n";
  str << "\t  offset: " << model.point_offset() << "\n";
  str << "\tDistortion Model: " << model.lens_distortion()->name() << "\n";
  str << "\t  " << *model.lens_distortion() << "\n";

  return str;
}
