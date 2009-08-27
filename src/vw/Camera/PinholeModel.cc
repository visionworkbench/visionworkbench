// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Camera/PinholeModel.h>
#include <vw/Math/EulerAngles.h>

// For std::setprecision
#include <iomanip>

// Reads in a file containing parameters of a pinhole model with
// a tsai lens distortion model. An example is provided at the end of this file.
void vw::camera::PinholeModel::read_file(std::string const& filename) {

  char line[2048];
  double fu, fv, cu, cv;
  Vector3 u_direction, v_direction, w_direction;
  Vector3 C;
  Matrix3x3 R;
  Vector4 distortion_params(0,0,0,0);


  FILE *cam_file = fopen(filename.c_str(), "r");
  if (cam_file == 0) vw_throw( IOErr() << "PinholeModel::read_file: Could not open file\n" );

  // Read intrinsic parameters
  fgets(line, sizeof(line), cam_file);
  if (sscanf(line,"fu = %lf", &fu) != 1) {
    fclose(cam_file);
    vw_throw( IOErr() << "PinholeModel::read_file(): Could not read x focal length\n" );
  }

  fgets(line, sizeof(line), cam_file);
  if (sscanf(line,"fv = %lf", &fv) != 1) {
    fclose(cam_file);
    vw_throw( IOErr() << "PinholeModel::read_file(): Could not read y focal length\n" );
  }

  fgets(line, sizeof(line), cam_file);
  if (sscanf(line,"cu = %lf", &cu) != 1) {
    fclose(cam_file);
    vw_throw( IOErr() << "PinholeModel::read_file(): Could not read x principal point\n" );
  }

  fgets(line, sizeof(line), cam_file);
  if (sscanf(line,"cv = %lf", &cv) != 1) {
    fclose(cam_file);
    vw_throw( IOErr() << "PinholeModel::read_file(): Could not read y principal point\n" );
  }

  fgets(line, sizeof(line), cam_file);
  if (sscanf(line,"u_direction = %lf %lf %lf", &u_direction(0), &u_direction(1), &u_direction(2)) != 3) {
    fclose(cam_file);
    vw_throw( IOErr() << "PinholeModel::read_file(): Could not read u direction vector\n" );
  }

  fgets(line, sizeof(line), cam_file);
  if (sscanf(line,"v_direction = %lf %lf %lf", &v_direction(0), &v_direction(1), &v_direction(2)) != 3) {
    fclose(cam_file);
    vw_throw( IOErr() << "PinholeModel::read_file(): Could not read v direction vector\n" );
  }

  fgets(line, sizeof(line), cam_file);
  if (sscanf(line,"w_direction = %lf %lf %lf", &w_direction(0), &w_direction(1), &w_direction(2)) != 3) {
    fclose(cam_file);
    vw_throw( IOErr() << "PinholeModel::read_file(): Could not read w direction vector\n" );
  }

  // Read extrinsic parameters
  fgets(line, sizeof(line), cam_file);
  if (sscanf(line,"C = %lf %lf %lf", &C(0), &C(1), &C(2)) != 3) {
    fclose(cam_file);
    vw_throw( IOErr() << "PinholeModel::read_file: Could not read C (camera center) vector\n" );
  }

  fgets(line, sizeof(line), cam_file);
  if ( sscanf(line, "R = %lf %lf %lf %lf %lf %lf %lf %lf %lf",
              &R(0,0), &R(0,1), &R(0,2),
              &R(1,0), &R(1,1), &R(1,2),
              &R(2,0), &R(2,1), &R(2,2)) != 9 ) {
      fclose(cam_file);
      vw_throw( IOErr() << "PinholeModel::read_file(): Could not read rotation matrix\n" );
  }

  // Read distortion parameters.
   fgets(line, sizeof(line), cam_file);
  if (sscanf(line,"k1 = %lf", &distortion_params[0] ) != 1) {
    fclose(cam_file);
    vw_throw( IOErr() << "PinholeModel::read_file(): Could not read tsai distortion parameter k1\n" );
  }

  fgets(line, sizeof(line), cam_file);
  if (sscanf(line,"k2 = %lf", &distortion_params[1] ) != 1) {
    fclose(cam_file);
    vw_throw( IOErr() << "PinholeModel::read_file(): Could not read tsai distortion parameter k2\n" );
  }

  fgets(line, sizeof(line), cam_file);
  if (sscanf(line,"p1 = %lf", &distortion_params[2] ) != 1) {
    fclose(cam_file);
    vw_throw( IOErr() << "PinholeModel::read_file(): Could not read tsai distortion parameter p1\n" );
  }

  fgets(line, sizeof(line), cam_file);
  if (sscanf(line,"p2 = %lf", &distortion_params[3] ) != 1) {
    fclose(cam_file);
    vw_throw( IOErr() << "PinholeModel::read_file(): Could not read tsai distortion parameter p2\n" );
  }

  fclose(cam_file);

  m_u_direction = u_direction;
  m_v_direction = v_direction;
  m_w_direction = w_direction;

  m_fu = fu;
  m_fv = fv;
  m_cu = cu;
  m_cv = cv;
  m_camera_center = C;

  m_rotation = R;
  this->rebuild_camera_matrix();

  if( distortion_params == Vector4(0,0,0,0))
    m_distortion.reset(new NullLensDistortion());
  else
    m_distortion.reset(new TsaiLensDistortion(distortion_params));
}


//   Write parameters of an exiting PinholeModel into a .tsai file for later use.
//   FIXME: does not output distortion parameters

void vw::camera::PinholeModel::write_file(std::string const& filename) const {
  std::ofstream cam_file(filename.c_str());
  if( !cam_file.is_open() ) vw_throw( IOErr() << "PinholeModel::write_file: Could not open file\n" );

  cam_file << std::setprecision(std::numeric_limits<double>::digits10);
  cam_file << "fu = " << m_fu << "\n";
  cam_file << "fv = " << m_fv << "\n";
  cam_file << "cu = " << m_cu << "\n";
  cam_file << "cv = " << m_cv << "\n";
  cam_file << "u_direction = " << m_u_direction[0] << " " << m_u_direction[1] << " " << m_u_direction[2] << "\n";
  cam_file << "v_direction = " << m_v_direction[0] << " " << m_v_direction[1] << " " << m_v_direction[2] << "\n";
  cam_file << "w_direction = " << m_w_direction[0] << " " << m_w_direction[1] << " " << m_w_direction[2] << "\n";
  cam_file << "C = " << m_camera_center[0] << " " << m_camera_center[1] << " " << m_camera_center[2] << "\n";
  cam_file << "R = " << m_rotation(0,0) << " " << m_rotation(0,1) << " " << m_rotation(0,2) << " " << m_rotation(1,0) << " " << m_rotation(1,1) << " " << m_rotation(1,2) << " " << m_rotation(2,0) << " " << m_rotation(2,1) << " " << m_rotation(2,2) << "\n";

  //  Distortion Parameters. This should be implemented by overloading the
  //  << operator for distortion models

  cam_file << *m_distortion << "\n";
  cam_file << "\n" << "\n" << " Parameters for a Pinhole camera model with tsai lens distortion model." << "\n";
  cam_file.close();
}


/* Contents of a sample .tsai file:

fu = 611.651
fv = 610.216
cu = 500.829
cv = 396.22
u_direction = 1 0 0
v_direction = 0 1 0
w_direction = 0 0 1
C = -0.328711 -0.0637059 -0.828905
R = 0.000412095 -0.99998 0.00624732 0.409245 0.00586886 0.912405 -0.912424 0.00218069 0.40924
k1 = 0
k2 = 0
p1 = 0
p2 = 0

 Parameters for a Pinhole camera model with tsai lens distortion model.
*/

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
