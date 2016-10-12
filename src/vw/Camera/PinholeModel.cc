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

#include <vw/Core/Log.h>
#include <vw/config.h>
#include <vw/Math/Matrix.h>
#include <vw/Math/Vector.h>
#include <vw/Math/Quaternion.h>
#include <vw/Camera/PinholeModel.h>
#include <vw/Camera/LensDistortion.h>

#if defined(VW_HAVE_PKG_LAPACK) && VW_HAVE_PKG_LAPACK==1
#include <vw/Math/LinearAlgebra.h>
#endif

#include <algorithm>
#include <sstream>
#include <iomanip>
#include <string>

// TODO: One day remove this protobuf dependency.
#if defined(VW_HAVE_PKG_PROTOBUF) && VW_HAVE_PKG_PROTOBUF==1
#include <vw/Camera/TsaiFile.pb.h>
using google::protobuf::RepeatedFieldBackInserter;
#endif

#include <boost/filesystem/convenience.hpp>
namespace fs = boost::filesystem;

using namespace vw;
using namespace camera;

PinholeModel::PinholeModel() : m_distortion(DistortPtr(new NullLensDistortion)),
                               m_camera_center(Vector3(0,0,0)),
                               m_fu(1), m_fv(1), m_cu(0), m_cv(0),
                               m_u_direction(Vector3(1,0,0)),
                               m_v_direction(Vector3(0,1,0)),
                               m_w_direction(Vector3(0,0,1)), m_pixel_pitch(1) {

  m_rotation.set_identity();
  this->rebuild_camera_matrix();
}

/// Initialize from a file on disk.
PinholeModel::PinholeModel(std::string const& filename) : m_distortion(DistortPtr(new NullLensDistortion)) {
  read(filename);
}

PinholeModel::PinholeModel(Vector3 camera_center, Matrix<double,3,3> rotation,
                           double f_u, double f_v, double c_u, double c_v,
                           Vector3 u_direction, Vector3 v_direction,
                           Vector3 w_direction,
                           LensDistortion const& distortion_model,
                           double pixel_pitch) : m_distortion(DistortPtr(distortion_model.copy())),
                                                 m_camera_center(camera_center),  
                                                 m_rotation(rotation),
                                                 m_fu(f_u), m_fv(f_v), m_cu(c_u), m_cv(c_v),
                                                 m_u_direction(u_direction),
                                                 m_v_direction(v_direction),
                                                 m_w_direction(w_direction),
                                                 m_pixel_pitch(pixel_pitch) {
  this->rebuild_camera_matrix();
}

PinholeModel::PinholeModel(Vector3 camera_center, Matrix<double,3,3> rotation,
                           double f_u, double f_v, double c_u, double c_v,
                           LensDistortion const& distortion_model,
                           double pixel_pitch) : m_distortion(DistortPtr(distortion_model.copy())),
                                                 m_camera_center(camera_center),
                                                 m_rotation(rotation),
                                                 m_fu(f_u), m_fv(f_v), m_cu(c_u), m_cv(c_v),
                                                 m_u_direction(Vector3(1,0,0)),
                                                 m_v_direction(Vector3(0,1,0)),
                                                 m_w_direction(Vector3(0,0,1)),
                                                 m_pixel_pitch(pixel_pitch) {
  rebuild_camera_matrix();
}


/// Construct a basic pinhole model with no lens distortion
PinholeModel::PinholeModel(Vector3 camera_center, Matrix<double,3,3> rotation,
                           double f_u, double f_v,
                           double c_u, double c_v,
                           double pixel_pitch) : m_distortion(DistortPtr(new NullLensDistortion)),
                                                 m_camera_center(camera_center),
                                                 m_rotation(rotation),
                                                 m_fu(f_u), m_fv(f_v), m_cu(c_u), m_cv(c_v),
                                                 m_u_direction(Vector3(1,0,0)),
                                                 m_v_direction(Vector3(0,1,0)),
                                                 m_w_direction(Vector3(0,0,1)),
                                                 m_pixel_pitch(pixel_pitch) {
  rebuild_camera_matrix();
}



void PinholeModel::read(std::string const& filename) {

  // Handle deprecated .pinhole protobuf files
  fs::path filename_path( filename );
  if ( filename_path.extension() == ".pinhole" )
    return read_protobuf_file(filename);

  // Open the input file
  std::ifstream cam_file;
  cam_file.open(filename.c_str());
  if (cam_file.fail())
    vw_throw( IOErr() << "PinholeModel::read_file: Could not open file: " << filename );

  // Check for version number on the first line
  int file_version = 1; // The default version written before the 2016 changes
  std::string line;
  std::getline(cam_file, line);
  if (line.find("VERSION") != std::string::npos) {
    sscanf(line.c_str(),"VERSION_%d", &file_version); // Parse the version of the input file
    
    // Right now there is only one version (VERSION_3) so if we find the version
    //  we just skip it and move on to the next line.  If the version is changed,
    //  handler logic needs to be implemented here.
    std::getline(cam_file, line);
  }
  
  // Start parsing all the parameters from the lines.
  if (!cam_file.good() || sscanf(line.c_str(),"fu = %lf", &m_fu) != 1) {
    cam_file.close();
    vw_throw( IOErr() << "PinholeModel::read_file(): Could not read x focal length\n" );
  }

  std::getline(cam_file, line);
  if (!cam_file.good() || sscanf(line.c_str(),"fv = %lf", &m_fv) != 1) {
    cam_file.close();
    vw_throw( IOErr() << "PinholeModel::read_file(): Could not read y focal length\n" );
  }

  std::getline(cam_file, line);
  if (!cam_file.good() || sscanf(line.c_str(),"cu = %lf", &m_cu) != 1) {
    cam_file.close();
    vw_throw( IOErr() << "PinholeModel::read_file(): Could not read x principal point\n" );
  }

  std::getline(cam_file, line);
  if (!cam_file.good() || sscanf(line.c_str(),"cv = %lf", &m_cv) != 1) {
    cam_file.close();
    vw_throw( IOErr() << "PinholeModel::read_file(): Could not read y principal point\n" );
  }

  std::getline(cam_file, line);
  if (!cam_file.good() || sscanf(line.c_str(),"u_direction = %lf %lf %lf", 
        &m_u_direction(0), &m_u_direction(1), &m_u_direction(2)) != 3) {
    cam_file.close();
    vw_throw( IOErr() << "PinholeModel::read_file(): Could not read u direction vector\n" );
  }

  std::getline(cam_file, line);
  if (!cam_file.good() || sscanf(line.c_str(),"v_direction = %lf %lf %lf", 
        &m_v_direction(0), &m_v_direction(1), &m_v_direction(2)) != 3) {
    cam_file.close();
    vw_throw( IOErr() << "PinholeModel::read_file(): Could not read v direction vector\n" );
  }

  std::getline(cam_file, line);
  if (!cam_file.good() || sscanf(line.c_str(),"w_direction = %lf %lf %lf", 
        &m_w_direction(0), &m_w_direction(1), &m_w_direction(2)) != 3) {
    cam_file.close();
    vw_throw( IOErr() << "PinholeModel::read_file(): Could not read w direction vector\n" );
  }

  // Read extrinsic parameters
  std::getline(cam_file, line);
  if (!cam_file.good() || sscanf(line.c_str(),"C = %lf %lf %lf", 
        &m_camera_center(0), &m_camera_center(1), &m_camera_center(2)) != 3) {
    cam_file.close();
    vw_throw( IOErr() << "PinholeModel::read_file: Could not read C (camera center) vector\n" );
  }

  std::getline(cam_file, line);
  if ( !cam_file.good() ||
       sscanf(line.c_str(), "R = %lf %lf %lf %lf %lf %lf %lf %lf %lf",
              &m_rotation(0,0), &m_rotation(0,1), &m_rotation(0,2),
              &m_rotation(1,0), &m_rotation(1,1), &m_rotation(1,2),
              &m_rotation(2,0), &m_rotation(2,1), &m_rotation(2,2)) != 9 ) {
    cam_file.close();
    vw_throw( IOErr() << "PinholeModel::read_file(): Could not read rotation matrix\n" );
  }

  // The pitch line does not exist in older files, so don't fail if it is missing.
  // - At this point we need to hang on to the position immediately before the 
  int lens_start = cam_file.tellg();
  std::getline(cam_file, line);
  if (line.find("pitch") != std::string::npos) {
    if (!cam_file.good() || sscanf(line.c_str(),"pitch = %lf", &m_pixel_pitch) != 1) {
      cam_file.close();
      vw_throw( IOErr() << "PinholeModel::read_file(): Could not read pixel pitch\n" );
    }
    lens_start = cam_file.tellg();
    std::getline(cam_file, line); // After reading the pitch, read the next line.
  }
  else // Pixel pitch not specified, use 1.0 as the default.
  {
    if (file_version > 2){
      cam_file.close();
      vw_throw( IOErr() << "PinholeModel::read_file(): Pitch value required in this file version!\n" );
    }
    m_pixel_pitch = 1.0;
  }

  // Now that we have loaded all the parameters, update the dependend class members.
  this->rebuild_camera_matrix();

  // This creates m_distortion but we still need to read the parameters.
  bool found_name = construct_lens_distortion(line);

  if (!found_name && (file_version > 2)){
    cam_file.close();
    vw_throw( IOErr() << "PinholeModel::read_file(): Distortion name required in this file version!\n" );
  }
  
  // If there was no line containing the distortion model name (true for old files)
  //  then we need to back up to before the distortion parameters begin in the file.
  if (!found_name)
    cam_file.seekg(lens_start ,std::ios_base::beg);    

  // The lens distortion class knows how to parse the rest of the input stream.
  m_distortion->read(cam_file);

  cam_file.close();    
}

bool PinholeModel::construct_lens_distortion(std::string const& config_line) {

  // Check if the passed in string contains the string for any of the
  //  recognized lens distortion models.
  if (config_line.find(NullLensDistortion::class_name()) != std::string::npos) {
    m_distortion.reset(new NullLensDistortion());
    return true;
  }
  if (config_line.find(BrownConradyDistortion::class_name()) != std::string::npos) {
    m_distortion.reset(new BrownConradyDistortion());
    return true;
  }
  if (config_line.find(AdjustableTsaiLensDistortion::class_name()) != std::string::npos) {
    m_distortion.reset(new AdjustableTsaiLensDistortion());
    return true;
  }
  if (config_line.find(PhotometrixLensDistortion::class_name()) != std::string::npos) {
    m_distortion.reset(new PhotometrixLensDistortion());
    return true;
  }

  // TSAI is the default model.  Older files which do not have a specifier string
  //  contain TSAI parameters.
  m_distortion.reset(new TsaiLensDistortion());
  
  if (config_line.find(TsaiLensDistortion::class_name()) != std::string::npos)
    return true;
  else
    return false;
}


// This file type is DEPRECATED!
void PinholeModel::read_protobuf_file(std::string const& filename) {

  fs::path filename_path( filename );
  if ( filename_path.extension() != ".pinhole" )
    vw_throw( IOErr() << "Not a protobuf file extension: " << filename_path.extension() );
  
#if defined(VW_HAVE_PKG_PROTOBUF) && VW_HAVE_PKG_PROTOBUF==1
  std::fstream input( filename.c_str(), std::ios::in | std::ios::binary );
  if ( !input )
    vw_throw( IOErr() << "Pinhole::read_file: Could not open " << filename << "\n" );
  TsaiFile file;
  if ( !file.ParseFromIstream( &input ) )
    vw_throw( IOErr() << "Pinhole::read_file: Protocol buffer failed to parse \"" << filename << "\"\n" );
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
  // If you hit this point, you need to install Google Protobuffers to to read this file type.
  vw_throw( IOErr() << "Pinhole::write_file: Camera IO not supported without Google Protobuffers" );
#endif

}

void PinholeModel::write(std::string const& filename) const {

  // Update this field whenever there is a significant change to the file.
  // - It can be used to keep backwards compatibility with future plain text changes.
  const std::string PINHOLE_VERSION = "VERSION_3";

  // Set the path an open the output file for writing
  std::string file_path = fs::path(filename).replace_extension(".tsai").string();
  std::ofstream cam_file(file_path.c_str());
  if( !cam_file.is_open() ) 
    vw_throw( IOErr() << "PinholeModel::write: Could not open file: " << file_path );

  // Write the pinhole camera model parts
  //   # digits to survive double->text->double conversion
  const size_t ACCURATE_DIGITS = 17; // = std::numeric_limits<double>::max_digits10
  cam_file << std::setprecision(ACCURATE_DIGITS); 
  cam_file << PINHOLE_VERSION << "\n";
  cam_file << "fu = " << m_fu << "\n";
  cam_file << "fv = " << m_fv << "\n";
  cam_file << "cu = " << m_cu << "\n";
  cam_file << "cv = " << m_cv << "\n";
  cam_file << "u_direction = " << m_u_direction[0] << " " << m_u_direction[1] << " " << m_u_direction[2] << "\n";
  cam_file << "v_direction = " << m_v_direction[0] << " " << m_v_direction[1] << " " << m_v_direction[2] << "\n";
  cam_file << "w_direction = " << m_w_direction[0] << " " << m_w_direction[1] << " " << m_w_direction[2] << "\n";
  cam_file << "C = " << m_camera_center[0] << " " << m_camera_center[1] << " " << m_camera_center[2] << "\n";
  cam_file << "R = " << m_rotation(0,0) << " " << m_rotation(0,1) << " " << m_rotation(0,2) << " " << m_rotation(1,0) << " " << m_rotation(1,1) << " " << m_rotation(1,2) << " " << m_rotation(2,0) << " " << m_rotation(2,1) << " " << m_rotation(2,2) << "\n";
  cam_file << "pitch = " << m_pixel_pitch << "\n";

  // Write the name of the distortion model, then use the distortion model
  //  << overload to write it to the file. 
  cam_file << m_distortion->name() << std::endl;
  cam_file << *m_distortion;
  cam_file.close();
}


Vector2 PinholeModel::point_to_pixel(Vector3 const& point) const {

  // Multiply the pixel location by the 3x4 camera matrix.
  // - The pixel coordinate is de-homogenized by dividing by the denominator.
  double denominator = m_camera_matrix(2,0)*point(0) + m_camera_matrix(2,1)*point(1) +
                       m_camera_matrix(2,2)*point(2) + m_camera_matrix(2,3);
  Vector2 pixel = Vector2( (m_camera_matrix(0,0)*point(0) + m_camera_matrix(0,1)*point(1) +
                            m_camera_matrix(0,2)*point(2) + m_camera_matrix(0,3)           ) / denominator,
                           (m_camera_matrix(1,0)*point(0) + m_camera_matrix(1,1)*point(1) +
                            m_camera_matrix(1,2)*point(2) + m_camera_matrix(1,3)           ) / denominator);

  // Apply the lens distortion model
  // - Divide by pixel pitch to convert from metric units to pixels if the intrinsic
  //   values were not specified in pixel units (in that case m_pixel_pitch == 1.0)
  return m_distortion->distorted_coordinates(*this, pixel)/m_pixel_pitch;
}

Vector2 PinholeModel::point_to_pixel_no_distortion(Vector3 const& point) const {

  // Multiply the pixel location by the 3x4 camera matrix.
  // - The pixel coordinate is de-homogenized by dividing by the denominator.
  double denominator = m_camera_matrix(2,0)*point(0) + m_camera_matrix(2,1)*point(1) +
                       m_camera_matrix(2,2)*point(2) + m_camera_matrix(2,3);
  Vector2 pixel = Vector2( (m_camera_matrix(0,0)*point(0) + m_camera_matrix(0,1)*point(1) +
                            m_camera_matrix(0,2)*point(2) + m_camera_matrix(0,3)           ) / denominator,
                           (m_camera_matrix(1,0)*point(0) + m_camera_matrix(1,1)*point(1) +
                            m_camera_matrix(1,2)*point(2) + m_camera_matrix(1,3)           ) / denominator);

  // Divide by pixel pitch to convert from metric units to pixels if the intrinsic
  //   values were not specified in pixel units (in that case m_pixel_pitch == 1.0)
  return pixel/m_pixel_pitch;
}

bool PinholeModel::projection_valid(Vector3 const& point) const {
  // z coordinate after extrinsic transformation
  double z = m_extrinsics(2, 0)*point(0) + m_extrinsics(2, 1)*point(1) +
    m_extrinsics(2, 2)*point(2) + m_extrinsics(2,3);
  return z > 0;
}

Vector3 PinholeModel::pixel_to_vector (Vector2 const& pix) const {
  // Apply the inverse lens distortion model
  Vector2 undistorted_pix = m_distortion->undistorted_coordinates(*this, pix*m_pixel_pitch);

  // Compute the direction of the ray emanating from the camera center.
  Vector3 p(0,0,1);
  subvector(p,0,2) = undistorted_pix;
  return normalize( m_inv_camera_transform * p);
}

Vector3 PinholeModel::camera_center(Vector2 const& /*pix*/ ) const {
  return m_camera_center;
};

void PinholeModel::set_camera_center(Vector3 const& position) {
  m_camera_center = position; 
  rebuild_camera_matrix();
}

Quaternion<double> PinholeModel::camera_pose(Vector2 const& /*pix*/ ) const {
  return Quaternion<double>(m_rotation);
}

void PinholeModel::set_camera_pose(Quaternion<double> const& pose) {
  m_rotation = pose.rotation_matrix(); 
  rebuild_camera_matrix();
}

void PinholeModel::set_camera_pose(Matrix<double,3,3> const& pose) {
  m_rotation = pose; 
  rebuild_camera_matrix();
}

void PinholeModel::coordinate_frame(Vector3 &u_vec, Vector3 &v_vec, Vector3 &w_vec) const {
  u_vec = m_u_direction;
  v_vec = m_v_direction;
  w_vec = m_w_direction;
}

void PinholeModel::set_coordinate_frame(Vector3 u_vec, Vector3 v_vec, Vector3 w_vec) {
  m_u_direction = u_vec;
  m_v_direction = v_vec;
  m_w_direction = w_vec;

  rebuild_camera_matrix();
}

Vector3 PinholeModel::coordinate_frame_u_direction() const { return m_u_direction; }
Vector3 PinholeModel::coordinate_frame_v_direction() const { return m_v_direction; }
Vector3 PinholeModel::coordinate_frame_w_direction() const { return m_w_direction; }

const LensDistortion* PinholeModel::lens_distortion() const { return m_distortion.get(); };
void PinholeModel::set_lens_distortion(LensDistortion const& distortion) {
  m_distortion = distortion.copy();
}

void PinholeModel::intrinsic_parameters(double& f_u, double& f_v,
                                        double& c_u, double& c_v) const {
  f_u = m_fu;  f_v = m_fv;  c_u = m_cu;  c_v = m_cv;
}

void PinholeModel::set_intrinsic_parameters(double f_u, double f_v,
                                            double c_u, double c_v) {
  m_fu = f_u;  m_fv = f_v;  m_cu = c_u;  m_cv = c_v;
  rebuild_camera_matrix();
}

Vector2 PinholeModel::focal_length() const { return Vector2(m_fu,m_fv); }
void PinholeModel::set_focal_length(Vector2 const& f, bool rebuild ) {
  m_fu = f[0]; m_fv = f[1];
  if (rebuild) rebuild_camera_matrix();
}
Vector2 PinholeModel::point_offset() const { return Vector2(m_cu,m_cv); }
void PinholeModel::set_point_offset(Vector2 const& c, bool rebuild ) {
  m_cu = c[0]; m_cv = c[1];
  if (rebuild) rebuild_camera_matrix();
}
double PinholeModel::pixel_pitch() const { return m_pixel_pitch; }
void PinholeModel::set_pixel_pitch( double pitch ) { m_pixel_pitch = pitch; }


void PinholeModel::set_camera_matrix( Matrix<double,3,4> const& p ) {
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
    vw_out(WarningMessage,"camera") << "Significant skew not modelled by pinhole camera\n";

  // Rebuild
  rebuild_camera_matrix();
#else
  vw_throw( NoImplErr() << "PinholeModel::set_Camera_Matrix is unavailable without LAPACK" );
#endif
}

Matrix<double,3,4> PinholeModel::camera_matrix() const {
  return m_camera_matrix;
}

void PinholeModel::rebuild_camera_matrix() {

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

  Matrix<double,3,3> rotation_inverse = transpose(m_rotation);
  submatrix(m_extrinsics,0,0,3,3) = uvwRotation * rotation_inverse;
  select_col(m_extrinsics,3) = uvwRotation * -rotation_inverse * m_camera_center;

  m_camera_matrix = m_intrinsics * m_extrinsics;
  m_inv_camera_transform = inverse(uvwRotation*rotation_inverse) * inverse(m_intrinsics);
}

// Apply a given rotation + translation + scale transform to a pinhole camera
void PinholeModel::apply_transform(vw::Matrix3x3 const & rotation,
                                   vw::Vector3   const & translation,
                                   double                scale) {

  // Extract current parameters
  vw::Vector3 position = this->camera_center();
  vw::Quat    pose     = this->camera_pose();
  
  vw::Quat rotation_quaternion(rotation);
  
  // New position and rotation
  position = scale*rotation*position + translation;
  pose     = rotation_quaternion*pose;
  this->set_camera_center(position);
  this->set_camera_pose  (pose);
}

PinholeModel
camera::scale_camera(PinholeModel const& camera_model, float scale) {
  if (scale == 0)
    vw_throw( ArgumentErr() << "PinholeModel::scale_camera cannot have zero scale value!" );
  // Scaling the camera is easy, just update the pixel pitch to account for the new image size.
  Vector2 focal  = camera_model.focal_length();
  Vector2 offset = camera_model.point_offset();
  boost::shared_ptr<LensDistortion> lens = camera_model.lens_distortion()->copy();
  return PinholeModel( camera_model.camera_center(),
                       camera_model.camera_pose().rotation_matrix(),
                       focal[0], focal[1], offset[0], offset[1],
                       camera_model.coordinate_frame_u_direction(),
                       camera_model.coordinate_frame_v_direction(),
                       camera_model.coordinate_frame_w_direction(),
                       *lens,
                       camera_model.pixel_pitch()/scale);
}


PinholeModel
camera::strip_lens_distortion(PinholeModel const& camera_model) {
  Vector2 focal  = camera_model.focal_length();
  Vector2 offset = camera_model.point_offset();
  NullLensDistortion distortion;
  return PinholeModel(camera_model.camera_center(),
                      camera_model.camera_pose().rotation_matrix(),
                      focal[0], focal[1], offset[0], offset[1],
                      camera_model.coordinate_frame_u_direction(),
                      camera_model.coordinate_frame_v_direction(),
                      camera_model.coordinate_frame_w_direction(),
                      distortion, camera_model.pixel_pitch());
}

std::ostream& camera::operator<<(std::ostream& str,
                                 PinholeModel const& model) {
  str << "Pinhole camera: \n";
  str << "\tCamera Center: " << model.camera_center() << "\n";
  str << "\tRotation Matrix: " << model.camera_pose() << "\n";
  str << "\tIntrinsics:\n";
  str << "\t  focal: "       << model.focal_length() << "\n";
  str << "\t  offset: "      << model.point_offset() << "\n";
  str << "\t  pixel pitch: " << model.pixel_pitch()  << "\n";
  str << "\tu direction: " << model.coordinate_frame_u_direction() << "\n";
  str << "\tv direction: " << model.coordinate_frame_v_direction() << "\n";
  str << "\tw direction: " << model.coordinate_frame_w_direction() << "\n";
  str << "\tDistortion Model: " << model.lens_distortion()->name() << "\n";
  str << *(model.lens_distortion()); // this will be multiple lines

  return str;
}
