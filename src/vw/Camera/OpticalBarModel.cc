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

#include <vw/Math/EulerAngles.h>
#include <vw/Camera/CameraModel.h>
#include <vw/Camera/CameraSolve.h>
#include <vw/Camera/OpticalBarModel.h>
#include <vw/Camera/OrbitalCorrections.h>
#include <vw/Math/NewtonRaphson.h>

#include <vw/Camera/LinescanErr.h>

#include <iomanip>

namespace vw {
namespace camera {

OpticalBarModel::OpticalBarModel(): 
  m_motion_compensation(1.0),
  m_have_velocity_vec(false)
  {}

OpticalBarModel::OpticalBarModel(std::string const& path):
  m_motion_compensation(1.0),
  m_have_velocity_vec(false)
   {
  // Create from file. This will read m_mean_earth_radius and m_mean_surface_elevation.    
  read(path);
} 

OpticalBarModel::OpticalBarModel(vw::Vector2i image_size,
                vw::Vector2  center_offset_pixels,
                double   pixel_size,
                double   focal_length,
                double   scan_time,
                bool     scan_left_to_right,
                double   forward_tilt_radians,
                vw::Vector3  initial_position,
                vw::Vector3  initial_orientation,
                double   speed,
                double   motion_compensation_factor,
                vw::Vector3 const& velocity):
    m_image_size          (image_size),
    m_center_loc_pixels   (center_offset_pixels),
    m_pixel_size          (pixel_size),
    m_focal_length        (focal_length),
    m_scan_time           (scan_time),
    m_scan_left_to_right  (scan_left_to_right),
    m_forward_tilt_radians(forward_tilt_radians),
    m_initial_position    (initial_position),
    m_initial_orientation (initial_orientation),
    m_speed               (speed),
    m_motion_compensation(motion_compensation_factor),
    m_mean_earth_radius(DEFAULT_EARTH_RADIUS),
    m_mean_surface_elevation(DEFAULT_SURFACE_ELEVATION),
    m_have_velocity_vec(false) { 
  compute_scan_rate();
}
vw::Vector2i OpticalBarModel::get_image_size() const { 
  return m_image_size;          
}
vw::Vector2 OpticalBarModel::get_optical_center() const { 
  return m_center_loc_pixels;   
}
double OpticalBarModel::get_focal_length() const { 
  return m_focal_length;        
}
double OpticalBarModel::get_scan_rate() const { 
  return m_scan_rate_radians;   
}
double OpticalBarModel::get_speed() const { 
  return m_speed;               
}
double OpticalBarModel::get_pixel_size() const { 
  return m_pixel_size;          
}
double OpticalBarModel::get_scan_time() const { 
  return m_scan_time;  
}
bool OpticalBarModel::get_scan_dir() const { 
  return m_scan_left_to_right;  
}
double OpticalBarModel::get_forward_tilt() const { 
  return m_forward_tilt_radians;
}

double OpticalBarModel::get_motion_compensation() const { 
  return m_motion_compensation; 
}

bool OpticalBarModel::get_have_velocity_vec() const {
  return m_have_velocity_vec;
}

void OpticalBarModel::set_velocity(vw::Vector3 const& velocity) {
  m_velocity = velocity;
  m_have_velocity_vec = (m_velocity != vw::Vector3(0, 0, 0));
}

void OpticalBarModel::set_camera_center(vw::Vector3 const& position) {
  m_initial_position  = position;
}

void OpticalBarModel::set_camera_pose(vw::Vector3 const& orientation) {
  m_initial_orientation = orientation;
}

void OpticalBarModel::set_camera_pose(vw::Quaternion<double> const& pose) {
  set_camera_pose(pose.axis_angle());
}

void OpticalBarModel::set_image_size(vw::Vector2i image_size) { 
  m_image_size = image_size;
  compute_scan_rate();
}

void OpticalBarModel::set_optical_center(vw::Vector2  optical_center) { 
  m_center_loc_pixels = optical_center;
  compute_scan_rate();
}

void OpticalBarModel::set_focal_length(double focal_length) {
  m_focal_length = focal_length;
  compute_scan_rate();
}

void OpticalBarModel::set_speed(double speed) { 
  m_speed = speed;
  compute_scan_rate();
}

void OpticalBarModel::set_pixel_size(double pixel_size) { 
  m_pixel_size = pixel_size;
  compute_scan_rate();
}
void OpticalBarModel::set_scan_time(double scan_time) { 
  m_scan_time   = scan_time;
  compute_scan_rate();
}

void OpticalBarModel::set_scan_dir(bool scan_l_to_r) { 
  m_scan_left_to_right = scan_l_to_r;
}

void OpticalBarModel::set_forward_tilt(double tilt_angle) { 
  m_forward_tilt_radians = tilt_angle;
}

void OpticalBarModel::set_motion_compensation(double mc_factor) { 
  m_motion_compensation = mc_factor;
}

Vector2 OpticalBarModel::pixel_to_sensor_plane(Vector2 const& pixel) const {
  Vector2 result = (pixel - m_center_loc_pixels) * m_pixel_size;
  return result;
}

double OpticalBarModel::sensor_to_alpha(vw::Vector2 const& sensor_loc) const {
  // This calculation comes from the focal point projected on to a circular surface.
  return sensor_loc[0] / m_focal_length;
}

void OpticalBarModel::compute_scan_rate() {

  // This is only needed when not modeling the velocity as a 3D vector. 
  if (m_have_velocity_vec)
    return;

  // Compute the scan angle using pixel information.
  Vector2 p1 = pixel_to_sensor_plane(Vector2(0,0));
  Vector2 p2 = pixel_to_sensor_plane(Vector2(m_image_size - Vector2i(1,1)));
  
  double alpha1     = sensor_to_alpha(p1);
  double alpha2     = sensor_to_alpha(p2);
  double scan_angle = alpha2 - alpha1;
  
  m_scan_rate_radians = scan_angle / m_scan_time;
}

// Return time at the given pixel. In the latest approach the time at last row is 1.
double OpticalBarModel::pixel_to_time_delta(Vector2 const& pix) const {

  // Since the camera sweeps a scan through columns, use that to
  //  determine the fraction of the way it is through the image.
  const int max_col = m_image_size[0]-1;
  double scan_fraction = 0;
  if (m_scan_left_to_right)
    scan_fraction = pix[0] / max_col; // TODO: Add 0.5 pixels to calculation?
  else // Right to left scan direction
    scan_fraction = (max_col - pix[0]) / max_col;
    
  // Normalized scan time in the newer approach
  if (m_have_velocity_vec)
    return scan_fraction;
  
  // In the prior approach, an actual scan time was used, but this is redundant
  // given other parameters.
  double time_delta = scan_fraction * m_scan_time;
  return time_delta;
}

Vector3 OpticalBarModel::get_velocity(vw::Vector2 const& pixel) const {

  // Independent velocity vector in the new approach
  if (m_have_velocity_vec)
    return m_velocity;

  // TODO: For speed, store the pose*velocity vector.
  // Convert the velocity from sensor coords to GCC coords
  Matrix3x3 pose = camera_pose(pixel).rotation_matrix();

  // Recover the satellite attitude relative to the tilted camera position
  Matrix3x3 m = pose*vw::math::rotation_x_axis(m_forward_tilt_radians);
  
  return m * Vector3(0, m_speed, 0);
}

Vector3 OpticalBarModel::camera_center(Vector2 const& pix) const {

  // We model with a constant velocity
  double dt = pixel_to_time_delta(pix);

  return m_initial_position + dt * get_velocity(pix);
}

// Camera pose is treated as constant for the duration of a scan.
Quat OpticalBarModel::camera_pose(Vector2 const& pix) const {
  return axis_angle_to_quaternion(m_initial_orientation);
}

Vector3 OpticalBarModel::pixel_to_vector(Vector2 const& pixel) const {
 
  Vector2 sensor_plane_pos = pixel_to_sensor_plane(pixel);
  Vector3 cam_center       = camera_center(pixel);
  Quat    cam_pose         = camera_pose  (pixel);

  // This is the horizontal angle away from the center point (from straight out of the camera)
  double alpha = sensor_to_alpha(sensor_plane_pos);

  // Distortion caused by compensation for the satellite's forward motion during the image.
  // The film was actually translated underneath the lens to compensate for the motion.
  double image_motion_compensation = 0.0;

  if (m_have_velocity_vec) {
    // Newer approach. Many params that are correlated are combined.
    image_motion_compensation = m_focal_length * sin(alpha) * m_motion_compensation;
  } else {
  
    // Older approach
    // Distance from the camera center to the ground.
    double H = norm_2(cam_center) - (m_mean_surface_elevation + m_mean_earth_radius);
    image_motion_compensation = ((m_focal_length * m_speed) / (H*m_scan_rate_radians))
                                        * sin(alpha) * m_motion_compensation;
                                      
  }
                                        
  if (!m_scan_left_to_right) // Sync alpha with motion compensation.
    image_motion_compensation *= -1.0;

  // This vector is ESD format, consistent with the linescan model. It is in 
  // camera coordinates.
  Vector3 r(m_focal_length * sin(alpha),
            sensor_plane_pos[1] + image_motion_compensation,
            m_focal_length * cos(alpha));
  r = normalize(r);

  // Convert the ray vector into ECEF coordinates.
  Vector3 result = cam_pose.rotate(r);

  return result;
}

// TODO(oalexan1): This could be sped up further, as done in the usgscsm linescan class,
// where an initial affine transform for ground-to-image is found.
Vector2 OpticalBarModel::point_to_pixel(Vector3 const& point) const {

  // Use the image center as the initial guess for the pixel
  Vector2 guess = m_image_size / 2.0;

  // The step size for numerical differentiation in pixel units. This is for
  // finding the Jacobian, so it need not be perfect. A small step size can
  // result in numerical error.   
  double step = 1e-2;

  // The error for this tol is measured at the ground level in LinescanErr.
  double tol = 1e-10;
  
  // Find pix so that err_func(pix) = ans, where guess is the pix initial guess
  vw::Vector2 ans(0, 0);
  LinescanErr err_func(this, point, guess);
  vw::math::NewtonRaphson nr(err_func);
  Vector2 pix = nr.solve(guess, ans, step, tol);
  
  return pix;
}

void OpticalBarModel::apply_transform(vw::Matrix3x3 const & rotation,
                                      vw::Vector3   const & translation,
                                      double                scale) {
  // Extract current parameters
  vw::Vector3 position = this->camera_center();
  vw::Quat    pose     = this->camera_pose();

  vw::Quat rotation_quaternion(rotation);
  
  // New position and rotation
  position = scale * rotation * position + translation;
  pose     = rotation_quaternion * pose;

  this->set_camera_center(position);
  this->set_camera_pose  (pose.axis_angle());
}

void OpticalBarModel::read(std::string const& filename) {

  // Open the input file
  std::ifstream cam_file;
  cam_file.open(filename.c_str());
  if (cam_file.fail())
    vw_throw( IOErr() << "OpticalBarModel::read_file: Could not open file: " << filename );

  // Check for version number on the first line
  std::string line;
  std::getline(cam_file, line);
  if (line.find("VERSION") == std::string::npos) {
    vw_throw( IOErr() << "OpticalBarModel::read_file(): Version missing!\n" );
  }

  int file_version = 1;
  sscanf(line.c_str(),"VERSION_%d", &file_version); // Parse the version of the input file
  if (file_version < 4)
    vw_throw( ArgumentErr() << "OpticalBarModel::read_file(): Versions prior to 4 are not supported!\n" );

  // Read the camera type
  std::getline(cam_file, line);
  if (line.find("OPTICAL_BAR") == std::string::npos)
        vw_throw( ArgumentErr() 
                 << "OpticalBarModel::read_file: Expected OPTICAL_BAR type, but got type "
                 << line );

  // Start parsing all the parameters from the lines.
  std::getline(cam_file, line);
  if (!cam_file.good() || sscanf(line.c_str(),"image_size = %d %d",
      &m_image_size[0], &m_image_size[1]) != 2) {
    cam_file.close();
    vw_throw( IOErr() << "OpticalBarModel::read_file(): Could not read the image size\n" );
  }

  std::getline(cam_file, line);
  if (!cam_file.good() || sscanf(line.c_str(),"image_center = %lf %lf",
      &m_center_loc_pixels[0], &m_center_loc_pixels[1]) != 2) {
    cam_file.close();
    vw_throw( IOErr() << "OpticalBarModel::read_file(): Could not read the image center\n" );
  }

  std::getline(cam_file, line);
  if (!cam_file.good() || sscanf(line.c_str(),"pitch = %lf", &m_pixel_size) != 1) {
    cam_file.close();
    vw_throw( IOErr() << "OpticalBarModel::read_file(): Could not read the pixel pitch\n" );
  }

  std::getline(cam_file, line);
  if (!cam_file.good() || sscanf(line.c_str(),"f = %lf", &m_focal_length) != 1) {
    cam_file.close();
    vw_throw( IOErr() << "OpticalBarModel::read_file(): Could not read the focal_length\n" );
  }

  std::getline(cam_file, line);
  if (!cam_file.good() || sscanf(line.c_str(),"scan_time = %lf", &m_scan_time) != 1) {
    cam_file.close();
    vw_throw( IOErr() << "OpticalBarModel::read_file(): Could not read the scan time\n" );
  }

  std::getline(cam_file, line);
  if (!cam_file.good() || sscanf(line.c_str(),"forward_tilt = %lf", &m_forward_tilt_radians) != 1) {
    cam_file.close();
    vw_throw( IOErr() << "OpticalBarModel::read_file(): Could not read the forward tilt angle\n" );
  }
  
  std::getline(cam_file, line);
  if (!cam_file.good() || sscanf(line.c_str(),"iC = %lf %lf %lf", 
        &m_initial_position(0), &m_initial_position(1), &m_initial_position(2)) != 3) {
    cam_file.close();
    vw_throw( IOErr() << "OpticalBarModel::read_file(): Could not read the initial position\n" );
  }

  // Read and convert the rotation matrix.
  Matrix3x3 rot_mat;
  std::getline(cam_file, line);
  if ( !cam_file.good() ||
       sscanf(line.c_str(), "iR = %lf %lf %lf %lf %lf %lf %lf %lf %lf",
              &rot_mat(0,0), &rot_mat(0,1), &rot_mat(0,2),
              &rot_mat(1,0), &rot_mat(1,1), &rot_mat(1,2),
              &rot_mat(2,0), &rot_mat(2,1), &rot_mat(2,2)) != 9 ) {
    cam_file.close();
    vw_throw( IOErr() << "OpticalBarModel::read_file(): Could not read the rotation matrix\n" );
  }
  Quat q(rot_mat);
  m_initial_orientation = q.axis_angle();

  std::getline(cam_file, line);
  if (!cam_file.good() || sscanf(line.c_str(),"speed = %lf", &m_speed) != 1) {
    cam_file.close();
    vw_throw( IOErr() << "OpticalBarModel::read_file(): Could not read the speed\n" );
  }

  std::getline(cam_file, line);
  if (!cam_file.good() || sscanf(line.c_str(),"mean_earth_radius = %lf", &m_mean_earth_radius) != 1) {
    cam_file.close();
    vw_throw( IOErr() << "OpticalBarModel::read_file(): Could not read the mean earth radius\n" );
  }
  
  std::getline(cam_file, line);
  if (!cam_file.good() || sscanf(line.c_str(),"mean_surface_elevation = %lf", &m_mean_surface_elevation) != 1) {
    cam_file.close();
    vw_throw( IOErr() << "OpticalBarModel::read_file(): Could not read the mean surface elevation\n" );
  }

  std::getline(cam_file, line);
  if (!cam_file.good() || sscanf(line.c_str(),"motion_compensation_factor = %lf",
                                 &m_motion_compensation) != 1) {
    cam_file.close();
    vw_throw( IOErr() << "OpticalBarModel::read_file(): Could not read motion compensation factor\n" );
  }

  std::getline(cam_file, line);
  m_scan_left_to_right = line.find("scan_dir = left") == std::string::npos;

  // Get the line with velocity. This is optional. If not set, use 0.
  std::getline(cam_file, line);
  if (line.find("velocity = ") != std::string::npos) {
    if (sscanf(line.c_str(),"velocity = %lf %lf %lf",
               &m_velocity[0], &m_velocity[1], &m_velocity[2]) != 3) {
      cam_file.close();
      vw_throw( IOErr() << "OpticalBarModel::read_file(): Could not read the velocity\n" );
    }
  } else {
    m_velocity = Vector3(0, 0, 0);
  }
  
  m_have_velocity_vec = (m_velocity != vw::Vector3(0, 0, 0));
  
  compute_scan_rate();
  
  cam_file.close();
}

void OpticalBarModel::write(std::string const& filename) const {

  // TODO: Make compatible with .tsai files!

  // Open the output file for writing
  std::ofstream cam_file(filename.c_str());
  if( !cam_file.is_open() ) 
    vw_throw( IOErr() << "OpticalBarModel::write: Could not open file: " << filename );

  // Write the pinhole camera model parts
  //   # digits to survive double->text->double conversion
  const size_t ACCURATE_DIGITS = 17; // = std::numeric_limits<double>::max_digits10
  cam_file << std::setprecision(ACCURATE_DIGITS); 
  cam_file << "VERSION_4\n";
  cam_file << "OPTICAL_BAR\n";
  cam_file << "image_size = "   << m_image_size[0] << " " 
                                << m_image_size[1]<< "\n";
  cam_file << "image_center = " << m_center_loc_pixels[0] << " "
                                << m_center_loc_pixels[1] << "\n";
  cam_file << "pitch = "        << m_pixel_size             << "\n";
  cam_file << "f = "            << m_focal_length           << "\n";
  cam_file << "scan_time = "   << m_scan_time     << "\n";
  cam_file << "forward_tilt = " << m_forward_tilt_radians   << "\n";
  cam_file << "iC = " << m_initial_position[0] << " "
                      << m_initial_position[1] << " "
                      << m_initial_position[2] << "\n";
  // Store in the same format as the pinhole camera model.
  Matrix3x3 rot_mat = camera_pose(Vector2(0,0)).rotation_matrix();
  cam_file << "iR = " << rot_mat(0,0) << " " << rot_mat(0,1) << " " << rot_mat(0,2) << " "
                      << rot_mat(1,0) << " " << rot_mat(1,1) << " " << rot_mat(1,2) << " "
                      << rot_mat(2,0) << " " << rot_mat(2,1) << " " << rot_mat(2,2) << "\n";
  cam_file << "speed = " << m_speed << "\n";
  cam_file << "mean_earth_radius = "          << m_mean_earth_radius      << "\n";
  cam_file << "mean_surface_elevation = "     << m_mean_surface_elevation << "\n";
  cam_file << "motion_compensation_factor = " << m_motion_compensation    << "\n";
  if (m_scan_left_to_right)
    cam_file << "scan_dir = right\n";
  else
    cam_file << "scan_dir = left\n";
  
  // Write the 3 values of m_velocity
  if (m_have_velocity_vec)
    cam_file << "velocity = " << m_velocity[0] << " "
                              << m_velocity[1] << " "
                              << m_velocity[2] << "\n";
                              
  cam_file.close();
}

std::ostream& operator<<( std::ostream& os, OpticalBarModel const& camera_model) {
  os << "\n------------------------ Optical Bar Model -----------------------\n\n";
  os << " Image size:             " << camera_model.m_image_size             << "\n";
  os << " Center loc (pixels):    " << camera_model.m_center_loc_pixels      << "\n";
  os << " Pixel size (m):         " << camera_model.m_pixel_size             << "\n";
  os << " Focal length (m):       " << camera_model.m_focal_length           << "\n";
  os << " Scan time (s):          " << camera_model.m_scan_time              << "\n";
  os << " Forward tilt (rad):     " << camera_model.m_forward_tilt_radians   << "\n";
  os << " Initial position:       " << camera_model.m_initial_position       << "\n";
  os << " Initial pose:           " << camera_model.m_initial_orientation    << "\n";
  os << " Speed:                  " << camera_model.m_speed                  << "\n";
  os << " Mean earth radius:      " << camera_model.m_mean_earth_radius      << "\n";
  os << " Mean surface elevation: " << camera_model.m_mean_surface_elevation << "\n";
  os << " Motion comp factor:     " << camera_model.m_motion_compensation    << "\n";
  os << " Left to right scan:     " << camera_model.m_scan_left_to_right     << "\n";
  os << " Velocity (m/s):         " << camera_model.m_velocity               << "\n";
  os << "\n------------------------------------------------------------------------\n\n";
  return os;
}

}} // namespace asp::camera
