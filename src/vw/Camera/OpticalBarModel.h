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


/// \file OpticalBarModel.h
///
/// A generic linescan camera model object  TODO REMOVE DUPLICATE
///
///
#ifndef _VW_CAMERA_OPTICALBAR_MODEL_H_
#define _VW_CAMERA_OPTICALBAR_MODEL_H_

#include <vw/Math/Quaternion.h>
#include <vw/Math/LevenbergMarquardt.h>
#include <vw/Camera/CameraModel.h>

namespace vw {
namespace camera {

  // A camera model to approximate the type of optical bar cameras
  // that were used in the Corona and Hexagon satellites.
  
  // This implements the logic in: 
  // Rigorous Panoramic Camera Model for DISP Imagery
  // Tony Schenk, Beata Csatho, Sung Woong Shin
  // - Proceedings of Joint Workshop of ISPRS Working Groups, 2003
  
  // With the motion compensation logic taken from:
  // Mathematical modelling of historical reconnaissance CORONA KH-4B imagery
  // HG Sohn, GH Kim, JH Yom - The Photogrammetric Record, 2004
  
  // Other good papers:
  
  // An evaluation of the stereoscopic capabilities of CORONA
  // declassified spy satellite image data, Nikolaos Galiatsatos,
  // Daniel N. M. Donoghue, and Graham Philip.

  // Stereo analysis, DEM extraction and orthorectification of CORONA
  // satellite imagery: archaeological applications from the Near East
  // Jesse Casana1 & Jackson Cothren

  // Data is at https://earthexplorer.usgs.gov/

  class OpticalBarModel : public vw::camera::CameraModel {

  public:
    //------------------------------------------------------------------
    // Constructors / Destructors
    //------------------------------------------------------------------

    OpticalBarModel();
    OpticalBarModel(std::string const& path);
    OpticalBarModel(vw::Vector2i image_size,
                    vw::Vector2  center_offset_pixels,
                    double   pixel_size,
                    double   focal_length,
                    double   scan_time,
                    bool     scan_left_to_right,
                    double   forward_tilt_radians,
                    vw::Vector3  initial_position,
                    vw::Vector3  initial_orientation,
                    double   speed,
                    double   motion_compensation_factor);
    virtual ~OpticalBarModel() {}

    virtual std::string type() const { return "OpticalBar"; }

    // TODO: Make compatible with .tsai files!
    /// Read / Write a from a file on disk.
    void read (std::string const& filename);
    void write(std::string const& filename) const;

    //------------------------------------------------------------------
    // Interface
    //------------------------------------------------------------------

    // -- This set of functions implements virtual functions from CameraModel.h --

    /// Get the pixel that observes a point in world coordinates.
    virtual vw::Vector2 point_to_pixel (vw::Vector3 const& point) const;

    /// Gives a pointing vector in the world coordinates.
    virtual vw::Vector3 pixel_to_vector(vw::Vector2 const& pix) const;

    /// Gives the camera position in world coordinates.
    virtual vw::Vector3 camera_center(vw::Vector2 const& pix) const;

    /// Gives a pose vector which represents the rotation from camera to world units
    virtual vw::Quat camera_pose(vw::Vector2 const& pix) const;

    // -- These are new functions --

    // These return the initial center/pose at time=0.
    vw::Vector3 camera_center() const {return m_initial_position;}
    vw::Quat    camera_pose  () const {return camera_pose(vw::Vector2(0,0));} // Constant

    /// Apply a given rotation + translation + scale transform to the camera.
    void apply_transform(vw::Matrix3x3 const & rotation,
                         vw::Vector3   const & translation,
                         double                scale);

    // Parameter accessors

    void set_camera_center(vw::Vector3 const& position   )     {m_initial_position    = position;}
    void set_camera_pose  (vw::Vector3 const& orientation)     {m_initial_orientation = orientation;}
    void set_camera_pose  (vw::Quaternion<double> const& pose) {set_camera_pose(pose.axis_angle());}

    /// Returns the image size in pixels
    vw::Vector2i get_image_size    () const { return m_image_size;          }
    vw::Vector2  get_optical_center() const { return m_center_loc_pixels;   }
    double       get_focal_length  () const { return m_focal_length;        }
    double       get_scan_rate     () const { return m_scan_rate_radians;   }
    double       get_speed         () const { return m_speed;               }
    double       get_pixel_size    () const { return m_pixel_size;          }
    double       get_scan_time     () const { return m_scan_time;  }
    bool         get_scan_dir      () const { return m_scan_left_to_right;  }
    double       get_forward_tilt  () const { return m_forward_tilt_radians;}

    // TODO(oalexan1): Move this to .cc

    void set_image_size(vw::Vector2i image_size) { 
      m_image_size = image_size;
    }
    void set_optical_center(vw::Vector2  optical_center) { 
      m_center_loc_pixels = optical_center;
    }
    void set_focal_length(double focal_length) {
      m_focal_length         = focal_length;
    }
    
    // TODO(oalexan1): Wipe most of these
    void set_speed(double speed) { 
      m_speed = speed;
    }
    
    void set_pixel_size(double pixel_size) { 
      m_pixel_size = pixel_size;
    }
     
    void set_scan_time(double scan_time) { 
      m_scan_time   = scan_time;
    }
         
    void set_scan_dir(bool scan_l_to_r) { 
      m_scan_left_to_right = scan_l_to_r;
    }
    
    void set_forward_tilt(double tilt_angle) { 
      m_forward_tilt_radians = tilt_angle;
    }
    
    // Set the constant velocity in ECEF
    void set_velocity(vw::Vector3 const& velocity) {
      m_velocity = velocity;
    }

    /// Returns the constant velocity in ECEF
    vw::Vector3 get_velocity() const {
      return m_velocity;
    }
    
    double get_motion_compensation() const { 
      return m_motion_compensation; 
    }
    
    void set_motion_compensation(double mc_factor) { 
      m_motion_compensation = mc_factor;
    }

    friend std::ostream& operator<<(std::ostream&, OpticalBarModel const&);

  private:

    /// Get the alpha (scanner rotation) angle from a location on the sensor (film) plane.
    double sensor_to_alpha(vw::Vector2 const& sensor_loc) const;

    /// Get position on the (flattened) sensor (film) plane in meters.
    vw::Vector2 pixel_to_sensor_plane(vw::Vector2 const& pixel) const;

    /// Compute the normalized time since start of scan for a given pixel.
    double pixel_to_time_delta(vw::Vector2 const& pixel) const;

  protected:

    /// Image size in pixels: [columns, rows]
    vw::Vector2i m_image_size;

    /// Offset from the top left of the image (film) to the camera center.
    vw::Vector2  m_center_loc_pixels;

    /// The physical size of each pixel in meters (scanner resolution).
    /// - Should be 7 or 14 microns.
    double m_pixel_size;

    /// The focal length in meters
    // Nominal focal length for Corona is 609.602mm
    // Nominal focal length for Hexagon is 1.5m inches
    double m_focal_length;

    /// The time it takes to complete one scan.
    /// - The Corona scan time is about 0.5 seconds.
    double m_scan_time;

    /// The angular velocity of the scanner.
    /// - The Corona scan rate is nominally 192 degrees/second
    double m_scan_rate_radians;

    /// The tilt of the camera forward relative to nadir.
    /// - Some optical bar systems have a forward and backward aimed camera.
    /// - This tilt angle is needed to properly account for the camera's motion.
    double m_forward_tilt_radians;

    /// If true, the first column corresponds to the inital position,
    ///  otherwise the last column is at the initial position.
    bool m_scan_left_to_right;
    
    vw::Vector3 m_initial_position;
    vw::Vector3 m_initial_orientation; // TODO: Record as matrix
    double      m_speed; /// Velocity in the sensor Y axis only.
    vw::Vector3 m_velocity; // velocity in ECEF, will supersede m_speed
    
    // These are used for ray corrections.
    double m_mean_earth_radius;
    double m_mean_surface_elevation;

    /// Apply this fraction of the nominal motion compensation.
    double m_motion_compensation;

  protected:

    /// Returns the radius of the Earth under the current camera position.
    double get_earth_radius() const;

  }; // End class OpticalBarModel


  /// Output stream method.
  std::ostream& operator<<( std::ostream& os, OpticalBarModel const& camera_model);


}}      // namespace asp::camera

#endif  //_ASP_CAMERA_OPTICALBAR_MODEL_H_
