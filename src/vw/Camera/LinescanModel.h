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


/// \file LinescanModel.h
///
/// A generic linescan camera model object
///
///
#ifndef _VW_CAMERA_LINESCAN_MODEL_H_
#define _VW_CAMERA_LINESCAN_MODEL_H_

#include <vw/Math/Quaternion.h>
#include <vw/Math/LevenbergMarquardt.h>
#include <vw/Camera/CameraModel.h>

namespace vw {
namespace camera {

  /// This is a generic line scan camera model that can be derived
  /// from to help implement specific cameras.  Some parts (velocity
  /// and atmospheric correction) currently only work for Earth.
  ///
  /// This expects the pose to be a rotation from the camera frame to
  /// the world frame. The position is a the camera's location in the
  /// world frame.
  ///
  /// The intrinisic model expects +Z to point out the camera. +X is
  /// the column direction of the image and is perpendicular to
  /// direction of flight. +Y is the row direction of the image (down
  /// the image); it is also the flight direction.  If this is not 
  /// accurate for your camera you can apply a rotation in PoseFuncT.

  class LinescanModel : public vw::camera::CameraModel {

  public:
    //------------------------------------------------------------------
    // Constructors / Destructors
    //------------------------------------------------------------------
    LinescanModel(Vector2i const& image_size,
		              bool            correct_velocity_aberration,
		              bool            correct_atmospheric_refraction) :
      m_image_size(image_size), 
      m_correct_velocity_aberration(correct_velocity_aberration),
      m_correct_atmospheric_refraction(correct_atmospheric_refraction){
      
      // Set default values for these constants which can be overridden later on.
      const double DEFAULT_EARTH_RADIUS      = 6371000.0;  // In meters.
      const double DEFAULT_SURFACE_ELEVATION = 0.0;
      m_mean_earth_radius      = DEFAULT_EARTH_RADIUS;
      m_mean_surface_elevation = DEFAULT_SURFACE_ELEVATION;
    }

    virtual ~LinescanModel() {}
    virtual std::string type() const { return "Linescan"; }

    //------------------------------------------------------------------
    // Interface
    //------------------------------------------------------------------
    
    // -- This set of functions implements virtual functions from CameraModel.h --
    
    /// Get the pixel that observes a point in world coordinates.
    virtual Vector2 point_to_pixel (Vector3 const& point) const;

    /// Gives a pointing vector in the world coordinates.
    virtual Vector3 pixel_to_vector(Vector2 const& pix) const;

    /// Gives the camera position in world coordinates.
    virtual Vector3 camera_center(Vector2 const& pix) const {
      return get_camera_center_at_time(get_time_at_line(pix.y()));
    }

    /// Gives a pose vector which represents the rotation from camera to world units
    virtual Quat camera_pose(Vector2 const& pix) const {
      return get_camera_pose_at_time(get_time_at_line(pix.y()));
    }

    // -- These are new functions --

    /// Returns the image size in pixels
    Vector2i get_image_size() const { return m_image_size; }
    
    int samples_per_line() const { return m_image_size[0]; }
    int number_of_lines () const { return m_image_size[1]; }
    
    /// Gives the camera velocity in world coordinates.
    Vector3 camera_velocity(vw::Vector2 const& pix) const {
      return get_camera_velocity_at_time(get_time_at_line(pix.y()));
    }
    
    // New functions for Linescan derived classes.
    // - Most of these deal with the fact that the camera is moving while
    //   acquiring the image.
    // - Derived classes need to implement these so that the functions inherited from
    //   CameraModel will work.  Consider using the functors in Extrinsics.h.

    /// Gives the camera center at a time.
    virtual Vector3 get_camera_center_at_time  (double time) const = 0;
    /// Gives the camera velocity at a time.
    virtual Vector3 get_camera_velocity_at_time(double time) const = 0;
    /// Get the pose at a time.
    virtual Quat    get_camera_pose_at_time    (double time) const = 0;
    /// Return the computed time for a given line.
    virtual double  get_time_at_line           (double line) const = 0;

    /// As pixel_to_vector, but in the local camera frame.
    virtual Vector3 get_local_pixel_vector(vw::Vector2 const& pix) const = 0;


    // Here we use an initial guess for the line number
    // - This class provides a generic implementation but specific
    //   linescan cameras may be able to use more specific implementation.
    virtual Vector2 point_to_pixel(vw::Vector3 const& point, double starty) const;

  protected:

    /// Image size in pixels: [num lines, num samples]
    Vector2i m_image_size;      

    double m_mean_earth_radius;
    double m_mean_surface_elevation;

    /// Set this flag to enable velocity aberration correction.
    /// - For satellites this makes a big difference, make sure it is set!
    bool m_correct_velocity_aberration;
    
    /// Set this flag to enable atmospheric refraction correction.
    bool m_correct_atmospheric_refraction;

  protected:

    /// Returns the radius of the Earth under the current camera position.
    double get_earth_radius() const;

  }; // End class LinescanModel
  
/*
  /// Output stream method for printing a summary of the linear
  /// pushbroom camera model parameters.
  std::ostream& operator<<( std::ostream& os, LinescanModel const& camera_model);
*/

}}      // namespace vw::camera

#endif  //_VW_CAMERA_LINEARPUSHBROOM_MODEL_H_
