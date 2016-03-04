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

  // This is a generic line scan camera model that should be
  // possible to use for most linescan cameras.
  //
  // This expects the pose to be a rotation from the camera frame to
  // the world frame. The position is a the camera's location in the
  // world frame.
  //
  // The intrinisic model expects +Z to point out the camera. +X is
  // the column direction of the image and is perpendicular to
  // direction of flight. +Y is the row direction of the image (down
  // the image); it is also the flight direction.  If this is not 
  // accurate for your camera you can apply a rotation in PoseFuncT.

  template <class PositionFuncT, class VelocityFuncT, class PoseFuncT, class TimeFuncT>
  class LinescanModel : public vw::camera::CameraModel {

  public:
    //------------------------------------------------------------------
    // Constructors / Destructors
    //------------------------------------------------------------------
    LinescanModel(PositionFuncT const& position,
		              VelocityFuncT const& velocity,
		              PoseFuncT     const& pose,
		              TimeFuncT     const& time,
		              vw::Vector2i  const& image_size,
		              vw::Vector2   const& detector_origin,
		              double focal_length,
		              bool   correct_velocity_aberration
		              ) :
      m_position_func(position), m_velocity_func(velocity),
      m_pose_func(pose), m_time_func(time),
      m_detector_origin(detector_origin),
      m_focal_length(focal_length),
      m_image_size(image_size), 
      m_correct_velocity_aberration(correct_velocity_aberration){}

    virtual ~LinescanModel() {}
    virtual std::string type() const { return "Linescan"; }

    //------------------------------------------------------------------
    // Interface
    //------------------------------------------------------------------
    
    
    vw::Vector2i get_image_size     () const {return m_image_size;     } ///< Returns the image size in pixels
    double       get_focal_length   () const {return m_focal_length;   } ///< Returns the focal length in pixels
    vw::Vector2  get_detector_origin() const {return m_detector_origin;} ///< Returns the detector origin in pixels
    
    virtual vw::Vector2 point_to_pixel    (vw::Vector3 const& point) const;

    // Here we use an initial guess for the line number
    vw::Vector2 point_to_pixel            (vw::Vector3 const& point, double starty) const;
    vw::Vector2 point_to_pixel_uncorrected(vw::Vector3 const& point, double starty) const; ///< Never  corrects velocity aberration
    vw::Vector2 point_to_pixel_corrected  (vw::Vector3 const& point, double starty) const; ///< Always corrects velocity aberration

    // Make these virtual as needed.  Templates should avoid that problem.

    /// Gives a pointing vector in the world coordinates.
    vw::Vector3 pixel_to_vector(vw::Vector2 const& pix) const;

    /// Gives the camera position in world coordinates.
    vw::Vector3 camera_center(vw::Vector2 const& pix) const {
      return m_position_func( m_time_func( pix.y() ) );
    }

    /// Gives the camera center at the specified time.
    vw::Vector3 camera_center(double time) const {
      return m_position_func( time );
    }

    /// Gives the camera velocity in world coordinates.
    vw::Vector3 camera_velocity(vw::Vector2 const& pix) const {
      return m_velocity_func( m_time_func( pix.y() ) );
    }
    /// Gives the camera velocity at a specified time.
    vw::Vector3 camera_velocity(double time) const {
      return m_velocity_func(time);
    }

    /// Gives a pose vector which represents the rotation from camera to world units
    vw::Quat camera_pose(vw::Vector2 const& pix) const {
      return m_pose_func( m_time_func( pix.y() ) );
    }
    /// Get the pose at the specified time.
    vw::Quat camera_pose(double time) const {
      return m_pose_func(time);
    }

    /// Return the computed time for a given line.
    double get_line_time(double line) const { return m_time_func(line); }

    PositionFuncT const& get_position_func() const {return m_position_func;} ///< Access the position function
    VelocityFuncT const& get_velocity_func() const {return m_velocity_func;} ///< Access the velocity function
    PoseFuncT     const& get_pose_func    () const {return m_pose_func;    } ///< Access the pose     function
    TimeFuncT     const& get_time_func    () const {return m_time_func;    } ///< Access the time     function

  protected:
    // Extrinsics
    PositionFuncT m_position_func; ///< Yields position at time T
    VelocityFuncT m_velocity_func; ///< Yields velocity at time T
    PoseFuncT     m_pose_func;     ///< Yields pose     at time T
    TimeFuncT     m_time_func;     ///< Yields time at a given line.

    // Intrinsics

    /// Location of (0,0) coordinate of the detector relative to the center of
    ///  the origin of the camera coordinate system.
    /// - Stored internally in pixels.
    vw::Vector2  m_detector_origin; 
    
    double       m_focal_length;    ///< The focal length, also stored in pixels.
    vw::Vector2i m_image_size;      ///< Image size in pixels: [num lines, num samples]

    /// Set this flag to enable velocity aberration correction.
    /// - For satellites this makes a big difference, make sure it is set!
    bool m_correct_velocity_aberration;


  protected:
    // Levenberg Marquardt solver for linescan number
    //
    // We solve for the line number of the image that position the
    // camera so that the projection into the camera model actually
    // hits the detector. The detector is normally offset in the y
    // direction on the optical plane.
    class LinescanLMA : public vw::math::LeastSquaresModelBase<LinescanLMA> {
      const LinescanModel* m_model;
      vw::Vector3 m_point;
    public:
      typedef vw::Vector<double> result_type;   // 1D error on the optical plane.
      typedef result_type        domain_type;   // 1D linescan number
      typedef vw::Matrix<double> jacobian_type;

      LinescanLMA( const LinescanModel* model, const vw::Vector3& pt ) :
        m_model(model), m_point(pt) {}

      inline result_type operator()( domain_type const& y ) const;
    };

    // Levenberg Marquardt solver for linescan number (y) and pixel
    // number (x) for the given point in space. The obtained solution
    // pixel (x, y) must be such that the vector from this camera
    // pixel goes through the given point. The extra complication as
    // compared to LinescanLMA is the non-linear velocity aberration
    // correction in pixel_to_vector. This makes it for a more complex
    // equation and we need to solve for both x and y, rather than
    // just for y and getting x for free.
    class LinescanCorrLMA : public vw::math::LeastSquaresModelBase<LinescanCorrLMA> {
      const LinescanModel* m_model;
      vw::Vector3 m_point;
    public:
      typedef vw::Vector2 domain_type;     // 2D pixel, input to cost function vector
      typedef vw::Vector3 result_type;     // 3D error, output of cost function vector
      typedef vw::Matrix<double, 3, 2> jacobian_type;

      LinescanCorrLMA( const LinescanModel* model, const vw::Vector3& pt ) :
        m_model(model), m_point(pt) {}

      inline result_type operator()( domain_type const& pix ) const {
        return m_model->pixel_to_vector(pix) - normalize(m_point - m_model->camera_center(pix));
      }

    };

  }; // End class LinescanModel
  
  /// Output stream method for printing a summary of the linear
  /// pushbroom camera model parameters.
  template <class PositionFuncT, class VelocityFuncT, class PoseFuncT, class TimeFuncT>
  std::ostream& operator<<( std::ostream& os, LinescanModel<PositionFuncT, VelocityFuncT, PoseFuncT, TimeFuncT> const& camera_model);

#include <vw/Camera/LinescanModel.tcc>

}}      // namespace vw::camera

#endif  //_VW_CAMERA_LINEARPUSHBROOM_MODEL_H_
