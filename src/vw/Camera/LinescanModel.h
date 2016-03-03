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
#include <vw/Camera/CameraModel.h>

namespace vw {
namespace camera {

  template <class PositionFuncT, class PoseFuncT>
  class LinescanModel : public CameraModel {

    // Extrinsics
    // - These functions need to be one of the classes from Extrinsics.h
    PositionFuncT m_position_func; ///< Function to compute the position at time T
    PoseFuncT     m_pose_func;     ///< Function to compute the pose     at time T

    // Intrinsics
    int     m_number_of_lines;  ///< "height" of the image
    int     m_samples_per_line; ///< "width" of the image
    int     m_sample_offset;    ///< Number of skipped pixels at the start of the scan line
    double  m_focal_length;
    double  m_across_scan_pixel_size; ///< Pixel size in "height" direction
    double  m_along_scan_pixel_size;  ///< Pixel size in "width"  direction
    std::vector<double> m_line_times; ///< The time of each line scan, size == m_number_of_lines
    Vector3 m_pointing_vec; ///< The sensor pointing direction relative to sensor pose.
    Vector3 m_u_vec;        ///< The column direction vector   relative to sensor pose.

  public:
    //------------------------------------------------------------------
    // Constructors / Destructors
    //------------------------------------------------------------------

    /// This version of the constructor assumes that the line
    /// integration varies one each scanline of the image.  The
    /// line_integration_times vector must be equal to the number_of_lines
    LinescanModel( int     number_of_lines,
                   int     samples_per_line,
                   int     sample_offset,
                   double  focal_length,
                   double  along_scan_pixel_size,
                   double  across_scan_pixel_size,
                   std::vector<double> const& line_times,
                   Vector3 pointing_vec,
                   Vector3 u_vec,
                   PositionFuncT const& position_func,
                   PoseFuncT     const& pose_func);

    /// This version of the constructor assumes that the line
    /// integration time is consistent throughout the image.
    LinescanModel( int     number_of_lines,
                   int     samples_per_line,
                   int     sample_offset,
                   double  focal_length,
                   double  along_scan_pixel_size,
                   double  across_scan_pixel_size,
                   double  line_integration_time,
                   Vector3 pointing_vec,
                   Vector3 u_vec,
                   PositionFuncT const& position_func,
                   PoseFuncT     const& pose_func);

    virtual ~LinescanModel() {}
    virtual std::string type() const { return "Linescan"; }

    //------------------------------------------------------------------
    // Interface
    //------------------------------------------------------------------
    virtual Vector2 point_to_pixel(Vector3 const& /*point*/) const {
      vw_throw( vw::NoImplErr() << "LinescanModel::point_to_pixel is not yet implemented." );
      return Vector2(); // never reached
    }

    /// Given a pixel in image coordinates, what is the pointing
    /// vector in 3-space if you apply the camera model.
    ///
    /// Important Note: For linear pushbroom sensors, the orientation
    /// of the image is important.  In the Camera module we adopt the
    /// following convention for pushbroom imagers:
    ///
    /// u is the crosstrack (i.e., perspective) direction
    /// v is the downtrack (i.e., orthographic) direction
    virtual Vector3 pixel_to_vector(Vector2 const& pix) const;

    virtual Vector3 camera_center(Vector2 const& pix = Vector2() ) const;

    /// Returns the pose (as a quaternion) of the camera for a given
    /// pixel. Pose is the rotation from camera coordinates to world.
    virtual Quaternion<double> camera_pose(Vector2 const& pix) const;

    //------------------------------------------------------------------
    // Public Methods
    //------------------------------------------------------------------

    // Accessors
    virtual double              number_of_lines       () const { return m_number_of_lines;        }
    virtual double              samples_per_line      () const { return m_samples_per_line;       }
    virtual double              sample_offset         () const { return m_sample_offset;          }
    virtual double              focal_length          () const { return m_focal_length;           }
    virtual double              along_scan_pixel_size () const { return m_along_scan_pixel_size;  }
    virtual double              across_scan_pixel_size() const { return m_across_scan_pixel_size; }
    virtual std::vector<double> line_times            () const { return m_line_times;             }
    
    // Retrieve the position and pose at the given time
    virtual Vector3             camera_position(double t) const { return m_position_func(t); }
    virtual Quaternion<double>  camera_pose    (double t) const { return m_pose_func    (t); }

    virtual void set_number_of_lines       (int    val) { m_number_of_lines       = val; }
    virtual void set_samples_per_line      (int    val) { m_samples_per_line      = val; }
    virtual void set_sample_offset         (int    val) { m_sample_offset         = val; }
    virtual void set_along_scan_pixel_size (double val) { m_along_scan_pixel_size = val; }
    virtual void set_across_scan_pixel_size(double val) {m_across_scan_pixel_size = val; }
    virtual void set_focal_length          (double val) {m_focal_length           = val; }
    virtual void set_line_times            (std::vector<double> val) { m_line_times = val; }
    
  protected:
    /// Estimates the line time for a given subpixel row location that may
    ///  fall between the known line times.
    inline double interp_line_time(double row) const;
    
  }; // End class LinescanModel

  /// Output stream method for printing a summary of the linear
  /// pushbroom camera model parameters.
  template <class PositionFuncT, class PoseFuncT>
  std::ostream& operator<<( std::ostream& os, LinescanModel<PositionFuncT, PoseFuncT> const& camera_model);


#include "LinescanModel.tcc"

}}      // namespace vw::camera

#endif  //_VW_CAMERA_LINEARPUSHBROOM_MODEL_H_
