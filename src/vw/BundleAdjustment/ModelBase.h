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


/// \file ModelBase.h
///
/// Classes used to define the bundle adjustment problem. Every
/// program implementing BA, must implement one of these.

#ifndef __VW_BUNDLEADJUSTMENT_MODEL_BASE_H__
#define __VW_BUNDLEADJUSTMENT_MODEL_BASE_H__

// Standard
#include <string>

// Vision Workbench
#include <vw/Math/Matrix.h>
#include <vw/Math/Vector.h>
#include <vw/BundleAdjustment/ControlNetwork.h>
#include <vw/Core/Log.h>
#include <vw/Camera/CameraModel.h>

// Boost
#include <boost/smart_ptr.hpp>

namespace vw {
namespace ba {

  // CRTP Base class for Bundle Adjustment functors.
  //---------------------------------------------------------
  template <class ImplT, size_t CameraParamsN, size_t PointParamsN>
  class ModelBase {
  public:
    static const size_t camera_params_n         = CameraParamsN;
    static const size_t point_params_n          = PointParamsN;
    static const size_t focal_length_params_n   = 1;
    static const size_t optical_center_params_n = 2;
    static const size_t nonlens_intrinsics_n    = focal_length_params_n + optical_center_params_n;

    /// \cond INTERNAL
    // Methods to access the derived type
    inline ImplT      & impl()       { return static_cast<ImplT      &>(*this); }
    inline ImplT const& impl() const { return static_cast<ImplT const&>(*this); }
    /// \endcond

    // Required access to camera
    Vector2 cam_pixel ( size_t i, size_t j,
                         Vector<double,camera_params_n> const& cam_j,
                         Vector<double,point_params_n> const& point_i ) {
      return impl().cam_pixel(i,j,cam_j,point_i);
    }

    // Approximate the jacobian for small variations in the cam_j
    // parameters (camera parameters).
    inline Matrix<double, 2, CameraParamsN> cam_jacobian ( size_t i, size_t j,
                                                         Vector<double, CameraParamsN> const& cam_j,
                                                         Vector<double, PointParamsN> const& point_i ) {

      // Jacobian is #outputs x #params
      Matrix<double, 2, CameraParamsN> J;

      Vector2 h0;
      try {
        // Get nominal function value
        h0 = impl().cam_pixel(i,j,cam_j,point_i);
      } catch (const camera::PointToPixelErr& e) {
        // Unable to project this point into the camera, so abort!
        return J;
      }

      // For each param dimension, add epsilon and re-evaluate h() to
      // get numerical derivative w.r.t. that parameter
      for ( size_t n=0; n < CameraParamsN; ++n ){
        Vector<double, CameraParamsN> cam_j_prime = cam_j;

        // Variable step size, depending on parameter value
        double epsilon = 1e-7 + fabs(cam_j(n)*1e-7);
        cam_j_prime(n) += epsilon;

        // Evaluate function with this step and compute the derivative
        // w.r.t. parameter i
        try {
          Vector2 hi = impl().cam_pixel(i,j,cam_j_prime, point_i);
          select_col(J,n) = (hi-h0)/epsilon;
        } catch (const camera::PointToPixelErr& e) {
          select_col(J,n) = Vector2();
        }
      }

      return J;
    }

    // Approximate the jacobian for small variations in the point_i
    // parameters (3d point locations).
    inline Matrix<double, 2, PointParamsN> point_jacobian ( size_t i, size_t j,
                                                        Vector<double, CameraParamsN> const& cam_j,
                                                        Vector<double, PointParamsN> const& point_i ) {

      // Jacobian is #outputs x #params
      Matrix<double, 2, PointParamsN> J;

      Vector2 h0;
      try {
        // Get nominal function value
        h0 = impl().cam_pixel(i,j,cam_j, point_i);
      } catch (const camera::PointToPixelErr& e) {
        // Unable to project this point into the camera so abort!
        return J;
      }

      // For each param dimension, add epsilon and re-evaluate h() to
      // get numerical derivative w.r.t. that parameter
      for ( size_t n=0; n < PointParamsN; ++n ){
        Vector<double, PointParamsN> point_i_prime = point_i;

        // Variable step size, depending on parameter value
        double epsilon = 1e-7 + fabs(point_i(n)*1e-7);
        point_i_prime(n) += epsilon;

        // Evaluate function with this step and compute the derivative
        // w.r.t. parameter i
        try {
          Vector2 hi = impl().cam_pixel(i,j,cam_j,point_i_prime);
          select_col(J,n) = (hi-h0)/epsilon;
        } catch (const camera::PointToPixelErr& e) {
          select_col(J,n) = Vector2();
        }
      }
      return J;
    }

    // -- Report Functions -------------------------------------------

    std::string image_unit          () { return "pixels";  }
    std::string camera_position_unit() { return "meters";  }
    std::string camera_pose_unit    () { return "radians"; }
    std::string gcp_unit            () { return "meters";  }
    
    // Forcing the user to define
    inline double image_compare( Vector2 const& meas,
                                 Vector2 const& obj ) {
      return impl().image_compare( meas, obj );
    }
    inline double position_compare( Vector<double, CameraParamsN> const& meas,
                                    Vector<double, CameraParamsN> const& obj ) {
      return impl().position_compare( meas, obj );
    }
    inline double pose_compare( Vector<double, CameraParamsN> const& meas,
                                Vector<double, CameraParamsN> const& obj) {
      return impl().pose_compare( meas, obj );
    }
    inline double gcp_compare( Vector<double, PointParamsN> const& meas,
                               Vector<double, PointParamsN> const& obj ) {
      return impl().gcp_compare( meas, obj );
    }
    boost::shared_ptr<ControlNetwork> control_network() {
      vw_throw( vw::NoImplErr() << "Programmer needs to implement ModelBase::control_network()\n" );
    }
  };

}} // namespace vw::ba

#endif//__VW_BUNDLEADJUSTMENT_MODEL_BASE_H__
