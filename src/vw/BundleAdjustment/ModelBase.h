// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file BundleAdjustModelBase.h
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
    static const size_t camera_params_n = CameraParamsN;
    static const size_t point_params_n = PointParamsN;

    /// \cond INTERNAL
    // Methods to access the derived type
    inline ImplT& impl() { return static_cast<ImplT&>(*this); }
    inline ImplT const& impl() const { return static_cast<ImplT const&>(*this); }
    /// \endcond

    // Required access to camera
    Vector2 operator() ( size_t i, size_t j,
                         Vector<double,camera_params_n> const& a_j,
                         Vector<double,point_params_n> const& b_i ) {
      return impl()(i,j,a_j,b_i);
    }

    // Approximate the jacobian for small variations in the a_j
    // parameters (camera parameters).
    inline Matrix<double, 2, CameraParamsN> A_jacobian ( size_t i, size_t j,
                                                         Vector<double, CameraParamsN> const& a_j,
                                                         Vector<double, PointParamsN> const& b_i ) {

      // Jacobian is #outputs x #params
      Matrix<double, 2, CameraParamsN> J;

      Vector2 h0;
      try {
        // Get nominal function value
        h0 = impl()(i,j,a_j,b_i);
      } catch ( camera::PixelToRayErr &e ) {
        // Unable to project this point into the camera, so abort!
        return J;
      }

      // For each param dimension, add epsilon and re-evaluate h() to
      // get numerical derivative w.r.t. that parameter
      for ( size_t n=0; n < CameraParamsN; ++n ){
        Vector<double, CameraParamsN> a_j_prime = a_j;

        // Variable step size, depending on parameter value
        double epsilon = 1e-7 + fabs(a_j(n)*1e-7);
        a_j_prime(n) += epsilon;

        // Evaluate function with this step and compute the derivative
        // w.r.t. parameter i
        try {
          Vector2 hi = impl()(i,j,a_j_prime, b_i);
          select_col(J,n) = (hi-h0)/epsilon;
        } catch ( camera::PixelToRayErr &e ) {
          select_col(J,n) = Vector2();
        }
      }

      return J;
    }

    // Approximate the jacobian for small variations in the b_i
    // parameters (3d point locations).
    inline Matrix<double, 2, PointParamsN> B_jacobian ( size_t i, size_t j,
                                                        Vector<double, CameraParamsN> const& a_j,
                                                        Vector<double, PointParamsN> const& b_i ) {

      // Jacobian is #outputs x #params
      Matrix<double, 2, PointParamsN> J;

      Vector2 h0;
      try {
        // Get nominal function value
        h0 = impl()(i,j,a_j, b_i);
      } catch ( camera::PixelToRayErr &e ) {
        // Unable to project this point into the camera so abort!
        return J;
      }

      // For each param dimension, add epsilon and re-evaluate h() to
      // get numerical derivative w.r.t. that parameter
      for ( size_t n=0; n < PointParamsN; ++n ){
        Vector<double, PointParamsN> b_i_prime = b_i;

        // Variable step size, depending on parameter value
        double epsilon = 1e-7 + fabs(b_i(n)*1e-7);
        b_i_prime(n) += epsilon;

        // Evaluate function with this step and compute the derivative
        // w.r.t. parameter i
        try {
          Vector2 hi = impl()(i,j,a_j,b_i_prime);
          select_col(J,n) = (hi-h0)/epsilon;
        } catch ( camera::PixelToRayErr &e ) {
          select_col(J,n) = Vector2();
        }
      }
      return J;
    }

    // -- Report Functions -------------------------------------------

    std::string image_unit() { return "pixels"; }
    std::string camera_position_unit() { return "meters"; }
    std::string camera_pose_unit() { return "radians"; }
    std::string gcp_unit() { return "meters"; }
    // Forced on the user to define
    void image_errors( std::vector<double>& x ) {
      impl().image_errors(x); }
    void camera_position_errors( std::vector<double>& x ) {
      impl().camera_position_errors(x); }
    void camera_pose_errors( std::vector<double>& x) {
      impl().camera_pose_errors(x); }
    void gcp_errors( std::vector<double>& x ) {
      impl().gcp_errors(x); }
    boost::shared_ptr<ControlNetwork> control_network() {
      vw_throw( vw::NoImplErr() << "Programmer needs to implement ModelBase::control_network()\n" );
    }
  };

}} // namespace vw::ba

#endif//__VW_BUNDLEADJUSTMENT_MODEL_BASE_H__
