// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

/// \file BundleAdjustModelBase.h
/// 
/// Classes used to define the bundle adjustment problem. Every
/// program implementing BA, must implement one of these.

#ifndef __VW_CAMERA_BUNDLE_ADJUST_MODEL_BASE_H__
#define __VW_CAMERA_BUNDLE_ADJUST_MODEL_BASE_H__

// Standard
#include <string>

// Vision Workbench
#include <vw/Math/Matrix.h>
#include <vw/Math/Vector.h>
#include <vw/Camera/ControlNetwork.h>
#include <vw/Core/Log.h>

// Boost
#include <boost/smart_ptr.hpp>

namespace vw {
namespace camera {

  // CRTP Base class for Bundle Adjustment functors.
  //--------------------------------------------------------- 
  template <class ImplT, unsigned CameraParamsN, unsigned PointParamsN>
  struct BundleAdjustmentModelBase {

    static const unsigned camera_params_n = CameraParamsN;
    static const unsigned point_params_n = PointParamsN;

    /// \cond INTERNAL
    // Methods to access the derived type
    inline ImplT& impl() { return static_cast<ImplT&>(*this); }
    inline ImplT const& impl() const { return static_cast<ImplT const&>(*this); }
    /// \endcond

    virtual ~BundleAdjustmentModelBase() {}
    
    // Required access to camera
    virtual Vector2 operator() ( unsigned i, unsigned j,
				 Vector<double,camera_params_n> const& a_j,
				 Vector<double,point_params_n> const& b_j ) const = 0;

    // Approximate the jacobian for small variations in the a_j
    // parameters (camera parameters). 
    inline Matrix<double, 2, CameraParamsN> A_jacobian ( unsigned i, unsigned j,
                                                         Vector<double, CameraParamsN> const& a_j, 
                                                         Vector<double, PointParamsN> const& b_i ) {

      // Get nominal function value
      Vector2 h0 = impl()(i,j,a_j,b_i);

      // Jacobian is #outputs x #params
      Matrix<double, 2, CameraParamsN> J;

      // For each param dimension, add epsilon and re-evaluate h() to
      // get numerical derivative w.r.t. that parameter
      for ( unsigned n=0; n < CameraParamsN; ++n ){
        Vector<double, CameraParamsN> a_j_prime = a_j;

        // Variable step size, depending on parameter value
        double epsilon = 1e-7 + fabs(a_j(n)*1e-7);
        a_j_prime(n) += epsilon;

        // Evaluate function with this step and compute the derivative
        // w.r.t. parameter i
        Vector2 hi = impl()(i,j,a_j_prime, b_i);
        select_col(J,n) = (hi-h0)/epsilon;
      }

      return J;
    }

    // Approximate the jacobian for small variations in the b_i
    // parameters (3d point locations). 
    inline Matrix<double, 2, PointParamsN> B_jacobian ( unsigned i, unsigned j,
                                                        Vector<double, CameraParamsN> const& a_j, 
                                                        Vector<double, PointParamsN> const& b_i ) {
      // Get nominal function value
      Vector2 h0 = impl()(i,j,a_j, b_i);

      // Jacobian is #outputs x #params
      Matrix<double, 2, PointParamsN> J;

      // For each param dimension, add epsilon and re-evaluate h() to
      // get numerical derivative w.r.t. that parameter
      for ( unsigned n=0; n < PointParamsN; ++n ){
        Vector<double, PointParamsN> b_i_prime = b_i;

        // Variable step size, depending on parameter value
        double epsilon = 1e-7 + fabs(b_i(n)*1e-7);
        b_i_prime(n) += epsilon;

        // Evaluate function with this step and compute the derivative
        // w.r.t. parameter i
        Vector2 hi = impl()(i,j,a_j,b_i_prime);
        select_col(J,n) = (hi-h0)/epsilon;
      }
      return J;
    }

    // Report functions
    virtual std::string image_unit( void ) { return "pixels"; }
    virtual std::string camera_position_unit( void ) { return "meters"; }
    virtual std::string camera_pose_unit( void ) { return "radians"; }
    virtual std::string gcp_unit( void ) { return "meters"; }
    
    // Forced on the user to define
    virtual void image_errors(std::vector<double>&) = 0;
    virtual void camera_position_errors(std::vector<double>&) = 0;
    virtual void camera_pose_errors(std::vector<double>&) = 0;
    virtual void gcp_errors(std::vector<double>&) = 0;
    virtual boost::shared_ptr<ControlNetwork> control_network(void) = 0;
  };

}} // namespace vw::camera

#endif//__VW_CAMERA_BUNDLE_ADJUST_MODEL_BASE_H__
