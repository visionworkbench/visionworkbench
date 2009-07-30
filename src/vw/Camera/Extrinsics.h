// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file Extrinsics.h
/// 
/// Utilities for describing extrinsic camera parameters (position and
/// pose).
///
/// 
#ifndef __VW_CAMERA_EXTRINSICS_H__
#define __VW_CAMERA_EXTRINSICS_H__

#include <vw/Math/Quaternion.h>
#include <vw/Math/LinearAlgebra.h>

#include <iostream>

namespace vw { 
namespace camera {

  // --------------------------------------------------------------------------
  // POSITION INTERPOLATION
  // --------------------------------------------------------------------------

  class LinearPositionInterpolation {
    Vector3 m_initial_position, m_velocity_vector; 

  public:
    LinearPositionInterpolation(Vector3 initial_position, Vector3 velocity_vector) :
      m_initial_position(initial_position), m_velocity_vector(velocity_vector) {}

    Vector3 operator()(double t) const {
      return m_initial_position + t * m_velocity_vector;
    }

  };

  class Curve3DPositionInterpolation {
    Matrix<double,3,3> m_cached_fit;

  public:
    /// The constructor fits a 3D curve to the provided points and
    /// caches the fit so that it can be quickly evaluated using the
    /// call operator.
    ///
    /// This routine will fit a 3d parametric curve (2nd degree
    /// polynomial) to a serios of N points in the Nx3 matrix points.
    /// 
    /// The curve is specified by the coefficient
    /// matrix, which has the form:
    /// 
    ///          | a b c |
    ///   A  =   | d e f |
    ///          | g h i |
    /// 
    /// and the curve is is evaluatied using the equation:
    /// 
    ///              |  1  |
    ///   x(t) = A * |  t  |
    ///              | t^2 |
    ///
    Curve3DPositionInterpolation(std::vector<Vector3> const& position_samples, double t0, double dt) {
      Matrix<double> Z(position_samples.size()*3, 9);
      //    fill(Z, 0.0);
      
      Vector<double> p(position_samples.size() * 3);
      // Reshape the position_samples matrix into a column vector 
      for (unsigned int i = 0; i < position_samples.size(); i++) {
        p(3*i) = position_samples[i][0];
        p(3*i+1) = position_samples[i][1];
        p(3*i+2) = position_samples[i][2];
      }
      
      Vector<double> t(position_samples.size());
      for (unsigned int i = 0; i < t.size(); i++) {
        t(i) = t0 + dt*i;
      }
      
      // Populate the Z matrix 
      for (unsigned int i = 0; i < position_samples.size(); i++) {
        Z(3*i  , 0) = 1.0;
        Z(3*i  , 1) = t(i);
        Z(3*i  , 2) = t(i)*t(i);
        Z(3*i+1, 3) = 1.0;
        Z(3*i+1, 4) = t(i);
        Z(3*i+1, 5) = t(i)*t(i);
        Z(3*i+2, 6) = 1.0;
        Z(3*i+2, 7) = t(i);
        Z(3*i+2, 8) = t(i)*t(i);
      }
      
      Vector<double> x = least_squares(Z,p);
      //    std::cout << Z << "    " << "\n" << p << "\n" << x << "\n";
      Matrix<double,3,3> coeff;
      coeff(0,0) = x(0);  coeff(0,1) = x(1); coeff(0,2) = x(2);
      coeff(1,0) = x(3);  coeff(1,1) = x(4); coeff(1,2) = x(5);
      coeff(2,0) = x(6);  coeff(2,1) = x(7); coeff(2,2) = x(8);
      
      m_cached_fit = coeff;
    }

    /// The call operator evaluates the curve at the given time t.
    /// See the constructor documentation for a summary of the
    /// mathematical operations that are carried out in this function.
    Vector3 operator()(double t) const {
      Vector3 T(1, t, t*t);
      return m_cached_fit * T;
    }
  };


  // --------------------------------------------------------------------------
  // POSE INTERPOLATION
  // --------------------------------------------------------------------------
  
  class ConstantPoseInterpolation {
    Quaternion<double> m_pose;
  public:
    ConstantPoseInterpolation(Quaternion<double> const& pose) : m_pose(pose) {}
    
    inline Quaternion<double> operator()(double /*t*/) const {
      return m_pose;
    }
  };

  
  /// Performs interpolation between sparse pose data points using the
  /// spherical linear interpolation algorithm.
  class SLERPPoseInterpolation {
    std::vector<Quaternion<double> > m_pose_samples;
    double m_t0, m_dt;

    static Quaternion<double> slerp(double alpha, Quaternion<double> const& a, Quaternion<double> const& b, int spin) {
      const double SLERP_EPSILON = 1.0E-6; 	        // a tiny number
      double beta;			// complementary interp parameter 
      double theta;			// angle between A and B 
      double sin_t, cos_t;		// sine, cosine of theta 
      double phi;			// theta plus spins 
      int bflip;			// use negation of B? 
      
      // cosine theta = dot product of A and B 
      cos_t = a(1)*b(1) + a(2)*b(2) + a(3)*b(3) + a(0)*b(0);
      
      // if B is on opposite hemisphere from A, use -B instead 
      if (cos_t < 0.0) {
        cos_t = -cos_t;
        bflip = true;
      } else {
        bflip = false;
      }
      
      // if B is (within precision limits) the same as A,
      // just linear interpolate between A and B.
      // Can't do spins, since we don't know what direction to spin.
      if (1.0 - cos_t < SLERP_EPSILON) {
        beta = 1.0 - alpha;
      } else {				/* normal case */
        theta = acos(cos_t);
        phi = theta + spin * M_PI;
        sin_t = sin(theta);
        beta = sin(theta - alpha*phi) / sin_t;
        alpha = sin(alpha*phi) / sin_t;
      }
      
      if (bflip)
        alpha = -alpha;
      
      // interpolate 
      Quaternion<double> q;
      q(1) = beta*a(1) + alpha*b(1);
      q(2) = beta*a(2) + alpha*b(2);
      q(3) = beta*a(3) + alpha*b(3);
      q(0) = beta*a(0) + alpha*b(0);
      return q;
    }

  public:
    SLERPPoseInterpolation(std::vector<Quaternion<double> > const& pose_samples, double t0, double dt) :
      m_pose_samples(pose_samples), m_t0(t0), m_dt(dt) {}
    
    /// Compute the pose at a given time t.  The pose will be an interpolated value
    Quaternion<double> operator()(double t) const {

      // Make sure that t lies within the range [t0, t0+dt*length(points)] 
      if ((t < m_t0) || (t > m_t0+m_dt*m_pose_samples.size())) {
        std::cout << "Time: " << t << "   min: " << m_t0 << "   max: " << (m_t0+m_dt*m_pose_samples.size()) <<"\n";
        vw_throw( ArgumentErr() << "Cannot extrapolate point for time " << t << ". Out of valid range." );
      }

      unsigned int low_ind = (unsigned int)floor( (t-m_t0) / m_dt );
      unsigned int high_ind = (unsigned int)ceil( (t-m_t0) / m_dt );
      
      // If there are not enough points to interpolate at the end, we
      // will limit the high_ind here.
      if (high_ind > m_pose_samples.size()) {
        vw_throw( ArgumentErr() << "Attempted to interpolate a quaternion past the last available control point." );
      } else if (high_ind == m_pose_samples.size()) {
        high_ind = m_pose_samples.size() - 1;
      }
      
      double low_t =  m_t0 + m_dt * low_ind; 
      double norm_t = (t - low_t)/m_dt;
      
      Quaternion<double> a = m_pose_samples[low_ind];
      Quaternion<double> b = m_pose_samples[high_ind];
      
      return this->slerp(norm_t,a,b,0);
    }

  };

}} // namespace vw::camera

#endif // __VW_CAMERA_EXTRINSICS_H__
