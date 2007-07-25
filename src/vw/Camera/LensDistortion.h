// 
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
// 
// Copyright 2006 Carnegie Mellon University. All rights reserved.
// 
// This software is distributed under the NASA Open Source Agreement
// (NOSA), version 1.3.  The NOSA has been approved by the Open Source
// Initiative.  See the file COPYING at the top of the distribution
// directory tree for the complete NOSA document.
// 
// THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
// KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
// LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
// SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
// A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
// THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
// DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
// 
// __END_LICENSE__

/// \file PinholeModel.h
/// 
/// This file contains the pinhole camera model.
///
#ifndef __VW_CAMERAMODEL_LENSDISTORTION_H__
#define __VW_CAMERAMODEL_LENSDISTORTION_H__

#include <vw/Camera/CameraModel.h>
#include <vw/Math/Vector.h>
#include <vw/Math/LevenbergMarquardt.h>

#include <boost/shared_ptr.hpp>

namespace vw { 
namespace camera {

  /// Base class for lens distortion models.  Do not derive subclasses
  /// from this base class -- instead, use the CRTP base class below.
  /// This extra layer of indirection allows us to work with pointers
  /// to this base class in the algorithms while still benefitting
  /// from certain features of templates that are provided by the
  /// second base class below.  Take a look at the TsaiLensDistortion
  /// class in PinholeModel.h for an example.
  class LensDistortion {

  public:
    // This method provides the lens distortion model with access to
    // the "parent" camera model class (e.g. if it needs access
    // intrinsic parameters when making the lens distortion
    // calculations, etc.).  Derived lens distortion classes can
    // access the "parent" camera model by calling
    // this->camera_model(), which is defined in the CRTP base class
    // below.  Subclasses should not override this method -- it is
    // defined in the CRTP base class but exists as a virtual function
    // here.
    virtual void set_parent_camera_model(void *model) = 0;

    // The camera model class can use this method to obtain a smart
    // pointer to a copy of the _derived_ class.  This is done through
    // the magic of the CRTP base class below.
    virtual boost::shared_ptr<LensDistortion> copy() const = 0;

    // Override these methods in the subclass to provide basic
    // functionality.
    virtual ~LensDistortion() {}
    virtual Vector2 get_distorted_coordinates(Vector2 const& v) const = 0;

    // This is the default method for computing the undistorted
    // coordinates.  It uses levenberg marquardt to iteratively
    // estimate the undistorted coordinate that are the inverse of
    // calling get_distorted_coordinates().  This is implemented in
    // the CRTP base class below.  Subclasses can override this
    // behavior by providing their own implementation of
    // get_undistorted_coordinates() which will perform faster than
    // the generic implementation based on iterative methods.
    virtual Vector2 get_undistorted_coordinates(Vector2 const& v) const = 0;
    
    
    // This is used in conjuction with a non-member function overloading
    // the << operator
    virtual void write(std::ostream & os) const =0;

  };

  // Why does removing "inline" cause compilation problems?    
  inline std::ostream & operator<<(std::ostream & os, const LensDistortion * ld){
    ld->write(os);
    return os;
  }
  




  /// \cond INTERNAL
  // Optimization functor for computing the undistorted coordinates
  // using levenberg marquardt.
  struct UndistortOptimizeFunctor : public vw::math::LeastSquaresModelBase<UndistortOptimizeFunctor> {
    typedef Vector2 result_type;
    typedef Vector2 domain_type;
    typedef Matrix<double> jacobian_type;
    
    const LensDistortion *m_model_ptr;
    UndistortOptimizeFunctor(const LensDistortion& m) : m_model_ptr(&m) {}
    
    inline result_type operator()( domain_type const& x ) const { 
      return m_model_ptr->get_distorted_coordinates(x);
    }
  };
  
  // This CRTP base class is inserted here in the inheritance chain to
  // provide a point where we can keep track of some templatized
  // values, such as the type of the implementation class and the
  // parent camera model class.  The subclass must specify (via the
  // ParentCameraT template argument), which classes of cameras the
  // lens distortion model applies to.
  template <class ImplT, class ParentCameraT>
  class LensDistortionBase : public LensDistortion {
    ParentCameraT *m_camera_model_ptr;

  public:
    LensDistortionBase() { m_camera_model_ptr = NULL; }
    virtual ~LensDistortionBase() {}
    //    virtual std::ostream& operator <<(std::ostream& os) const = 0;
    virtual void write(std::ostream & os) const {
      vw_throw( NoImplErr() << "LensDistortionBase: write has not been implemented." );
    }

  protected:
    // Subclasses can call this method to gain access to the parent
    // camera model object.
    ParentCameraT& camera_model() const {
      if (!m_camera_model_ptr) { throw LogicErr() << "LensDistortionBase: parent camera model has not yet been specified."; }
      else { return *m_camera_model_ptr; }
    }

    // The camera model class is responsible for calling this method
    // to set the lens distortion model's parent.
    virtual void set_parent_camera_model(void *model) { 
      m_camera_model_ptr = static_cast<ParentCameraT*>(model); 
    }

    // Methods to access the derived type
    inline ImplT& impl() { return static_cast<ImplT&>(*this); }
    inline ImplT const& impl() const { return static_cast<ImplT const&>(*this); }

    // Method to return a shared pointer of a copy of the derived type
    virtual boost::shared_ptr<LensDistortion> copy() const { 
      return boost::shared_ptr<LensDistortion>(new ImplT(impl()));
    }

    // This is the default method for computing the undistorted
    // coordinates.  It uses levenberg marquardt to iteratively
    // estimate the undistorted coordinate that are the inverse of
    // calling get_distorted_coordinates().  Subclasses can override
    // this behavior by providing their own implementation of
    // get_undistorted_coordinates() which will perform faster than
    // this generic iterative implementation.
    virtual Vector2 get_undistorted_coordinates(Vector2 const& p) const {
      UndistortOptimizeFunctor model(*this);
      int status;
      Vector2 solution = levenberg_marquardt( model, p, p, status, 0.1, 0.1 ); // tol = 0.1 pixels
      VW_DEBUG_ASSERT( status != vw::math::optimization::eConvergedRelTolerance,
                       PixelToRayErr() << "get_undistorted_coordinates: failed to converge." );
      return solution;
    }

  };
  /// \endcond

  
  /// A NULL lens distortion model.  
  struct NullLensDistortion : public LensDistortionBase<NullLensDistortion, CameraModel> {
    virtual ~NullLensDistortion() {}
    virtual Vector2 get_distorted_coordinates(Vector2 const& v) const { return v; }
    virtual Vector2 get_undistorted_coordinates(Vector2 const& v) const { return v; }

    virtual void write(std::ostream & os) const {
      os << "k1 = " << 0 << "\n";
      os << "k2 = " << 0 << "\n";
      os << "p1 = " << 0 << "\n";
      os << "p2 = " << 0 << "\n";
    }

    friend std::ostream & operator<<(std::ostream & os, const NullLensDistortion nld) {
      nld.write(os);
      return os;
    }

  };

}} // namespace vw::camera

#endif // __VW_CAMERAMODEL_LENSDISTORTION_H__
