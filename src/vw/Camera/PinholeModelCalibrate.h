// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file PinholeModelCalibrate.h
/// 
/// This file contains the methods for calibrating a pinhole camera model based on a set of known 3D points a
/// and corresponding pixel coordinates.
///
#ifndef __VW_CAMERAMODEL_PINHOLE_CALIBRATE_H__
#define __VW_CAMERAMODEL_PINHOLE_CALIBRATE_H__

#include <vw/config.h>

#include <vw/Math/RANSAC.h>
#include <vw/Math/Vector.h>
#include <vw/Camera/PinholeModel.h>

namespace vw { 
namespace camera {

namespace {
/// \brief Intrinsic PinholeModel parameter extraction class
///
/// This class converts from/to a PinholeModel to/from a Vector of doubles
/// containing the focal lengths f_u, f_v and camera center coordinates c_u, c_v
///
/// Pass this class as template parameter to serialize_pinholemodel() or deserialize_pinholemodel()
class PinholeModelSerializeIntrinsic {
public:
  static const unsigned int size = 4; /// number of parameters this class will extract

  /// serialize the passed PinholeModel into a Vector of doubles
  static inline const vw::Vector<double, size> get(const PinholeModel& m)  {
    vw::Vector<double, size> data;
    m.intrinsic_parameters(data(0), data(1), data(2), data(3));
    return data;
  }
  
  /// Set the parameters of the PinholeModel to the values contained in the Vector
  static inline void set(PinholeModel& m, const vw::Vector<double, size>& v) {
    m.set_intrinsic_parameters(v(0), v(1), v(2), v(3));
  }
};

/// \brief Extrinsic Rotation parameter extraction class
///
/// This class converts from/to a PinholeModel to/from a Vector of doubles
/// containing the rotation vector components r_x, r_y and r_z
///
/// Pass this class as template parameter to serialize_pinholemodel() or deserialize_pinholemodel()
class PinholeModelSerializeRotation {
public:
  static const unsigned int size = 3; /// number of parameters this class will extract

  /// Serialize the passed PinholeModel into a Vector of doubles
  static inline const vw::Vector<double, size> get(const PinholeModel& m) {
    vw::Vector<double, size> data = m.camera_pose().axis_angle();
    return data;
  }

  /// Set the parameters of the PinholeModel to the values contained in the Vector
  static inline void set(PinholeModel& m, const vw::Vector<double, size>& v) {
    m.set_camera_pose( axis_angle_to_quaternion( v ));
  }
};

/// \brief Extrinsic Translation parameter extraction class
///
/// This class converts from/to a PinholeModel to/from a Vector of doubles
/// containing the translations t_x, t_y and t_z
///
/// Pass this class as template parameter to serialize_pinholemodel() or deserialize_pinholemodel()
class PinholeModelSerializeTranslation {
public:
  static const unsigned int size = 3; /// number of parameters this class will extract

  /// Serialize the passed PinholeModel into a Vector of doubles
  static inline const vw::Vector<double, size> get(const PinholeModel& m) {
    vw::Vector<double, size> data = m.camera_center(); 
    return data;
  }

  /// Set the parameters of the PinholeModel to the values contained in the Vector
  static inline void set(PinholeModel& m, const vw::Vector<double, size>& v) {
    m.set_camera_center(v);
  }
};

/// \brief Tsai lens distortion parameter extraction class
///
/// This class converts from/to a PinholeModel to/from a Vector of four doubles
/// containing the TSAI distortion parameters k_1, k_2, p_1 and p_2
///
/// Pass this class as template parameter to serialize_pinholemodel() or deserialize_pinholemodel()
class PinholeModelSerializeTSAI {
public:
  static const unsigned int size = 4; /// number of parameters this class will extract

  /// Serialize the passed PinholeModel into a Vector of doubles
  static inline const vw::Vector<double, size> get(const PinholeModel& m) {
    try {
      return dynamic_cast<vw::camera::TsaiLensDistortion&>(*m.lens_distortion()).distortion_parameters();
    }
    catch(const std::bad_cast& e) {
      vw::vw_out(vw::DebugMessage, "Camera") <<  "Bad lens distortion model cast in PinholeModelSerializeTSAI::get()" << std::endl;
    }

    return vw::Vector<double, size>();
  }

  /// Set the parameters of the PinholeModel to the values contained in the Vector
  static inline void set(PinholeModel& m, const vw::Vector<double, size>& v) {
    vw::camera::TsaiLensDistortion t(v);
    m.set_lens_distortion(t); // set_lens_distortion creates a copy of t on the heap, no memory problem here...
  }
};

/// \brief PinholeModel parameter extraction class doing nothing
///
/// This class is for internal use only in order to avoid code quadruplication
class PinholeModelSerializeNull {
public:
  static const unsigned int size = 0; /// number of parameters this class will extract

  /// serialize the passed PinholeModel into a Vector of doubles
  static inline const vw::Vector<double, size> get(const PinholeModel& m)  {
    return vw::Vector<double, size>();
  }
  
  /// Set the parameters of the PinholeModel to the values contained in the Vector
  static inline void set(PinholeModel& m, const vw::Vector<double, size>& v) {
  }
};
} /* anonymous namespace */

/// \brief Serialize PinholeModel using the specified serialization classes
///
/// Specialize on one or several of PinholeModelSerializeIntrinsic, PinholeModelSerializeRotation, PinholeModelSerializeTranslation,
/// PinholeModelSerializeTSAI or a custom serialization class
template<class Ex1T, class Ex2T, class Ex3T, class Ex4T>
inline const vw::math::Vector<double, Ex1T::size + Ex2T::size + Ex3T::size + Ex4T::size> serialize_pinholemodel(const PinholeModel& m) {
  vw::Vector<double, Ex1T::size + Ex2T::size + Ex3T::size + Ex4T::size> res;
  vw::math::subvector(res, 0, Ex1T::size) = Ex1T::get(m);
  vw::math::subvector(res, Ex1T::size, Ex2T::size) = Ex2T::get(m);
  vw::math::subvector(res, Ex1T::size + Ex2T::size, Ex3T::size) = Ex3T::get(m);
  vw::math::subvector(res, Ex1T::size + Ex2T::size + Ex3T::size, Ex4T::size) = Ex4T::get(m);
  return res;
}

/// \brief Serialize PinholeModel using the specified serialization classes
///
/// Specialize on one or several of PinholeModelSerializeIntrinsic, PinholeModelSerializeRotation, PinholeModelSerializeTranslation,
/// PinholeModelSerializeTSAI or a custom serialization class
template<class Ex1T, class Ex2T, class Ex3T>
inline const vw::math::Vector<double, Ex1T::size + Ex2T::size + Ex3T::size> serialize_pinholemodel(const PinholeModel& m) {
  vw::Vector<double, Ex1T::size + Ex2T::size + Ex3T::size> res;
  vw::math::subvector(res, 0, Ex1T::size) = Ex1T::get(m);
  vw::math::subvector(res, Ex1T::size, Ex2T::size) = Ex2T::get(m);
  vw::math::subvector(res, Ex1T::size + Ex2T::size, Ex3T::size) = Ex3T::get(m);
  return res;
}

/// \brief Serialize PinholeModel using the specified serialization classes
///
/// Specialize on one or several of PinholeModelSerializeIntrinsic, PinholeModelSerializeRotation, PinholeModelSerializeTranslation,
/// PinholeModelSerializeTSAI or a custom serialization class
template<class Ex1T, class Ex2T>
inline const vw::math::Vector<double, Ex1T::size + Ex2T::size> serialize_pinholemodel(const PinholeModel& m) {
  vw::Vector<double, Ex1T::size + Ex2T::size> res;
  vw::math::subvector(res, 0, Ex1T::size) = Ex1T::get(m);
  vw::math::subvector(res, Ex1T::size, Ex2T::size) = Ex2T::get(m);
  return res;
}

/// \brief Serialize PinholeModel using the specified serialization classes
///
/// Specialize on one or several of PinholeModelSerializeIntrinsic, PinholeModelSerializeRotation, PinholeModelSerializeTranslation,
/// PinholeModelSerializeTSAI or a custom serialization class
template<class Ex1T>
inline const vw::math::Vector<double, Ex1T::size> serialize_pinholemodel(const PinholeModel& m) {
  return Ex1T::get(m);
}

/// \brief Deserialize PinholeModel using the specified serialization classes
///
/// Specialize on one or several of PinholeModelSerializeIntrinsic, PinholeModelSerializeRotation, PinholeModelSerializeTranslation,
/// PinholeModelSerializeTSAI or a custom serialization class
template<class Ex1T, class Ex2T, class Ex3T, class Ex4T>
inline void deserialize_pinholemodel(PinholeModel& m, const vw::math::Vector<double, Ex1T::size + Ex2T::size + Ex3T::size + Ex4T::size>& in) {
  Ex1T::set(m, vw::math::subvector(in, 0, Ex1T::size));
  Ex2T::set(m, vw::math::subvector(in, Ex1T::size, Ex2T::size));
  Ex3T::set(m, vw::math::subvector(in, Ex1T::size + Ex2T::size, Ex3T::size));
  Ex4T::set(m, vw::math::subvector(in, Ex1T::size + Ex2T::size + Ex3T::size, Ex4T::size));
}

/// \brief Deserialize PinholeModel using the specified serialization classes
///
/// Specialize on one or several of PinholeModelSerializeIntrinsic, PinholeModelSerializeRotation, PinholeModelSerializeTranslation,
/// PinholeModelSerializeTSAI or a custom serialization class
template<class Ex1T, class Ex2T, class Ex3T>
inline void deserialize_pinholemodel(PinholeModel& m, const vw::math::Vector<double, Ex1T::size + Ex2T::size + Ex3T::size>& in) {
  Ex1T::set(m, vw::math::subvector(in, 0, Ex1T::size));
  Ex2T::set(m, vw::math::subvector(in, Ex1T::size, Ex2T::size));
  Ex3T::set(m, vw::math::subvector(in, Ex1T::size + Ex2T::size, Ex3T::size));
}

/// \brief Deserialize PinholeModel using the specified serialization classes
///
/// Specialize on one or several of PinholeModelSerializeIntrinsic, PinholeModelSerializeRotation, PinholeModelSerializeTranslation,
/// PinholeModelSerializeTSAI or a custom serialization class
template<class Ex1T, class Ex2T>
inline void deserialize_pinholemodel(PinholeModel& m, const vw::math::Vector<double, Ex1T::size + Ex2T::size>& in) {
  Ex1T::set(m, vw::math::subvector(in, 0, Ex1T::size));
  Ex2T::set(m, vw::math::subvector(in, Ex1T::size, Ex2T::size));
}

/// \brief Deserialize PinholeModel using the specified serialization classes
///
/// Specialize on one or several of PinholeModelSerializeIntrinsic, PinholeModelSerializeRotation, PinholeModelSerializeTranslation,
/// PinholeModelSerializeTSAI or a custom serialization class
template<class Ex1T>
inline void deserialize_pinholemodel(PinholeModel& m, const vw::math::Vector<double, Ex1T::size>& in) {
  Ex1T::set(m, vw::math::subvector(in, 0, Ex1T::size));
}

namespace {
/// The LeastSquaresModel for the LevenbergMarquardt method used for calibrating the potentially 
/// nonlinear (tangential and radial distortion-corrected) camera model
///
/// The template arguments allow customization of the camera (de)serialization process in the optimization
/// and thus specification of which parameters to optimize and which parameters to leave untouched
template<class Ex1T, class Ex2T, class Ex3T, class Ex4T>
class PinholeModelLeastSquares : public vw::math::LeastSquaresModelBase<PinholeModelLeastSquares<Ex1T, Ex2T, Ex3T, Ex4T> > 
{
  PinholeModel& m_camera;
  const std::vector<vw::Vector2>& m_pixels;
  const std::vector<vw::Vector3>& m_points;
public:
  typedef vw::math::Vector<double> result_type;
  typedef vw::math::Vector<double> domain_type;
  typedef vw::math::Matrix<double> jacobian_type;

  inline PinholeModelLeastSquares(PinholeModel& camera, const std::vector<vw::Vector2>& pixels, const std::vector<vw::Vector3>& points) 
    : vw::math::LeastSquaresModelBase<PinholeModelLeastSquares<Ex1T, Ex2T, Ex3T, Ex4T> >(), 
      m_camera(camera), m_pixels(pixels), m_points(points) 
  {}

  inline result_type operator()( domain_type const& x )  const {
    deserialize_pinholemodel<Ex1T, Ex2T, Ex3T, Ex4T>(m_camera, x);

    result_type diff(m_points.size()*2);
    for (vw::uint32 i = 0; i < m_points.size(); i++) {
      vw::Vector2 proj(m_camera.point_to_pixel(m_points[i]));	
      diff[2*i] = proj.x();
      diff[2*i+1] = proj.y();
    }      
    return diff;
  }
};
} /* anonymous namespace */

/// Calibrate using the given camera model as starting point and using the given (de)serialization functors - thus some camera parameters
/// can be hold fixed 
///
/// \param Ex1T Serialization one of PinholeModelSerializeIntrinsic, PinholeModelSerializeRotation, PinholeModelSerializeTranslation,
///         PinholeModelSerializeTSAI or a custom serialization class
/// \param m Initial camera model guess
/// \param pixels Pixel locations in image plane
/// \param points Points in 3D euclidian space
/// \param niter Maximum number of Levenberg Marquardt iterations
template<class Ex1T, class Ex2T, class Ex3T, class Ex4T> 
void pinholemodel_calibrate(PinholeModel& m, const std::vector<vw::Vector2>& pixels, const std::vector<vw::Vector3>& points, unsigned int niter = 1000) {  
  typedef typename PinholeModelLeastSquares<Ex1T, Ex2T, Ex3T, Ex4T>::domain_type domain_type;
  typedef typename PinholeModelLeastSquares<Ex1T, Ex2T, Ex3T, Ex4T>::result_type result_type;

  PinholeModelLeastSquares<Ex1T, Ex2T, Ex3T, Ex4T> model(m, pixels, points);

  // reference pixel coordinates - function optimization goal
  result_type observations(pixels.size()*2);
  for (vw::uint32 i = 0; i < pixels.size(); i++) {
    observations[2*i] = pixels[i].x();
    observations[2*i+1] = pixels[i].y();
  }
    
  domain_type seed(serialize_pinholemodel<Ex1T, Ex2T, Ex3T, Ex4T>(m));
 
  int status = 0;
  domain_type result = vw::math::levenberg_marquardt(model, serialize_pinholemodel<Ex1T, Ex2T, Ex3T, Ex4T>(m), 
						     observations, status, 1e-16, 1e-16, niter);
    
  vw::vw_out(vw::DebugMessage, "Camera") << "Nonlinear calibration found: camera_parameter_vector " << result << std::endl;
}

/// Convinience overload with less template parameters
template<class Ex1T, class Ex2T, class Ex3T> 
inline void pinholemodel_calibrate(PinholeModel& m, const std::vector<vw::Vector2>& pixels, const std::vector<vw::Vector3>& points, unsigned int niter = 1000) {
  pinholemodel_calibrate<Ex1T, Ex2T, Ex3T, PinholeModelSerializeNull>(m, pixels, points);
}
  
/// Convinience overload with less template parameters
template<class Ex1T, class Ex2T> 
inline void pinholemodel_calibrate(PinholeModel& m, const std::vector<vw::Vector2>& pixels, const std::vector<vw::Vector3>& points, unsigned int niter = 1000) {
  pinholemodel_calibrate<Ex1T, Ex2T, PinholeModelSerializeNull>(m, pixels, points);
}

/// Convinience overload with less template parameters
template<class Ex1T> 
inline void pinholemodel_calibrate(PinholeModel& m, const std::vector<vw::Vector2>& pixels, const std::vector<vw::Vector3>& points, unsigned int niter = 1000) {
  pinholemodel_calibrate<Ex1T, PinholeModelSerializeNull>(m, pixels, points);
}

namespace { 
/// The RANSAC fitting functor
template<class Ex1T, class Ex2T, class Ex3T, class Ex4T> 
class PinholeModelRansacFitting {
  PinholeModel m_cam;
  int m_levenberg_marquardt_iter;
public:
  typedef PinholeModel result_type;

  PinholeModelRansacFitting(const PinholeModel& cam, int levenberg_marquardt_iter) : m_cam(cam), m_levenberg_marquardt_iter(levenberg_marquardt_iter)
  {}

  PinholeModel operator()(const std::vector<vw::Vector2>& pixels, const std::vector<vw::Vector3>& points) {
    PinholeModel tmp_cam(m_cam);
    pinholemodel_calibrate<Ex1T, Ex2T, Ex3T, Ex4T>(tmp_cam, pixels, points, m_levenberg_marquardt_iter);
    return tmp_cam;
  }

  template <class DummyT>
  unsigned int min_elements_needed_for_fit(const DummyT& d) {
    return (Ex1T::size + Ex2T::size + Ex3T::size + Ex4T::size + 1)/2;
  }
};

/// The RANSAC error functor
class PinholeModelRansacError {
public:
  double operator()(const PinholeModel& cam, const vw::Vector2& pixel, const vw::Vector3& point) const {
    return vw::math::norm_2( pixel - cam.point_to_pixel(point));
  }
};
} /* anonymous namespace */
  
/// Calibrate using the given camera model as starting point and using the given (de)serialization functors - thus some camera parameters
/// can be held fixed 
///
/// This method uses RANSAC and runs pinholemodel_calibrate at every iteration to find the best camera model for the 
/// set of points selected by RANSAC
///
/// \param Ex1T Serialization one of PinholeModelSerializeIntrinsic, PinholeModelSerializeRotation, PinholeModelSerializeTranslation,
///         PinholeModelSerializeTSAI or a custom serialization class
/// \param m Initial camera model guess
/// \param pixels Pixel locations in image plane
/// \param points Points in 3D euclidian space
/// \parm inlier_threshold Distance in image space below which pixels are considered inliers
/// \parm ransac_iter number of RANSAC iterations, 0 uses the default value specified in the RANSAC implementation (currently twice the pixel/point vector size)
/// \parm lm_iter number of Levenberg Marquardt iterations to run at every RANSAC iteration
/// \return indices of points and pixels RANSAC considers to be inliers
template<class Ex1T, class Ex2T, class Ex3T, class Ex4T> 
std::vector<int> pinholemodel_calibrate_ransac(PinholeModel& m, const std::vector<vw::Vector2>& pixels, const std::vector<vw::Vector3>& points, double inlier_threshold = 10, int ransac_iter = 0, double lm_iter = 100) {
  PinholeModelRansacFitting<Ex1T, Ex2T, Ex3T, Ex4T> fitting(m, lm_iter);
  PinholeModelRansacError error;
  
  vw::math::RandomSampleConsensus<PinholeModelRansacFitting<Ex1T, Ex2T, Ex3T, Ex4T>, PinholeModelRansacError> ransac(fitting, error, inlier_threshold);
  m = ransac(pixels, points, ransac_iter);
  
  std::vector<int> inlier_indices(ransac.inlier_indices(m, pixels, points));
  vw::vw_out(vw::DebugMessage, "Camera") << "RANSAC classified as inliers: " << inlier_indices.size() << '/' << pixels.size() << std::endl;
  return inlier_indices;
}

/// Convinience overload with less template parameters
template<class Ex1T, class Ex2T, class Ex3T> 
inline std::vector<int> pinholemodel_calibrate_ransac(PinholeModel& m, const std::vector<vw::Vector2>& pixels, const std::vector<vw::Vector3>& points, double inlier_threshold = 10, int ransac_iter = 0, double lm_iter = 100) {
  return pinholemodel_calibrate_ransac<Ex1T, Ex2T, Ex3T, PinholeModelSerializeNull>(m, pixels, points, inlier_threshold, ransac_iter, lm_iter);
}
  
/// Convinience overload with less template parameters
template<class Ex1T, class Ex2T> 
inline std::vector<int> pinholemodel_calibrate_ransac(PinholeModel& m, const std::vector<vw::Vector2>& pixels, const std::vector<vw::Vector3>& points, double inlier_threshold = 10, int ransac_iter = 0, double lm_iter = 100) {  
  return pinholemodel_calibrate_ransac<Ex1T, Ex2T, PinholeModelSerializeNull>(m, pixels, points, inlier_threshold, ransac_iter, lm_iter);
}

/// Convinience overload with less template parameters
template<class Ex1T> 
inline std::vector<int> pinholemodel_calibrate_ransac(PinholeModel& m, const std::vector<vw::Vector2>& pixels, const std::vector<vw::Vector3>& points, double inlier_threshold = 10, int ransac_iter = 0, double lm_iter = 100) {
  return pinholemodel_calibrate_ransac<Ex1T, PinholeModelSerializeNull>(m, pixels, points, inlier_threshold, ransac_iter, lm_iter);
}

}}

#endif /*  __VW_CAMERAMODEL_PINHOLE_CALIBRATE_H__ */
