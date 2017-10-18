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


/// \file LensDistortion.h
///
/// This file contains the pinhole camera model.
///
#ifndef __VW_CAMERA_LENSDISTORTION_H__
#define __VW_CAMERA_LENSDISTORTION_H__

#include <vw/Math/Vector.h>

#include <iosfwd>
#include <string>

#include <boost/smart_ptr/shared_ptr.hpp>

namespace vw {
namespace camera {

  // Forward declaration
  class PinholeModel;

  /// Base class which all distortion models inherit from.
  /// - Trivia: Back in 2009 this was implemented using CRTP.  See commit e2749b36d3db37f3176acd8907434dbf4ab29096.
  class VW_API LensDistortion {
  protected:
    std::vector<std::string> m_distortion_param_names;
  public:
    LensDistortion() {}

    virtual ~LensDistortion() {}

    // For the two functions below, default implementations are provided in which a
    //  solver attempts to use the *other* function to find the answer.
    
    /// From an undistorted input coordinate, compute the distorted coordinate.
    /// - The input location is in the same units as the focal length that was provided to 
    ///   the PinholeModel class.
    virtual Vector2 distorted_coordinates  (const PinholeModel&, Vector2 const&) const;
   
    /// From a distorted input coordinate, compute the undistorted coordinate.
    /// - The input location is in the same units as the focal length that was provided to 
    ///   the PinholeModel class.
    virtual Vector2 undistorted_coordinates(const PinholeModel&, Vector2 const&) const;
    
    /// Write all the distortion parameters to the stream
    virtual void write(std::ostream & os) const = 0;
    
    /// Read all the distortion parameters from the stream
    virtual void read(std::istream & os) = 0;
    
    /// Return a pointer to a copy of this distortion object
    virtual boost::shared_ptr<LensDistortion> copy() const = 0;
    
    /// Return a vector containing all the distortion parameters.
    virtual Vector<double> distortion_parameters() const;

    /// Initialize the object from a set of distortion parameters.
    virtual void set_distortion_parameters(Vector<double> const& params);

    /// Each derived model needs to have a string name.
    virtual std::string name() const = 0;
    
    /// Used to scale distortion w/ image size
    virtual void scale(float scale) = 0;
    
    /// Used to scale distortion w/ image size
    std::vector<std::string> distortion_param_names() const { return m_distortion_param_names; }
    
  }; // End class LensDistortion

  /// Write any derived lens distortion class to the stream.
  VW_API std::ostream& operator<<(std::ostream& os, const LensDistortion& ld);


  // ------------------------------------------------------------------------------
  // -- Derived classes section


  /// A NULL lens distortion model.
  struct VW_API NullLensDistortion : public LensDistortion {
    virtual Vector2 distorted_coordinates  (const PinholeModel&, Vector2 const& v) const { return v; }
    virtual Vector2 undistorted_coordinates(const PinholeModel&, Vector2 const& v) const { return v; }

    virtual boost::shared_ptr<LensDistortion> copy() const;
    virtual void write(std::ostream& os) const;
    virtual void read (std::istream& os);
    static  std::string class_name()       { return "NULL";       }
    virtual std::string name      () const { return class_name(); }
    virtual void        scale(float /*scale*/);
  };

  /// TSAI Lens Distortion Model
  ///
  /// For a given set of observed (distorted) pixel coordinates, return the
  /// location where the pixel would have appeared if there were no lens distortion.
  ///
  /// The equations which produce these are:
  /// (u, v) = undistorted coordinates
  /// (u', v') = observed (distorted) coordinates
  /// (x, y) = object coordinates of projected point
  /// r2 = x * x + y * y   -- square of distance from object to primary vector
  /// k1, k2 are radial distortion parameters; p1, p2 are tangential distortion
  /// parameters. principal point is at (cx, cy).
  ///
  /// u' = u + (u - cx) * (k1 * r2 + k2 * r4 + 2 * p1 * y + p2 * (r2/x + 2x))
  /// v' = v + (v - cy) * (k1 * r2 + k2 * r4 + 2 * p2 * x + p1 * (r2/y + 2y))
  ///
  /// k1 is distortion[0], k2 is distortion[1],  p1 is distortion[2], p2 is distortion[3]
  /// The input pixel value is normalivzed before the calculations are performed.
  ///
  /// References: Roger Tsai, A Versatile Camera Calibration Technique for a High-Accuracy 3D
  /// Machine Vision Metrology Using Off-the-shelf TV Cameras and Lenses
  ///
  /// Be careful when you find a camera calibration result with K1/K2, even though the names
  ///  are the same they could be used in a different way than the TSAI model!!!

  class VW_API TsaiLensDistortion : public LensDistortion {
    Vector4 m_distortion;
  public:
    TsaiLensDistortion();
    TsaiLensDistortion(Vector4 const& params);
    virtual Vector<double> distortion_parameters() const;
    virtual void set_distortion_parameters(Vector<double> const& params);
    virtual boost::shared_ptr<LensDistortion> copy() const;

    virtual Vector2 distorted_coordinates(const PinholeModel& cam, Vector2 const& p) const;
    
    virtual void write(std::ostream& os) const;
    virtual void read (std::istream& os);

    static  std::string class_name()       { return "TSAI";       }
    virtual std::string name      () const { return class_name(); }
    virtual void        scale( float scale );
    void init_distortion_param_names();
  };

  /// Brown Conrady Distortion
  ///
  /// Also known as the plumb-bob distortion model as that was how it could be
  /// solved for. This is used for 'old' camera models.
  ///
  /// References:
  /// Decentering Distortion of Lenses - D.C. Brown,
  ///   Photometric Engineering, pages 444-462, Vol. 32, No. 3, 1966
  /// Close-Range Camera Calibration - D.C. Brown,
  ///   Photogrammetric Engineering, pages 855-866, Vol. 37, No. 8, 1971
  class VW_API BrownConradyDistortion : public LensDistortion {
    Vector2 m_principal_point;      // xp, yp
    Vector3 m_radial_distortion;    // K1, K2, K3
    Vector2 m_centering_distortion; // P1, P2
    double  m_centering_angle;      // phi
  public:
    BrownConradyDistortion();
    BrownConradyDistortion( Vector<double> const& params );
    BrownConradyDistortion( Vector<double> const& principal,
                            Vector<double> const& radial,
                            Vector<double> const& centering,
                            double const& angle );

    virtual Vector<double> distortion_parameters() const;
    virtual void set_distortion_parameters(Vector<double> const& params);
    virtual boost::shared_ptr<LensDistortion> copy() const;

    virtual Vector2 undistorted_coordinates(const PinholeModel&, Vector2 const&) const;

    virtual void write(std::ostream& os) const;
    virtual void read (std::istream& os);
    static  std::string class_name()       { return "BrownConrady"; }
    virtual std::string name      () const { return class_name();   }
    virtual void        scale( float /*scale*/ );
    void init_distortion_param_names();
  };

  /// Adjustable Tsai Distortion
  ///
  /// This is another implementation of TSAI but it supports arbitrary
  /// number of radial coefficients ( but only on the even terms ). This
  /// model is also different in that it's math follows what is available in
  /// the Matlab Camera Calibration Tool Box.
  ///
  /// Coefficients are [r2,r4,r6, ...., t1, t2, alpha]. The last 3
  /// elements tangential and alpha(=skew) are always supplied.
  ///
  /// References:
  /// http://www.vision.caltech.edu/bouguetj/calib_doc/htmls/parameters.html
  class VW_API AdjustableTsaiLensDistortion : public LensDistortion {
    Vector<double> m_distortion;
  public:
    AdjustableTsaiLensDistortion() {}
    AdjustableTsaiLensDistortion(Vector<double> params);
    virtual Vector<double> distortion_parameters() const;
    virtual void set_distortion_parameters(Vector<double> const& params);
    virtual boost::shared_ptr<LensDistortion> copy() const;

    virtual Vector2 distorted_coordinates(PinholeModel const&, Vector2 const&) const;

    virtual void write(std::ostream& os) const;
    virtual void read (std::istream& os);
    
    static  std::string class_name()       { return "AdjustableTSAI"; }
    virtual std::string name      () const { return class_name(); }
    virtual void scale( float /*scale*/ );
  };
  
  
  
  /// Photometrix Lens Distortion Model
  ///
  /// This model is similar to the TSAI model above but it differs slightly
  /// to match the conventions used by the Australis software from Photometrix.
  /// This type of calibration was originally seen for the NASA IceBridge cameras.
  ///
  /// Parameters used: c, xp, yp, K1, K2, K3, P1, P2, B1, B2
  ///  - c (focal length) comes from the base class so 
  ///    the parameters stored here are [xp, py, K1, K2, K3, P1, P2, B1, B2]
  ///
  /// As copied from a sample output calibration file:
  ///
  /// x = x(meas) - xp
  /// y = y(meas) - yp
  /// x and y are now with respect to the principle point.
  ///
  /// r2 = x * x + y * y
  /// dr = K1*r3 + K2*r5 + K3*r7
  /// 
  /// x(corr) = x(meas) - xp + x*dr/r + P1*(r2 +2x^2) + 2*P2*x*y
  /// y(corr) = y(meas) - yp + y*dr/r + P2*(r2 +2y^2) + 2*P1*x*y
  ///
  /// k1, k2 are radial distortion parameters; p1, p2 are tangential distortion
  /// parameters. principal point is at (xp, yp). B1 and B2 are not used yet.
  ///
  class VW_API PhotometrixLensDistortion : public LensDistortion {
    Vector<float64,9> m_distortion;
  public:
    PhotometrixLensDistortion();
    PhotometrixLensDistortion(Vector<float64,9> const& params);
    virtual Vector<double> distortion_parameters() const;
    virtual void set_distortion_parameters(Vector<double> const& params);
    virtual boost::shared_ptr<LensDistortion> copy() const;

    virtual Vector2 undistorted_coordinates(const PinholeModel& cam, Vector2 const& p) const;
    
    virtual void write(std::ostream& os) const;
    virtual void read (std::istream& os);

    static  std::string class_name()       { return "Photometrix"; }
    virtual std::string name      () const { return class_name();  }

    virtual void scale( float scale );
    void init_distortion_param_names();
  };
  
  
  

}} // namespace vw::camera

#endif // __VW_CAMERA_LENSDISTORTION_H__
