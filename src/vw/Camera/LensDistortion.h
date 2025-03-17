// __BEGIN_LICENSE__
//  Copyright (c) 2006-2025, United States Government as represented by the
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
#ifndef __VW_CAMERA_LENS_DISTORTION_H__
#define __VW_CAMERA_LENS_DISTORTION_H__

#include <vw/Math/Vector.h>

#include <iosfwd>
#include <string>
#include <boost/smart_ptr/shared_ptr.hpp>

namespace vw {
namespace camera {

  // Forward declaration
  class PinholeModel;

  /// Base class which all distortion models inherit from.
  /// Trivia: Back in 2009 this was implemented using CRTP. See commit
  /// e2749b36d3db37f3176acd8907434dbf4ab29096.
  class LensDistortion {
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
    /// Many derived classes override this with a faster implementation.
    virtual Vector2 distorted_coordinates(const PinholeModel& cam, Vector2 const&) const;

    /// From a distorted input coordinate, compute the undistorted coordinate.
    /// - The input location is in the same units as the focal length that was provided to
    ///   the PinholeModel class.
    /// A derived class must reimplement this.
    virtual Vector2 undistorted_coordinates(const PinholeModel& cam, 
                                            Vector2 const& p) const = 0;

    /// Return true if the distorted_coordinates() implementation does not use a solver.
    virtual bool has_fast_distort  () const {return false;}

    /// Return true if the undistorted_coordinates() implementation does not use a solver.
    virtual bool has_fast_undistort() const {return false;}

    /// Write all the distortion parameters to the stream
    virtual void write(std::ostream & os) const = 0;

    /// Read all the distortion parameters from the stream
    virtual void read(std::istream & is) = 0;

    /// Return a pointer to a copy of this distortion object
    virtual boost::shared_ptr<LensDistortion> copy() const = 0;

    /// Return a vector containing all the distortion parameters.
    virtual Vector<double> distortion_parameters() const;

    /// Initialize the object from a set of distortion parameters.
    virtual void set_distortion_parameters(Vector<double> const& params);

    // Number of distortion parameters
    virtual int num_dist_params() const = 0;

    /// Each derived model needs to have a string name.
    virtual std::string name() const = 0;

    /// Used to scale distortion with image size
    virtual void scale(double scale) = 0;

    /// Used to scale distortion with image size
    std::vector<std::string> distortion_param_names() const { return m_distortion_param_names; }
  }; // End class LensDistortion

  /// Write any derived lens distortion class to the stream.
  std::ostream& operator<<(std::ostream& os, const LensDistortion& ld);


  // ------------------------------------------------------------------------------
  // -- Derived classes section


  /// A NULL lens distortion model.
  struct NullLensDistortion: public LensDistortion {
    virtual Vector2 distorted_coordinates  (const PinholeModel& cam, Vector2 const& v) const { return v; }
    virtual Vector2 undistorted_coordinates(const PinholeModel& cam, Vector2 const& v) const { return v; }
    virtual int num_dist_params() const { return 0; };

    virtual bool has_fast_distort  () const {return true;}
    virtual bool has_fast_undistort() const {return true;}

    virtual boost::shared_ptr<LensDistortion> copy() const;
    virtual void write(std::ostream& os) const;
    virtual void read (std::istream& is);
    static  std::string class_name() { return "NULL";       }
    virtual std::string name      () const { return class_name(); }
    virtual void        scale(double /*scale*/);
  };

  /// TSAI Lens Distortion Model
  /// This was validated to be in perfect agreement with the OpenCV implementation.
  /// https://docs.opencv.org/4.x/d9/d0c/group__calib3d.html

  /// The function cv::projectPoints() was used for validation, with no rotation
  /// or translation.

  /// References: Roger Tsai, A Versatile Camera Calibration Technique for a High-Accuracy 3D
  /// Machine Vision Metrology Using Off-the-shelf TV Cameras and Lenses
  ///
  /// Be careful when you find a camera calibration result with K1/K2, even though the names
  /// are the same they could be used in a different way than the TSAI model!!!

  class TsaiLensDistortion: public LensDistortion {
  public:
    TsaiLensDistortion();
    TsaiLensDistortion(Vector<double> const& params);
    virtual Vector<double> distortion_parameters() const;
    virtual void set_distortion_parameters(Vector<double> const& params);
    virtual int num_dist_params() const { return m_distortion.size(); }
    virtual boost::shared_ptr<LensDistortion> copy() const;

    // Apply the distortion to a normalized pixel as function object. To be used in
    // Newton-Raphson.
    vw::Vector2 operator()(vw::Vector2 const& p) const;

    virtual Vector2 distorted_coordinates(const PinholeModel& cam, Vector2 const& p) const;
    virtual Vector2 undistorted_coordinates(const PinholeModel& cam, Vector2 const& p) const;

    virtual bool has_fast_distort  () const {return true;}

    virtual void write(std::ostream& os) const;
    virtual void read (std::istream& is);

    static  std::string class_name() { return "TSAI";       }
    virtual std::string name      () const { return class_name(); }
    virtual void        scale(double scale);
    void init_distortion_param_names();

  private:
    Vector<double> m_distortion;
  };

  // This agrees with rig_calibrator
  class FovLensDistortion: public LensDistortion {
  public:
    FovLensDistortion();
    FovLensDistortion(Vector<double> const& params);
    virtual Vector<double> distortion_parameters() const;
    virtual void set_distortion_parameters(Vector<double> const& params);
    virtual int num_dist_params() const { return m_distortion.size(); }
    virtual boost::shared_ptr<LensDistortion> copy() const;

    virtual Vector2 distorted_coordinates(const PinholeModel& cam, Vector2 const& p) const;
    virtual Vector2 undistorted_coordinates(const PinholeModel& cam, Vector2 const& p) const;

    virtual bool has_fast_distort  () const {return true;}

    virtual void write(std::ostream& os) const;
    virtual void read (std::istream& is);

    static  std::string class_name() { return "FOV";       }
    virtual std::string name      () const { return class_name(); }
    virtual void        scale(double scale);
    void init_distortion_param_names();

  private:
    Vector<double> m_distortion;
    double m_distortion_precalc1, m_distortion_precalc2;
  };

  // Fisheye lens distortion. Agrees with OpenCV and rig_calibrator. When
  // comparing with OpenCV's fisheye distortion function, need to ensure that
  // function is passed in normalized and centered pixels.
  class FisheyeLensDistortion: public LensDistortion {
  public:
    FisheyeLensDistortion();
    FisheyeLensDistortion(Vector<double> const& params);
    virtual Vector<double> distortion_parameters() const;
    virtual void set_distortion_parameters(Vector<double> const& params);
    virtual int num_dist_params() const { return m_distortion.size(); }
    virtual boost::shared_ptr<LensDistortion> copy() const;

    virtual Vector2 distorted_coordinates(const PinholeModel& cam, Vector2 const& p) const;
    virtual Vector2 undistorted_coordinates(const PinholeModel& cam, Vector2 const& p) const;

    // Apply the distortion to a normalized pixel a function object. To be used in
    // Newton-Raphson.
    vw::Vector2 operator()(vw::Vector2 const& p) const;

    virtual bool has_fast_distort  () const {return true;}

    virtual void write(std::ostream& os) const;
    virtual void read (std::istream& is);

    static  std::string class_name() { return "FISHEYE"; }
    virtual std::string name      () const { return class_name(); }
    virtual void        scale(double scale);
    void init_distortion_param_names();

  private:
    Vector<double> m_distortion;
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
  class BrownConradyDistortion: public LensDistortion {
  public:
    static const size_t num_distortion_params = 8;
    BrownConradyDistortion();
    BrownConradyDistortion(Vector<double> const& params);
    BrownConradyDistortion(Vector<double> const& principal,
                            Vector<double> const& radial,
                            Vector<double> const& centering,
                            double const& angle);

    virtual Vector<double> distortion_parameters() const;
    virtual void set_distortion_parameters(Vector<double> const& params);
    virtual int num_dist_params() const { return num_distortion_params; }
    virtual boost::shared_ptr<LensDistortion> copy() const;

    virtual Vector2 undistorted_coordinates(const PinholeModel& cam, Vector2 const& p) const;

    virtual bool has_fast_undistort() const {return true;}

    virtual void write(std::ostream& os) const;
    virtual void read (std::istream& is);
    static  std::string class_name() { return "BrownConrady"; }
    virtual std::string name      () const { return class_name();   }
    virtual void        scale(double /*scale*/);
    void init_distortion_param_names();
  private:
    Vector2 m_principal_point;      // xp, yp
    Vector3 m_radial_distortion;    // K1, K2, K3
    Vector2 m_centering_distortion; // P1, P2
    double  m_centering_angle;      // phi
  };

  /// Adjustable Tsai Distortion
  ///
  /// This is another implementation of TSAI but it supports arbitrary
  /// number of radial coefficients (but only on the even terms). This
  /// model is also different in that it's math follows what is available in
  /// the Matlab Camera Calibration Tool Box.
  ///
  /// Coefficients are [r2,r4,r6, ...., t1, t2, alpha]. The last 3
  /// elements tangential and alpha(=skew) are always supplied.
  ///
  /// References:
  /// http://www.vision.caltech.edu/bouguetj/calib_doc/htmls/parameters.html
  class AdjustableTsaiLensDistortion: public LensDistortion {
    Vector<double> m_distortion;
  public:
    AdjustableTsaiLensDistortion() {}
    AdjustableTsaiLensDistortion(Vector<double> params);
    virtual Vector<double> distortion_parameters() const;
    virtual void set_distortion_parameters(Vector<double> const& params);
    virtual int num_dist_params() const { return m_distortion.size(); }
    virtual boost::shared_ptr<LensDistortion> copy() const;

    // Apply the distortion to a normalized pixel as function object. To be used in
    // Newton-Raphson.
    vw::Vector2 operator()(vw::Vector2 const& p) const;

    virtual Vector2 distorted_coordinates(PinholeModel const&, Vector2 const&) const;
    virtual Vector2 undistorted_coordinates(const PinholeModel& cam, Vector2 const& p) const;

    virtual bool has_fast_distort  () const {return true;}

    virtual void write(std::ostream& os) const;
    virtual void read (std::istream& is);

    static  std::string class_name() { return "AdjustableTSAI"; }
    virtual std::string name      () const { return class_name(); }
    virtual void scale(double /*scale*/);
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
  class PhotometrixLensDistortion: public LensDistortion {
    Vector<double> m_distortion;
  public:
    PhotometrixLensDistortion();
    PhotometrixLensDistortion(Vector<double> const& params);
    virtual Vector<double> distortion_parameters() const;
    virtual void set_distortion_parameters(Vector<double> const& params);
    virtual int num_dist_params() const { return m_distortion.size(); }
    virtual boost::shared_ptr<LensDistortion> copy() const;

    virtual Vector2 undistorted_coordinates(const PinholeModel& cam, Vector2 const& p) const;

    virtual bool has_fast_undistort() const {return true;}

    virtual void write(std::ostream& os) const;
    virtual void read (std::istream& is);

    static  std::string class_name() { return "Photometrix"; }
    virtual std::string name      () const { return class_name();  }

    virtual void scale(double scale);
    void init_distortion_param_names();
  };

  // RPC lens distortion of arbitrary degree.  Undistortion is done
  // analogously using a second set of coefficients.  undistortion
  // parameters are computed.
  // TODO: Make undistortion computation a member of this class.
  class RPCLensDistortion: public LensDistortion {
    int m_rpc_degree;
    Vector2i m_image_size;
    Vector<double> m_distortion;

  public:

    RPCLensDistortion();
    RPCLensDistortion(Vector<double> const& params);
    void reset(int rpc_degree);  // Form the identity transform
    static int rpc_degree(int num_dist_params) {
      return int(round(sqrt(2.0*num_dist_params + 5.0)/2.0 - 1.5));
    }
    int rpc_degree() const { return m_rpc_degree; }
    virtual Vector<double> distortion_parameters() const;
    Vector<double> undistortion_parameters() const;
    void set_image_size(Vector2i const& image_size);
    Vector2i image_size() const { return m_image_size; }
    virtual void set_distortion_parameters(Vector<double> const& params);
    void set_undistortion_parameters(Vector<double> const& params);
    virtual int num_dist_params() const { return m_distortion.size(); }
    static int num_dist_params(int rpc_degree) { return 2*(rpc_degree+1)*(rpc_degree+2)-2; }
    virtual boost::shared_ptr<LensDistortion> copy() const;

    virtual Vector2 distorted_coordinates(const PinholeModel& cam, Vector2 const& p) const;
    virtual Vector2 undistorted_coordinates(const PinholeModel& cam, Vector2 const& p) const;

    // Apply the distortion to a normalized pixel a function object. To be used in
    // Newton-Raphson.
    vw::Vector2 operator()(vw::Vector2 const& p) const;

    virtual bool has_fast_distort  () const {return true;}
    virtual bool has_fast_undistort() const {return true;}

    virtual void write(std::ostream& os) const;
    virtual void read (std::istream& is);

    static  std::string class_name() { return "RPC"; }
    virtual std::string name      () const { return class_name();  }

    virtual void scale(double scale);

    static void init_as_identity(Vector<double> & params);
    static void increment_degree(Vector<double> & params);
  private:
    static void validate_distortion_params(Vector<double> const& params);
    static void unpack_params(Vector<double> const& params,
                              Vector<double> & num_x, Vector<double> & den_x,
                              Vector<double> & num_y, Vector<double> & den_y);
    static void pack_params(Vector<double> & params,
                            Vector<double> const& num_x, Vector<double> const& den_x,
                            Vector<double> const& num_y, Vector<double> const& den_y);
  };

}} // namespace vw::camera

#endif // __VW_CAMERA_LENS_DISTORTION_H__
