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

  // forward decl
  class PinholeModel;

  class LensDistortion {
  public:
    LensDistortion();

    virtual ~LensDistortion();
    virtual Vector2 distorted_coordinates(const PinholeModel&, Vector2 const&) const;
    virtual Vector2 undistorted_coordinates(const PinholeModel&, Vector2 const&) const;
    virtual void write(std::ostream & os) const = 0;
    virtual boost::shared_ptr<LensDistortion> copy() const = 0;
    virtual Vector<double> distortion_parameters() const;

    virtual std::string name() const = 0;
    virtual void scale(float scale) = 0; // Used to scale distortion w/ image size
  };

  std::ostream& operator<<(std::ostream& os, const LensDistortion& ld);

  /// A NULL lens distortion model.
  struct NullLensDistortion : public LensDistortion {
    inline Vector2 distorted_coordinates(const PinholeModel&, Vector2 const& v) const { return v; }
    inline Vector2 undistorted_coordinates(const PinholeModel&, Vector2 const& v) const { return v; }

    boost::shared_ptr<LensDistortion> copy() const;
    void write(std::ostream & os) const;
    std::string name() const;
    void scale(float /*scale*/);
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
  ///
  /// References: Roger Tsai, A Versatile Camera Calibration Technique for a High-Accuracy 3D
  /// Machine Vision Metrology Using Off-the-shelf TV Cameras and Lenses

  class TsaiLensDistortion : public LensDistortion {
    Vector4 m_distortion;
  public:
    TsaiLensDistortion(Vector4 const& params);
    Vector<double> distortion_parameters() const;
    boost::shared_ptr<LensDistortion> copy() const;

    //  Location where the given pixel would have appeared if there were no lens distortion.
    Vector2 distorted_coordinates(const PinholeModel&, Vector2 const&) const;
    void write(std::ostream & os) const;

    std::string name() const;

    void scale( float scale );
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
  class BrownConradyDistortion : public LensDistortion {
    Vector2 m_principal_point;
    Vector3 m_radial_distortion;
    Vector2 m_centering_distortion;
    double m_centering_angle;
  public:
    BrownConradyDistortion( Vector<double> const& params );
    BrownConradyDistortion( Vector<double> const& principal,
                            Vector<double> const& radial,
                            Vector<double> const& centering,
                            double const& angle );

    Vector<double> distortion_parameters() const;
    boost::shared_ptr<LensDistortion> copy() const;

    Vector2 undistorted_coordinates(const PinholeModel&, Vector2 const&) const;

    void write(std::ostream& os) const;
    std::string name() const;
    void scale( float /*scale*/ );
  };

  /// Adjustable Tsai Distortion
  ///
  /// This is another implementation of TSAI but it supports arbitrary
  /// number of radial coefficients ( but only on the even terms ). This
  /// model is also different in that it's math follows what is available in
  /// the Matlab Camera Calibration Tool Box.
  ///
  /// Coefficients are [r2,r4,r6, ...., t1, t2, alpha]. The last 3
  /// elements tangential and alpha are always supplied.
  ///
  /// References:
  /// http://www.vision.caltech.edu/bouguetj/calib_doc/htmls/parameters.html
  class AdjustableTsaiLensDistortion : public LensDistortion {
    Vector<double> m_distortion;
  public:
    AdjustableTsaiLensDistortion(Vector<double> params);
    Vector<double> distortion_parameters() const;
    boost::shared_ptr<LensDistortion> copy() const;

    //  Location where the given pixel would have appeared if there were no lens distortion.
    Vector2 distorted_coordinates(PinholeModel const&, Vector2 const&) const;

    void write(std::ostream & os) const;
    std::string name() const;
    void scale( float /*scale*/ );
  };

}} // namespace vw::camera

#endif // __VW_CAMERA_LENSDISTORTION_H__
