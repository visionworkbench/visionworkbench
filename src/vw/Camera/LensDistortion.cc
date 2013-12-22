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


#include <vw/Camera/LensDistortion.h>
#include <vw/Camera/PinholeModel.h>
#include <vw/Math/LevenbergMarquardt.h>

using namespace vw;
using namespace camera;

// Special LMA Models to figure out foward and backward ---------

// Optimization functor for computing the undistorted coordinates
// using levenberg marquardt.
struct UndistortOptimizeFunctor : public math::LeastSquaresModelBase<UndistortOptimizeFunctor> {
  typedef Vector2 result_type;
  typedef Vector2 domain_type;
  typedef Matrix<double> jacobian_type;

  const camera::PinholeModel& m_cam;
  const camera::LensDistortion &m_distort;
  UndistortOptimizeFunctor(const camera::PinholeModel& cam, const camera::LensDistortion& d) : m_cam(cam), m_distort(d) {}

  inline result_type operator()( domain_type const& x ) const {
    return m_distort.distorted_coordinates(m_cam, x);
  }
};

struct DistortOptimizeFunctor :  public math::LeastSquaresModelBase<DistortOptimizeFunctor> {
  typedef Vector2 result_type;
  typedef Vector2 domain_type;
  typedef Matrix<double> jacobian_type;

  const camera::PinholeModel& m_cam;
  const camera::LensDistortion &m_distort;
  DistortOptimizeFunctor(const camera::PinholeModel& cam, const camera::LensDistortion& d) : m_cam(cam), m_distort(d) {}
  inline result_type operator()( domain_type const& x ) const {
    return m_distort.undistorted_coordinates(m_cam, x);
  }
};

// Backup implemenations for Lens Distortion -------------------

LensDistortion::LensDistortion() {}

LensDistortion::~LensDistortion() {}

Vector<double>
LensDistortion::distortion_parameters() const { return Vector<double>(); }

Vector2
LensDistortion::undistorted_coordinates(const camera::PinholeModel& cam, Vector2 const& v) const {
  UndistortOptimizeFunctor model(cam, *this);
  int status;
  Vector2 solution =
    math::levenberg_marquardt( model, v, v, status, 1e-6, 1e-6, 50 );
  VW_DEBUG_ASSERT( status != math::optimization::eConvergedRelTolerance, PixelToRayErr() << "undistorted_coordinates: failed to converge." );
  return solution;
}

vw::Vector2
LensDistortion::distorted_coordinates(const camera::PinholeModel& cam, Vector2 const& v) const {
  DistortOptimizeFunctor model(cam, *this);
  int status;
  Vector2 solution =
    math::levenberg_marquardt( model, v, v, status, 1e-6, 1e-6, 50 );
  VW_DEBUG_ASSERT( status != math::optimization::eConvergedRelTolerance, PixelToRayErr() << "distorted_coordinates: failed to converge." );
  return solution;
}

// Specific Implementations -------------------------------------

boost::shared_ptr<LensDistortion> NullLensDistortion::copy() const {
  return boost::shared_ptr<NullLensDistortion>(new NullLensDistortion(*this));
}

void NullLensDistortion::write(std::ostream & os) const {
  os << "No distortion applied.\n";
}

std::string NullLensDistortion::name() const { return "NULL"; }

void NullLensDistortion::scale(float scale) { }


TsaiLensDistortion::TsaiLensDistortion(Vector4 const& params) : m_distortion(params) {}

Vector<double>
TsaiLensDistortion::distortion_parameters() const { return m_distortion; }

boost::shared_ptr<LensDistortion>
TsaiLensDistortion::copy() const {
  return boost::shared_ptr<TsaiLensDistortion>(new TsaiLensDistortion(*this));
}

Vector2
TsaiLensDistortion::distorted_coordinates(const camera::PinholeModel& cam, Vector2 const& p) const {

  Vector2 focal = cam.focal_length();
  Vector2 offset = cam.point_offset();

  if (focal[0] < 1e-300 || focal[1] < 1e-300)
    return Vector2(HUGE_VAL, HUGE_VAL);

  Vector2 p_0 = elem_quot(p - offset, focal); // represents x and y
  double r2 = norm_2_sqr( p_0 );
  Vector2 distortion( m_distortion[3], m_distortion[2] );
  Vector2 p_1 = elem_quot(distortion, p_0);
  Vector2 p_3 = 2*elem_prod(distortion,p_0);

  Vector2 b =  elem_prod(r2,p_1);
  b = elem_sum(b,r2*(m_distortion[0] + r2 * m_distortion[1]) + sum(p_3));

  // Prevent divide by zero at the origin or along the x and y center line
  Vector2 result = p + elem_prod(b,p - offset);
  if (p[0] == offset[0])
    result[0] = p[0];
  if (p[1] == offset[1])
    result[1] = p[1];

  return result;
}

void TsaiLensDistortion::write(std::ostream & os) const {
  os << "k1 = " << m_distortion[0] << "\n";
  os << "k2 = " << m_distortion[1] << "\n";
  os << "p1 = " << m_distortion[2] << "\n";
  os << "p2 = " << m_distortion[3] << "\n";
}

std::string TsaiLensDistortion::name() const { return "TSAI"; }

void TsaiLensDistortion::scale( float scale ) {
  m_distortion *= scale;
}

BrownConradyDistortion::BrownConradyDistortion( Vector<double> const& params ) {
  VW_ASSERT( params.size() == 8,
             ArgumentErr() << "BrownConradyDistortion: requires constructor input of size 8.");
  m_principal_point = subvector(params,0,2);
  m_radial_distortion = subvector(params,2,3);
  m_centering_distortion = subvector(params,5,2);
  m_centering_angle = params[7];
}

BrownConradyDistortion::BrownConradyDistortion( Vector<double> const& principal,
                                                Vector<double> const& radial,
                                                Vector<double> const& centering,
                                                double const& angle ) :
  m_principal_point(principal), m_radial_distortion(radial),
  m_centering_distortion(centering), m_centering_angle( angle ) {}

boost::shared_ptr<LensDistortion>
BrownConradyDistortion::copy() const {
  return boost::shared_ptr<BrownConradyDistortion>(new BrownConradyDistortion(*this));
}

Vector<double> BrownConradyDistortion::distortion_parameters() const {
  Vector<double,8> output;
  subvector(output,0,2) = m_principal_point;
  subvector(output,2,3) = m_radial_distortion;
  subvector(output,5,2) = m_centering_distortion;
  output[7] = m_centering_angle;
  return output;
}

Vector2
BrownConradyDistortion::undistorted_coordinates(const camera::PinholeModel& cam, Vector2 const& p) const {
  Vector2 offset = cam.point_offset();
  Vector2 intermediate = p - m_principal_point - offset;
  double r2 = norm_2_sqr(intermediate);
  double radial = 1 + m_radial_distortion[0]*r2 +
    m_radial_distortion[1]*r2*r2 + m_radial_distortion[2]*r2*r2*r2;
  double tangental = m_centering_distortion[0]*r2 + m_centering_distortion[1]*r2*r2;
  intermediate *= radial;
  intermediate[0] -= tangental*sin(m_centering_angle);
  intermediate[1] += tangental*cos(m_centering_angle);
  return intermediate+offset;
}

void BrownConradyDistortion::write(std::ostream& os) const {
  os << distortion_parameters() << "\n";
}

std::string BrownConradyDistortion::name() const { return "BROWNCONRADY"; }

void BrownConradyDistortion::scale( float scale ) {
  vw_throw( NoImplErr() << "BrownConradyDistortion doesn't support scaling" );
}


AdjustableTsaiLensDistortion::AdjustableTsaiLensDistortion(Vector<double> params) : m_distortion(params) {
  VW_ASSERT( params.size() > 3, ArgumentErr() << "Requires at least 4 coefficients for distortion. Last 3 are always the distortion coefficients and alpha. All leading elements are even radial distortion coefficients." );
}

Vector<double>
AdjustableTsaiLensDistortion::distortion_parameters() const {
  return m_distortion;
}

boost::shared_ptr<LensDistortion>
AdjustableTsaiLensDistortion::copy() const {
  return boost::shared_ptr<AdjustableTsaiLensDistortion>(new AdjustableTsaiLensDistortion(*this));
}


Vector2
AdjustableTsaiLensDistortion::distorted_coordinates(const camera::PinholeModel& cam, Vector2 const& p )  const {
  Vector2 focal = cam.focal_length();
  Vector2 offset = cam.point_offset();

  if (focal[0] < 1e-300 || focal[1] < 1e-300)
    return Vector2(HUGE_VAL, HUGE_VAL);

  // Create normalized coordinates
  Vector2 p_0 = elem_quot(p - offset, focal); // represents x and y
  double r2 = norm_2_sqr( p_0 );

  // Calculating Radial effects
  double r_n = 1, radial = 0;
  for ( unsigned i = 0; i < m_distortion.size()-3; i++ ) {
    r_n *= r2;
    radial += m_distortion[i]*r_n;
  }

  // Calculating Tangential effects
  Vector2 tangent;
  Vector2 swap_coeff(m_distortion[m_distortion.size()-2],
                     m_distortion[m_distortion.size()-3]);
  tangent += elem_prod(swap_coeff,elem_sum(r2,2*elem_prod(p_0,p_0)));
  tangent += 2*prod(p_0)*subvector(m_distortion,m_distortion.size()-3,2);

  // Final normalized result
  Vector2 result = p_0 + tangent + radial*p_0;

  // Running back through intrinsic matrix (with alpha or skew)
  return elem_prod(result+Vector2(m_distortion[m_distortion.size()-1]*result.y(),0),focal)+offset;
}

void AdjustableTsaiLensDistortion::write(std::ostream & os) const {
  os << "Radial Coeff: " << subvector(m_distortion,0,m_distortion.size()-3) << "\n";
  os << "Tangental Coeff: " << subvector(m_distortion,m_distortion.size()-3,2) << "\n";
  os << "Alpha: " << m_distortion[m_distortion.size()-1] << "\n";
}

std::string AdjustableTsaiLensDistortion::name() const { return "AdjustableTSAI"; }

void AdjustableTsaiLensDistortion::scale( float /*scale*/ ) {
  vw_throw( NoImplErr() << "AdjustableTsai doesn't support scaling." );
}


std::ostream& camera::operator<<(std::ostream & os,
                                 const camera::LensDistortion& ld) {
  ld.write(os);
  return os;
}
