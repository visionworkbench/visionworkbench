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

// Default implemenations for Lens Distortion -------------------


Vector<double>
LensDistortion::distortion_parameters() const { return Vector<double>(); }

void LensDistortion::set_distortion_parameters(Vector<double> const& params) {}

Vector2
LensDistortion::undistorted_coordinates(const camera::PinholeModel& cam, Vector2 const& v) const {
  UndistortOptimizeFunctor model(cam, *this);
  int status;
  Vector2 solution = math::levenberg_marquardt( model, v, v, status, 1e-6, 1e-6, 50 );
  //VW_DEBUG_ASSERT( status != math::optimization::eConvergedRelTolerance, 
  //                 PixelToRayErr() << "undistorted_coordinates: failed to converge." );
  return solution;
}

vw::Vector2
LensDistortion::distorted_coordinates(const camera::PinholeModel& cam, Vector2 const& v) const {
  DistortOptimizeFunctor model(cam, *this);
  int status;
  Vector2 solution = math::levenberg_marquardt( model, v, v, status, 1e-6, 1e-6, 50 );
  //VW_DEBUG_ASSERT( status != math::optimization::eConvergedRelTolerance, 
  //                 PixelToRayErr() << "distorted_coordinates: failed to converge." );
  return solution;
}


std::ostream& camera::operator<<(std::ostream & os,
                                 const camera::LensDistortion& ld) {
  ld.write(os);
  return os;
}

// Specific Implementations -------------------------------------

// ======== NullLensDistortion ========

boost::shared_ptr<LensDistortion> NullLensDistortion::copy() const {
  return boost::shared_ptr<NullLensDistortion>(new NullLensDistortion(*this));
}

void NullLensDistortion::write(std::ostream & os) const {
  os << "NULL\n";
}

void NullLensDistortion::read(std::istream & is) {
}

void NullLensDistortion::scale(float scale) { }

// ======== TsaiLensDistortion ========

TsaiLensDistortion::TsaiLensDistortion(Vector4 const& params) : m_distortion(params) {}

Vector<double>
TsaiLensDistortion::distortion_parameters() const { return m_distortion; }

void TsaiLensDistortion::set_distortion_parameters(Vector<double> const& params) {
  m_distortion = params;
}

boost::shared_ptr<LensDistortion>
TsaiLensDistortion::copy() const {
  return boost::shared_ptr<TsaiLensDistortion>(new TsaiLensDistortion(*this));
}

Vector2
TsaiLensDistortion::distorted_coordinates(const camera::PinholeModel& cam, Vector2 const& p) const {

  Vector2 focal  = cam.focal_length(); // = [fu, fv] 
  Vector2 offset = cam.point_offset(); // = [cu, cv]

  if (focal[0] < 1e-300 || focal[1] < 1e-300)
    return Vector2(HUGE_VAL, HUGE_VAL);

  Vector2 dudv(p - offset); // = [u-cx, v-cy]
  Vector2 p_0 = elem_quot(dudv, focal); // = dudv / f = [x, y] // Normalized pixel coordinates (1 == f)
  double r2 = norm_2_sqr( p_0 ); // = x^2 + y^2
  Vector2 distortion( m_distortion[3], m_distortion[2] ); // [p2, p1]
  Vector2 p_1 =   elem_quot(distortion, p_0); // = [  p2/x,   p1/y]
  Vector2 p_3 = 2*elem_prod(distortion, p_0); // = [2*p2*x, 2*p1*y]

  Vector2 b =  elem_prod(r2,p_1); // = [r2*p2/x, r2*p1/y]
  b = elem_sum(b,r2*(m_distortion[0] + r2 * m_distortion[1]) + sum(p_3));
  // = elem_sum(b, k1*r2 + k2*r4 + 2*p2*x + 2*p1*y);
  // = [ k1*r2 + k2*r4 + 2*p2*x + 2*p1*y + r2*p2/x, 
  //     k1*r2 + k2*r4 + 2*p2*x + 2*p1*y + r2*p1/y ]

  // Note that after the multiplication step below, this matches the commonly seen equations:
  // = [ x(k1*r2 + k2*r4) + 2*p2*x^2 + 2*p1*y*x + r2*p2, 
  //     y(k1*r2 + k2*r4) + 2*p2*x*y + 2*p1*y^2 + r2*p1 ]
  // = [ x(k1*r2 + k2*r4) + 2*p1*x*y + p2(r2 + 2x^2), 
  //     y(k1*r2 + k2*r4) + 2*p2*x*y + p1(r2 + 2y^2) ]

  // Prevent divide by zero at the origin or along the x and y center line
  Vector2 result = p + elem_prod(b, dudv); // = p + [du, dv]*(b)
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

void TsaiLensDistortion::read(std::istream & is) {
  std::string line;
  std::getline(is, line);
  if (!is.good() || sscanf(line.c_str(),"k1 = %lf", &m_distortion[0]) != 1)
    vw_throw( IOErr() << "TsaiLensDistortion::read(): Could not read k1\n" );
  std::getline(is, line);
  if (!is.good() || sscanf(line.c_str(),"k2 = %lf", &m_distortion[1]) != 1)
    vw_throw( IOErr() << "TsaiLensDistortion::read(): Could not read k2\n" );
  std::getline(is, line);
  if (!is.good() || sscanf(line.c_str(),"p1 = %lf", &m_distortion[2]) != 1)
    vw_throw( IOErr() << "TsaiLensDistortion::read(): Could not read p1\n" );
  std::getline(is, line);
  if (!is.good() || sscanf(line.c_str(),"p2 = %lf", &m_distortion[3]) != 1)
    vw_throw( IOErr() << "TsaiLensDistortion::read(): Could not read p2\n" );
}

void TsaiLensDistortion::scale( float scale ) {
  m_distortion *= scale;
}

// ======== BrownConradyDistortion ========

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

void BrownConradyDistortion::set_distortion_parameters(Vector<double> const& params) {
  m_principal_point      = subvector(params,0,2);
  m_radial_distortion    = subvector(params,2,3);
  m_centering_distortion = subvector(params,5,2);
  m_centering_angle      = params[7];
}

Vector2
BrownConradyDistortion::undistorted_coordinates(const camera::PinholeModel& cam, Vector2 const& p) const {
  Vector2 offset       = cam.point_offset();
  Vector2 intermediate = p - m_principal_point - offset;
  double r2        = norm_2_sqr(intermediate);
  double radial    = 1 + m_radial_distortion   [0]*r2 + m_radial_distortion   [1]*r2*r2
                                                      + m_radial_distortion   [2]*r2*r2*r2;
  double tangental =     m_centering_distortion[0]*r2 + m_centering_distortion[1]*r2*r2;
  intermediate *= radial;
  intermediate[0] -= tangental*sin(m_centering_angle);
  intermediate[1] += tangental*cos(m_centering_angle);
  return intermediate+offset;
}

void BrownConradyDistortion::write(std::ostream& os) const {
  os << "xp  = " << m_principal_point[0]      << "\n";
  os << "yp  = " << m_principal_point[1]      << "\n";
  os << "k1  = " << m_radial_distortion[0]    << "\n";
  os << "k2  = " << m_radial_distortion[1]    << "\n";
  os << "k3  = " << m_radial_distortion[2]    << "\n";
  os << "p1  = " << m_centering_distortion[0] << "\n";
  os << "p2  = " << m_centering_distortion[1] << "\n";
  os << "phi = " << m_centering_angle         << "\n";
}

void BrownConradyDistortion::read(std::istream & is) {
  std::string line;
  std::getline(is, line);
  if (!is.good() || sscanf(line.c_str(),"xp = %lf", &m_principal_point[0]) != 1)
    vw_throw( IOErr() << "BrownConradyDistortion::read(): Could not read xp\n" );
  std::getline(is, line);
  if (!is.good() || sscanf(line.c_str(),"yp = %lf", &m_principal_point[1]) != 1)
    vw_throw( IOErr() << "BrownConradyDistortion::read(): Could not read yp\n" );
  std::getline(is, line);
  if (!is.good() || sscanf(line.c_str(),"k1 = %lf", &m_radial_distortion[0]) != 1)
    vw_throw( IOErr() << "BrownConradyDistortion::read(): Could not read k1\n" );
  std::getline(is, line);
  if (!is.good() || sscanf(line.c_str(),"k2 = %lf", &m_radial_distortion[1]) != 1)
    vw_throw( IOErr() << "BrownConradyDistortion::read(): Could not read k2\n" );
  std::getline(is, line);
  if (!is.good() || sscanf(line.c_str(),"k3 = %lf", &m_radial_distortion[2]) != 1)
    vw_throw( IOErr() << "BrownConradyDistortion::read(): Could not read k3\n" );
  std::getline(is, line);
  if (!is.good() || sscanf(line.c_str(),"p1 = %lf", &m_centering_distortion[0]) != 1)
    vw_throw( IOErr() << "BrownConradyDistortion::read(): Could not read p1\n" );
  std::getline(is, line);
  if (!is.good() || sscanf(line.c_str(),"p2 = %lf", &m_centering_distortion[1]) != 1)
    vw_throw( IOErr() << "BrownConradyDistortion::read(): Could not read p2\n" );
  std::getline(is, line);
  if (!is.good() || sscanf(line.c_str(),"phi = %lf", &m_centering_angle) != 1)
    vw_throw( IOErr() << "BrownConradyDistortion::read(): Could not read phi\n" );
}


void BrownConradyDistortion::scale( float scale ) {
  vw_throw( NoImplErr() << "BrownConradyDistortion doesn't support scaling" );
}

// ======== AdjustableTsaiLensDistortion ========

AdjustableTsaiLensDistortion::AdjustableTsaiLensDistortion(Vector<double> params) : m_distortion(params) {
  VW_ASSERT( params.size() > 3, ArgumentErr() << "Requires at least 4 coefficients for distortion. Last 3 are always the distortion coefficients and alpha. All leading elements are even radial distortion coefficients." );
}

Vector<double>
AdjustableTsaiLensDistortion::distortion_parameters() const {
  return m_distortion;
}

void AdjustableTsaiLensDistortion::set_distortion_parameters(Vector<double> const& params) {
  m_distortion = params;
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
  // Since we don't know how many parameters, print them out in three rows.
  os << "Radial Coeff: "    << subvector(m_distortion,0,m_distortion.size()-3) << "\n";
  os << "Tangental Coeff: " << subvector(m_distortion,m_distortion.size()-3,2) << "\n";
  os << "Alpha: "           << m_distortion[m_distortion.size()-1] << "\n";
}

void AdjustableTsaiLensDistortion::read(std::istream & is) {
  Vector<float64,0> radial_vec, tangential_vec;
  double alpha;
  
  std::string label, line;
  try {
    is >> label >> label >> radial_vec;
    is >> label >> label >> tangential_vec;
  }
  catch(...) {
    vw_throw( IOErr() << "AdjustableTsaiLensDistortion::read(): Could not read vector params!\n" );
  }
  while (is.peek() == '\n') // Get to the start of the next line
     is.ignore(1,' ');
  std::getline(is, line);
  if (!is.good() || (sscanf(line.c_str(),"Alpha: %lf", &alpha) != 1))
    vw_throw( IOErr() << "AdjustableTsaiLensDistortion::read(): Could not read alpha\n" );

  // Pack in the values we read
  size_t total_size = radial_vec.size() + tangential_vec.size() + 1;
  m_distortion.set_size(total_size);
  subvector(m_distortion,0,total_size-3) = radial_vec;
  subvector(m_distortion,total_size-3,2) = tangential_vec;
  m_distortion[m_distortion.size()-1]    = alpha;
}

void AdjustableTsaiLensDistortion::scale( float /*scale*/ ) {
  vw_throw( NoImplErr() << "AdjustableTsai doesn't support scaling." );
}


// ======== PhotometrixLensDistortion ========

PhotometrixLensDistortion::PhotometrixLensDistortion(Vector<float64,9> const& params) 
  : m_distortion(params) {
}

Vector<double>
PhotometrixLensDistortion::distortion_parameters() const { 
  return m_distortion; 
}

void PhotometrixLensDistortion::set_distortion_parameters(Vector<double> const& params) {
  m_distortion = params;
}

boost::shared_ptr<LensDistortion>
PhotometrixLensDistortion::copy() const {
  return boost::shared_ptr<PhotometrixLensDistortion>(new PhotometrixLensDistortion(*this));
}

Vector2
PhotometrixLensDistortion::undistorted_coordinates(const camera::PinholeModel& cam, Vector2 const& p) const {

  double x_meas = p[0];
  double y_meas = p[1];
  
  // PinholeModel contains the nominal center and we add the offsets from calibration here
  // - These two values need to be specified in the same units!
  Vector2 center_point = cam.point_offset(); // = [cu, cv]
  double xp  = m_distortion[0] + center_point[0];
  double yp  = m_distortion[1] + center_point[1];

  double x   = x_meas - xp;
  double y   = y_meas - yp;

  double x2  = x*x;
  double y2  = y*y;
  double r2  = x2 + y2;
  
  double K1 = m_distortion[2];
  double K2 = m_distortion[3];
  double K3 = m_distortion[4];
  
  double drr = K1*r2 + K2*r2*r2 + K3*r2*r2*r2; // This is dr/r, not dr

  //std::cout << "x = " << x << std::endl;
  //std::cout << "y = " << y << std::endl;
  //std::cout << "drr = " << drr << std::endl;
  
  double P1 = m_distortion[5];
  double P2 = m_distortion[6];
  
  //std::cout << "dx_t = " << P1*(r2 +2.0*x2) + 2.0*P2*x*y << std::endl;
  //std::cout << "dy_t = " << P2*(r2 +2.0*y2) + 2.0*P1*x*y << std::endl;
  
  // The cal document includes -xp/yp in these lines, but they are removed
  //  so that the results are relative to 0,0 instead of relative to the principle point.
  double x_corr = x_meas + x*drr + P1*(r2 +2.0*x2) + 2.0*P2*x*y;
  double y_corr = y_meas + y*drr + P2*(r2 +2.0*y2) + 2.0*P1*x*y;

  // Note that parameters B1 and B2 are not used.  The software output provides them
  // but did not specify their use since they were zero.  If you see an example that 
  // includes them, update the calculations above!

  return Vector2(x_corr, y_corr);
}

void PhotometrixLensDistortion::write(std::ostream & os) const {
  os << "xp = " << m_distortion[0] << "\n";
  os << "yp = " << m_distortion[1] << "\n";
  os << "k1 = " << m_distortion[2] << "\n";
  os << "k2 = " << m_distortion[3] << "\n";
  os << "k3 = " << m_distortion[4] << "\n";
  os << "p1 = " << m_distortion[5] << "\n";
  os << "p2 = " << m_distortion[6] << "\n";
  os << "b1 = " << m_distortion[7] << "\n";
  os << "b2 = " << m_distortion[8] << "\n";
}

void PhotometrixLensDistortion::read(std::istream & is) {
  std::string line;
  std::getline(is, line);
  if (!is.good() || sscanf(line.c_str(),"xp = %lf", &m_distortion[0]) != 1)
    vw_throw( IOErr() << "PhotometrixLensDistortion::read(): Could not read xp\n" );
  std::getline(is, line);
  if (!is.good() || sscanf(line.c_str(),"yp = %lf", &m_distortion[1]) != 1)
    vw_throw( IOErr() << "PhotometrixLensDistortion::read(): Could not read yp\n" );
  std::getline(is, line);
  if (!is.good() || sscanf(line.c_str(),"k1 = %lf", &m_distortion[2]) != 1)
    vw_throw( IOErr() << "PhotometrixLensDistortion::read(): Could not read k1\n" );
  std::getline(is, line);
  if (!is.good() || sscanf(line.c_str(),"k2 = %lf", &m_distortion[3]) != 1)
    vw_throw( IOErr() << "PhotometrixLensDistortion::read(): Could not read k2\n" );
  std::getline(is, line);
  if (!is.good() || sscanf(line.c_str(),"k3 = %lf", &m_distortion[4]) != 1)
    vw_throw( IOErr() << "PhotometrixLensDistortion::read(): Could not read k3\n" );
  std::getline(is, line);
  if (!is.good() || sscanf(line.c_str(),"p1 = %lf", &m_distortion[5]) != 1)
    vw_throw( IOErr() << "PhotometrixLensDistortion::read(): Could not read p1\n" );
  std::getline(is, line);
  if (!is.good() || sscanf(line.c_str(),"p2 = %lf", &m_distortion[6]) != 1)
    vw_throw( IOErr() << "PhotometrixLensDistortion::read(): Could not read p2\n" );
  std::getline(is, line);
  if (!is.good() || sscanf(line.c_str(),"b1 = %lf", &m_distortion[7]) != 1)
    vw_throw( IOErr() << "PhotometrixLensDistortion::read(): Could not read b1\n" );
  std::getline(is, line);
  if (!is.good() || sscanf(line.c_str(),"b2 = %lf", &m_distortion[8]) != 1)
    vw_throw( IOErr() << "PhotometrixLensDistortion::read(): Could not read b2\n" );
    
}

void PhotometrixLensDistortion::scale( float scale ) {
  m_distortion *= scale;
}





