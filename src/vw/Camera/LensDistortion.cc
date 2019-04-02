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

// Isolate these local utilities to this file.
namespace {
  
// Pull all lines from the stream. Search for "name = val". Store the
// values in the order given in "names". Complain if some fields were
// not populated. This has the advantage that the order of lines in
// the stream is not important, and it won't complain if there are
// extraneous fields.
template<class VectorT>
void read_fields_in_vec(std::vector<std::string> const& names, VectorT & vals, std::istream & is){

  std::map<std::string, double> name2val;
  std::string line;
  char name[1024];
  double val;
  while (is.good()){
    std::getline(is, line);
    if ( sscanf(line.c_str(),"%s = %lf", name, &val) != 2 )
      continue;
    
    name2val[name] = val;
  }

  // Populate the output
  for (size_t i = 0; i < names.size(); i++) {
    if (i >= vals.size())
      vw_throw( IOErr() << "Not enough room allocated for output.\n" );

    std::map<std::string, double>::iterator it = name2val.find(names[i]);
    if (it == name2val.end()) 
      vw_throw( IOErr() << "Could not read a value for " << names[i] << ".\n" );
    
    vals[i] = it->second;
  }
}

// Read a line of the form: name = a b c. The number of values is not known.
template<class VectorT>
void read_param_vec(std::string const& param_name, std::istream & is, VectorT & vals){

  if (!is.good())
    vw_throw( IOErr() << "Could not read LensDistortion.\n" );

  std::string line;
  std::getline(is, line);
  std::istringstream iss(line);
  
  std::string val;
  if ( !(iss >> val) || val != param_name )
    vw_throw( IOErr() << "In LensDistortion, got " << val << " but expected " << param_name );
  
  if ( !(iss >> val)  || val != "=" )
    vw_throw( IOErr() << "In LensDistortion, got " << val << " but expected the equal sign.");

  // First read in a vector of doubles
  double dval = 0.0;
  std::vector<double> local_vals;
  while (iss >> dval) 
    local_vals.push_back(dval);
  
  // Copy to the vector vals.
  int num_vals = local_vals.size();
  vals.set_size(num_vals);
  for (int it = 0; it < num_vals; it++) 
    vals[it] = local_vals[it];
}

// Read a line of the form: name = a b c. The number of parameters must be param_len.
template<class VectorT>
void read_param_vec(std::string const& param_name, size_t param_len,
                    std::istream & is, VectorT & vals){

  read_param_vec(param_name, is, vals);
  if (vals.size() != param_len) {
    vw_throw( IOErr() << "Book-keeping failure in reading LensDistortion. Expecting "
              << vals.size() << " == " << param_len << "\n" );
  }
}

// Read a line of the form: name = a.
template<class T>
void read_param(std::string const& param_name, std::istream & is, T & val){
  Vector<T> vals;
  read_param_vec(param_name, 1, is, vals);
  val = vals[0];
}
  
// Write a line of the form: name = a b c
void write_param_vec(std::string const& param_name, std::ostream & os, Vector<double> const& vals){
  os << param_name << " = ";
  for (size_t p = 0; p < vals.size(); p++){
    os << vals[p];
    if (p + 1 < vals.size()) 
      os << " "; // write a whitespace after each number except the last
  }
  os << "\n";
}

// Write a line of the form: name = val
void write_param(std::string const& param_name, std::ostream & os, double val){
  Vector<double> vals;
  vals.set_size(1);
  vals[0] = val;
  write_param_vec(param_name, os, vals);
}

} // end of anonymous namespace

// Default implemenations for Lens Distortion -------------------


Vector<double>
LensDistortion::distortion_parameters() const { return Vector<double>(); }

void LensDistortion::set_distortion_parameters(Vector<double> const& params) {}

Vector2
LensDistortion::undistorted_coordinates(const camera::PinholeModel& cam, Vector2 const& v) const {
  UndistortOptimizeFunctor model(cam, *this);
  int status;

  // Must push the solver really hard, to make sure bundle adjust gets accurate values
  // to play with.
  Vector2 solution = math::levenberg_marquardt( model, v, v, status, 1e-16, 1e-16, 100 );

  Vector2 dist = this->distorted_coordinates(cam, solution);
  double err = norm_2(dist - v)/std::max(norm_2(v), 0.1); // don't make this way too strict
  double tol = 1e-10;
  if (err > tol)
    vw_throw( PointToPixelErr() << "LensDistortion: Did not converge.\n" );
  
  //VW_ASSERT((status == math::optimization::eConvergedAbsTolerance), 
  //                 PixelToRayErr() << "distorted_coordinates: failed to converge." );
  //double error = norm_2(model(solution) - v);
  //std::cout << "status = " << status << ", input = " << v 
  //          << ", pixel = " << solution << ", error = " << error << std::endl;
  return solution;
}

vw::Vector2
LensDistortion::distorted_coordinates(const camera::PinholeModel& cam, Vector2 const& v) const {
  DistortOptimizeFunctor model(cam, *this);
  int status;
  // Must push the solver really hard, to make sure bundle adjust gets accurate values
  // to play with.
  Vector2 solution = math::levenberg_marquardt(model, v, v, status, 1e-16, 1e-16, 100);
  //VW_DEBUG_ASSERT( status != math::optimization::eConvergedRelTolerance, 
  //                 PixelToRayErr() << "distorted_coordinates: failed to converge." );

  // Check if it failed badly to converge. That it did not converge is not on its own
  // unreasonable, sometimes the inputs are bad. But then the user must know about it.
  Vector2 undist = this->undistorted_coordinates(cam, solution);
  double err = norm_2(undist - v)/std::max(norm_2(v), 0.1); // don't make this way too strict
  double tol = 1e-10;
  if (err > tol)
    vw_throw( PointToPixelErr() << "LensDistortion: Did not converge.\n" );

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
  // Nothing to write. The "NULL" string is already
  // written before this, when we saved the distortion type.
  // os << "NULL\n";
}

void NullLensDistortion::read(std::istream & is) {
  // Nothing to read
}

void NullLensDistortion::scale(double scale) { }

// ======== TsaiLensDistortion ========

TsaiLensDistortion::TsaiLensDistortion(){
  TsaiLensDistortion::init_distortion_param_names();
  m_distortion.set_size(m_distortion_param_names.size());
}

TsaiLensDistortion::TsaiLensDistortion(Vector<double> const& params) : m_distortion(params) {
  TsaiLensDistortion::init_distortion_param_names();
  if (m_distortion.size() != m_distortion_param_names.size())
    vw_throw( IOErr() << "TsaiLensDistortion: Incorrect number of parameters was passed in.");
}

Vector<double>
TsaiLensDistortion::distortion_parameters() const { return m_distortion; }

void TsaiLensDistortion::init_distortion_param_names(){
  std::string names[] = {"k1", "k2", "p1", "p2"};
  size_t num_names = sizeof(names)/sizeof(std::string);
  m_distortion_param_names.resize(num_names);
  for (size_t p = 0; p < num_names; p++) 
    m_distortion_param_names[p] = names[p];
}

void TsaiLensDistortion::set_distortion_parameters(Vector<double> const& params) {
  TsaiLensDistortion::init_distortion_param_names();
  m_distortion = params;
  if (m_distortion.size() != m_distortion_param_names.size())
    vw_throw( IOErr() << "TsaiLensDistortion: Incorrect number of parameters was passed in.");
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
  for (size_t p = 0; p < m_distortion_param_names.size(); p++) 
    os << m_distortion_param_names[p] << " = " << m_distortion[p] << "\n";
}

void TsaiLensDistortion::read(std::istream & is) {
  m_distortion.set_size(m_distortion_param_names.size());
  read_fields_in_vec(m_distortion_param_names, m_distortion, is);
}

void TsaiLensDistortion::scale( double scale ) {
  m_distortion *= scale;
}

// ======== BrownConradyDistortion ========

BrownConradyDistortion::BrownConradyDistortion() {
  BrownConradyDistortion::init_distortion_param_names();
}

void BrownConradyDistortion::init_distortion_param_names(){
  std::string names[] = {"xp", "yp", "k1", "k2", "k3", "p1", "p2", "phi"};
  size_t num_names = sizeof(names)/sizeof(std::string);
  m_distortion_param_names.resize(num_names);
  for (size_t p = 0; p < num_names; p++)
    m_distortion_param_names[p] = names[p];
}

BrownConradyDistortion::BrownConradyDistortion( Vector<double> const& params ) {
  VW_ASSERT( params.size() == 8,
             ArgumentErr() << "BrownConradyDistortion: requires constructor input of size 8.");
  m_principal_point      = subvector(params,0,2);
  m_radial_distortion    = subvector(params,2,3);
  m_centering_distortion = subvector(params,5,2);
  m_centering_angle      = params[7];
  BrownConradyDistortion::init_distortion_param_names();
}

BrownConradyDistortion::BrownConradyDistortion( Vector<double> const& principal,
                                                Vector<double> const& radial,
                                                Vector<double> const& centering,
                                                double const& angle ) :
  m_principal_point(principal), m_radial_distortion(radial),
  m_centering_distortion(centering), m_centering_angle( angle ) {
  BrownConradyDistortion::init_distortion_param_names();
}

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
  double r2            = norm_2_sqr(intermediate);
  double radial        = 1 + m_radial_distortion   [0]*r2 + m_radial_distortion   [1]*r2*r2
                                                      + m_radial_distortion   [2]*r2*r2*r2;
  double tangental     = m_centering_distortion[0]*r2 + m_centering_distortion[1]*r2*r2;
  intermediate *= radial;
  intermediate[0] -= tangental*sin(m_centering_angle);
  intermediate[1] += tangental*cos(m_centering_angle);
  return intermediate+offset;
}

void BrownConradyDistortion::write(std::ostream& os) const {
  int p = 0;
  os << m_distortion_param_names[p] << "  = " << m_principal_point[0]      << "\n"; p++;
  os << m_distortion_param_names[p] << "  = " << m_principal_point[1]      << "\n"; p++;
  os << m_distortion_param_names[p] << "  = " << m_radial_distortion[0]    << "\n"; p++;
  os << m_distortion_param_names[p] << "  = " << m_radial_distortion[1]    << "\n"; p++;
  os << m_distortion_param_names[p] << "  = " << m_radial_distortion[2]    << "\n"; p++;
  os << m_distortion_param_names[p] << "  = " << m_centering_distortion[0] << "\n"; p++;
  os << m_distortion_param_names[p] << "  = " << m_centering_distortion[1] << "\n"; p++;
  os << m_distortion_param_names[p] << " = "  << m_centering_angle         << "\n"; p++; // phi
}

void BrownConradyDistortion::read(std::istream & is) {
  Vector<double, num_distortion_params> distortion;
  read_fields_in_vec(m_distortion_param_names, distortion, is);
  int p = 0;
  m_principal_point[0]      = distortion[p]; p++;
  m_principal_point[1]      = distortion[p]; p++;
  m_radial_distortion[0]    = distortion[p]; p++;
  m_radial_distortion[1]    = distortion[p]; p++;
  m_radial_distortion[2]    = distortion[p]; p++;
  m_centering_distortion[0] = distortion[p]; p++;
  m_centering_distortion[1] = distortion[p]; p++;
  m_centering_angle         = distortion[p]; p++;
}


void BrownConradyDistortion::scale( double scale ) {
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
  Vector<double,0> radial_vec, tangential_vec;
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

void AdjustableTsaiLensDistortion::scale( double /*scale*/ ) {
  vw_throw( NoImplErr() << "AdjustableTsai doesn't support scaling." );
}


// ======== PhotometrixLensDistortion ========

PhotometrixLensDistortion::PhotometrixLensDistortion(){
  PhotometrixLensDistortion::init_distortion_param_names();
  m_distortion.set_size(m_distortion_param_names.size());
}

PhotometrixLensDistortion::PhotometrixLensDistortion(Vector<double> const& params) 
  : m_distortion(params) {
  PhotometrixLensDistortion::init_distortion_param_names();
  if (m_distortion.size() != m_distortion_param_names.size())
    vw_throw( IOErr() << "PhotometrixLensDistortion: Incorrect number of parameters was passed in.");
}

void PhotometrixLensDistortion::init_distortion_param_names(){
  std::string names[] = {"xp", "yp",
                         "k1", "k2", "k3", "p1", "p2",
                         "b1", "b2"};
  
  size_t num_names = sizeof(names)/sizeof(std::string);
  m_distortion_param_names.resize(num_names);
  for (size_t p = 0; p < num_names; p++) 
    m_distortion_param_names[p] = names[p];
}

Vector<double>
PhotometrixLensDistortion::distortion_parameters() const { 
  return m_distortion; 
}

void PhotometrixLensDistortion::set_distortion_parameters(Vector<double> const& params) {
  m_distortion = params;
  if (m_distortion.size() != m_distortion_param_names.size())
    vw_throw( IOErr() << "PhotometrixLensDistortion: Incorrect number of parameters was passed in.");
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
  double P1 = m_distortion[5];
  double P2 = m_distortion[6];

  double drr = K1*r2 + K2*r2*r2 + K3*r2*r2*r2; // This is dr/r, not dr

  // The cal document includes -xp/yp in these lines, but they are removed
  //  so that the results are relative to 0,0 instead of relative to the principal point.
  double x_corr = x_meas + x*drr + P1*(r2 + 2.0*x2) + 2.0*P2*x*y;
  double y_corr = y_meas + y*drr + P2*(r2 + 2.0*y2) + 2.0*P1*x*y;

  // Note that parameters B1 and B2 are not used. The software output provides them
  // but did not specify their use since they were zero.  If you see an example that 
  // includes them, update the calculations above!

  return Vector2(x_corr, y_corr);
}

void PhotometrixLensDistortion::write(std::ostream & os) const {
  for (size_t p = 0; p < m_distortion_param_names.size(); p++) 
    os << m_distortion_param_names[p] << " = " << m_distortion[p] << "\n";
}

void PhotometrixLensDistortion::read(std::istream & is) {
  m_distortion.set_size(m_distortion_param_names.size());
  read_fields_in_vec(m_distortion_param_names, m_distortion, is);
}

void PhotometrixLensDistortion::scale( double scale ) {
  m_distortion *= scale;
}


// ======== RPCLensDistortion ========
// This class is not fully formed until both distortion and
// undistortion parameters are computed.
// One must always call set_undistortion_parameters()
// only after set_distortion_parameters().
RPCLensDistortion::RPCLensDistortion(){
  m_rpc_degree = 0;
  m_can_undistort = false;
}

RPCLensDistortion::RPCLensDistortion(Vector<double> const& params): m_distortion(params) {
  validate_distortion_params(params);
  m_rpc_degree = rpc_degree(params.size());
  m_can_undistort = false;
}

Vector<double>
RPCLensDistortion::distortion_parameters() const { 
  return m_distortion; 
}

Vector<double>
RPCLensDistortion::undistortion_parameters() const { 
  return m_undistortion; 
}

void RPCLensDistortion::set_distortion_parameters(Vector<double> const& params) {
  validate_distortion_params(params);

  // If the distortion parameters changed, one cannot undistort until the undistortion
  // coefficients are computed.
  if (params.size() != m_distortion.size()) {
    m_can_undistort = false;
  } else{
    for (size_t it = 0; it < params.size(); it++) {
      if (m_distortion[it] != params[it]) {
        m_can_undistort = false;
        break;
      }
    }
  }
  m_distortion = params;
  m_rpc_degree = rpc_degree(params.size());
}

void RPCLensDistortion::set_undistortion_parameters(Vector<double> const& params) {
  if (params.size() != num_dist_params()) 
    vw_throw( IOErr() << class_name()
              << ": The number of distortion and undistortion parameters must agree.");
  m_undistortion = params;
  m_can_undistort = true;
}

void RPCLensDistortion::set_image_size(Vector2i const& image_size){
  m_image_size = image_size;
}

// Form the identity transform
void RPCLensDistortion::reset(int rpc_degree){
  if (rpc_degree <= 0) 
    vw_throw( IOErr() << class_name() << ": The RPC degree must be positive.");

  m_rpc_degree = rpc_degree;
  int num_params = num_dist_params(rpc_degree);
  m_distortion.set_size(num_params);
  m_undistortion.set_size(num_params);
  init_as_identity(m_distortion);
  init_as_identity(m_undistortion);
  m_can_undistort = true;
}

boost::shared_ptr<LensDistortion>
RPCLensDistortion::copy() const {
  return boost::shared_ptr<RPCLensDistortion>(new RPCLensDistortion(*this));
}

Vector2
RPCLensDistortion::distorted_coordinates(const camera::PinholeModel& cam, Vector2 const& p) const {
  // Here, unlike in the older RPCLensDistortion4 and RPCLensDistortion5, we
  // put the origin at the optical center.
  Vector2 offset = cam.point_offset(); // = [cu, cv]
  return RPCLensDistortion::compute_rpc(p - offset, m_distortion) + offset;
}

Vector2
RPCLensDistortion::undistorted_coordinates(const camera::PinholeModel& cam, Vector2 const& p) const {
  if (!m_can_undistort) 
    vw_throw( IOErr() << class_name() << ": Undistorted coefficients are not up to date.\n" );
  Vector2 offset = cam.point_offset(); // = [cu, cv]
  return RPCLensDistortion::compute_rpc(p - offset, m_undistortion) + offset;
}

// Evaluate the RPC with given coefficients.
Vector2
RPCLensDistortion::compute_rpc(Vector2 const& p, Vector<double> const& coeffs) const {
  
  if (num_dist_params() != (int)coeffs.size()) 
    vw_throw( IOErr() << "Book-keeping failure in RPCLensDistortion.\n" );

  int rpc_deg = rpc_degree(coeffs.size());
  double x = p[0];
  double y = p[1];
  
  // Precompute x^n and y^m values
  std::vector<double> powx(rpc_deg + 1), powy(rpc_deg + 1);
  double valx = 1.0, valy = 1.0;
  for (int deg = 0; deg <= rpc_deg; deg++) {
    powx[deg] = valx; valx *= x; 
    powy[deg] = valy; valy *= y; 
  }

  // Evaluate the RPC expression. The denominator always
  // has a 1 as the 0th-coefficient.
  int rpc_count = 0;
  double vals[] = {0.0, 1.0, 0.0, 1.0};
  for (int it = 0; it < 4; it++) {
    int start = 0; // starting degree for numerator
    if (it == 1 || it == 3)
      start = 1; // starting degree for denominator

    for (int deg = start; deg <= rpc_deg; deg++) { // start at 0
      for (int degy = 0; degy <= deg; degy++) {
        int degx = deg - degy;
        vals[it] += coeffs[rpc_count] * powx[degx] * powy[degy];
        rpc_count++;
      }
    }
  }

  if (rpc_count != (int)coeffs.size()) 
    vw_throw( IOErr() << "Book-keeping failure in RPCLensDistortion.\n" );

  return Vector2(vals[0]/vals[1], vals[2]/vals[3]);
}

void RPCLensDistortion::scale(double scale) {
  m_distortion *= scale;
  m_undistortion *= scale;
}

void RPCLensDistortion::validate_distortion_params(Vector<double> const& params) {
  int num_params = params.size();
  int deg = rpc_degree(num_params);
  if (num_dist_params(deg) != num_params || deg <= 0 || deg != deg)
    vw_throw( IOErr() << class_name() << ": Incorrect number of parameters was passed in: "
              << num_params << ".");
}

// Make RPC coefficients so that the RPC transform is the identity.
// The vector params must already have the right size.
void RPCLensDistortion::init_as_identity(Vector<double> & params){
  RPCLensDistortion::validate_distortion_params(params);

  params.set_all(0);
  Vector<double> num_x, den_x, num_y, den_y;
  unpack_params(params, num_x, den_x,  num_y, den_y);

  // Initialize the transform (x, y) -> (x, y), which is 
  // ( (0 + 1*x + 0*y)/(1 + 0*x + 0*y), (0 + 0*x + 1*y)/(1 + 0*x + 0*y) )
  // hence set num_x and num_y accordingly. As always, we do not
  // store the 1 values in the denominator.
  num_x[1] = 1; num_y[2] = 1;
  pack_params(params, num_x, den_x, num_y, den_y);
}

void RPCLensDistortion::unpack_params(Vector<double> const& params,
                                       Vector<double> & num_x, Vector<double> & den_x,
                                       Vector<double> & num_y, Vector<double> & den_y){
  RPCLensDistortion::validate_distortion_params(params);
  
  int num_params = params.size();
  int num_len = (num_params + 2)/4;
  int den_len = num_len - 1; // because the denominator always starts with 1.
  num_x = subvector(params, 0,                           num_len);
  den_x = subvector(params, num_len,                     den_len);
  num_y = subvector(params, num_len + den_len,           num_len);
  den_y = subvector(params, num_len + den_len + num_len, den_len);
}

void RPCLensDistortion::pack_params(Vector<double> & params,
                                     Vector<double> const& num_x, Vector<double> const& den_x,
                                     Vector<double> const& num_y, Vector<double> const& den_y){

  int num_len = num_x.size();
  int den_len = den_x.size();

  if (num_len != den_len + 1 || num_len != int(num_y.size()) || den_len != int(den_y.size()))
    vw_throw( IOErr() << "Book-keeping failure in RPCLensDistortion.\n" );

  params.set_size(2*num_len + 2*den_len);
  
  subvector(params, 0,                           num_len) = num_x; 
  subvector(params, num_len,                     den_len) = den_x;
  subvector(params, num_len + den_len,           num_len) = num_y;
  subvector(params, num_len + den_len + num_len, den_len) = den_y;
  RPCLensDistortion::validate_distortion_params(params);
}

namespace {
  // Prepend a 1 to a vector
  void prepend_1(Vector<double> & params) {
    int len = params.size();
    Vector<double> tmp_params = params;
    params.set_size(len + 1);
    params[0] = 1;
    subvector(params, 1, len) = tmp_params;
  }
  // Remove 1 from first position in a vector
  void remove_1(Vector<double> & params) {
    int len = params.size();
    if (len <= 0)
      vw_throw( IOErr() << "Found an unexpected empty vector.\n" );
    Vector<double> tmp_params = params;
    params = subvector(tmp_params, 1, len - 1);
  }
}
  
void RPCLensDistortion::write(std::ostream & os) const {

  if (!m_can_undistort) 
    vw_throw( IOErr() << class_name() << ": Undistorted coefficients are not up to date.\n" );

  // TODO: Add domain of validity for the distored and undistorted pixels. This is needed
  // because otherwise the RPC model can return wrong results which confuses
  // bundle adjustment. 
  write_param("rpc_degree", os, m_rpc_degree);
  write_param_vec("image_size", os, m_image_size);

  // Write the distortion
  Vector<double> num_x, den_x, num_y, den_y;
  unpack_params(m_distortion, num_x, den_x, num_y, den_y);
  prepend_1(den_x);
  prepend_1(den_y);
  // Leave spaces below so that it will display nicely
  write_param_vec("distortion_num_x  ", os, num_x);
  write_param_vec("distortion_den_x  ", os, den_x);
  write_param_vec("distortion_num_y  ", os, num_y);
  write_param_vec("distortion_den_y  ", os, den_y);
  
  // Write the undistortion
  unpack_params(m_undistortion, num_x, den_x, num_y, den_y);
  prepend_1(den_x);
  prepend_1(den_y);
  write_param_vec("undistortion_num_x", os, num_x);
  write_param_vec("undistortion_den_x", os, den_x);
  write_param_vec("undistortion_num_y", os, num_y);
  write_param_vec("undistortion_den_y", os, den_y);
}

void RPCLensDistortion::read(std::istream & is) {

  read_param("rpc_degree", is, m_rpc_degree);
  read_param_vec("image_size", 2, is, m_image_size);

  int num_params = num_dist_params(m_rpc_degree);
  int quarter = (num_params+2)/4; // to account for the two 1's in the denominator

  // Read the distortion
  Vector<double> num_x, den_x, num_y, den_y;
  read_param_vec("distortion_num_x", quarter, is, num_x);
  read_param_vec("distortion_den_x", quarter, is, den_x);
  read_param_vec("distortion_num_y", quarter, is, num_y);
  read_param_vec("distortion_den_y", quarter, is, den_y);
  remove_1(den_x);
  remove_1(den_y);
  pack_params(m_distortion, num_x, den_x, num_y, den_y);

  // Read the undistortion
  read_param_vec("undistortion_num_x", quarter, is, num_x);
  read_param_vec("undistortion_den_x", quarter, is, den_x);
  read_param_vec("undistortion_num_y", quarter, is, num_y);
  read_param_vec("undistortion_den_y", quarter, is, den_y);
  remove_1(den_x);
  remove_1(den_y);
  pack_params(m_undistortion, num_x, den_x, num_y, den_y);

  m_can_undistort = true; 
}

// Old RPC lens distortion of degree 4.

RPCLensDistortion4::RPCLensDistortion4(){
  m_distortion.set_size(num_distortion_params);
  m_undistortion.set_size(num_distortion_params);
  m_can_undistort = false;
}

RPCLensDistortion4::RPCLensDistortion4(Vector<double> const& params) 
  : m_distortion(params) {
  if (m_distortion.size() != num_distortion_params)
    vw_throw( IOErr() << class_name() << ": Incorrect number of parameters was passed in.");
  m_can_undistort = false;
}

Vector<double>
RPCLensDistortion4::distortion_parameters() const { 
  return m_distortion; 
}

Vector<double>
RPCLensDistortion4::undistortion_parameters() const { 
  return m_undistortion; 
}

void RPCLensDistortion4::set_distortion_parameters(Vector<double> const& params) {
  if (params.size() != num_distortion_params) 
    vw_throw( IOErr() << class_name() << ": Incorrect number of parameters was passed in.");
  // If the distortion parameters changed, one cannot undistort until the undistortion
  // coefficients are computed.
  for (size_t it = 0; it < num_distortion_params; it++) {
    if (m_distortion[it] != params[it]) {
      m_can_undistort = false;
      break;
    }
  }
  m_distortion = params;
}

void RPCLensDistortion4::set_undistortion_parameters(Vector<double> const& params) {
  m_undistortion = params;
  m_can_undistort = true;
}

void RPCLensDistortion4::set_image_size(Vector2i const& image_size){
  m_image_size = image_size;
}

boost::shared_ptr<LensDistortion>
RPCLensDistortion4::copy() const {
  return boost::shared_ptr<RPCLensDistortion4>(new RPCLensDistortion4(*this));
}

Vector2
RPCLensDistortion4::distorted_coordinates(const camera::PinholeModel& cam, Vector2 const& p) const {
  return RPCLensDistortion4::compute_rpc(p, m_distortion);
}

Vector2
RPCLensDistortion4::undistorted_coordinates(const camera::PinholeModel& cam, Vector2 const& v) const {
  if (!m_can_undistort) 
    vw_throw( IOErr() << class_name() << ": Undistorted coefficients are not up to date.\n" );
  
  return RPCLensDistortion4::compute_rpc(v, m_undistortion);
}

// Compute the RPC transform. Note that if the RPC coefficients are all zero,
// the obtained transform is the identity.
Vector2
RPCLensDistortion4::compute_rpc(Vector2 const& p, Vector<double> const& coeffs) const {
  
  if (num_distortion_params != coeffs.size()) 
    vw_throw( IOErr() << "Book-keeping failure in RPCLensDistortion4.\n" );

  int i = 0;
  double a00 = coeffs[i]; i++;
  double a10 = coeffs[i]; i++;
  double a01 = coeffs[i]; i++;
  double a20 = coeffs[i]; i++;
  double a11 = coeffs[i]; i++;
  double a02 = coeffs[i]; i++;
  double a30 = coeffs[i]; i++;
  double a21 = coeffs[i]; i++;
  double a12 = coeffs[i]; i++;
  double a03 = coeffs[i]; i++;
  double a40 = coeffs[i]; i++;
  double a31 = coeffs[i]; i++;
  double a22 = coeffs[i]; i++;
  double a13 = coeffs[i]; i++;
  double a04 = coeffs[i]; i++;

  double b10 = coeffs[i]; i++;
  double b01 = coeffs[i]; i++;
  double b20 = coeffs[i]; i++;
  double b11 = coeffs[i]; i++;
  double b02 = coeffs[i]; i++;
  double b30 = coeffs[i]; i++;
  double b21 = coeffs[i]; i++;
  double b12 = coeffs[i]; i++;
  double b03 = coeffs[i]; i++;
  double b40 = coeffs[i]; i++;
  double b31 = coeffs[i]; i++;
  double b22 = coeffs[i]; i++;
  double b13 = coeffs[i]; i++;
  double b04 = coeffs[i]; i++;

  double c00 = coeffs[i]; i++;
  double c10 = coeffs[i]; i++;
  double c01 = coeffs[i]; i++;
  double c20 = coeffs[i]; i++;
  double c11 = coeffs[i]; i++;
  double c02 = coeffs[i]; i++;
  double c30 = coeffs[i]; i++;
  double c21 = coeffs[i]; i++;
  double c12 = coeffs[i]; i++;
  double c03 = coeffs[i]; i++;
  double c40 = coeffs[i]; i++;
  double c31 = coeffs[i]; i++;
  double c22 = coeffs[i]; i++;
  double c13 = coeffs[i]; i++;
  double c04 = coeffs[i]; i++;

  double d10 = coeffs[i]; i++;
  double d01 = coeffs[i]; i++;
  double d20 = coeffs[i]; i++;
  double d11 = coeffs[i]; i++;
  double d02 = coeffs[i]; i++;
  double d30 = coeffs[i]; i++;
  double d21 = coeffs[i]; i++;
  double d12 = coeffs[i]; i++;
  double d03 = coeffs[i]; i++;
  double d40 = coeffs[i]; i++;
  double d31 = coeffs[i]; i++;
  double d22 = coeffs[i]; i++;
  double d13 = coeffs[i]; i++;
  double d04 = coeffs[i]; i++;

  if (size_t(i) != coeffs.size()) 
    vw_throw( IOErr() << "Book-keeping failure in RPCLensDistortion4.\n" );
  
  double x = p[0];
  double y = p[1];

  // Note how we add 1.0 to a10 and a01.
  
  double x_corr =
    (a00 + (1.0+a10)*x + a01*y + a20*x*x + a11*x*y + a02*y*y +
     a30*x*x*x + a21*x*x*y + a12*x*y*y + a03*y*y*y +
     (a40*x*x*x*x + a31*x*x*x*y + a22*x*x*y*y + a13*x*y*y*y + + a04*y*y*y*y)/10000 ) /
      (1 + b10*x + b01*y + b20*x*x + b11*x*y + b02*y*y +
       b30*x*x*x + b21*x*x*y + b12*x*y*y + b03*y*y*y +
       (b40*x*x*x*x + b31*x*x*x*y + b22*x*x*y*y + b13*x*y*y*y + + b04*y*y*y*y)/10000 
       );
  
  double y_corr =
    (c00 + c10*x + (1.0 + c01)*y + c20*x*x + c11*x*y + c02*y*y +
     c30*x*x*x + c21*x*x*y + c12*x*y*y + c03*y*y*y +
     (c40*x*x*x*x + c31*x*x*x*y + c22*x*x*y*y + c13*x*y*y*y + + c04*y*y*y*y)/10000 
     ) /
      (1 + d10*x + d01*y + d20*x*x + d11*x*y + d02*y*y +
       d30*x*x*x + d21*x*x*y + d12*x*y*y + d03*y*y*y+
       (d40*x*x*x*x + d31*x*x*x*y + d22*x*x*y*y + d13*x*y*y*y + + d04*y*y*y*y)/10000 
       );

  return Vector2(x_corr, y_corr);
}

void RPCLensDistortion4::write(std::ostream & os) const {

  if (!m_can_undistort) 
    vw_throw( IOErr() << class_name() << ": Undistorted coefficients are not up to date.\n" );

  // TODO: Add domain of validity for the distored and undistorted pixels. This is needed
  // because otherwise the RPC model can return wrong results which confuse
  // bundle adjustment. 
  write_param_vec("image_size", os, m_image_size);
  write_param_vec("distortion_coeffs", os, m_distortion);
  write_param_vec("undistortion_coeffs", os, m_undistortion);
}

void RPCLensDistortion4::read(std::istream & is) {

  m_image_size.set_size(2);
  m_distortion.set_size(num_distortion_params);
  m_undistortion.set_size(num_distortion_params);

  read_param_vec("image_size", m_image_size.size(), is, m_image_size);
  read_param_vec("distortion_coeffs", num_distortion_params, is, m_distortion);
  read_param_vec("undistortion_coeffs", num_distortion_params, is, m_undistortion);

  m_can_undistort = true; 
}

void RPCLensDistortion4::scale( double scale ) {
  m_distortion *= scale;
  m_undistortion *= scale;
}

// Old RPC lens distortion of degree 5.

RPCLensDistortion5::RPCLensDistortion5(){
  m_distortion.set_size(num_distortion_params);
  m_undistortion.set_size(num_distortion_params);
  m_can_undistort = false;
}

RPCLensDistortion5::RPCLensDistortion5(Vector<double> const& params) 
  : m_distortion(params) {
  if (m_distortion.size() != num_distortion_params)
    vw_throw( IOErr() << class_name() << ": Incorrect number of parameters was passed in.");
  m_can_undistort = false;
}

Vector<double>
RPCLensDistortion5::distortion_parameters() const { 
  return m_distortion; 
}

Vector<double>
RPCLensDistortion5::undistortion_parameters() const { 
  return m_undistortion; 
}

void RPCLensDistortion5::set_distortion_parameters(Vector<double> const& params) {
  if (params.size() != num_distortion_params) 
    vw_throw( IOErr() << class_name() << ": Incorrect number of parameters was passed in.");
  // If the distortion parameters changed, one cannot undistort until the undistortion
  // coefficients are computed.
  for (size_t it = 0; it < num_distortion_params; it++) {
    if (m_distortion[it] != params[it]) {
      m_can_undistort = false;
      break;
    }
  }
  m_distortion = params;
}

void RPCLensDistortion5::set_undistortion_parameters(Vector<double> const& params) {
  m_undistortion = params;
  m_can_undistort = true;
}

void RPCLensDistortion5::set_image_size(Vector2i const& image_size){
  m_image_size = image_size;
}

boost::shared_ptr<LensDistortion>
RPCLensDistortion5::copy() const {
  return boost::shared_ptr<RPCLensDistortion5>(new RPCLensDistortion5(*this));
}

Vector2
RPCLensDistortion5::distorted_coordinates(const camera::PinholeModel& cam, Vector2 const& p) const {
  return RPCLensDistortion5::compute_rpc(p, m_distortion);
}

Vector2
RPCLensDistortion5::undistorted_coordinates(const camera::PinholeModel& cam, Vector2 const& v) const {
  if (!m_can_undistort) 
    vw_throw( IOErr() << class_name() << ": Undistorted coefficients are not up to date.\n" );
  
  return RPCLensDistortion5::compute_rpc(v, m_undistortion);
}

// Compute the RPC transform. Note that if the RPC coefficients are all zero,
// the obtained transform is the identity.
Vector2
RPCLensDistortion5::compute_rpc(Vector2 const& p, Vector<double> const& coeffs) const {
  
  if (num_distortion_params != coeffs.size()) 
    vw_throw( IOErr() << "Book-keeping failure in RPCLensDistortion5.\n" );

  int i = 0;
  double a00 = coeffs[i]; i++;
  double a10 = coeffs[i]; i++;
  double a01 = coeffs[i]; i++;
  double a20 = coeffs[i]; i++;
  double a11 = coeffs[i]; i++;
  double a02 = coeffs[i]; i++;
  double a30 = coeffs[i]; i++;
  double a21 = coeffs[i]; i++;
  double a12 = coeffs[i]; i++;
  double a03 = coeffs[i]; i++;
  double a40 = coeffs[i]; i++;
  double a31 = coeffs[i]; i++;
  double a22 = coeffs[i]; i++;
  double a13 = coeffs[i]; i++;
  double a04 = coeffs[i]; i++;

  double b10 = coeffs[i]; i++;
  double b01 = coeffs[i]; i++;
  double b20 = coeffs[i]; i++;
  double b11 = coeffs[i]; i++;
  double b02 = coeffs[i]; i++;
  double b30 = coeffs[i]; i++;
  double b21 = coeffs[i]; i++;
  double b12 = coeffs[i]; i++;
  double b03 = coeffs[i]; i++;
  double b40 = coeffs[i]; i++;
  double b31 = coeffs[i]; i++;
  double b22 = coeffs[i]; i++;
  double b13 = coeffs[i]; i++;
  double b04 = coeffs[i]; i++;

  double c00 = coeffs[i]; i++;
  double c10 = coeffs[i]; i++;
  double c01 = coeffs[i]; i++;
  double c20 = coeffs[i]; i++;
  double c11 = coeffs[i]; i++;
  double c02 = coeffs[i]; i++;
  double c30 = coeffs[i]; i++;
  double c21 = coeffs[i]; i++;
  double c12 = coeffs[i]; i++;
  double c03 = coeffs[i]; i++;
  double c40 = coeffs[i]; i++;
  double c31 = coeffs[i]; i++;
  double c22 = coeffs[i]; i++;
  double c13 = coeffs[i]; i++;
  double c04 = coeffs[i]; i++;

  double d10 = coeffs[i]; i++;
  double d01 = coeffs[i]; i++;
  double d20 = coeffs[i]; i++;
  double d11 = coeffs[i]; i++;
  double d02 = coeffs[i]; i++;
  double d30 = coeffs[i]; i++;
  double d21 = coeffs[i]; i++;
  double d12 = coeffs[i]; i++;
  double d03 = coeffs[i]; i++;
  double d40 = coeffs[i]; i++;
  double d31 = coeffs[i]; i++;
  double d22 = coeffs[i]; i++;
  double d13 = coeffs[i]; i++;
  double d04 = coeffs[i]; i++;

  double a50 = coeffs[i]; i++;
  double a41 = coeffs[i]; i++;
  double a32 = coeffs[i]; i++;
  double a23 = coeffs[i]; i++;
  double a14 = coeffs[i]; i++;
  double a05 = coeffs[i]; i++;

  double b50 = coeffs[i]; i++;
  double b41 = coeffs[i]; i++;
  double b32 = coeffs[i]; i++;
  double b23 = coeffs[i]; i++;
  double b14 = coeffs[i]; i++;
  double b05 = coeffs[i]; i++;

  double c50 = coeffs[i]; i++;
  double c41 = coeffs[i]; i++;
  double c32 = coeffs[i]; i++;
  double c23 = coeffs[i]; i++;
  double c14 = coeffs[i]; i++;
  double c05 = coeffs[i]; i++;

  double d50 = coeffs[i]; i++;
  double d41 = coeffs[i]; i++;
  double d32 = coeffs[i]; i++;
  double d23 = coeffs[i]; i++;
  double d14 = coeffs[i]; i++;
  double d05 = coeffs[i]; i++;

  if (size_t(i) != coeffs.size()) 
    vw_throw( IOErr() << "Book-keeping failure in RPCLensDistortion5.\n" );
  
  double x = p[0];
  double y = p[1];

  // Note how we add 1.0 to a10 and a01.
  
  double x_corr =
    (a00 + (1.0+a10)*x + a01*y + a20*x*x + a11*x*y + a02*y*y +
     a30*x*x*x + a21*x*x*y + a12*x*y*y + a03*y*y*y +
     (a40*x*x*x*x + a31*x*x*x*y + a22*x*x*y*y + a13*x*y*y*y + + a04*y*y*y*y)/10000 +
     (a50*x*x*x*x*x + a41*x*x*x*x*y + a32*x*x*x*y*y + a23*x*x*y*y*y + + a14*x*y*y*y*y + a05*y*y*y*y*y)/1000000 
     ) /
      (1 + b10*x + b01*y + b20*x*x + b11*x*y + b02*y*y +
       b30*x*x*x + b21*x*x*y + b12*x*y*y + b03*y*y*y +
       (b40*x*x*x*x + b31*x*x*x*y + b22*x*x*y*y + b13*x*y*y*y + + b04*y*y*y*y)/10000 +
       (b50*x*x*x*x*x + b41*x*x*x*x*y + b32*x*x*x*y*y + b23*x*x*y*y*y + + b14*x*y*y*y*y + b05*y*y*y*y*y)/1000000 
       );
  
  double y_corr =
    (c00 + c10*x + (1.0 + c01)*y + c20*x*x + c11*x*y + c02*y*y +
     c30*x*x*x + c21*x*x*y + c12*x*y*y + c03*y*y*y +
     (c40*x*x*x*x + c31*x*x*x*y + c22*x*x*y*y + c13*x*y*y*y + + c04*y*y*y*y)/10000 +
     (c50*x*x*x*x*x + c41*x*x*x*x*y + c32*x*x*x*y*y + c23*x*x*y*y*y + + c14*x*y*y*y*y + c05*y*y*y*y*y)/1000000 
     ) /
      (1 + d10*x + d01*y + d20*x*x + d11*x*y + d02*y*y +
       d30*x*x*x + d21*x*x*y + d12*x*y*y + d03*y*y*y+
       (d40*x*x*x*x + d31*x*x*x*y + d22*x*x*y*y + d13*x*y*y*y + + d04*y*y*y*y)/10000 +
       (d50*x*x*x*x*x + d41*x*x*x*x*y + d32*x*x*x*y*y + d23*x*x*y*y*y + + d14*x*y*y*y*y + d05*y*y*y*y*y)/1000000 
       );

  return Vector2(x_corr, y_corr);
}

void RPCLensDistortion5::write(std::ostream & os) const {

  if (!m_can_undistort) 
    vw_throw( IOErr() << class_name() << ": Undistorted coefficients are not up to date.\n" );

  // TODO: Add domain of validity for the distored and undistorted pixels. This is needed
  // because otherwise the RPC model can return wrong results which confuse
  // bundle adjustment. 
  write_param_vec("image_size", os, m_image_size);
  write_param_vec("distortion_coeffs", os, m_distortion);
  write_param_vec("undistortion_coeffs", os, m_undistortion);
}

void RPCLensDistortion5::read(std::istream & is) {
  
  m_image_size.set_size(2);
  m_distortion.set_size(num_distortion_params);
  m_undistortion.set_size(num_distortion_params);

  read_param_vec("image_size", m_image_size.size(), is, m_image_size);
  read_param_vec("distortion_coeffs", num_distortion_params, is, m_distortion);
  read_param_vec("undistortion_coeffs", num_distortion_params, is, m_undistortion);

  m_can_undistort = true;
}

void RPCLensDistortion5::scale( double scale ) {
  m_distortion *= scale;
  m_undistortion *= scale;
}
