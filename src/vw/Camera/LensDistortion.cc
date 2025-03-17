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

#include <vw/Camera/LensDistortion.h>
#include <vw/Camera/PinholeModel.h>
#include <vw/Math/LevenbergMarquardt.h>
#include <vw/Math/NewtonRaphson.h>

#include <set>
using namespace vw;
using namespace camera;

// Special LMA Models to reverse the distortion.
// TODO(oalexan1): Standardize all to using Newton-Raphson.
struct DistortOptimizeFunctor:  public math::LeastSquaresModelBaseFixed<DistortOptimizeFunctor, 2, 2> {
  typedef Vector2 result_type;
  typedef Vector2 domain_type;
  typedef Matrix<double, 2, 2> jacobian_type;

  const camera::PinholeModel   &m_cam;
  const camera::LensDistortion &m_distort;
  DistortOptimizeFunctor(const camera::PinholeModel& cam, const camera::LensDistortion& d):
    m_cam(cam), m_distort(d) {}
  inline result_type operator()(domain_type const& x) const {
    return m_distort.undistorted_coordinates(m_cam, x);
  }
};

// Isolate these local utilities to this file.
namespace {

// Pull all lines from the stream. Search for "name = val". Store the values in
// the order given in "names". Complain if some fields were not populated,
// unless in the set named missing_ok, in which case they are set to zero.
// This has the advantage that the order of lines in
// the stream is not important, and it won't complain if there are
// extraneous fields.
template<class VectorT>
void read_fields_in_vec(std::vector<std::string> const& names, VectorT & vals,
                        std::istream & is,
                        std::set<std::string> const& missing_ok = std::set<std::string>()) {

  std::map<std::string, double> name2val;
  std::string line;
  char name[1024];
  double val;
  while (is.good()) {
    std::getline(is, line);
    // Replace "=" with a space. This way need not worry about spaces around the equal sign.
    for (size_t i = 0; i < line.size(); i++) {
      if (line[i] == '=')
        line[i] = ' ';
    }
    if (sscanf(line.c_str(),"%s %lf", name, &val) != 2)
      continue;

    name2val[name] = val;
  }

  // Populate the output
  for (size_t i = 0; i < names.size(); i++) {
    if (i >= vals.size() + missing_ok.size())
      vw_throw(IOErr() << "Not enough room allocated for output.\n");

    std::map<std::string, double>::iterator it = name2val.find(names[i]);
    if (it == name2val.end()) {
      // See if this is in the missing_ok set, then set to zero
      if (missing_ok.find(names[i]) != missing_ok.end()) {
        vals[i] = 0;
        continue;
      }
      vw_throw(IOErr() << "LensDistortion: Could not read a value for "
                       << names[i] << ".\n");
    }

    vals[i] = it->second;
  }
}

// Read a line of the form: name = a b c. The number of values is not known.
template<class VectorT>
void read_param_vec(std::string const& param_name, std::istream & is, VectorT & vals) {

  if (!is.good())
    vw_throw(IOErr() << "Could not read LensDistortion.\n");

  std::string line;
  std::getline(is, line);
  std::istringstream iss(line);

  std::string val;
  if (!(iss >> val) || val != param_name)
    vw_throw(IOErr() << "In LensDistortion, got " << val << " but expected " << param_name);

  if (!(iss >> val)  || val != "=")
    vw_throw(IOErr() << "In LensDistortion, got " << val << " but expected the equal sign.");

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
                    std::istream & is, VectorT & vals) {

  read_param_vec(param_name, is, vals);
  if (vals.size() != param_len) {
    vw_throw(IOErr() << "Book-keeping failure in reading LensDistortion. Expecting "
              << vals.size() << " == " << param_len << "\n");
  }
}

// Read a line of the form: name = a.
template<class T>
void read_param(std::string const& param_name, std::istream & is, T & val) {
  Vector<T> vals;
  read_param_vec(param_name, 1, is, vals);
  val = vals[0];
}

// Write a line of the form: name = a b c
void write_param_vec(std::string const& param_name, std::ostream & os, Vector<double> const& vals) {
  os << param_name << " = ";
  for (size_t p = 0; p < vals.size(); p++) {
    os << vals[p];
    if (p + 1 < vals.size())
      os << " "; // write a whitespace after each number except the last
  }
  os << "\n";
}

// Write a line of the form: name = val
void write_param(std::string const& param_name, std::ostream & os, double val) {
  Vector<double> vals;
  vals.set_size(1);
  vals[0] = val;
  write_param_vec(param_name, os, vals);
}

} // end of anonymous namespace

// Default implementations for Lens Distortion

Vector<double>
LensDistortion::distortion_parameters() const { return Vector<double>(); }

void LensDistortion::set_distortion_parameters(Vector<double> const& params) {}

// Base class distorted_coordinates implementation. Many lens distortion
// models override this with a faster implementation.
vw::Vector2
LensDistortion::distorted_coordinates(const PinholeModel& cam, Vector2 const& v) const {
  DistortOptimizeFunctor model(cam, *this);
  int status;
  // Must push the solver really hard, to make sure bundle adjust gets accurate values
  // to play with.
  Vector2 solution = math::levenberg_marquardtFixed(model, v, v, status, 1e-16, 1e-16, 100);
  //VW_DEBUG_ASSERT(status != math::optimization::eConvergedRelTolerance,
  //                 PixelToRayErr() << "distorted_coordinates: failed to converge.");

  // Check if it failed badly to converge. That it did not converge is not on its own
  // unreasonable, sometimes the inputs are bad. But then the user must know about it.
  // TODO(oalexan1): Look at this
  Vector2 undist = this->undistorted_coordinates(cam, solution);
  double err = norm_2(undist - v)/std::max(norm_2(v), 0.1); // don't make this way too strict
  double tol = 1e-10;
  if (err > tol)
    vw_throw(PointToPixelErr() << "LensDistortion: Did not converge.\n");

  return solution;
}

std::ostream& camera::operator<<(std::ostream & os,
                                 const camera::LensDistortion& ld) {
  ld.write(os);
  return os;
}

// Specific Implementations

// NullLensDistortion 

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

// TsaiLensDistortion 

TsaiLensDistortion::TsaiLensDistortion() {
  TsaiLensDistortion::init_distortion_param_names();
  m_distortion.set_size(m_distortion_param_names.size());
}

TsaiLensDistortion::TsaiLensDistortion(Vector<double> const& params): m_distortion(params) {
  TsaiLensDistortion::init_distortion_param_names();
  if (m_distortion.size() != m_distortion_param_names.size())
    vw_throw(IOErr() << "TsaiLensDistortion: Incorrect number of parameters was passed in.");
}

Vector<double>
TsaiLensDistortion::distortion_parameters() const { return m_distortion; }

void TsaiLensDistortion::init_distortion_param_names() {
  std::string names[] = {"k1", "k2", "p1", "p2", "k3"}; // k3 must be last
  size_t num_names = sizeof(names)/sizeof(std::string);
  m_distortion_param_names.resize(num_names);
  for (size_t p = 0; p < num_names; p++)
    m_distortion_param_names[p] = names[p];
}

void TsaiLensDistortion::set_distortion_parameters(Vector<double> const& params) {
  TsaiLensDistortion::init_distortion_param_names();
  m_distortion = params;
  if (m_distortion.size() != m_distortion_param_names.size())
    vw_throw(IOErr() << "TsaiLensDistortion: Incorrect number of parameters was passed in.");
}

boost::shared_ptr<LensDistortion>
TsaiLensDistortion::copy() const {
  return boost::shared_ptr<TsaiLensDistortion>(new TsaiLensDistortion(*this));
}

// Tsai distortion model, after normalizing the point to the unit focal plane
vw::Vector2 TsaiDistortionNorm(vw::Vector2 const& P, vw::Vector<double> const& distortion) {
  
  double x = P[0];
  double y = P[1];
  double k1 = distortion[0];
  double k2 = distortion[1];
  double p1 = distortion[2];
  double p2 = distortion[3];
  double k3 = distortion[4]; // k3 must be last

  double r2 = x * x + y * y;
  double rdist = 1.0 + k1 * r2 + k2 * r2 * r2 + k3 * r2 * r2 * r2;
  double x_out = x * rdist + (2.0 * p1 * x * y + p2 * (r2 + 2.0 * x * x));
  double y_out = y * rdist + (p1 * (r2 + 2.0 * y * y) + 2.0 * p2 * x * y);
  
  return vw::Vector2(x_out, y_out);
}

// Apply the distortion to a normalized pixel a function object. To be used in
// Newton-Raphson.
vw::Vector2 TsaiLensDistortion::operator()(vw::Vector2 const& P) const {
  return TsaiDistortionNorm(P, m_distortion);
}

// Compute the Jacobian for the Tsai distortion model. The step
// size is not used but is part of the interface.
vw::Vector<double> TsaiDistortionJacobian(vw::Vector2 const& P, double step,
                                          vw::Vector<double> const& distortion) {

  // The Jacobian has 4 elements
  vw::Vector<double> jacobian(4);

  double x = P[0];
  double y = P[1];
  double k1 = distortion[0];
  double k2 = distortion[1];
  double p1 = distortion[2];
  double p2 = distortion[3];
  double k3 = distortion[4]; // k3 must be last

  double r2 = x * x + y * y;
  double dr2dx = 2.0 * x;
  double dr2dy = 2.0 * y;
  double rdist = 1.0 + k1 * r2 + k2 * r2 * r2 + k3 * r2 * r2 * r2;

  // dfx / dx
  jacobian[0] = rdist
              + x * (k1 * dr2dx + k2 * dr2dx * 2.0 * r2 + k3 * dr2dx * 3.0 * r2 * r2)
              + 2.0 * p1 * y + p2 * (dr2dx + 4.0 * x);

  // dfx / dy
  jacobian[1] = x * (k1 * dr2dy + k2 * dr2dy * 2.0 * r2 + k3 * dr2dy * 3.0 * r2 * r2)
              + 2.0 * p1 * x  + p2 * dr2dy;

  // dfy / dx
  jacobian[2] = y * (k1 * dr2dx + k2 * dr2dx * 2.0 * r2 + k3 * dr2dx * 3.0 * r2 * r2)
              + (p1 * dr2dx + 2.0 * p2 * y);

  // dfy / dy
  jacobian[3] = rdist
              + y * (k1 * dr2dy + k2 * dr2dy * 2.0 * r2 + k3 * dr2dy * 3.0 * r2 * r2)
              + p1 * (dr2dy + 4.0 * y) + 2.0 * p2 * x;
              
  return jacobian;
}

// Function object around the the Jacobian for the Tsai distortion model.
// Needed by the Newton-Raphson solver.
struct TsaiDistortionJacFun {
  
  TsaiDistortionJacFun(vw::Vector<double> const& distortion): m_distortion(distortion) {}
  
  vw::Vector<double> operator()(vw::Vector2 const& P, double step) {
    return TsaiDistortionJacobian(P, step, m_distortion);
  }
  
  vw::Vector<double> m_distortion;
};

// This was validated to be in perfect agreement with the OpenCV implementation.
// https://docs.opencv.org/4.x/d9/d0c/group__calib3d.html
// The function cv::projectPoints() was used for validation, with no rotation
// or translation.
// For comparing with rig_calibrator, see the convention at:
// Vector2 FisheyeLensDistortion::distorted_coordinates().
Vector2 TsaiLensDistortion::distorted_coordinates(const PinholeModel& cam,
                                                  Vector2 const& p) const {

  Vector2 focal  = cam.focal_length(); // = [fu, fv]
  Vector2 offset = cam.point_offset(); // = [cu, cv]

  if (focal[0] < 1e-300 || focal[1] < 1e-300)
    return Vector2(HUGE_VAL, HUGE_VAL);

  // Normalize the pixel
  Vector2 dudv = p - offset; // Subtract the offset
  Vector2 p_0 = elem_quot(dudv, focal); // Divide by focal length
  double x = p_0[0];
  double y = p_0[1];

  // Normalized distorted coordinates
  vw::Vector2 dist_norm_pix = TsaiDistortionNorm(vw::Vector2(x, y), m_distortion); 
  double dx = dist_norm_pix[0];
  double dy = dist_norm_pix[1];

  // Multiply by focal length and add the offset
  dx = dx * focal[0] + offset[0];
  dy = dy * focal[1] + offset[1];
  return Vector2(dx, dy);
}

Vector2 TsaiLensDistortion::undistorted_coordinates(const PinholeModel& cam,
                                                    Vector2 const& p) const {

  Vector2 focal  = cam.focal_length(); // = [fu, fv]
  Vector2 offset = cam.point_offset(); // = [cu, cv]

  if (focal[0] < 1e-300 || focal[1] < 1e-300)
    return Vector2(HUGE_VAL, HUGE_VAL);

  // Normalize the distorted pixel
  Vector2 dudv = p - offset; // Subtract the offset
  Vector2 p_0 = elem_quot(dudv, focal); // Divide by focal length
  double dx = p_0[0];
  double dy = p_0[1];
  
  double step = 1e-6; // not used, part of the interface
  double tol = 1e-9; // stop when the change is less than this
  
  // Newton-Raphson solver with the analytical Jacobian
  vw::math::NewtonRaphson nr(*this, TsaiDistortionJacFun(m_distortion));
  vw::Vector2 undist_guess_pix(dx, dy), distorted_pix(dx, dy);
  Vector2 U = nr.solve(undist_guess_pix, distorted_pix, step, tol);

  // Undo the normalization
  double ux = U[0], uy = U[1];
  ux = ux * focal[0] + offset[0];
  uy = uy * focal[1] + offset[1];

  return Vector2(ux, uy);
}

void TsaiLensDistortion::write(std::ostream & os) const {
  for (size_t p = 0; p < m_distortion_param_names.size(); p++)
    os << m_distortion_param_names[p] << " = " << m_distortion[p] << "\n";
}

void TsaiLensDistortion::read(std::istream & is) {
  m_distortion.set_size(m_distortion_param_names.size());
  std::set<std::string> missing_ok;
  missing_ok.insert("k3"); // this may be missing, then set to zero
  read_fields_in_vec(m_distortion_param_names, m_distortion, is, missing_ok);
}

void TsaiLensDistortion::scale(double scale) {
  m_distortion *= scale;
}

// FovLensDistortion 
// Single-parameter wide-angle lens distortion model.
FovLensDistortion::FovLensDistortion() {
  FovLensDistortion::init_distortion_param_names();
  m_distortion.set_size(m_distortion_param_names.size());
}

FovLensDistortion::FovLensDistortion(Vector<double> const& params) {
  FovLensDistortion::init_distortion_param_names();
  set_distortion_parameters(params);
}

Vector<double>
FovLensDistortion::distortion_parameters() const { return m_distortion; }

void FovLensDistortion::init_distortion_param_names() {
  std::string names[] = {"k1"};
  size_t num_names = sizeof(names)/sizeof(std::string);
  m_distortion_param_names.resize(num_names);
  for (size_t p = 0; p < num_names; p++)
    m_distortion_param_names[p] = names[p];
}

void FovLensDistortion::set_distortion_parameters(Vector<double> const& params) {
  FovLensDistortion::init_distortion_param_names();
  m_distortion = params;
  if (m_distortion.size() != m_distortion_param_names.size())
    vw_throw(IOErr() << "FovLensDistortion: Incorrect number of parameters was passed in.");

   // Distortion coefficient must be positive
   if (m_distortion[0] <= 0)
     vw_throw(IOErr() << "FovLensDistortion: The distortion coefficient "
                        << "must be positive.\n");

    // Precompute some values
    m_distortion_precalc1 = 1 / m_distortion[0];
    m_distortion_precalc2 = 2 * tan(m_distortion[0] / 2);
}

boost::shared_ptr<LensDistortion>
FovLensDistortion::copy() const {
  return boost::shared_ptr<FovLensDistortion>(new FovLensDistortion(*this));
}

Vector2 FovLensDistortion::distorted_coordinates(const PinholeModel& cam,
                                                  Vector2 const& p) const {

  Vector2 focal  = cam.focal_length(); // = [fu, fv]
  Vector2 offset = cam.point_offset(); // = [cu, cv]

  if (focal[0] < 1e-300 || focal[1] < 1e-300)
    return Vector2(HUGE_VAL, HUGE_VAL);

  // Normalize the pixel
  Vector2 dudv = p - offset; // Subtract the offset
  Vector2 p_0 = elem_quot(dudv, focal); // Divide by focal length

  double ru = norm_2(p_0);
  double rd = atan(ru * m_distortion_precalc2) * m_distortion_precalc1;
  double conv = 1.0;
   if (ru > 1e-8)
    conv = rd / ru; // Avoid division by a small number

  // Apply the fov distortion model to the normalized pixel
  Vector2 p_norm  = conv * p_0;

  // Multiply by focal length and add the offset
  Vector2 p_dist = elem_prod(p_norm, focal) + offset;

  return p_dist;
}

Vector2 FovLensDistortion::undistorted_coordinates(const PinholeModel& cam,
                                                   Vector2 const& p) const {

  Vector2 focal  = cam.focal_length(); // = [fu, fv]
  Vector2 offset = cam.point_offset(); // = [cu, cv]

  if (focal[0] < 1e-300 || focal[1] < 1e-300)
    return Vector2(HUGE_VAL, HUGE_VAL);

  Vector2 dudv = p - offset; // Subtract the offset
  Vector2 p_0 = elem_quot(dudv, focal); // Divide by focal length

  double rd = norm_2(p_0);
  double ru = tan(rd * m_distortion[0]) / m_distortion_precalc2;
  double conv = 1.0;
  if (rd > 1e-8)
    conv = ru / rd; // Avoid division by a small number

  // Apply the fov distortion model to the normalized pixel
  Vector2 p_norm = conv * p_0;

  // Multiply by focal length and add the offset
  Vector2 p_undist = elem_prod(p_norm, focal) + offset;

  return p_undist;
}

void FovLensDistortion::write(std::ostream & os) const {
  for (size_t p = 0; p < m_distortion_param_names.size(); p++)
    os << m_distortion_param_names[p] << " = " << m_distortion[p] << "\n";
}

void FovLensDistortion::read(std::istream & is) {
  m_distortion.set_size(m_distortion_param_names.size());
  read_fields_in_vec(m_distortion_param_names, m_distortion, is);
  set_distortion_parameters(m_distortion); // This does some pre-computations and checks
}

void FovLensDistortion::scale(double scale) {
  vw::vw_throw(vw::NoImplErr() << "FovLensDistortion::scale() is not implemented.");
}

// FisheyeLensDistortion 
// Four-parameter wide-angle lens distortion model.
FisheyeLensDistortion::FisheyeLensDistortion() {
  FisheyeLensDistortion::init_distortion_param_names();
  m_distortion.set_size(m_distortion_param_names.size());
}

FisheyeLensDistortion::FisheyeLensDistortion(Vector<double> const& params):
  m_distortion(params) {
  FisheyeLensDistortion::init_distortion_param_names();
  if (m_distortion.size() != m_distortion_param_names.size())
    vw_throw(IOErr() << "FisheyeLensDistortion: Incorrect number of parameters "
                     << "was passed in.");
}

Vector<double> FisheyeLensDistortion::distortion_parameters() const {
  return m_distortion;
}

void FisheyeLensDistortion::init_distortion_param_names() {
  std::string names[] = {"k1", "k2", "k3", "k4"};
  size_t num_names = sizeof(names)/sizeof(std::string);
  m_distortion_param_names.resize(num_names);
  for (size_t p = 0; p < num_names; p++)
    m_distortion_param_names[p] = names[p];
}

void FisheyeLensDistortion::set_distortion_parameters(Vector<double> const& params) {
  FisheyeLensDistortion::init_distortion_param_names();
  m_distortion = params;
  if (m_distortion.size() != m_distortion_param_names.size())
    vw_throw(IOErr() << "FisheyeLensDistortion: Incorrect number of parameters was "
                     << "passed in.");
}

boost::shared_ptr<LensDistortion>
FisheyeLensDistortion::copy() const {
  return boost::shared_ptr<FisheyeLensDistortion>(new FisheyeLensDistortion(*this));
}

// Apply the fisheye distortion model. Input and output are normalized pixels.
vw::Vector2 fishEyeDistortionNorm(vw::Vector2 const& P, vw::Vector<double> const& dist) {

  double k1 = dist[0];
  double k2 = dist[1];
  double k3 = dist[2];
  double k4 = dist[3];

  double x = P[0];
  double y = P[1];
  double r2 = x*x + y*y;
  double r = sqrt(r2);
  double theta = atan(r);

  double theta1 = theta*theta;   // theta^2
  double theta2 = theta1*theta1; // theta^4
  double theta3 = theta2*theta1; // theta^6
  double theta4 = theta2*theta2; // theta^8
  double theta_d = theta*(1 + k1*theta1 + k2*theta2 + k3*theta3 + k4*theta4);

  // Careful with the case where r is very small
  double scale = 1.0;
  if (r > 1e-8)
    scale = theta_d / r;

  return vw::Vector2(x*scale, y*scale);
}

// Apply the distortion to a normalized pixel a function object. To be used in
// Newton-Raphson.
vw::Vector2 FisheyeLensDistortion::operator()(vw::Vector2 const& P) const {
  return fishEyeDistortionNorm(P, m_distortion);
}

// This was validated to be in perfect agreement with rig_calibrator Fisheye
// model, as long as the approach is as follows. (a) Take the undistorted pixel
// and distort it with WV. (b) Take the undistorted pixel, subtract optical
// center, add half the undistorted image size, and then apply the
// rig_calibrator distortion model. For distortion the process is in reverse.
// This is a minor convention difference.
Vector2 FisheyeLensDistortion::distorted_coordinates(const PinholeModel& cam,
                                                     Vector2 const& p) const {

  Vector2 focal  = cam.focal_length(); // = [fu, fv]
  Vector2 offset = cam.point_offset(); // = [cu, cv]

  if (focal[0] < 1e-300 || focal[1] < 1e-300)
    return Vector2(HUGE_VAL, HUGE_VAL);

  // Normalize the pixel
  Vector2 dudv = p - offset; // Subtract the offset
  Vector2 p_0 = elem_quot(dudv, focal); // Divide by focal length

  // Apply the fisheye distortion model to the normalized pixel
  Vector2 p_norm = fishEyeDistortionNorm(p_0, m_distortion);

  // Multiply by focal length and add the offset
  Vector2 p_dist = elem_prod(p_norm, focal) + offset;

  return p_dist;
}

Vector2 FisheyeLensDistortion::undistorted_coordinates(const PinholeModel& cam,
                                                       Vector2 const& p) const {

  Vector2 focal  = cam.focal_length(); // = [fu, fv]
  Vector2 offset = cam.point_offset(); // = [cu, cv]

  if (focal[0] < 1e-300 || focal[1] < 1e-300)
    return Vector2(HUGE_VAL, HUGE_VAL);

  Vector2 dudv = p - offset; // Subtract the offset
  Vector2 p_0 = elem_quot(dudv, focal); // Divide by focal length

  // Find the normalized undistorted pixel using Newton-Raphson and the
  // numerical Jacobian. The step size for numerical differentiation 
  // better not be too small. A tolerance of 1e-9 is likely adequate.
  vw::Vector2 guess = p_0;
  double step = 1e-6, tol = 1e-9;
  vw::math::NewtonRaphson nr(*this);
  Vector2 U = nr.solve(guess, p_0, step, tol);

  // Multiply by focal length and add the offset
  Vector2 p_undist = elem_prod(U, focal) + offset;

  return p_undist;
}

void FisheyeLensDistortion::write(std::ostream & os) const {
  for (size_t p = 0; p < m_distortion_param_names.size(); p++)
    os << m_distortion_param_names[p] << " = " << m_distortion[p] << "\n";
}

void FisheyeLensDistortion::read(std::istream & is) {
  m_distortion.set_size(m_distortion_param_names.size());
  read_fields_in_vec(m_distortion_param_names, m_distortion, is);
  set_distortion_parameters(m_distortion); // This does some checks
}

void FisheyeLensDistortion::scale(double scale) {
  vw::vw_throw(vw::NoImplErr() << "FisheyeLensDistortion::scale() is not implemented.");
}

// BrownConradyDistortion 

BrownConradyDistortion::BrownConradyDistortion() {
  BrownConradyDistortion::init_distortion_param_names();
}

void BrownConradyDistortion::init_distortion_param_names() {
  std::string names[] = {"xp", "yp", "k1", "k2", "k3", "p1", "p2", "phi"};
  size_t num_names = sizeof(names)/sizeof(std::string);
  m_distortion_param_names.resize(num_names);
  for (size_t p = 0; p < num_names; p++)
    m_distortion_param_names[p] = names[p];
}

BrownConradyDistortion::BrownConradyDistortion(Vector<double> const& params) {
  VW_ASSERT(params.size() == 8,
             ArgumentErr() << "BrownConradyDistortion: requires constructor input of size 8.");
  m_principal_point      = subvector(params,0,2);
  m_radial_distortion    = subvector(params,2,3);
  m_centering_distortion = subvector(params,5,2);
  m_centering_angle      = params[7];
  BrownConradyDistortion::init_distortion_param_names();
}

BrownConradyDistortion::BrownConradyDistortion(Vector<double> const& principal,
                                               Vector<double> const& radial,
                                               Vector<double> const& centering,
                                               double const& angle):
  m_principal_point(principal), m_radial_distortion(radial),
  m_centering_distortion(centering), m_centering_angle(angle) {
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

// BrownConrady uses the base class distorted_coordinates implementation.

Vector2
BrownConradyDistortion::undistorted_coordinates(const PinholeModel& cam,
                                                Vector2 const& p) const {
  Vector2 offset       = cam.point_offset();
  Vector2 intermediate = p - m_principal_point - offset;
  double r2            = norm_2_sqr(intermediate);
  double radial        = 1 + m_radial_distortion[0]*r2 + m_radial_distortion[1]*r2*r2 +
                         m_radial_distortion[2]*r2*r2*r2;
  double tangential = m_centering_distortion[0]*r2 + m_centering_distortion[1]*r2*r2;
  intermediate *= radial;
  intermediate[0] -= tangential * sin(m_centering_angle);
  intermediate[1] += tangential * cos(m_centering_angle);
  return intermediate+offset;
}

void BrownConradyDistortion::write(std::ostream & os) const {
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


void BrownConradyDistortion::scale(double scale) {
  vw_throw(NoImplErr() << "BrownConradyDistortion doesn't support scaling");
}

// AdjustableTsaiLensDistortion

AdjustableTsaiLensDistortion::AdjustableTsaiLensDistortion(Vector<double> params):
  m_distortion(params) {
  VW_ASSERT(params.size() > 3, ArgumentErr() << "Requires at least 4 coefficients for distortion. Last 3 are always the distortion coefficients and alpha. All leading elements are even radial distortion coefficients.");
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

// Adjustable Tsai distortion model, after normalizing the point to the unit focal plane
vw::Vector2 AdjTsaiDistortionNorm(vw::Vector2 const& p_0, 
                                  vw::Vector<double> const& distortion) {

  // Calculate radial effects
  double r2 = norm_2_sqr(p_0);
  double r_n = 1, radial = 0;
  for (unsigned i = 0; i < distortion.size()-3; i++) {
    r_n *= r2;
    radial += distortion[i] * r_n;
  }

  // Calculate tangential effects
  Vector2 tangent;
  Vector2 swap_coeff(distortion[distortion.size()-2],
                     distortion[distortion.size()-3]);
  tangent += elem_prod(swap_coeff,elem_sum(r2,2*elem_prod(p_0, p_0)));
  tangent += 2*prod(p_0) * subvector(distortion, distortion.size()-3, 2);

  // Add the contributions
  Vector2 result = p_0 + tangent + radial * p_0;

  // Run back through intrinsic matrix (with alpha or skew)
  result += Vector2(distortion[distortion.size()-1] * result.y(), 0);

  return result;
}

// Apply the distortion to a normalized pixel a function object. To be used in
// Newton-Raphson.
vw::Vector2 AdjustableTsaiLensDistortion::operator()(vw::Vector2 const& P) const {
  return AdjTsaiDistortionNorm(P, m_distortion);
}

Vector2
AdjustableTsaiLensDistortion::distorted_coordinates(const PinholeModel& cam,
                                                    Vector2 const& p) const {
  Vector2 focal = cam.focal_length();
  Vector2 offset = cam.point_offset();

  if (focal[0] < 1e-300 || focal[1] < 1e-300)
    return Vector2(HUGE_VAL, HUGE_VAL);

  // Create normalized undistorted coordinates
  Vector2 p_0 = elem_quot(p - offset, focal); // represents x and y

  // Compute the normalized distorted coordinates
  vw::Vector2 result = AdjTsaiDistortionNorm(p_0, m_distortion);
  
  // Undo the normalization
  return elem_prod(result, focal) + offset;
}

Vector2 AdjustableTsaiLensDistortion::undistorted_coordinates(const PinholeModel& cam,
                                                              Vector2 const& p) const {

  Vector2 focal  = cam.focal_length(); // = [fu, fv]
  Vector2 offset = cam.point_offset(); // = [cu, cv]

  if (focal[0] < 1e-300 || focal[1] < 1e-300)
    return Vector2(HUGE_VAL, HUGE_VAL);

  // Normalize the distorted pixel
  Vector2 dudv = p - offset; // Subtract the offset
  Vector2 p_0 = elem_quot(dudv, focal); // Divide by focal length
  double dx = p_0[0];
  double dy = p_0[1];
  
  double step = 1e-6; // Numerical Jacobian differentiation step
  double tol = 1e-9; // stop when the change is less than this
  
  // Newton-Raphson solver with the numerical Jacobian
  vw::math::NewtonRaphson nr(*this);
  vw::Vector2 undist_guess_pix(dx, dy), distorted_pix(dx, dy);
  Vector2 U = nr.solve(undist_guess_pix, distorted_pix, step, tol);

  // Undo the normalization
  double ux = U[0], uy = U[1];
  ux = ux * focal[0] + offset[0];
  uy = uy * focal[1] + offset[1];

  return Vector2(ux, uy);
}

void AdjustableTsaiLensDistortion::write(std::ostream & os) const {
  // Since we don't know how many parameters, print them out in three rows.
  os << "Radial Coeff: "     << subvector(m_distortion,0,m_distortion.size()-3) << "\n";
  os << "Tangential Coeff: " << subvector(m_distortion,m_distortion.size()-3,2) << "\n";
  os << "Alpha: "            << m_distortion[m_distortion.size()-1] << "\n";
}

void AdjustableTsaiLensDistortion::read(std::istream & is) {
  Vector<double,0> radial_vec, tangential_vec;
  double alpha;

  std::string label, line;
  try {
    is >> label >> label >> radial_vec;
    is >> label >> label >> tangential_vec;
  } catch(...) {
    vw_throw(IOErr()
             << "AdjustableTsaiLensDistortion::read(): Could not read vector params.\n");
  }
  while (is.peek() == '\n') // Get to the start of the next line
     is.ignore(1, ' ');
  std::getline(is, line);
  if (!is.good() || (sscanf(line.c_str(),"Alpha: %lf", &alpha) != 1))
    vw_throw(IOErr() << "AdjustableTsaiLensDistortion::read(): Could not read alpha\n");

  // Pack in the values we read
  size_t total_size = radial_vec.size() + tangential_vec.size() + 1;
  m_distortion.set_size(total_size);
  subvector(m_distortion,0,total_size-3) = radial_vec;
  subvector(m_distortion,total_size-3,2) = tangential_vec;
  m_distortion[m_distortion.size()-1]    = alpha;
}

void AdjustableTsaiLensDistortion::scale(double scale) {
  vw_throw(NoImplErr() << "AdjustableTsai doesn't support scaling.");
}

// PhotometrixLensDistortion 

PhotometrixLensDistortion::PhotometrixLensDistortion() {
  PhotometrixLensDistortion::init_distortion_param_names();
  m_distortion.set_size(m_distortion_param_names.size());
}

PhotometrixLensDistortion::PhotometrixLensDistortion(Vector<double> const& params):
  m_distortion(params) {
  PhotometrixLensDistortion::init_distortion_param_names();
  if (m_distortion.size() != m_distortion_param_names.size())
    vw_throw(IOErr()
             << "PhotometrixLensDistortion: Incorrect number of parameters was passed in.");
}

void PhotometrixLensDistortion::init_distortion_param_names() {
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
    vw_throw(IOErr() << "PhotometrixLensDistortion: Incorrect number of parameters was passed in.");
}

boost::shared_ptr<LensDistortion>
PhotometrixLensDistortion::copy() const {
  return boost::shared_ptr<PhotometrixLensDistortion>(new PhotometrixLensDistortion(*this));
}

Vector2
PhotometrixLensDistortion::undistorted_coordinates(const PinholeModel& cam,
                                                   Vector2 const& p) const {

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
  // so that the results are relative to 0,0 instead of relative to the principal point.
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

void PhotometrixLensDistortion::scale(double scale) {
  m_distortion *= scale;
}

// RPCLensDistortion 
// This class is not fully formed until both distortion and
// undistortion parameters are computed.
// One must always call set_undistortion_parameters()
// only after set_distortion_parameters().
RPCLensDistortion::RPCLensDistortion() {
  m_rpc_degree = 0;
}

RPCLensDistortion::RPCLensDistortion(Vector<double> const& params): m_distortion(params) {
  validate_distortion_params(params);
  m_rpc_degree = rpc_degree(params.size());
}

Vector<double>
RPCLensDistortion::distortion_parameters() const {
  return m_distortion;
}

void RPCLensDistortion::set_distortion_parameters(Vector<double> const& params) {
  validate_distortion_params(params);
  m_distortion = params;
  m_rpc_degree = rpc_degree(params.size());
}

// Form the identity transform
void RPCLensDistortion::reset(int rpc_degree) {
  if (rpc_degree <= 0)
    vw_throw(IOErr() << class_name() << ": The RPC degree must be positive.");

  m_rpc_degree = rpc_degree;
  int num_params = num_dist_params(rpc_degree);
  m_distortion.set_size(num_params);
  init_as_identity(m_distortion);
}

boost::shared_ptr<LensDistortion>
RPCLensDistortion::copy() const {
  return boost::shared_ptr<RPCLensDistortion>(new RPCLensDistortion(*this));
}

// RPC distortion model, after shifting relative to the principal point. So,
// the pixel is normalized and the output is also normalized.
vw::Vector2 RpcDistortion(vw::Vector2 const& P, Vector<double> const& distortion) {

  int rpc_deg = RPCLensDistortion::rpc_degree(distortion.size());

  // Precompute x^n and y^m values
  double x = P[0];
  double y = P[1];
  std::vector<double> powx(rpc_deg + 1), powy(rpc_deg + 1);
  double valx = 1.0, valy = 1.0;
  for (int deg = 0; deg <= rpc_deg; deg++) {
    powx[deg] = valx; valx *= x;
    powy[deg] = valy; valy *= y;
  }

  // Evaluate the RPC expression. The denominator always
  // has a 1 as the 0th-coefficient.

  // Loop four times, for output first coordinate numerator and
  // denominator, then for output second coordinate numerator and
  // denominator.

  int coeff_index = 0;
  double vals[] = {0.0, 1.0, 0.0, 1.0};
  for (int count = 0; count < 4; count++) {
    int start = 0; // starting degree for numerator
    if (count == 1 || count == 3)
      start = 1; // starting degree for denominator

    for (int deg = start; deg <= rpc_deg; deg++) { // start at 0
      for (int i = 0; i <= deg; i++) {
        // Add coeff * x^(deg-i) * y^i
        vals[count] += distortion[coeff_index] * powx[deg - i] * powy[i];
        coeff_index++;
      }
    }
  }

  if (coeff_index != (int)distortion.size())
    vw_throw(IOErr() << "Book-keeping failure in RPCLensDistortion.\n");

  return vw::Vector2(vals[0]/vals[1], vals[2]/vals[3]);
}

// Apply the distortion to a normalized pixel a function object. To be used in
// Newton-Raphson.
vw::Vector2 RPCLensDistortion::operator()(vw::Vector2 const& p) const {
  return RpcDistortion(p, m_distortion);
}

// A function object that applies the RPC distortion model.
Vector2
RPCLensDistortion::distorted_coordinates(const PinholeModel& cam,
                                          Vector2 const& p) const {
  Vector2 focal  = cam.focal_length(); // = [fu, fv]
  Vector2 offset = cam.point_offset(); // = [cu, cv]

  if (focal[0] < 1e-300 || focal[1] < 1e-300)
    return Vector2(HUGE_VAL, HUGE_VAL);

  // Normalize the pixel
  Vector2 dudv = p - offset; // Subtract the offset
  Vector2 p_0 = elem_quot(dudv, focal); // Divide by focal length

  // Apply the fisheye distortion model to the normalized pixel
  Vector2 p_norm = RpcDistortion(p_0, m_distortion);

  // Multiply by focal length and add the offset
  Vector2 p_dist = elem_prod(p_norm, focal) + offset;

  return p_dist;
}

// Use the numerical Jacobian to undistort the points.
Vector2 RPCLensDistortion::undistorted_coordinates(const PinholeModel& cam,
                                                    Vector2 const& p) const {

  Vector2 focal  = cam.focal_length(); // = [fu, fv]
  Vector2 offset = cam.point_offset(); // = [cu, cv]

  if (focal[0] < 1e-300 || focal[1] < 1e-300)
    return Vector2(HUGE_VAL, HUGE_VAL);

  Vector2 dudv = p - offset; // Subtract the offset
  Vector2 p_0 = elem_quot(dudv, focal); // Divide by focal length

  // Find the normalized undistorted pixel using Newton-Raphson and the
  // numerical Jacobian. The step size for numerical differentiation 
  // better not be too small. A tolerance of 1e-9 is likely adequate.
  vw::Vector2 guess = p_0;
  double step = 1e-6, tol = 1e-9;
  vw::math::NewtonRaphson nr(*this);
  Vector2 U = nr.solve(guess, p_0, step, tol);

  // Multiply by focal length and add the offset
  Vector2 p_undist = elem_prod(U, focal) + offset;

  return p_undist;
}

void RPCLensDistortion::scale(double scale) {
  // Throw an error. This is not a well-defined operation.
  vw_throw(NoImplErr() << "RPCLensDistortion::scale() is not implemented.");
  m_distortion *= scale;
}

void RPCLensDistortion::validate_distortion_params(Vector<double> const& params) {
  int num_params = params.size();
  int deg = rpc_degree(num_params);
  if (num_dist_params(deg) != num_params || deg <= 0 || std::isnan(deg))
    vw_throw(IOErr() << class_name() << ": Incorrect number of parameters was passed in: "
              << num_params << ".");
}

// Make RPC coefficients so that the RPC transform is the identity.
// The vector params must already have the right size.
void RPCLensDistortion::init_as_identity(Vector<double> & params) {
  RPCLensDistortion::validate_distortion_params(params);

  params.set_all(0);
  Vector<double> num_x, den_x, num_y, den_y;
  unpack_params(params, num_x, den_x, num_y, den_y);

  // Initialize the transform (x, y) -> (x, y), which is
  // ((0 + 1*x + 0*y)/(1 + 0*x + 0*y), (0 + 0*x + 1*y)/(1 + 0*x + 0*y))
  // hence set num_x and num_y accordingly. As always, we do not
  // store the 1 values in the denominator.
  num_x[1] = 1; num_y[2] = 1;
  pack_params(params, num_x, den_x, num_y, den_y);
}

namespace {

  // A little function to append zeros to a Vector.
  inline void append_zeros_to_vector(Vector<double> & vec, int num) {

    int len = vec.size();

    Vector<double> out_vec;
    out_vec.set_size(len + num);
    for (int it = 0; it < len; it++)
      out_vec[it] = vec[it];

    for (int it = len; it < len + num; it++) {
      out_vec[it] = 0.0;
    }
    vec = out_vec;
  }

}

// Given the RPC coefficients corresponding to the four polynomials,
// increase the degree of each polynomial by 1 and set the new
// coefficients to 0.
void RPCLensDistortion::increment_degree(Vector<double> & params) {
  RPCLensDistortion::validate_distortion_params(params);

  Vector<double> num_x, den_x, num_y, den_y;
  unpack_params(params, num_x, den_x, num_y, den_y);

  int r = rpc_degree(params.size());

  // The next monomials to add will be
  // x^(r+1), x^r*y, ..., x*y^r, y^(r+1)
  // and there are r + 2 of them.
  // Set their coefficients to zero.
  int num = r + 2;

  append_zeros_to_vector(num_x, num);
  append_zeros_to_vector(den_x, num);
  append_zeros_to_vector(num_y, num);
  append_zeros_to_vector(den_y, num);

  pack_params(params, num_x, den_x, num_y, den_y);

}

void RPCLensDistortion::unpack_params(Vector<double> const& params,
                                       Vector<double> & num_x, Vector<double> & den_x,
                                       Vector<double> & num_y, Vector<double> & den_y) {
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
                                    Vector<double> const& num_x,
                                    Vector<double> const& den_x,
                                    Vector<double> const& num_y,
                                    Vector<double> const& den_y) {

  int num_len = num_x.size();
  int den_len = den_x.size();

  if (num_len != den_len + 1 || num_len != int(num_y.size()) || den_len != int(den_y.size()))
    vw::vw_throw(vw::IOErr() << "Book-keeping failure in RPCLensDistortion.\n");

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
      vw_throw(IOErr() << "Found an unexpected empty vector.\n");
    Vector<double> tmp_params = params;
    params = subvector(tmp_params, 1, len - 1);
  }
}

void RPCLensDistortion::write(std::ostream & os) const {

  // TODO: Add domain of validity for the distorted and undistorted pixels. This is needed
  // because otherwise the RPC model can return wrong results which confuses
  // bundle adjustment.
  write_param("rpc_degree", os, m_rpc_degree);

  // Write the distortion
  Vector<double> num_x, den_x, num_y, den_y;
  unpack_params(m_distortion, num_x, den_x, num_y, den_y);
  prepend_1(den_x);
  prepend_1(den_y);
  write_param_vec("distortion_num_x", os, num_x);
  write_param_vec("distortion_den_x", os, den_x);
  write_param_vec("distortion_num_y", os, num_y);
  write_param_vec("distortion_den_y", os, den_y);
}

void RPCLensDistortion::read(std::istream & is) {

  read_param("rpc_degree", is, m_rpc_degree);

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
}
