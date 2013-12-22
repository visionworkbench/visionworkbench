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


#include <vw/Core/Log.h>
#include <vw/Camera/CAHVORModel.h>
#include <vw/Camera/CAHVModel.h>
#include <fstream>
#include <boost/foreach.hpp>

using namespace vw;
using namespace camera;

// Overloaded constructor - this one reads in the file name
// where the CAVHOR camera model is saved.
CAHVORModel::CAHVORModel(std::string const& filename) {

  try {
    std::ifstream input(filename.c_str(), std::ifstream::in);
    input.exceptions(std::ifstream::failbit | std::ifstream::badbit);

    vw_out(InfoMessage, "camera") << "Reading CAHVOR file: "
                                  << filename << ".\n";

    char r1,r2;

    while (true) {
      input.ignore(1024, 'C');
      input >> r1;
      if (r1 == '=')
        break;
    }

    input >> C(0) >> C(1) >> C(2);

    input >> r1 >> r2;
    if (r1 != 'A' || r2 != '=')
      vw_throw( IOErr() << "CAHVORModel: Could not read A vector\n" );
    input >> A(0) >> A(1) >> A(2);

    input >> r1 >> r2;
    if (r1 != 'H' || r2 != '=')
      vw_throw( IOErr() << "CAHVORModel: Could not read H vector\n" );
    input >> H(0) >> H(1) >> H(2);

    input >> r1 >> r2;
    if (r1 != 'V' || r2 != '=')
      vw_throw( IOErr() << "CAHVORModel: Could not read V vector\n" );
    input >> V(0) >> V(1) >> V(2);

    input >> r1 >> r2;
    if (r1 != 'O' || r2 != '=')
      vw_throw( IOErr() << "CAHVORModel: Could not read O vector\n" );
    input >> O(0) >> O(1) >> O(2);

    input >> r1 >> r2;
    if (r1 != 'R' || r2 != '=')
      vw_throw( IOErr() << "CAHVORModel: Could not read R vector\n" );
    input >> R(0) >> R(1) >> R(2);

  } catch (const std::ifstream::failure& e) {
    vw_throw( IOErr() << "CAHVORModel: Could not read file: " << filename << " (" << e.what() << ")" );
  }
}

/// Initialize the CAHVOR vectors directly in the native CAHVOR format.
CAHVORModel::CAHVORModel(Vector3 const& C_vec, Vector3 const& A_vec,
                         Vector3 const& H_vec, Vector3 const& V_vec,
                         Vector3 const& O_vec, Vector3 const& R_vec) :
  C(C_vec), A(A_vec), H(H_vec), V(V_vec), O(O_vec), R(R_vec) {}


CAHVORModel::~CAHVORModel() {}

std::string
CAHVORModel::type() const { return "CAHVOR"; }

// Write CAHVOR model to file.
void CAHVORModel::write(std::string const& filename) {

  try {
    std::ofstream output(filename.c_str(), std::ofstream::out);
    output.exceptions(std::ofstream::failbit | std::ofstream::badbit);
    output.precision(20);

    vw_out(InfoMessage, "camera") << "Writing CAHVOR file: " << filename << "\n";

    output << "C = " << C(0) << " " << C(1) << " " << C(2) << "\n"
           << "A = " << A(0) << " " << A(1) << " " << A(2) << "\n"
           << "H = " << H(0) << " " << H(1) << " " << H(2) << "\n"
           << "V = " << V(0) << " " << V(1) << " " << V(2) << "\n"
           << "O = " << O(0) << " " << O(1) << " " << O(2) << "\n"
           << "R = " << R(0) << " " << R(1) << " " << R(2) << "\n";

  } catch (const std::ofstream::failure& e) {
    vw_throw( IOErr() << "CAHVORModel: Could not write file: " << filename << " (" << e.what() << ")" );
  }
}

// Set iteration and convergence constants
#define VW_CAHVOR_MAXITER  20     // maximum number of iterations allowed
#define VW_CAHVOR_CONV   1.0e-6   // covergence tolerance - check adequacy for application

// CAHVOR pixel_to_vector with partial_derivative output
Vector3 CAHVORModel::pixel_to_vector(Vector2 const& pix,
                                     Matrix<double> &partial_derivatives) const {

  // Note, vec is actually the output vector
  // Based on JPL_CMOD_CAHVOR_2D_TO_3D
  Vector3 vec;
  int i, j;
  double omega, omega_2, tau, mu, u, u_2, du, k1, k3, k5, poly, deriv;

  Vector3 f, g, rr, pp, wo, lambda;

  Matrix3x3 m33, n33;
  double sgn, magv, magi;
  Vector3 t, w, v3, u3;
  Matrix3x3 irrt;

  double dudt;
  Vector3 drpdx, drpdy;
  Matrix3x3 dldr;

  Vector3 drdx, drdy;
  Matrix3x3 drpdr, drpdri;

  //  The projection point is merely the C of the camera model.
  //  vec = C; This isn't needed, only vec is returned for our
  //  usage.

  // Calculate the projection ray assuming normal vector directions,
  // neglecting distortion.

  f = pix.y() * A; // pos2[1] in JPL code is our y, and pos2[0] is x...
  f = V - f;

  g = pix.x() * A;
  g = H - g;
  rr = cross_prod(f,g);
  magi = 1.0/norm_2(rr);
  rr = magi * rr;

  // Check and optionally correct for vector directions.
  sgn = 1;
  t = cross_prod(V,H);

  if (dot_prod(t, A) < 0) {
    rr = -1.0 * rr;
    sgn = -1;
  }

  // Optionally compute partial for non-linear part of the model
  irrt.set_identity();

  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      irrt(i,j) -= (rr(i) * rr(j));
      // used to be A.as_ref() in vxl
      t = cross_prod(f, A);
      w = irrt * t;
      drpdx(0) = partial_derivatives(0,0) = -sgn * w(0) * magi;
      drpdx(1) = partial_derivatives(1,0) = -sgn * w(1) * magi;
      drpdx(2) = partial_derivatives(2,0) = -sgn * w(2) * magi;

      // used to be A.as_ref() in vxl
      t = cross_prod(g, A);
      w = irrt * t;
      drpdy(0) = partial_derivatives(0,1) = sgn * w(0) * magi;
      drpdy(1) = partial_derivatives(1,1) = sgn * w(1) * magi;
      drpdy(2) = partial_derivatives(2,1) = sgn * w(2) * magi;
    }
  }

  // Remove the radial lens distortion.  Preliminary values of
  // omega, lambda, and tau are computed from the rr vector
  // including distortion, in order to obtain the coefficients of
  // the equation k5*u^5 + k3*u^3 + k1*u = 1, which is solved for u
  // by means of Newton's method.  This value is used to compute the
  // corrected rr.
  omega = dot_prod(rr, O);
  omega_2 = omega * omega;
  wo = omega * O;
  lambda = rr - wo;
  tau = dot_prod(lambda, lambda) / omega_2;

  k1 = 1 + R(0);                //  1 + rho0
  k3 = R(1) * tau;              //  rho1*tau
  k5 = R(2) * tau*tau;  //  rho2*tau^2

  mu = R(0) + k3 + k5;
  u = 1.0 - mu; // initial approximation for iterations

  for (i=0; i<VW_CAHVOR_MAXITER; i++) {
    u_2 = u*u;
    poly  =  ((k5*u_2  +  k3)*u_2 + k1)*u - 1;
    deriv = (5*k5*u_2 + 3*k3)*u_2 + k1;

    if (deriv <= 0) {
      vw_out(InfoMessage, "camera") << "CAHVORModel.pixel_to_vector(): Distortion is too negative\n";
      break;
    } else {
      du = poly/deriv;
      u -= du;
      if (fabs(du) < VW_CAHVOR_CONV)
        break;
    }
  }

  if (i >= VW_CAHVOR_MAXITER) {
    vw_out(InfoMessage, "camera") << "CAHVORModel.pixel_to_vector(): Too many iterations (" << i << ")\n";
  }

  mu = 1 - u;
  pp = mu * lambda;
  vec = rr - pp;
  magv = norm_2(vec);
  vec = 1.0/magv * vec;

  // Note:  If partial derivatives are to be computed, corrected values
  // of omega, lambda, tau, and mu must be computed.

  // Recompute omega, lambda, tau, and mu
  omega = dot_prod(vec, O);
  omega_2 = omega * omega;
  wo = omega * O;
  lambda = vec - wo;
  tau = dot_prod(lambda, lambda) / omega_2;
  mu = R(0) + R(1)*tau + R(2)*tau*tau;

  // Compute the partial derivatives for distortion
  dldr.set_identity();

  Matrix<double, 3, 1> Ocol;
  Matrix<double, 1, 3> Orow;

  select_col(Ocol, 0) = O;
  select_row(Orow, 0) = O;

  m33 = Ocol*Orow; // outer product, hopefully
  dldr = dldr - m33;

  dudt = R(1) + (2 * R(2) * tau);


  v3 = transpose(transpose(lambda) * dldr);
  v3 = (2/omega_2) * v3;
  u3 = (2 * tau / omega) * O;
  v3 = v3 - u3;

  Matrix<double, 3, 1> lambdacol;
  Matrix<double, 1, 3> v3row;

  select_col(lambdacol,0) = lambda;
  select_row(v3row,0) = v3;

  m33 = lambdacol * v3row; // outer product, hopefully
  m33 = dudt * m33;
  n33 = mu * dldr;
  drpdr = m33 + n33;

  m33.set_identity();
  drpdr = m33 + drpdr;
  drpdr = magv * drpdr;

  // Apply these partials to get the final result
  drpdri = inverse(drpdr);
  drdx = drpdri * drpdx;
  drdy = drpdri * drpdy;

  partial_derivatives(0,0) = drdx(0);
  partial_derivatives(1,0) = drdx(1);
  partial_derivatives(2,0) = drdx(2);
  partial_derivatives(0,1) = drdy(0);
  partial_derivatives(1,1) = drdy(1);
  partial_derivatives(2,1) = drdy(2);

  return vec;
}


// pixel_to_vector (no returned partial matrix)
Vector3 CAHVORModel::pixel_to_vector(Vector2 const& pix) const {
  // Based on JPL_CMOD_CAHVOR_2D_TO_3D

  // Calculate the projection ray assuming normal vector directions,
  // neglecting distortion.
  Vector3 rr = normalize(cross_prod(V - pix.y() * A,
                                    H - pix.x() * A));

  // Check and optionally correct for vector directions.
  if (dot_prod(cross_prod(V,H), A) < 0)
    rr = -1.0 * rr;

  // Remove the radial lens distortion.  Preliminary values of
  // omega, lambda, and tau are computed from the rr vector
  // including distortion, in order to obtain the coefficients of
  // the equation k5*u^5 + k3*u^3 + k1*u = 1, which is solved for u
  // by means of Newton's method.  This value is used to compute the
  // corrected rr.
  double omega = dot_prod(rr, O);
  Vector3 lambda = rr - omega * O;
  const double tau = dot_prod(lambda, lambda) / (omega*omega);

  const double k1 = 1 + R(0);        //  1 + rho0
  const double k3 = R(1) * tau;      //  rho1*tau
  const double k5 = R(2) * tau*tau;  //  rho2*tau^2

  double u = 1.0 - (R(0) + k3 + k5); // initial approximation for iterations

  double du, poly, deriv;
  for (int32 i=0;; i++) {
    if (i >= VW_CAHVOR_MAXITER) {
      vw_out(InfoMessage, "camera") << "CAHVORModel.pixel_to_vector(): Too many iterations (" << i << ")\n";
      break;
    }

    double u_2 = u*u;
    poly  =  ((k5*u_2  +  k3)*u_2 + k1)*u - 1;
    deriv = (5*k5*u_2 + 3*k3)*u_2 + k1;

    if (deriv <= 0) {
      vw_out(InfoMessage, "camera") << "CAHVORModel.pixel_to_vector(): Distortion is too negative\n";
      break;
    } else {
      du = poly/deriv;
      u -= du;
      if (fabs(du) < VW_CAHVOR_CONV)
        break;
    }
  }

  return normalize(rr - (1 - u)*lambda);
}

Vector3 CAHVORModel::camera_center( Vector2 const& pix ) const { return C; }

// vector_to_pixel with partial_derivatives
Vector2 CAHVORModel::point_to_pixel(Vector3 const& point,
                                    Matrix<double> &partial_derivatives) const {

  Vector3 vec = point - C;
  // Based on JPL 3D to 2D POINT (not the 3D to 2D function alone).
  Vector2 pix;
  double alpha, beta, gamma, xh, yh;
  double omega, omega_2, tau, mu;
  Vector3 pp_c, wo, lambda;
  Matrix3x3 dldp, dppdp, m33, n33;
  Vector3 dxhdpp, dyhdpp, v3, u3;
  double dudt;

  // Calculate necessary quantities
  omega = dot_prod(vec,O);
  omega_2 = omega * omega;
  wo = omega * O;
  lambda = vec - wo;
  tau = dot_prod(lambda, lambda) / omega_2;
  mu = R(0) + (R(1) * tau) + (R(2) * tau * tau);
  pp_c = mu * lambda;
  pp_c = vec + pp_c;

  // Calculate alpha, beta, gamma,
  // dotted with a, h, v, respectively
  alpha = dot_prod(pp_c, A);
  beta = dot_prod(pp_c, H);
  gamma = dot_prod(pp_c, V);

  // Calculate the projection
  pix.x() = xh = beta / alpha;
  pix.y() = yh = gamma /alpha;

  // Calculate the approximate partial derivatives
  v3 = xh * A;
  v3 = H - v3;
  dxhdpp = 1/alpha * v3;
  v3 = yh * A;
  v3 = V - v3;
  dyhdpp = 1/alpha * v3;

  // Complete the calculations for accuracy
  dldp.set_identity();

  Matrix<double, 3, 1> Ocol;
  Matrix<double, 1, 3> Orow;

  select_col(Ocol, 0) = O;
  select_row(Orow, 0) = O;

  m33 = Ocol*Orow; // outer product, hopefully

  dldp = dldp - m33;
  dudt = R(1) + (2 * R(2) * tau);
  lambda = dldp * v3;
  v3 = 2/omega_2 * v3;
  u3 = (2 * tau/omega) * O;
  v3 = v3 - u3;

  Matrix<double, 3, 1> lambdacol;
  Matrix<double, 1, 3> v3row;

  select_col(lambdacol,0) = lambda;
  select_row(v3row,0) = v3;

  m33 = lambdacol * v3row; // outer product, hopefully

  m33 = dudt * m33;
  n33 = mu * dldp;
  dppdp = m33 + n33;
  m33.set_identity();
  dppdp = dppdp + m33;
  select_row(partial_derivatives,0) = transpose(transpose(dxhdpp) * dppdp);
  select_row(partial_derivatives,1) = transpose(transpose(dyhdpp) * dppdp);
  return pix;
}

// vector_to_pixel without partial_derivatives
Vector2 CAHVORModel::point_to_pixel(Vector3 const& point) const {

  // Convert to directional
  Vector3 vec = point - C;

  // Calculate necessary quantities
  double omega = dot_prod(vec,O);
  Vector3 lambda = vec - omega * O;
  double tau = dot_prod(lambda,lambda) / (omega*omega);
  double mu = R(0) + (R(1) * tau) + (R(2) * tau * tau);
  Vector3 pp_c = vec + mu * lambda;

  double alpha = dot_prod(pp_c, A);

  // Calculate the projection
  return Vector2( dot_prod(pp_c,H) / alpha,
                  dot_prod(pp_c,V) / alpha );
}

// linearize_camera
//
// Takes CAHVOR camera --> CAHV camera
// Requires knowledge of size of image from CAHVOR camera
//
// This function warps a camera model so that it is purely
// linear. The parameter C will not change. The parameters O
// (identical to A) and R (all terms zero) will not be output. Note
// that image warping will be necessary in order to use the new
// models.
CAHVModel camera::linearize_camera( CAHVORModel const& camera_model,
                                    Vector2i const& cahvor_image_size,
                                    Vector2i const& cahv_image_size ) {

  CAHVModel output_camera;
  output_camera.C = camera_model.C;

  static const bool minfov = true; // set to 0 if you do not want to
                                   // minimize to a common field of
                                   // view

  // Record the landmark 2D coordinates around the perimeter of the image
  Vector2 hpts[6], vpts[6];
  hpts[0] = Vector2();
  hpts[1] = Vector2(0,(cahvor_image_size[1]-1)/2.0);
  hpts[2] = Vector2(0,cahvor_image_size[1]-1);
  hpts[3] = Vector2(cahvor_image_size[0]-1,0);
  hpts[4] = Vector2(cahvor_image_size[0]-1,(cahvor_image_size[1]-1)/2.0);
  hpts[5] = cahvor_image_size - Vector2(1,1);
  vpts[0] = Vector2();
  vpts[1] = Vector2((cahvor_image_size[0]-1)/2.0,0);
  vpts[2] = Vector2(cahvor_image_size[0]-1,0);
  vpts[3] = Vector2(0,cahvor_image_size[1]-1);
  vpts[4] = Vector2((cahvor_image_size[0]-1)/2,cahvor_image_size[1]-1);
  vpts[5] = cahvor_image_size - Vector2(1,1);

  BOOST_FOREACH( Vector2 const& local, vpts ) {
    output_camera.A += camera_model.pixel_to_vector(local);
  }
  BOOST_FOREACH( Vector2 const& local, hpts ) {
    output_camera.A += camera_model.pixel_to_vector(local);
  }
  output_camera.A = normalize(output_camera.A);

  // Compute the original right and down vectors
  Vector3 dn = cross_prod(camera_model.A, camera_model.H); // down vector
  Vector3 rt = normalize(cross_prod(dn, camera_model.A)); // right vector
  dn = normalize(dn);

  // Adjust the right and down vectors to be orthogonal to new axis
  rt = cross_prod(dn, output_camera.A);
  dn = normalize(cross_prod(output_camera.A, rt));
  rt = normalize(rt);

  // Find horizontal and vertical fields of view

  // Horizontal
  double hmin = 1, hmax = -1;
  BOOST_FOREACH( Vector2 const& loop, hpts ) {
    Vector3 u3 = camera_model.pixel_to_vector(loop);
    double sn = norm_2(cross_prod(output_camera.A,
                                  normalize(u3 - dot_prod(dn, u3) * dn)));
    if (hmin > sn) hmin = sn;
    if (hmax < sn) hmax = sn;
  }

  // Vertical
  double vmin = 1, vmax = -1;
  BOOST_FOREACH( Vector2 const& loop, vpts ) {
    Vector3 u3 = camera_model.pixel_to_vector(loop);
    double sn = norm_2(cross_prod(output_camera.A,
                                  normalize(u3 - dot_prod(rt, u3) * rt) ) );
    if (vmin > sn) vmin = sn;
    if (vmax < sn) vmax = sn;
  }

  // Compute the all-encompassing scale factors
  Vector2 scale_factors;
  Vector2 image_center = (cahv_image_size-Vector2(1,1))/2.0;
  Vector2 image_center_2 = elem_prod(image_center,image_center);
  if ( minfov ) {
    // Use min value
    scale_factors =
      sqrt( elem_quot(image_center_2,
                      Vector2(hmin*hmin,vmin*vmin)) - image_center_2 );
  } else {
    // Use max value
    scale_factors =
      sqrt( elem_quot(image_center_2,
                      Vector2(hmax*hmax,vmax*vmax)) - image_center_2 );
  }

  // Assign idealized image centers and coordinate angles
  output_camera.H = scale_factors[0] * rt + image_center[0] * output_camera.A;
  output_camera.V = scale_factors[1] * dn + image_center[1] * output_camera.A;

  return output_camera;
}

CAHVModel camera::linearize_camera( CAHVORModel const& camera_model,
                                    int32 ix, int32 iy, int32 ox, int32 oy ) {
  return vw::camera::linearize_camera( camera_model,
                                       Vector2i( ix, iy), Vector2i(ox, oy) );
}

// Note: the second derivatives with respect to motion of the point are the
// same as the second derivatives with respect to motion of the camera, and
// the first derivatives are opposite.
void CAHVORModel::get_point_derivatives( Vector3 const& P, double& u, double& v,
                                         Vector3& grad_u, Vector3& grad_v,
                                         Matrix3x3& hess_u, Matrix3x3& hess_v ) const {
  // Compute the image plane position
  double xi = dot_prod(P-C,O);
  Vector3 lambda = P-C-xi*O;
  double tau = dot_prod(lambda,lambda)/(xi*xi);
  double mu = R[0]+tau*(R[1]+tau*R[2]);
  Vector3 PP = P + mu*lambda;
  Vector3 PPC = PP - C;
  double denom = dot_prod(PPC,A);
  u = dot_prod(PPC,H)/denom;
  v = dot_prod(PPC,V)/denom;

  // Compute the gradients
  Vector3 grad_tau = 2*(lambda-tau*xi*O)/(xi*xi);
  Vector3 grad_mu = 2*(R[1]+2*R[2]*tau)*(lambda/xi-tau*O)/xi;
  Vector3 HuA = H - u*A;
  Vector3 VvA = V - v*A;
  grad_u = ( (1+mu)*HuA - mu*dot_prod(O,HuA)*O + dot_prod(lambda,HuA)*grad_mu ) / denom;
  grad_v = ( (1+mu)*VvA - mu*dot_prod(O,VvA)*O + dot_prod(lambda,VvA)*grad_mu ) / denom;

  // Compute the Hessians
  Matrix3x3 I; I.set_identity();
  Matrix3x3 hess_tau = 2/(xi*xi) * ( I + (3*tau-1)*outer_prod(O,O)
                                     - 2*(outer_prod(O,lambda)+outer_prod(lambda,O))/xi );
  Matrix3x3 hess_mu = (R[1]+2*R[2]*tau)*hess_tau + 2*R[2]*outer_prod(grad_tau,grad_tau);
  Matrix3x3 tmp_u = outer_prod( grad_mu, HuA - dot_prod(O,HuA)*O )
    + outer_prod( grad_u, mu*dot_prod(O,A)*O - dot_prod(lambda,A)*grad_mu - (1+mu)*A );
  hess_u = ( tmp_u + transpose(tmp_u) + dot_prod(lambda,HuA)*hess_mu ) / denom;
  Matrix3x3 tmp_v = outer_prod( grad_mu, VvA - dot_prod(O,VvA)*O )
    + outer_prod( grad_v, mu*dot_prod(O,A)*O - dot_prod(lambda,A)*grad_mu - (1+mu)*A );
  hess_v = ( tmp_v + transpose(tmp_v) + dot_prod(lambda,VvA)*hess_mu ) / denom;
}
