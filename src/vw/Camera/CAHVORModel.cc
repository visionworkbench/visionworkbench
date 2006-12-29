// __BEGIN_LICENSE__
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
//
#include <vw/Camera/CAHVORModel.h>

using namespace std;

// Overloaded constructor - this one reads in the file name
// where the CAVHOR camera model is saved.
vw::camera::CAHVORModel::CAHVORModel(std::string const& filename) {
    
  FILE *cahvorFP = fopen(filename.c_str(), "r");
  if (cahvorFP == 0)
    vw_throw( IOErr() << "CAHVORModel: Could not open file." );

  char line[4096];    
  printf("Reading CAHVOR file: %s\n", filename.c_str());
    
  // Scan through comments 
  fgets(line, sizeof(line), cahvorFP);
  while(line[0] == '#') 
    fgets(line, sizeof(line), cahvorFP);
  
  if (sscanf(line,"C = %lf %lf %lf", &C(0), &C(1), &C(2)) != 3) {
    vw_throw( IOErr() << "CAHVORModel: Could not read C vector\n" );
    fclose(cahvorFP);
  }
  
  fgets(line, sizeof(line), cahvorFP);
  if (sscanf(line,"A = %lf %lf %lf", &A(0),&A(1), &A(2)) != 3) {
    vw_throw( IOErr() << "CAHVORModel: Could not read A vector\n" );
    fclose(cahvorFP);
  }
  
  fgets(line, sizeof(line), cahvorFP);
  if (sscanf(line,"H = %lf %lf %lf", &H(0), &H(1), &H(2)) != 3) {
    vw_throw( IOErr() << "CAHVORModel: Could not read H vector\n" );
    fclose(cahvorFP);
  }
  
  fgets(line, sizeof(line), cahvorFP);
  if (sscanf(line,"V = %lf %lf %lf", &V(0), &V(1), &V(2)) != 3) {
    vw_throw( IOErr() << "CAHVORModel: Could not read V vector\n" );
    fclose(cahvorFP);
  }
  
  fgets(line, sizeof(line), cahvorFP);
  if (sscanf(line,"O = %lf %lf %lf", &O(0), &O(1), &O(2)) != 3) {
    vw_throw( IOErr() << "CAHVORModel: Could not read O vector\n" );
    fclose(cahvorFP);
  }
  
  fgets(line, sizeof(line), cahvorFP);
  if (sscanf(line,"R = %lf %lf %lf", &R(0), &R(1), &R(2)) != 3) {
    vw_throw( IOErr() << "CAHVORModel: Could not read R vector\n" );
    fclose(cahvorFP);
  }
  
  fclose(cahvorFP);
  
  // For debugging:
  //     cout << "CameraModel_CAHV: C vector:" << C << "\n";
  //     cout << "CameraModel_CAHV: A vector:" << A << "\n";
  //     cout << "CameraModel_CAHV: H vector:" << H << "\n";
  //     cout << "CameraModel_CAHV: V vector:" << V << "\n";
  //     cout << "CameraModel_CAHV: O vector:" << O << "\n";
  //     cout << "CameraModel_CAHV: R vector:" << R << "\n";
}

// Set iteration and convergence constants
#define VW_CAHVOR_MAXITER  20     // maximum number of iterations allowed 
#define VW_CAHVOR_CONV   1.0e-6   // covergence tolerance - check adequacy for application 
  
// CAHVOR pixel_to_vector with partial_derivative output
vw::Vector3 vw::camera::CAHVORModel::pixel_to_vector(vw::Vector2 const& pix, vw::Matrix<double> &partial_derivatives) const {

  // Note, vec is actually the output vector
  // Based on JPL_CMOD_CAHVOR_2D_TO_3D
  Vector3 vec;
  int i, j;
  double omega, omega_2, tau, mu, u, u_2, du, k1, k3, k5, poly, deriv;
  
  Vector3 f, g, rr, pp, wo, lambda; 
  
  Matrix<double, 3, 3> m33, n33;
  double sgn, magv, magi;
  Vector3 t, w, v3, u3;
  Matrix<double, 3, 3> irrt;
  
  double dudt;
  Vector3 drpdx, drpdy;
  Matrix<double, 3, 3> dldr;
  
  Vector3 drdx, drdy;
  Matrix<double, 3, 3> drpdr, drpdri;
  
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
  
  k1 = 1 + R(0);		//  1 + rho0 
  k3 = R(1) * tau;		//  rho1*tau  
  k5 = R(2) * tau*tau;	//  rho2*tau^2  
  
  mu = R(0) + k3 + k5;
  u = 1.0 - mu;	// initial approximation for iterations 
  
  for (i=0; i<VW_CAHVOR_MAXITER; i++) {
    u_2 = u*u;
    poly  =  ((k5*u_2  +  k3)*u_2 + k1)*u - 1;
    deriv = (5*k5*u_2 + 3*k3)*u_2 + k1;
    
    if (deriv <= 0) {
      printf("CAHVORModel.pixel_to_vector(): Distortion is too negative\n");
      break;
    } else {
      du = poly/deriv;
      u -= du;
      if (fabs(du) < VW_CAHVOR_CONV)
        break;
    } 
  } 
  
  if (i >= VW_CAHVOR_MAXITER) {
    printf("CAHVORModel.pixel_to_vector(): Too many iterations (%d)\n", i);
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
vw::Vector3 vw::camera::CAHVORModel::pixel_to_vector(vw::Vector2 const& pix) const {

  // vec is actually the output vector, vec is just the C vector of the camera
  // Based on JPL_CMOD_CAHVOR_2D_TO_3D
  // cout << "Starting CAHVOR PixelToVector w/o Partial... " << endl;
    
  Vector3 vec;
  int i;
  double omega, omega_2, tau, mu, u, u_2, du, k1, k3, k5, poly, deriv;
    
  Vector3 f, g, rr, pp, wo, lambda; 
  Matrix<double, 3, 3> m33, n33;
  
  double sgn, magv, magi;
  Vector3 t, w, v3, u3;
  Matrix<double, 3, 3> irrt;
    
  Vector3 drpdx, drpdy;
  Matrix<double, 3, 3> dldr;
  
  Vector3 drdx, drdy;
  Matrix<double, 3, 3> drpdr, drpdri;

  //  The projection point is merely the C of the camera model. 
  //  vec = C; not used here, as output is vec

  // Calculate the projection ray assuming normal vector directions, 
  // neglecting distortion.                                          
  f = pix.y() * A; // pos2[1] in JPL code is y, and pos2[0] is x...
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
    //cout << "CAHVOR PixelToVector changes rr sign" << endl;
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
    
  k1 = 1 + R(0);		//  1 + rho0 
  k3 = R(1) * tau;		//  rho1*tau  
  k5 = R(2) * tau*tau;	//  rho2*tau^2  
 
  mu = R(0) + k3 + k5;
  u = 1.0 - mu;	// initial approximation for iterations 
    
  for (i=0; i<VW_CAHVOR_MAXITER; i++) {

    u_2 = u*u;
    poly  =  ((k5*u_2  +  k3)*u_2 + k1)*u - 1;
    deriv = (5*k5*u_2 + 3*k3)*u_2 + k1;

    if (deriv <= 0) {
      printf("CAHVORModel.pixel_to_vector(): Distortion is too negative\n");
      break;
    } else {
      du = poly/deriv;
      u -= du;
      if (fabs(du) < VW_CAHVOR_CONV)
        break;
    } 
  } 
  

  if (i >= VW_CAHVOR_MAXITER) {
    printf("CAHVORModel.pixel_to_vector(): Too many iterations (%d)\n", i);
  } 

  mu = 1 - u;
  pp = mu * lambda;
  vec = rr - pp;
  magv = norm_2(vec);
  vec = 1.0/magv * vec;
  return vec;
}


// vector_to_pixel with partial_derivatives 
vw::Vector2 vw::camera::CAHVORModel::point_to_pixel(vw::Vector3 const& point, vw::Matrix<double> &partial_derivatives) const {

  Vector3 vec = point - C;
  // Based on JPL 3D to 2D POINT (not the 3D to 2D function alone).
  Vector2 pix;
  double alpha, beta, gamma, xh, yh;
  double omega, omega_2, tau, mu;
  Vector3 pp_c, wo, lambda;
  Matrix<double, 3, 3> dldp, dppdp, m33, n33;
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
vw::Vector2 vw::camera::CAHVORModel::point_to_pixel(vw::Vector3 const& point) const {

  Vector3 vec = point - C;
  Vector2 pix;
  // cout << "Starting CAHVOR VectorToPixel w/o partial... " << endl;
  double alpha, beta, gamma, xh, yh;
  double omega, omega_2, tau, mu;
  Vector3 pp_c, wo, lambda;
    
  // Calculate necessary quantities 
  omega = dot_prod(vec,O);
  omega_2 = omega * omega;
  wo = omega * O;
  lambda = vec - wo;
  tau = dot_prod(lambda, lambda) / omega_2;
  mu = R(0) + (R(1) * tau) + (R(2) * tau * tau);
  pp_c = mu * lambda;
  pp_c = vec + pp_c;
    
  // Calculate alpha, beta, gamma, which are 
  // dotted with a, h, v, respectively  
  alpha = dot_prod(pp_c, A);
  beta = dot_prod(pp_c, H);
  gamma = dot_prod(pp_c, V);
    
  // Calculate the projection 
  pix.x() = xh = beta / alpha;
  pix.y() = yh = gamma /alpha; 
  return pix;
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
vw::camera::CAHVModel vw::camera::linearize_camera(vw::camera::CAHVORModel &camera_model, 
                                                   unsigned cahvor_image_width, unsigned cahvor_image_height,
                                                   unsigned cahv_image_width, unsigned cahv_image_height) {
    
  unsigned int minfov = 1; // set to 0 if you do not want to minimize to a common field of view
  unsigned int i;

  Vector3 a2, h2, v2;
  Vector3 rt, dn, p3, u3, vec1, vec2;
  double sn, x, hmin, hmax, vmin, vmax;
  double hs, vs, hc, vc, theta; // these used to be output parameters
  Matrix<double> pts(12,2);
  Matrix<double> hpts(6,2);
  Matrix<double> vpts(6,2);

  Vector<unsigned int> idims(2);
  idims(0) = cahvor_image_width; 
  idims(1) = cahvor_image_height;

  Vector<unsigned int> odims(2);
  odims(0) = cahv_image_width;
  odims(1) = cahv_image_height;

  // Record the landmark 2D coordinates around the perimeter of the image 
  pts(0,0) = hpts(0,0) = 0;
  pts(0,1) = hpts(0,1) = 0;
  
  pts(1,0) = hpts(1,0) = 0;
  pts(1,1) = hpts(1,1) = (idims(1)-1)/2.0;
  
  pts(2,0) = hpts(2,0) = 0;
  pts(2,1) = hpts(2,1) = idims(1)-1;

  pts(3,0) = hpts(3,0) = idims(0)-1;
  pts(3,1) = hpts(3,1) = 0;

  pts(4,0) = hpts(4,0) = idims(0)-1;
  pts(4,1) = hpts(4,1) = (idims(1)-1)/2.0;

  pts(5,0) = hpts(5,0) = idims(0)-1;
  pts(5,1) = hpts(5,1) = idims(1)-1;

  pts(6,0) = vpts(0,0) = 0;
  pts(6,1) = vpts(0,1) = 0;

  pts(7,0) = vpts(1,0) = (idims(0)-1)/2.0;
  pts(7,1) = vpts(1,1) = 0;

  pts(8,0) = vpts(2,0) = idims(0)-1;
  pts(8,1) = vpts(2,1) = 0;

  pts(9,0) = vpts(3,0) = 0;
  pts(9,1) = vpts(3,1) = idims(1)-1;

  pts(10,0) = vpts(4,0) = (idims(0)-1)/2.0;
  pts(10,1) = vpts(4,1) = idims(1)-1;

  pts(11,0) = vpts(5,0) = idims(0)-1;
  pts(11,1) = vpts(5,1) = idims(1)-1;


  // Choose a camera axis in the middle of the perimeter 
    
  // Think VW autofills with zeros
  //a2.fill(0);
  //h2.fill(0);
  //v2.fill(0);

  Vector2 LoopPixel;
    
  for (i=0; i<12; i++) {
      
    LoopPixel.x() = pts(i,0);
    LoopPixel.y() = pts(i,1);
      
    u3 = camera_model.pixel_to_vector(LoopPixel);
    a2 = u3 + a2;
      
  } // end get a2 loop
    
  a2 = a2/norm_2(a2);
  // Compute the original right and down vectors 
  dn = cross_prod(camera_model.A, camera_model.H); // down vector 
  rt = cross_prod(dn, camera_model.A); // right vector 

  dn = dn/norm_2(dn);
  rt = rt/norm_2(rt);

  // Adjust the right and down vectors to be orthogonal to new axis 
  rt = cross_prod(dn, a2);
  dn = cross_prod(a2, rt);

  dn = dn/norm_2(dn);
  rt = rt/norm_2(rt);

  // Find horizontal and vertical fields of view 

  // Horizontal
  hmin =  1;
  hmax = -1;
  Vector2 HLoopPix;

  for (i=0; i<6; i++) {
    HLoopPix.x() = hpts(i,0);
    HLoopPix.y() = hpts(i,1);

    u3 = camera_model.pixel_to_vector(HLoopPix);

    x = dot_prod(dn, u3);
    vec1 = x * dn;
    vec2 = u3 - vec1;
    vec2 = vec2/norm_2(vec2);

    vec1 = cross_prod(a2, vec2);
    sn = norm_2(vec1);

    if (hmin > sn)
      hmin = sn;
    if (hmax < sn)
      hmax = sn;
  } // end HLoop

    // Vertical
  vmin =  1;
  vmax = -1;
  Vector2 VLoopPix;

  for (i=0; i<6; i++) {

    VLoopPix.x() = vpts(i,0);
    VLoopPix.y() = vpts(i,1);
      
    u3 = camera_model.pixel_to_vector(VLoopPix);

    x = dot_prod(rt, u3);
    vec1 = x * rt;
    vec2 = u3 - vec1;
    vec2 = vec2/norm_2(vec2);

    vec1 = cross_prod(a2, vec2);
    sn = norm_2(vec1);

    if (vmin > sn)
      vmin = sn;
    if (vmax < sn)
      vmax = sn;
  } // end VLoop
  
    // Compute the all-encompassing scale factors 
  sn = (minfov ? hmin : hmax);
  x = odims(0) / 2.0;
  hs = sqrt((x*x)/(sn*sn) - x*x);
  sn = (minfov ? vmin : vmax);
  x = odims(1) / 2.0;
  vs = sqrt((x*x)/(sn*sn) - x*x);
  
  // Assign idealized image centers and coordinate angles 
  hc = (odims[0] - 1) / 2.0;
  vc = (odims[1] - 1) / 2.0;
  theta = -M_PI / 2.0;
  
  /* Construct H and V */
  vec1 = hs * rt;
  vec2 = hc * a2;
  
  h2 = vec1 + vec2;

  vec1 = vs * dn;
  vec2 = vc * a2;

  v2 = vec1 + vec2;

  CAHVModel output_camera;
  output_camera.C = camera_model.C;
  output_camera.A = a2;
  output_camera.H = h2;
  output_camera.V = v2;

  // For debugging:
  //       cout << "CAHVOR linearize_camera completed: Converted to CAHV: " << endl;
  //       cout << "CAHVORModel: Linearized C vector:" << output_camera.C << "\n";
  //       cout << "CAHVORModel: Linearized A vector:" << output_camera.A << "\n";
  //       cout << "CAHVORModel: Linearized H vector:" << output_camera.H << "\n";
  //       cout << "CAHVORModel: Linearized V vector:" << output_camera.V << "\n";
  
  return output_camera;
} 
  
// Note: the second derivatives with respect to motion of the point are the 
// same as the second derivatives with respect to motion of the camera, and 
// the first derivatives are opposite.
void vw::camera::CAHVORModel::get_point_derivatives( Vector3 const& P, double& u, double& v,
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
  // hess_u(i,j) = ( dot_prod( grad_mu[j]*(I[i]-O[i]*O) + grad_mu[i]*(I[j]-O[j]*O) + lambda*hess_mu(i,j), HuA )
  //               - dot_prod( grad_u[j]*(I[i]+mu*(I[i]-O[i]*O)+lambda*grad_mu[i], A )
  //               - dot_prod( grad_u[i]*(I[j]+mu*(I[j]-O[j]*O)+lambda*grad_mu[j], A ) ) / denom;
  //             = ( grad_mu[j]*(HuA[i]-O[i]*dot_prod(O,HuA))
  //               + grad_mu[i]*(HuA[j]-O[j]*dot_prod(O,HuA))
  //               - grad_mu[j]*((1+mu)*A[i]-mu*O[i]*dot_prod(O,A)+dot_prod(lambda,A)*grad_mu[i])
  //               - grad_mu[i]*((1+mu)*A[j]-mu*O[j]*dot_prod(O,A)+dot_prod(lambda,A)*grad_mu[j])
  //               + dot_prod(lambda,HuA)*hess_mu ) / denom;
  Matrix3x3 tmp_u = outer_prod( grad_mu, HuA - dot_prod(O,HuA)*O )
    + outer_prod( grad_u, mu*dot_prod(O,A)*O - dot_prod(lambda,A)*grad_mu - (1+mu)*A );
  hess_u = ( tmp_u + transpose(tmp_u) + dot_prod(lambda,HuA)*hess_mu ) / denom;
  Matrix3x3 tmp_v = outer_prod( grad_mu, VvA - dot_prod(O,VvA)*O )
    + outer_prod( grad_v, mu*dot_prod(O,A)*O - dot_prod(lambda,A)*grad_mu - (1+mu)*A );
  hess_v = ( tmp_v + transpose(tmp_v) + dot_prod(lambda,VvA)*hess_mu ) / denom;
}
