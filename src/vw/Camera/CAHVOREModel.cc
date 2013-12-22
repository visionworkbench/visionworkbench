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
#include <vw/Camera/CAHVModel.h>
#include <vw/Camera/CAHVOREModel.h>
#include <fstream>
#include <boost/foreach.hpp>

using namespace vw;
using namespace vw::camera;

bool CAHVOREModel::check_line( std::istream& istream, char letter ) {
  char r1,r2;
  istream >> r1 >> r2;
  if ( r1 != letter || r2 != '=' )
    return false;
  return true;
}

// Overloaded constructor - this one reads in the file name
// where the CAHVORE camera model is saved.
CAHVOREModel::CAHVOREModel(std::string const& filename) {

  try {
    std::ifstream input(filename.c_str(), std::ifstream::in);
    input.exceptions(std::ifstream::failbit | std::ifstream::badbit);

    vw_out(InfoMessage, "camera") << "Reading CAHVORE file: "
                                  << filename << ".\n";

    char r1;

    while (true) {
      input.ignore(1024, 'C');
      input >> r1;
      if (r1 == '=')
        break;
    }
    input >> C(0) >> C(1) >> C(2);

    if ( !check_line( input, 'A') )
      vw_throw( IOErr() << "CAHVOREModel: Could not read A vector\n" );
    input >> A(0) >> A(1) >> A(2);

    if ( !check_line( input, 'H') )
      vw_throw( IOErr() << "CAHVOREModel: Could not read H vector\n" );
    input >> H(0) >> H(1) >> H(2);

    if ( !check_line( input, 'V') )
      vw_throw( IOErr() << "CAHVOREModel: Could not read V vector\n" );
    input >> V(0) >> V(1) >> V(2);

    if ( !check_line( input, 'O') )
      vw_throw( IOErr() << "CAHVOREModel: Could not read O vector\n" );
    input >> O(0) >> O(1) >> O(2);

    if ( !check_line( input, 'R') )
      vw_throw( IOErr() << "CAHVOREModel: Could not read R vector\n" );
    input >> R(0) >> R(1) >> R(2);

    if ( !check_line( input, 'E') )
      vw_throw( IOErr() << "CAHVOREModel: could not read R vector\n" );
    input >> E(0) >> E(1) >> E(2);

    if ( !check_line( input, 'T') )
      vw_throw( IOErr() << "CAHVOREModel: could not read T element\n" );
    int T;
    input >> T;

    if ( !check_line( input, 'P') )
      vw_throw( IOErr() << "CAHVOREModel: could not read P element\n" );
    input >> P;

    switch( T ) {
    case 1: P = 1.0; break;
    case 2: P = 0.0; break;
    case 3:
      if ( P < 0 || P > 1 ) vw_throw( ArgumentErr() << "Invalid P value: "
                                      << P << "\n" );
      break;
    default: vw_throw( ArgumentErr() << "Unknown CAHVORE type: " << T << "\n" );
    }

  } catch (const std::ifstream::failure& e) {
    vw_throw( IOErr() << "CAHVOREModel: Could not read file: " << filename << " (" << e.what() << ")" );
  }
}

CAHVOREModel::CAHVOREModel(Vector3 const& C_vec, Vector3 const& A_vec,
                           Vector3 const& H_vec, Vector3 const& V_vec,
                           Vector3 const& O_vec, Vector3 const& R_vec,
                           Vector3 const& E_Vec) :
  C(C_vec), A(A_vec), H(H_vec), V(V_vec), O(O_vec), R(R_vec), E(E_Vec), P(1.0) {}

CAHVOREModel::CAHVOREModel(Vector3 const& C_vec, Vector3 const& A_vec,
                           Vector3 const& H_vec, Vector3 const& V_vec,
                           Vector3 const& O_vec, Vector3 const& R_vec,
                           Vector3 const& E_Vec, int T, double P_v) :
  C(C_vec), A(A_vec), H(H_vec), V(V_vec), O(O_vec), R(R_vec), E(E_Vec), P(P_v) {
  switch( T ) {
  case 1: P = 1.0; break;
  case 2: P = 0.0; break;
  case 3:
    if ( P < 0 || P > 1 ) vw_throw( ArgumentErr() << "Invalid P value: " << P_v << "\n" );
    break;
  default: vw_throw( ArgumentErr() << "Unknown CAHVORE type: " << T << "\n" );
  }
}

CAHVOREModel::CAHVOREModel(Vector3 const& C_vec, Vector3 const& A_vec,
                           Vector3 const& H_vec, Vector3 const& V_vec,
                           Vector3 const& O_vec, Vector3 const& R_vec,
                           Vector3 const& E_Vec, double P_v) :
  C(C_vec), A(A_vec), H(H_vec), V(V_vec), O(O_vec), R(R_vec), E(E_Vec), P(P_v) {
  if ( P < 0 || P > 1 ) vw_throw( ArgumentErr() << "Invalid P value: " << P_v << "\n" );
}

CAHVOREModel::~CAHVOREModel() {}

std::string CAHVOREModel::type() const { return "CAHVORE"; }

// Write CAHVOR model to file.
void CAHVOREModel::write(std::string const& filename) {
  try {
    std::ofstream output(filename.c_str(), std::ofstream::out);
    output.exceptions(std::ofstream::failbit | std::ofstream::badbit);
    output.precision(20);

    vw_out(InfoMessage, "camera") << "Writing CAHVORE file: " << filename << "\n";

    output << "C = " << C(0) << " " << C(1) << " " << C(2) << "\n"
           << "A = " << A(0) << " " << A(1) << " " << A(2) << "\n"
           << "H = " << H(0) << " " << H(1) << " " << H(2) << "\n"
           << "V = " << V(0) << " " << V(1) << " " << V(2) << "\n"
           << "O = " << O(0) << " " << O(1) << " " << O(2) << "\n"
           << "R = " << R(0) << " " << R(1) << " " << R(2) << "\n"
           << "E = " << E(0) << " " << E(1) << " " << E(2) << "\n";
    if ( P == 0 ) {
      output << "T = 2\n";
    } else if ( P == 1 ) {
      output << "T = 1\n";
    } else {
      output << "T = 3\n";
    }
    output << "P = " << P << "\n";
  } catch ( const std::ofstream::failure& e ) {
    vw_throw( IOErr() << "CAHVOREModel: Could not write file: " << filename << " (" << e.what() << ")" );
  }
}

vw::Vector3 CAHVOREModel::pixel_to_vector(vw::Vector2 const& pix) const {
  // Based on JPL's cmod_cahvore_2d_to_3d
  Vector3 result;

  // Calculate initial terms
  Vector3 w3 = cross_prod(V - pix[1]*A,
                          H - pix[0]*A);
  Vector3 rp = (1/dot_prod(A,cross_prod(V,H))) * w3;

  double zetap = dot_prod(rp,O);
  Vector3 lambdap3 = rp - zetap*O;
  double lambdap = norm_2(lambdap3);
  double chip = lambdap / zetap;

  if (chip < 1e-8) {
    // Approximations for small angles
    result = O;
  } else {
    // Full calculations

    // Calculate chi using Newton's Method
    double chi = chip;
    double dchi = 1;
    for (int32 n = 0;; ++n) {
      // Checking exit conditions
      if (n > 100)
        vw_throw( PixelToRayErr() << "CAHVOREModel: Did not converge.\n" );
      if (fabs(dchi) < 1e-8)
        break;

      // Compute terms from the current value of chi
      double chi2 = chi * chi;
      double chi3 = chi * chi2;
      double chi4 = chi * chi3;
      double chi5 = chi * chi4;

      // Update chi
      double deriv = (1 + R[0]) + 3*R[1]*chi2 + 5*R[2]*chi4;
      dchi = ((1 + R[0])*chi + R[1]*chi3 + R[2]*chi5 - chip) / deriv;
      chi -= dchi;
    }

    // Compute the incoming ray's angle
    double linchi, theta;
    linchi = P * chi;
    if (P < -1e-15)
      theta = asin(linchi) / P;
    else if (P > 1e-15)
      theta = atan(linchi) / P;
    else
      theta = chi;

    double theta2 = theta*theta;
    double mu = E[2];
    mu = E[1] + theta2*mu;
    mu = E[0] + theta2*mu;

    result = sin(theta)*normalize(lambdap3) + cos(theta)*O;
  }

  return result;
}

Vector3 CAHVOREModel::camera_center(Vector2 const& pix ) const { return C; }

Vector2 CAHVOREModel::point_to_pixel(vw::Vector3 const& point) const {
  // Base on JPL's cmod_cahvore_3d_to_2d_general

  // Calculate initial terms
  Vector3 p_c = point - C;
  double zeta = dot_prod(p_c, O);
  Vector3 lambda3 = p_c - zeta * O;
  double lambda = norm_2(lambda3);

  // Calculate theta using Newton's Method
  double theta = atan2(lambda, zeta);
  double dtheta = 1;
  for (int32 n = 0;;++n) {

    // Checking exit conditions
    if (n > 100)
      vw_throw( PointToPixelErr() << "CAHVOREModel: Did not converge.\n" );
    if (fabs(dtheta) < 1e-8)
      break;

    // Compute terms from the current value of theta
    double costh = cos(theta);
    double sinth = sin(theta);
    double theta2 = theta * theta;
    double theta3 = theta * theta2;
    double theta4 = theta * theta3;
    double upsilon = zeta*costh + lambda*sinth
      - (1     - costh) * (E[0] +  E[1]*theta2 +   E[2]*theta4)
      - (theta - sinth) * (      2*E[1]*theta  + 4*E[2]*theta3);

    // Update theta
    dtheta = (
              zeta*sinth - lambda*costh
              - (theta - sinth) * (E[0] + E[1]*theta2 + E[2]*theta4)
              ) / upsilon;
    theta -= dtheta;
  }

  // Check the value of theta
  if ((theta * fabs(P)) > M_PI/2)
    vw_throw( PointToPixelErr() << "CAHVOREModel: Theta out of bounds.\n" );

  // Approximations for small theta
  Vector3 rp;
  if (theta < 1e-8) {
    rp = p_c;
  } else {
    // Full calculations
    double linth, chi;

    linth = P * theta;
    if (P < -1e-15)
      chi = sin(linth) / P;
    else if (P > 1e-15)
      chi = tan(linth) / P;
    else
      chi = theta;

    double chi2 = chi*chi;
    double mu = R[2];
    mu = R[1] + chi2 * mu;
    mu = R[0] + chi2 * mu;
    rp = (lambda / chi) * O + (1+mu)*lambda3;
  }

  // Calculate the projection
  double alpha  = dot_prod(rp, A);
  return Vector2( dot_prod(rp,H) / alpha,
                  dot_prod(rp,V) / alpha );
}

CAHVModel camera::linearize_camera( CAHVOREModel const& camera_model,
                                    Vector2i const& cahvore_image_size,
                                    Vector2i const& cahv_image_size ) {
  // Limit to field of view
  const static double limfov = M_PI * 3/4; // 135 degrees field of view.
  const bool minfov          = true;       // Yes, minimize the common field of view

  CAHVModel cahv_model;
  cahv_model.C = camera_model.C;

  // Record the landmark 2D coordinates around the perimeter of the image
  Vector2 hpts[6], vpts[6];
  hpts[0] = Vector2();
  hpts[1] = Vector2(0,(cahvore_image_size[1]-1)/2);
  hpts[2] = Vector2(0,cahvore_image_size[1]-1);
  hpts[3] = Vector2(cahvore_image_size[0]-1,0);
  hpts[4] = Vector2(cahvore_image_size[0]-1,(cahvore_image_size[1]-1)/2);
  hpts[5] = cahvore_image_size - Vector2(1,1);
  vpts[0] = Vector2();
  vpts[1] = Vector2((cahvore_image_size[0]-1)/2,0);
  vpts[2] = Vector2(cahvore_image_size[0]-1,0);
  vpts[3] = Vector2(0,cahvore_image_size[1]-1);
  vpts[4] = Vector2((cahvore_image_size[0]-1)/2,cahvore_image_size[1]-1);
  vpts[5] = cahvore_image_size - Vector2(1,1);

  // Choose a camera axis in the middle of the image
  Vector2 p2 = (cahvore_image_size-Vector2i(1,1))/2.0;
  cahv_model.A = camera_model.pixel_to_vector( p2 );

  // Compute the original right and down vectors
  Vector3 dn = cross_prod(camera_model.A, camera_model.H);
  Vector3 rt = normalize(cross_prod(dn,   camera_model.A));
  dn = normalize( dn );

  // Adjust the right and down vectors to be orthogonal to new axis
  rt = cross_prod(dn, cahv_model.A);
  dn = normalize(cross_prod(cahv_model.A, rt));
  rt = normalize( rt );

  // Find horizontal and vertical fields of view
  double hmin = 1, hmax = -1;
  BOOST_FOREACH( Vector2 const& loop, hpts ) {
    const Vector3 u3 = camera_model.pixel_to_vector(loop);
    double cs = dot_prod(cahv_model.A, normalize(u3 - dot_prod(dn, u3) * dn));
    if (hmin > cs) hmin = cs;
    if (hmax < cs) hmax = cs;
  }
  double vmin = 1, vmax = -1;
  BOOST_FOREACH( Vector2 const& loop, vpts ) {
    const Vector3 u3 = camera_model.pixel_to_vector(loop);
    double cs = dot_prod(cahv_model.A,normalize(u3 - dot_prod(rt, u3)*rt));
    if (vmin > cs) vmin = cs;
    if (vmax < cs) vmax = cs;
  }

  // Compute the all-encompassing scale factors
  Vector2 cosines;
  if ( minfov ) {
    // use max
    cosines = Vector2(hmax,vmax);
  } else {
    // use min
    cosines = Vector2(hmin,vmin);
  }
  if ( acos(cosines[0]) > limfov ) cosines[0] = cos(limfov);
  if ( acos(cosines[1]) > limfov ) cosines[1] = cos(limfov);
  Vector2 scalars =
    elem_quot(elem_prod(cahv_image_size/2.0,cosines),
              sqrt(Vector2(1,1) - elem_prod(cosines,cosines)));

  // Assign idealized image centers and coordinate angles
  Vector2 centers = (cahv_image_size - Vector2(1,1))/2.0;

  // Construct H and V
  cahv_model.H = scalars[0] * rt + centers[0] * cahv_model.A;
  cahv_model.V = scalars[1] * dn + centers[1] * cahv_model.A;

  return cahv_model;
}
