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

// Implementation of the Newton-Raphson method. See the .h file for details.

#include <vw/Math/NewtonRaphson.h>

namespace vw {
namespace math {

numericalJacobian::numericalJacobian(FuncType func): m_func(func) {}
    
vw::Vector<double> numericalJacobian::operator()(vw::Vector2 const& P, double step) {
  
  // The Jacobian has 4 elements
  vw::Vector<double> jacobian(4);
  
  // First column
  vw::Vector2 dx(step, 0);
  vw::Vector2 JX = (m_func(P + dx) - m_func(P - dx)) / (2*step);
 
  // Second column
  vw::Vector2 dy(0, step);
  vw::Vector2 JY = (m_func(P + dy) - m_func(P - dy)) / (2*step);
  
  // Form the Jacobian
  jacobian[0] = JX[0];
  jacobian[1] = JY[0];
  jacobian[2] = JX[1];
  jacobian[3] = JY[1];  
  
  return jacobian;
}

NewtonRaphson::NewtonRaphson(FuncType func, JacType jac): m_func(func) {
  
  // If a Jacobian function was passed in, use it. Otherwise use the numerical jacobian.
  if (jac)
    m_jac = jac;
   else
    m_jac = numericalJacobian(func);
}

vw::Vector2 NewtonRaphson::solve(vw::Vector2 const& guessX,  // initial guess
                                 vw::Vector2 const& outY,    // desired output
                                 double step_size, double tol) {

  // Start with initial guess
  vw::Vector2 X = guessX;

  // If an iteration fails, will fall back to best found so far.
  // This is unlikely, but cannot be ruled out.
  vw::Vector2 bestX = X;
  double best_err = std::numeric_limits<double>::max();
  
  // Normally 20 iterations are enough. If not, likely it will not converge.
  int count = 1, maxTries = 20;
  while (count < maxTries) {
    
    vw::Vector2 F;
    try {
      F = m_func(X) - outY;
      
      if (norm_2(F) < best_err || F != F) {
        best_err = norm_2(F);
        bestX = X;
      }
      
    } catch (...) {
      // Something went wrong. Cannot continue. Return most recent result.
      return bestX;
    }
    
    // Compute the Jacobian
    vw::Vector<double> J(4);
    try {
      J = m_jac(X, step_size);
    } catch (...) {
      // Something went wrong. Cannot continue. Return most recent result.
      return bestX;
    }
    
    // Find the determinant
    double det = J[0]*J[3] - J[1]*J[2];
    if (std::abs(det) < 1e-6 || std::isnan(det)) 
      return bestX; // Bad determinant. Cannot continue. Return most recent result.
     
    // Update X
    vw::Vector2 DX;
    DX[0] = (J[3]*F[0] - J[1]*F[1]) / det;
    DX[1] = (J[0]*F[1] - J[2]*F[0]) / det;
    X -= DX;
    
    // If DX is small enough, we are done
    if (norm_2(DX) < tol)
      return X;
    
    count++;
  }
  
  // Fallback result, if the loop did not finish
  return bestX;
}

// Newton-Raphson method with analytical Jacobian. 
// TODO(oalexan1): Integrate with the numerical jacobian version.
void newtonRaphson(double dx, double dy, double &ux, double &uy,
                    Vector<double> const& extraArgs,
                    const double tolerance,
                    FuncType2 func,
                    JacType2 jac) {

  const int maxTries = 20;

  double x, y, jacobian[4];

  // Initial guess for the root
  x = dx;
  y = dy;

  vw::Vector2 fval = func(vw::Vector2(x, y), extraArgs);
  double fx = fval[0];
  double fy = fval[1];

  for (int count = 1;
        ((fabs(fx) + fabs(fy)) > tolerance) && (count < maxTries); count++) {
  
    fval = func(vw::Vector2(x, y), extraArgs);
    fx = fval[0];
    fy = fval[1];

    fx = dx - fx;
    fy = dy - fy;

    double step = 1e-6; // not used, part of the interface
    vw::Vector<double> jacobian = jac(vw::Vector2(x, y), step, extraArgs);

    // Jxx * Jyy - Jxy * Jyx
    double determinant =
        jacobian[0] * jacobian[3] - jacobian[1] * jacobian[2];
    if (fabs(determinant) < 1e-6) {
      ux = x;
      uy = y;
      // Near-zero determinant. Cannot continue. Return most recent result.
      return;
    }

    x = x + (jacobian[3] * fx - jacobian[1] * fy) / determinant;
    y = y + (jacobian[0] * fy - jacobian[2] * fx) / determinant;
  }

  if ((fabs(fx) + fabs(fy)) <= tolerance) {
    // The method converged to a root.
    ux = x;
    uy = y;

    return;
  }
}

}} // end namespace vw::math
