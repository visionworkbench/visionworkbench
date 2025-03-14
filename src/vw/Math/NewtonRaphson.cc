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

// See the .h file for description of this function.
vw::Vector<double> NewtonRaphson::numericalJacobian(vw::Vector2 const& P,
                                                    double step) {

  // The Jacobian has 4 elements
  vw::Vector<double> jacobian(4);
  
  // First column
  vw::Vector2 JX = (m_func(P + vw::Vector2(step, 0)) - 
                    m_func(P - vw::Vector2(step, 0))) / (2*step);
  // Second column
  vw::Vector2 JY = (m_func(P + vw::Vector2(0, step)) - 
                    m_func(P - vw::Vector2(0, step))) / (2*step);
  
  // Put in the jacobian matrix
  jacobian[0] = JX[0];
  jacobian[1] = JY[0];
  jacobian[2] = JX[1];
  jacobian[3] = JY[1];  
  
  return jacobian;
}

// See the .h file for description of this function.
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
      J = numericalJacobian(X, step_size);
    } catch(...) {
      // Something went wrong. Cannot continue. Return most recent result.
      return bestX;
    }
    
    // Find the determinant
    double det = J[0]*J[3] - J[1]*J[2];
    if (std::abs(det) < 1e-6 || std::isnan(det)) 
      return bestX; // bad determinant. Cannot continue. Return most recent result.
     
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

}} // end namespace vw::math
