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


/// \file NewtonRaphson.h 
///
/// Class for inverting a 2D to 2D function via the Newton-Raphson method. 
///  - Must create an object that implements operator() that will be the
///    function.

///  - Need not use templates or virtual functions, as the function signature
///    is fixed. 

/// -  The function must be differentiable. The Jacobian is computed
///    numerically. 
/// -  The function and the Jacobian can throw exceptions on failure.

#ifndef __VW_MATH_NEWTON_RAPHSON_H__
#define __VW_MATH_NEWTON_RAPHSON_H__

// Vision Workbench
#include <vw/Core/Log.h>
#include <vw/Math/Matrix.h>
#include <vw/Math/Vector.h>
#include <vw/Math/LinearAlgebra.h>

namespace vw {
namespace math {

// TODO(oalexan1): Temporary
using FuncType = std::function<vw::Vector2(vw::Vector2)>;
using JacType = std::function<vw::Vector<double>(vw::Vector2 const& P, double step)>;

// TODO(oalexan1): Temporary
// Find the Jacobian of a function at a given point using numerical
// differentiation. A good value for the step is 1e-6 or larger.  
// It should be notably larger than the tolerance for finding the function
// value or the magnitude of any noise.
struct numericalJacobian {
  
  numericalJacobian(FuncType func);    
  vw::Vector<double> operator()(vw::Vector2 const& P, double step);
  FuncType m_func;       
};

// Newton-Raphson method to find the root of a function. Use either the provided
// Jacobian or the internal numerical Jacobian. The function and the jacobian
// to evaluate must be function objects implementing operator() with expected
// signatures.
class NewtonRaphson {
public:
    
    NewtonRaphson(FuncType func, JacType jac = {});
    
    // Newton-Raphson method to find X such that func(X) = outY.
    // Great care is needed for the step size and tol. See also the 
    // numericalJacobian function for an explanation of the step size.
    // The user must check how small func(X) - outY is.
    vw::Vector2 solve(vw::Vector2 const& guessX, vw::Vector2 const& outY,
                      double step_size, double tol);
  
private:
    FuncType m_func;
    JacType m_jac;
};

// TODO(oalexan1): Temporary
using FuncType2 = std::function<vw::Vector2(vw::Vector2, vw::Vector<double> const&)>;
using JacType2 = std::function<vw::Vector<double>(vw::Vector2 const& P, double step,
                                                  vw::Vector<double> const& extra)>;

// Newton-Raphson method with analytical Jacobian. 
// TODO(oalexan1): Integrate with the numerical jacobian version.
void newtonRaphson(double dx, double dy, double &ux, double &uy,
                    Vector<double> const& extraArgs,
                    const double tolerance,
                    FuncType2 func,
                    JacType2 jac);

}} // namespace vw::math

#endif // __VW_MATH_NEWTON_RAPHSON_H__
