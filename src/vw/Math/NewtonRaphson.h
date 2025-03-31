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

#ifndef __VW_MATH_NEWTON_RAPHSON_H__
#define __VW_MATH_NEWTON_RAPHSON_H__

#include <vw/Core/Log.h>
#include <vw/Math/Matrix.h>
#include <vw/Math/Vector.h>
#include <vw/Math/LinearAlgebra.h>

namespace vw {
namespace math {

// Newton-Raphson method to find X in 2D such that func(X) = outY. Use either the
// supplied Jacobian or the internal numerical Jacobian. The function and the
// Jacobian to evaluate must be function objects implementing operator() with
// signature of NewtFuncType and NewtJacType respectively.

// Great care is needed for the step size and tol. See also the
// numericalJacobian function for an explanation of the step size. The user must
// check how small func(X) - outY is.

using NewtFuncType = std::function<vw::Vector2(vw::Vector2 const&)>;
using NewtJacType = std::function<vw::Vector<double>(vw::Vector2 const& P, double step)>;

class NewtonRaphson {
    
public:
    
    NewtonRaphson(NewtFuncType func, NewtJacType jac = {});
    
    vw::Vector2 solve(vw::Vector2 const& guessX, vw::Vector2 const& outY,
                      double step_size, double tol);
  
private:
    NewtFuncType m_func;
    NewtJacType m_jac;
};

// Find the Jacobian of a function at a given point using numerical
// differentiation. A good value for the step is 1e-6 or larger. It should be
// notably larger than the tolerance for finding the function value or the
// magnitude of any noise.
struct numericalJacobian {
  numericalJacobian(NewtFuncType func);    
  vw::Vector<double> operator()(vw::Vector2 const& P, double step);
  NewtFuncType m_func;       
};

}} // namespace vw::math

#endif // __VW_MATH_NEWTON_RAPHSON_H__
