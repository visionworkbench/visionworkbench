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


/// \file IntegralInterestOperator.cc
///
/// Classes that detect interest by using integral images
///
#include <vw/InterestPoint/IntegralInterestOperator.h>

namespace vw {
namespace ip {

  const double OBALoGInterestOperator::SCALE_LOG_SIGMA[10] = {1.2, 1.5, 1.875, 2.34375, 2.9296875, 3.662109, 4.577636, 5.722046, 7.152557, 8.9406967};
  const int OBALoGInterestOperator::SCALE_BOX_WIDTH[10][6] =
    {{11,7,5,3,3,1},     // S 0
     {13,9,7,5,3,1},     // S 1
     {17,11,9,5,5,3},    // S 2
     {21,13,11,7,7,3},   // S 3
     {27,17,13,9,9,5},   // S 4
     {33,21,17,11,11,5}, // S 5
     {41,27,21,13,13,7}, // S 6
     {51,33,27,17,17,7}, // S 7
     {3,3,3,3,3,3},      // S 8
     {3,3,3,3,3,3}};     // S 9
  const int OBALoGInterestOperator::SCALE_BOX_HEIGHT[10][6] =
    {{11,5,7,3,1,3},     // S 0
     {13,7,9,1,3,5},     // S 1
     {17,9,11,5,3,5},    // S 2
     {21,11,13,3,7,7},   // S 3
     {27,13,17,5,9,9},   // S 4
     {33,17,21,5,11,11}, // S 5
     {41,21,27,7,13,13}, // S 6
     {51,27,33,7,17,17}, // S 7
     {3,3,3,3,3,3},      // S 8
     {3,3,3,3,3,3}};     // S 9
  const double OBALoGInterestOperator::SCALE_BOX_WEIGHT[10][6] =
    {{0.0006179999909,0.008283999749,0.008283999749,-0.03556599841,-0.05576400086,-0.05576400086}, // S 0
     {0.0002659999882,0.004889999982,0.004889999982,-0.0221420005,-0.0488499999,-0.0221420005}, // S 1
     {0.0002659424644,0.003207367443,0.003207367443,-0.001847865087,-0.022190776,-0.022190776}, // S 2
     {0.000218920823,0.002005755751,0.002005755751,-0.01499466244,-0.0008247125278,-0.01499466244}, // S 3
     {0.0001030362087,0.00137811946,0.00137811946,-0.009606413672,0.002226325845,-0.009606413672}, // S 4
     {7.198028372e-05,0.0008259398524,0.0008259398524,-0.006072390568,-1.177469443e-06,-0.006072390568}, // S 5
     {3.803569417e-05,0.0005251058097,0.0005251058097,-0.003885340628,0.0002823550512,-0.003885340628}, // S 6
     {2.892054582e-05,0.0003287680581,0.0003287680581,-0.002447736656,-0.0002717150913,-0.002447736656}, // S 7
     {1,-1,1,-1,1,-1}, // S 8
     {1,-1,1,-1,1,-1}};// S 9

}
}
