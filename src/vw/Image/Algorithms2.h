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


/// \file Algorithms2.h
///
/// Some very simple code which needs to be tightened. It is used only
/// in very few places and is not very generic, hence not in the main
/// Algorithms.h.
///
#ifndef __VW_IMAGE_ALGORITHMS2_H__
#define __VW_IMAGE_ALGORITHMS2_H__

#include <vw/Image/Algorithms.h>
#include <vw/Image/Manipulation.h>

namespace vw {

// Round a double value to int, and clamp to the min and max of this type of int.
template<class IntType>
IntType round_and_clamp(double val){

  // Round. Be careful to work with doubles before clamping, to not overflow.
  val = round(val);

  // Clamp. Note that this won't work right unless IntType is an int type.
  if (val < double(std::numeric_limits<IntType>::min()) )
    val = double(std::numeric_limits<IntType>::min());
  if (val > double(std::numeric_limits<IntType>::max()) )
    val = double(std::numeric_limits<IntType>::max());

  return static_cast<IntType>(val);
}

// Apply this to a an image pixel. If more than one channel, use the first. 
template<class IntType, class InputType>
class RoundAndClamp: public ReturnFixedType<IntType> {
public:
  IntType operator()(InputType const& v) const {
    
    // I could not figure out in reasonable time how to select a channel
    // from a pixel which can be single or compound. Hence use the
    // functionality for selecting a channel from an image.
    ImageView<InputType> A(1, 1);
    A(0, 0) = v;
    
    // First cast to double, as the values of A could be out of range
    // for the OutputPixelT data type.
    ImageView<double> B = select_channel(A, 0);
    
    return round_and_clamp<IntType>(B(0,0));
  }
};
  
} // namespace vw

#endif // __VW_IMAGE_ALGORITHMS2_H__
