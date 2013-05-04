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


#ifndef GPU_Interpolate_H
#define GPU_Interpolate_H


#include <vw/Image/Interpolation.h>


namespace vw {

// Default Interpolation Used By GPU Functions
typedef BilinearInterpolation DefaultInterpolation;

namespace GPU {
  /*
  enum EdgeExtensionType {
    INTERPOLATION_BILINEAR,
    INTERPOLATION_BICUBIC,
    INTERPOLATION_NEARESTPIXEL
  };
  */
   template <class InterpT>
 struct TraitsForInterpT {
   static const char* ShaderString() {
     static const char* string = "";
     return string;
   }
 };

   template <>
 struct TraitsForInterpT<NearestPixelInterpolation> {
   static const char* ShaderString() {
     static const char* string = "Interpolation/interpolation-nearest-pixel";
     return string;
   }
   static const int quality = 1;
 };


   template <>
 struct TraitsForInterpT<BilinearInterpolation> {
   static const char* ShaderString() {
     static const char* string = "Interpolation/interpolation-bilinear";
     return string;
   }
   static const int quality = 2;
 };

   template <>
 struct TraitsForInterpT<BicubicInterpolation> {
   static const char* ShaderString() {
     static const char* string = "Interpolation/interpolation-bicubic";
     return string;
   }
   static const int quality = 3;
 };



} // namespaces GPU

}


#endif
