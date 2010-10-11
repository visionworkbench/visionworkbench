// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
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
