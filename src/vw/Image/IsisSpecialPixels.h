// __BEGIN_LICENSE__
//  Copyright (c) 2006-2026, United States Government as represented by the
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

/// \file IsisSpecialPixels.h
///
/// This file encodes the ISIS special-pixel sentinel values, which are a fixed
/// IEEE bit-pattern and integer-sentinel file-format convention. 

/// This does not include or link the ISIS library. It is low-level pixel-type
/// logic.

#ifndef __VW_IMAGE_ISISSPECIALPIXELS_H__
#define __VW_IMAGE_ISISSPECIALPIXELS_H__

#include <vw/Image/PixelTypeInfo.h>
#include <vw/Image/PixelMask.h>
#include <vw/Image/ImageViewRef.h>

#include <cstdint>
#include <cstring>
#include <cfloat>

namespace vw {

// The ISIS special-pixel sentinel values, transcribed from ISIS SpecialPixel.h.
// The float and double values are IEEE bit patterns, so they are formed by a
// byte-wise reinterpret (std::memcpy, which is not constexpr in C++17, hence
// the inline functions and inline const values below).
namespace isis {

  // Reinterpret a 32-bit pattern as a float, byte-exact.
  inline float floatFromBits(uint32_t bits) {
    float f = 0.0f;
    std::memcpy(&f, &bits, sizeof(f));
    return f;
  }

  // Reinterpret a 64-bit pattern as a double, byte-exact. The hex values used
  // below are the little-endian in-memory representation of the ISIS double
  // sentinels (all supported platforms are little-endian).
  inline double doubleFromBits(uint64_t bits) {
    double d = 0.0;
    std::memcpy(&d, &bits, sizeof(d));
    return d;
  }

  // 8-byte (double) sentinels. In ISIS these are NaN bit patterns.
  inline const double VALID_MIN8      = doubleFromBits(0xFFFFFFFAFFEFFFFFull);
  inline const double ValidMinimum    = VALID_MIN8;
  inline const double NULL8           = doubleFromBits(0xFFFFFFFBFFEFFFFFull);
  inline const double LOW_REPR_SAT8   = doubleFromBits(0xFFFFFFFCFFEFFFFFull);
  inline const double LOW_INSTR_SAT8  = doubleFromBits(0xFFFFFFFDFFEFFFFFull);
  inline const double HIGH_INSTR_SAT8 = doubleFromBits(0xFFFFFFFEFFEFFFFFull);
  inline const double HIGH_REPR_SAT8  = doubleFromBits(0xFFFFFFFFFFEFFFFFull);
  inline const double ValidMaximum    = DBL_MAX;

  // 4-byte (float) sentinels, from the ISIS IEEE bit patterns.
  inline const float VALID_MIN4      = floatFromBits(0xFF7FFFFAu);
  inline const float NULL4           = floatFromBits(0xFF7FFFFBu);
  inline const float LOW_REPR_SAT4   = floatFromBits(0xFF7FFFFCu);
  inline const float LOW_INSTR_SAT4  = floatFromBits(0xFF7FFFFDu);
  inline const float HIGH_INSTR_SAT4 = floatFromBits(0xFF7FFFFEu);
  inline const float HIGH_REPR_SAT4  = floatFromBits(0xFF7FFFFFu);
  inline const float VALID_MAX4      = FLT_MAX;

  // 4-byte signed sentinels (the true ISIS SignedInteger specials).
  constexpr int VALID_MINI4      = -2147483632;
  constexpr int NULLI4           = -2147483648;
  constexpr int LOW_REPR_SATI4   = -2147483647;
  constexpr int LOW_INSTR_SATI4  = -2147483646;
  constexpr int HIGH_INSTR_SATI4 = -2147483645;
  constexpr int HIGH_REPR_SATI4  = -2147483644;
  constexpr int VALID_MAXI4      =  2147483647;

  // 2-byte signed sentinels.
  constexpr short VALID_MIN2      = -32752;
  constexpr short NULL2           = -32768;
  constexpr short LOW_REPR_SAT2   = -32767;
  constexpr short LOW_INSTR_SAT2  = -32766;
  constexpr short HIGH_INSTR_SAT2 = -32765;
  constexpr short HIGH_REPR_SAT2  = -32764;
  constexpr short VALID_MAX2      =  32767;

  // 2-byte unsigned sentinels.
  constexpr unsigned short VALID_MINU2      = 3;
  constexpr unsigned short NULLU2           = 0;
  constexpr unsigned short LOW_REPR_SATU2   = 1;
  constexpr unsigned short LOW_INSTR_SATU2  = 2;
  constexpr unsigned short HIGH_INSTR_SATU2 = 65534;
  constexpr unsigned short HIGH_REPR_SATU2  = 65535;
  constexpr unsigned short VALID_MAXU2      = 65522;

  // 1-byte unsigned sentinels.
  constexpr unsigned char VALID_MIN1      = 1;
  constexpr unsigned char NULL1           = 0;
  constexpr unsigned char LOW_REPR_SAT1   = 0;
  constexpr unsigned char LOW_INSTR_SAT1  = 0;
  constexpr unsigned char HIGH_INSTR_SAT1 = 255;
  constexpr unsigned char HIGH_REPR_SAT1  = 255;
  constexpr unsigned char VALID_MAX1      = 254;

} // namespace isis

// Valid data range for an ISIS cube of the given VW channel type. Mirrors the
// per-type switch that ASP's DiskImageResourceIsis valid_minimum/valid_maximum
// use, keyed on the GDAL-reported channel type.
void isisValidRange(ChannelTypeEnum channel_type,
                    double & valid_min, double & valid_max);

// Mask the special pixels of an already-read float image, given the cube's
// channel type and the file nodata (the ISIS Null that GDAL maps). Returns a
// masked image with the special and out-of-range pixels invalid.
ImageViewRef<PixelMask<float>>
maskIsisSpecialPixels(ImageViewRef<float> const& image,
                      ChannelTypeEnum channel_type, float nodata);

} // namespace vw

#endif // __VW_IMAGE_ISISSPECIALPIXELS_H__
