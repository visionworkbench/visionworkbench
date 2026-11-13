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

#include <gtest/gtest_VW.h>

#include <vw/Image/IsisSpecialPixels.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/PixelMask.h>
#include <test/Helpers.h>

#include <cstring>
#include <cstdint>
#include <cfloat>
#include <vector>

using namespace vw;

// Build a one-row float image from the given values, mask ISIS specials for the
// given channel type and nodata, and return the rasterized masked result.
static ImageView<PixelMask<float>> runMask(std::vector<float> const& vals,
                                           ChannelTypeEnum ct, float nodata) {
  ImageView<float> img(vals.size(), 1);
  for (size_t i = 0; i < vals.size(); i++)
    img(i, 0) = vals[i];
  ImageView<PixelMask<float>> masked = maskIsisSpecialPixels(img, ct, nodata);
  return masked;
}

// Assert the first pixel is valid and equals the valid value, and the remaining
// pixels (the five ISIS specials) are all invalid.
static void checkOneValidRestSpecial(ImageView<PixelMask<float>> const& m,
                                     float valid_value) {
  ASSERT_EQ(m.cols(), 6);
  EXPECT_TRUE(m(0, 0).valid());
  EXPECT_EQ(m(0, 0).child(), valid_value);
  for (int i = 1; i < 6; i++)
    EXPECT_FALSE(m(i, 0).valid());
}

TEST(IsisSpecialPixels, ValidRange) {
  double mn = 0.0, mx = 0.0;

  isisValidRange(VW_CHANNEL_FLOAT32, mn, mx);
  EXPECT_EQ(mn, double(isis::VALID_MIN4));
  EXPECT_EQ(mx, double(isis::VALID_MAX4));

  // The FLOAT64 sentinel is a NaN, so compare the bit patterns, not the values.
  isisValidRange(VW_CHANNEL_FLOAT64, mn, mx);
  uint64_t got = 0, want = 0;
  std::memcpy(&got,  &mn,                sizeof(got));
  std::memcpy(&want, &isis::ValidMinimum, sizeof(want));
  EXPECT_EQ(got, want);
  EXPECT_EQ(mx, DBL_MAX);

  // INT32 uses the true ISIS SignedInteger valid minimum.
  isisValidRange(VW_CHANNEL_INT32, mn, mx);
  EXPECT_EQ(mn, -2147483632.0);
  EXPECT_EQ(mx,  2147483647.0);

  isisValidRange(VW_CHANNEL_INT16, mn, mx);
  EXPECT_EQ(mn, -32752.0);
  EXPECT_EQ(mx,  32767.0);

  isisValidRange(VW_CHANNEL_UINT16, mn, mx);
  EXPECT_EQ(mn, 3.0);
  EXPECT_EQ(mx, 65522.0);

  // UINT8 and the default types (INT8, UINT32) all map to the 1-byte range.
  isisValidRange(VW_CHANNEL_UINT8, mn, mx);
  EXPECT_EQ(mn, 1.0);
  EXPECT_EQ(mx, 254.0);
  isisValidRange(VW_CHANNEL_INT8, mn, mx);
  EXPECT_EQ(mn, 1.0);
  EXPECT_EQ(mx, 254.0);
  isisValidRange(VW_CHANNEL_UINT32, mn, mx);
  EXPECT_EQ(mn, 1.0);
  EXPECT_EQ(mx, 254.0);
}

TEST(IsisSpecialPixels, Float32) {
  std::vector<float> v = {100.0f, isis::NULL4, isis::LOW_REPR_SAT4,
                          isis::LOW_INSTR_SAT4, isis::HIGH_INSTR_SAT4,
                          isis::HIGH_REPR_SAT4};
  ImageView<PixelMask<float>> m = runMask(v, VW_CHANNEL_FLOAT32, isis::NULL4);
  checkOneValidRestSpecial(m, 100.0f);
}

TEST(IsisSpecialPixels, Int16) {
  std::vector<float> v = {100.0f, float(isis::NULL2), float(isis::LOW_REPR_SAT2),
                          float(isis::LOW_INSTR_SAT2), float(isis::HIGH_INSTR_SAT2),
                          float(isis::HIGH_REPR_SAT2)};
  ImageView<PixelMask<float>> m = runMask(v, VW_CHANNEL_INT16, float(isis::NULL2));
  checkOneValidRestSpecial(m, 100.0f);
}

TEST(IsisSpecialPixels, UInt16) {
  std::vector<float> v = {100.0f, float(isis::NULLU2), float(isis::LOW_REPR_SATU2),
                          float(isis::LOW_INSTR_SATU2), float(isis::HIGH_INSTR_SATU2),
                          float(isis::HIGH_REPR_SATU2)};
  ImageView<PixelMask<float>> m = runMask(v, VW_CHANNEL_UINT16, float(isis::NULLU2));
  checkOneValidRestSpecial(m, 100.0f);
}

TEST(IsisSpecialPixels, UInt8) {
  std::vector<float> v = {100.0f, float(isis::NULL1), float(isis::LOW_REPR_SAT1),
                          float(isis::LOW_INSTR_SAT1), float(isis::HIGH_INSTR_SAT1),
                          float(isis::HIGH_REPR_SAT1)};
  ImageView<PixelMask<float>> m = runMask(v, VW_CHANNEL_UINT8, float(isis::NULL1));
  checkOneValidRestSpecial(m, 100.0f);
}

TEST(IsisSpecialPixels, Int32) {
  std::vector<float> v = {1000.0f, float(isis::NULLI4), float(isis::LOW_REPR_SATI4),
                          float(isis::LOW_INSTR_SATI4), float(isis::HIGH_INSTR_SATI4),
                          float(isis::HIGH_REPR_SATI4)};
  ImageView<PixelMask<float>> m = runMask(v, VW_CHANNEL_INT32, float(isis::NULLI4));
  checkOneValidRestSpecial(m, 1000.0f);
}
