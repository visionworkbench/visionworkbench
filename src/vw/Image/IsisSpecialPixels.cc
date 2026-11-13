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

#include <vw/Image/IsisSpecialPixels.h>
#include <vw/Core/Functors.h>
#include <vw/Image/MaskViews.h>
#include <vw/Image/Filter.h>

#include <cmath>
#include <limits>

namespace vw {

// Per-pixel functor that replaces the ISIS float special pixels. The image is
// already read and converted to float, so only the 4-byte ISIS sentinels apply.
class IsisSpecialPixelFunc: public UnaryReturnSameType {
  float m_replacement_low;
  float m_replacement_high;
  float m_replacement_null;
public:
  IsisSpecialPixelFunc(float pix_l, float pix_h, float pix_n):
    m_replacement_low(pix_l), m_replacement_high(pix_h), m_replacement_null(pix_n) {}

  float operator()(float const& pix) const {
    if (pix < isis::VALID_MIN4) {
      if (pix == isis::HIGH_INSTR_SAT4 || pix == isis::HIGH_REPR_SAT4)
        return m_replacement_high;
      else if (pix == isis::NULL4)
        return m_replacement_null;
      else
        return m_replacement_low;
    }
    return pix;
  }
};

// Replace the ISIS float special pixel values with given replacement values.
ImageViewRef<float>
removeIsisSpecialPixels(ImageViewRef<float> const& image,
                        float r_low, float r_high, float r_null) {
  IsisSpecialPixelFunc obj(r_low, r_high, r_null);
  return per_pixel_filter(image.impl(), obj);
}

void isisValidRange(ChannelTypeEnum channel_type,
                    double & valid_min, double & valid_max) {
  switch (channel_type) {
    case VW_CHANNEL_FLOAT64:
      valid_min = isis::ValidMinimum; valid_max = isis::ValidMaximum; break;
    case VW_CHANNEL_FLOAT32:
      valid_min = isis::VALID_MIN4;   valid_max = isis::VALID_MAX4;   break;
    case VW_CHANNEL_INT32:
      valid_min = isis::IVALID_MIN4;  valid_max = isis::IVALID_MAX4;  break;
    case VW_CHANNEL_INT16:
      valid_min = isis::VALID_MIN2;   valid_max = isis::VALID_MAX2;   break;
    case VW_CHANNEL_UINT16:
      valid_min = isis::VALID_MINU2;  valid_max = isis::VALID_MAXU2;  break;
    default:
      valid_min = isis::VALID_MIN1;   valid_max = isis::VALID_MAX1;   break;
  }
}

ImageViewRef<PixelMask<float>>
maskIsisSpecialPixels(ImageViewRef<float> const& image,
                      ChannelTypeEnum channel_type, float nodata) {

  double valid_min = 0.0, valid_max = 0.0;
  isisValidRange(channel_type, valid_min, valid_max);
  float isis_lo = float(valid_min);
  float isis_hi = float(valid_max);

  // Force the low value to be greater than the nodata value
  if (!std::isnan(nodata) && nodata >= isis_lo) {
    isis_lo = std::nextafter(nodata, std::numeric_limits<float>::max());
    if (isis_hi < isis_lo)
      isis_hi = isis_lo;
  }

  ImageViewRef<PixelMask<float>> masked_image = create_mask(image, nodata);
  ImageViewRef<float> processed_image = apply_mask(masked_image, isis_lo);
  // Invalidate pixels outside range
  masked_image = create_mask(processed_image, isis_lo, isis_hi);
  processed_image = apply_mask(masked_image, isis_lo);
  // ISIS-aware removal of special pixels
  processed_image = removeIsisSpecialPixels(processed_image, isis_lo, isis_hi, isis_lo);
  masked_image = create_mask(processed_image, isis_lo);

  return masked_image;
}

} // namespace vw
