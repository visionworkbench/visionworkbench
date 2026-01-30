// __BEGIN_LICENSE__
//  Copyright (c) 2009-2026, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NGT platform is licensed under the Apache License, Version 2.0 (the
//  "License"); you may not use this file except in compliance with the
//  License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__

// Explicit template instantiations for common ImageView types.
//
// This file compiles common ImageView template instantiations exactly once
// into libVwImage. This prevents the compiler from regenerating the same
// template code in every translation unit that uses these types.
//
// The corresponding extern declarations are in ImageView.h.
// This significantly reduces compilation time and binary size across projects
// that use VisionWorkbench (e.g., Ames Stereo Pipeline).
//
// Only add types here that are truly common across the codebase. Exotic
// pixel types should remain as header-only templates to avoid bloat.

#include <vw/Image/ImageView.h>

namespace vw {

  // Common scalar types
  template class ImageView<float>;
  template class ImageView<double>;
  template class ImageView<uint8>;
  template class ImageView<int16>;
  template class ImageView<uint16>;
  template class ImageView<int32>;
  template class ImageView<uint32>;

  // Common pixel types
  template class ImageView<PixelGray<float>>;
  template class ImageView<PixelGray<double>>;
  template class ImageView<PixelGray<uint8>>;
  template class ImageView<PixelRGB<uint8>>;
  template class ImageView<PixelRGB<float>>;
  template class ImageView<PixelRGB<double>>;
  template class ImageView<PixelRGBA<uint8>>;
  
  // Masked types
  template class ImageView<PixelMask<uint8>>;
  template class ImageView<PixelMask<float>>;
  template class ImageView<PixelMask<double>>;
  template class ImageView<PixelMask<PixelGray<uint8>>>;
  template class ImageView<PixelMask<PixelGray<float>>>;
  template class ImageView<PixelMask<Vector2f>>;
  template class ImageView<PixelMask<Vector2i>>;
  
  // Vector types
  template class ImageView<Vector2>;
  template class ImageView<Vector3>;

} // namespace vw
