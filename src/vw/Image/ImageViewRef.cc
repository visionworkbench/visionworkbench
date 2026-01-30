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

// Explicit template instantiations for common ImageViewRef types.
//
// This file compiles common ImageViewRef template instantiations exactly once
// into libVwImage. This prevents the compiler from regenerating the same
// template code in every translation unit that uses these types.
//
// The corresponding extern declarations are in ImageViewRef.h.
// This significantly reduces compilation time and binary size across projects
// that use VisionWorkbench (e.g., Ames Stereo Pipeline).
//
// Only add types here that are truly common across the codebase. Exotic
// pixel types should remain as header-only templates to avoid bloat.

#include <vw/Image/ImageViewRef.h>

namespace vw {

  // ImageViewRef types
  template class ImageViewRef<float>;
  template class ImageViewRef<double>;
  template class ImageViewRef<uint8>;
  template class ImageViewRef<PixelGray<float>>;
  template class ImageViewRef<PixelGray<uint8>>;
  template class ImageViewRef<PixelRGB<uint8>>;
  template class ImageViewRef<PixelMask<float>>;
  template class ImageViewRef<PixelMask<double>>;
  template class ImageViewRef<PixelMask<uint8>>;
  template class ImageViewRef<PixelMask<Vector2f>>;
  template class ImageViewRef<PixelMask<PixelRGB<uint8>>>;
  template class ImageViewRef<Vector2>;
  template class ImageViewRef<Vector3>;

} // namespace vw
