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


#ifndef ImageStatistics_H
#define ImageStatistics_H

#include <vw/GPU/GPUImage.h>
#include <vw/GPU/GPUProgram.h>

namespace vw { namespace GPU {

  /// Compute the minimum value stored in all of the channels of all of the planes of the images.
  float min_channel_value(const GPUImageBase& image);

  /// Compute the maximum value stored in all of the channels of all of the planes of the images.
  float max_channel_value(const GPUImageBase& image);

  /// Simultaneously compute the min and max value in all of the
  /// channels of all of the planes of the image.
  inline void min_max_channel_values(const GPUImageBase& image, float& min, float& max) {
    min = min_channel_value(image);
    max = max_channel_value(image);
  }

  /// Compute the sum of the pixels of the channels of all of the planes of the image.
  float sum_channel_value(const GPUImageBase& image);

  /// Compute the mean value stored in all of the channels of all of the planes of the image.
  inline float mean_channel_value(const GPUImageBase& image) {
        return sum_channel_value(image) / (float) (image.width() * image.height() * image.num_channels());
  }

  /// Compute the standard deviation of the values stored in all of the channels of all of the planes of the image.
  float stddev_channel_value(const GPUImageBase& image);

  } // namespace GPU
} // namespace vw

#endif
