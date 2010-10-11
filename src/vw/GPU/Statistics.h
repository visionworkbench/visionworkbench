// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
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
