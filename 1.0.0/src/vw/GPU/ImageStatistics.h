
#ifndef ImageStatistics_H
#define ImageStatistics_H

#include <vw/GPU/GPUImage.h>
#include <vw/GPU/GPUProgram.h>

namespace vw { namespace GPU {

  float min_channel_value(const GPUImageBase& image);

  float max_channel_value(const GPUImageBase& image);

  void min_max_channel_values(const GPUImageBase& image, float& min, float& max);

  float mean_channel_value(const GPUImageBase& image); 
 
  } // namespace GPU
} // namespace vw

#endif
