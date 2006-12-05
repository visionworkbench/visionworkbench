
#ifndef ImageAlgorithms_H
#define ImageAlgorithms_H

#include <vw/GPU/GPUImage.h>
#include <vw/GPU/GenericShaders.h>
#include <vw/GPU/GPUProgram.h>
#include <vw/GPU/ImageStatistics.h>

namespace vw { 
  namespace GPU {

    // fill

    void fill(GPUImageBase& image, float red, float green, float blue, float alpha);

    template <class ChannelT>
    inline void fill(GPUImage<PixelGray<ChannelT> >& image, const PixelGray<ChannelT>& color) {
      fill(image, color.v(), 0.0, 0.0, 0.0);
    }

    template <class ChannelT>
    inline void fill(GPUImage<PixelRGB<ChannelT> >& image, const PixelRGB<ChannelT>& color) {
      fill(image, color.r(), color.g(), color.b(), 0.0);
    }
    
    template <class ChannelT>
    inline void fill(GPUImage<PixelRGBA<ChannelT> >& image, const PixelRGBA<ChannelT>& color) {
      fill(image, color.r(), color.g(), color.b(), color.a());
    }
      
    // copy

    GENERIC_FRAGMENT_SHADER_FUNCTION_1i0f(copy, "ImageAlgorithms/copy-1i0f")
    
    // clamp

    GENERIC_FRAGMENT_SHADER_FUNCTION_1i2f(clamp, low, high, "ImageAlgorithms/clamp-1i2f")

    // normalize

    inline GPUImageBase normalize(const GPUImageBase& image, float low = 0.0, float high = 1.0) {
      static GenericFragmentShader_1i3f shader("ImageAlgorithms/normalize-1i3f");
      float old_min = 0.0;
      float old_max = 1.0;
      min_max_channel_values(image, old_min, old_max);
      return shader(image, old_min, low, (high - low) / (old_max - old_min));
    }

    template <class PixelT> 
    inline GPUImage<PixelT> normalize(const GPUImage<PixelT>& image, float low = 0.0, float high = 1.0) {
      return normalize((GPUImageBase&) image, low, high);
    }

    // threshold

    GENERIC_FRAGMENT_SHADER_FUNCTION_1i3f(threshold, thresh, low, high, "ImageAlgorithms/threshold-1i3f");
    

  } // namespace vw
} // namespace GPU

#endif
 
