
#ifndef ImageOperators_H
#define ImageOperators_H

#include <vw/GPU/GPUImage.h>
#include <vw/GPU/GenericShaders.h>
#include <vw/GPU/GPUProgram.h>

namespace vw { 
  namespace GPU {
    
 // operator+    
    GENERIC_FRAGMENT_SHADER_FUNCTION_2i0f(operator+, "ImageOperators/add-II-2i0f");
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i1f_IF(operator+, scalar, "ImageOperators/add-IF-1i1f");
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i1f_FI(operator+, scalar, "ImageOperators/add-IF-1i1f");

    template <class PixelT>
    inline GPUImage<PixelT>& operator+=(GPUImage<PixelT>& image1, const GPUImage<PixelT>& image2) {
      image1 = image1 + image2;
      return image1;
    }

    template <class PixelT>
    GPUImage<PixelT> operator+=(GPUImage<PixelT>& image, float scalar) {
      image = image + scalar;
      return image;
    }

// operator-  
    GENERIC_FRAGMENT_SHADER_FUNCTION_2i0f(operator-, "ImageOperators/subtract-II-2i0f");
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i1f_IF(operator-, scalar, "ImageOperators/subtract-IF-1i1f");
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i1f_FI(operator-, scalar, "ImageOperators/subtract-FI-1i1f");

    template <class PixelT>
    inline GPUImage<PixelT>& operator-=(GPUImage<PixelT>& image1, const GPUImage<PixelT>& image2) {
      image1 = image1 - image2;
      return image1;
    }

    template <class PixelT>
    GPUImage<PixelT> operator-=(GPUImage<PixelT>&  image, float scalar) {
      image = image - scalar;
      return image;
    }

// operator*
    GENERIC_FRAGMENT_SHADER_FUNCTION_2i0f(operator*, "ImageOperators/multiply-II-2i0f");
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i1f_IF(operator*, scalar, "ImageOperators/multiply-IF-1i1f");
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i1f_FI(operator*, scalar, "ImageOperators/multiply-IF-1i1f");

    template <class PixelT>
    inline GPUImage<PixelT>& operator*=(GPUImage<PixelT>& image1, const GPUImage<PixelT>& image2) {
      image1 = image1 * image2;
      return image1;
    }

    template <class PixelT>
    GPUImage<PixelT> operator*=(GPUImage<PixelT>& image, float scalar) {
      image = image * scalar;
      return image;
    }


// operator/
    GENERIC_FRAGMENT_SHADER_FUNCTION_2i0f(operator/, "ImageOperators/divide-II-2i0f");
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i1f_IF(operator/, scalar, "ImageOperators/divide-IF-1i1f");
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i1f_FI(operator/, scalar, "ImageOperators/divide-FI-1i1f");

    template <class PixelT>
    inline GPUImage<PixelT>& operator/=(GPUImage<PixelT>& image1, const GPUImage<PixelT>& image2) {
      image1 = image1 / image2;
      return image1;
    }

    template <class PixelT>
    GPUImage<PixelT> operator/=(GPUImage<PixelT>& image, float scalar) {
      image = image / scalar;
      return image;
    }

  } // namespace vw
} // namespace GPU

#endif
