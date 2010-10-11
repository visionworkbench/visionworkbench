// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef GenericShaders_H
#define GenericShaders_H

#include <vw/GPU/GPUImage.h>
#include <vw/GPU/GPUProgram.h>

namespace vw { namespace GPU {


  class GenericFragmentShader_1i0f {
  protected:
    std::string path;
  public:
    GenericFragmentShader_1i0f(const char* fragmentShaderBaseName);

    template <class PixelT>
    GPUImage<PixelT> operator()(const GPUImage<PixelT>& image1) {
      return operator()((GPUImageBase&) image1);
    }

    GPUImageBase operator()(const GPUImageBase& image1);
  };




  class GenericFragmentShader_2i0f {
  protected:
    std::string path;
  public:
    GenericFragmentShader_2i0f(const char* fragmentShaderBaseName);

    template <class PixelT>
    GPUImage<PixelT> operator()(const GPUImage<PixelT>& image1, const GPUImage<PixelT>& image2)  {
      return operator()((GPUImageBase&) image1, (GPUImageBase&) image2);
    }

    GPUImageBase operator()(const GPUImageBase& image1, const GPUImageBase& image2);
  };




  class GenericFragmentShader_1i1f {
  protected:
    std::string path;
  public:
    GenericFragmentShader_1i1f(const char* fragmentShaderBaseName);

    template <class PixelT>
    GPUImage<PixelT> operator()(const GPUImage<PixelT>& image1, float float1)  {
      return operator()((GPUImageBase&) image1, float1);
    }

    GPUImageBase operator()(const GPUImageBase& image1, float float1);
  };



  class GenericFragmentShader_1i2f {
  protected:
    std::string path;
  public:
    GenericFragmentShader_1i2f(const char* fragmentShaderBaseName);

    template <class PixelT>
    GPUImage<PixelT> operator()(const GPUImage<PixelT>& image1, float float1, float float2)  {
      return operator()((GPUImageBase&) image1, float1, float2);
    }

    GPUImageBase operator()(const GPUImageBase& image1, float float1, float float2);
  };



  class GenericFragmentShader_1i3f {
  protected:
    std::string path;
  public:
    GenericFragmentShader_1i3f(const char* fragmentShaderBaseName);

    template <class PixelT>
    GPUImage<PixelT> operator()(const GPUImage<PixelT>& image1, float float1, float float2, float float3)  {
      return operator()((GPUImageBase&) image1, float1, float2, float3);
    }

    GPUImageBase operator()(const GPUImageBase& image1, float float1, float float2, float float3);
  };



  class GenericFragmentShader_1i4f {
  protected:
    std::string path;
  public:
    GenericFragmentShader_1i4f(const char* fragmentShaderBaseName);

    template <class PixelT>
    GPUImage<PixelT> operator()(const GPUImage<PixelT>& image1, float float1, float float2, float float3, float float4)  {
      return operator()((GPUImageBase&) image1, float1, float2, float3, float4);
    }

    GPUImageBase operator()(const GPUImageBase& image1, float float1, float float2, float float3, float float4);
  };

#define DUMMY 1


#define GENERIC_FRAGMENT_SHADER_FUNCTION_1i0f(NAME, FILENAME)         \
  /* \cond INTERNAL */  \
  inline GPUImageBase NAME(const GPUImageBase& image1) {       \
    static GenericFragmentShader_1i0f shader(FILENAME);      \
    return shader(image1);      \
  }   \
  /* \endcond */  \
  template<class PixelT>                                     \
  inline GPUImage<PixelT> NAME(const GPUImage<PixelT>& image1) {  \
    return NAME((GPUImageBase&) image1);                        \
  }

#define GENERIC_FRAGMENT_SHADER_FUNCTION_2i0f(NAME, FILENAME)         \
  /* \cond INTERNAL */ \
  inline GPUImageBase NAME(const GPUImageBase& image1, const GPUImageBase& image2) {       \
    static GenericFragmentShader_2i0f shader(FILENAME);      \
    return shader(image1, image2);                                    \
  }  \
  /* \endcond */ \
  template<class PixelT>                                     \
  inline GPUImage<PixelT> NAME(const GPUImage<PixelT>& image1, const GPUImage<PixelT>& image2) {  \
    return NAME((GPUImageBase&) image1, (GPUImageBase&) image2);                                   \
  }

#define GENERIC_FRAGMENT_SHADER_FUNCTION_1i1f_IF(NAME, F1, FILENAME)      \
  /* \cond INTERNAL */  \
  inline GPUImageBase NAME(const GPUImageBase& image1, float N1) {       \
    static GenericFragmentShader_1i1f shader(FILENAME);      \
    return shader(image1, N1);                                    \
  }         \
  /* \endcond */  \
  template<class PixelT>                                     \
  inline GPUImage<PixelT> NAME(const GPUImage<PixelT>& image1, float N1) {  \
    return NAME((GPUImageBase&) image1, N1);                                   \
  }

#define GENERIC_FRAGMENT_SHADER_FUNCTION_1i1f_FI(NAME, F1, FILENAME)      \
  /* \cond INTERNAL */  \
  inline GPUImageBase NAME(float F1, const GPUImageBase& image1) {       \
    static GenericFragmentShader_1i1f shader(FILENAME);      \
    return shader(image1, F1);                                    \
  }   \
  /* \endcond */  \
  template<class PixelT>                                     \
  inline GPUImage<PixelT> NAME(float F1, const GPUImage<PixelT>& image1) {  \
    return NAME(F1, (GPUImageBase&) image1);                                     \
  }

#define GENERIC_FRAGMENT_SHADER_FUNCTION_1i2f(NAME, F1, F2, FILENAME)    \
  /* \cond INTERNAL */  \
  inline GPUImageBase NAME(const GPUImageBase& image1, float F1, float F2) {       \
    static GenericFragmentShader_1i2f shader(FILENAME);  \
    return shader(image1, F1, F2);                                    \
  }  \
  /* \endcond */   \
  template<class PixelT>         \
  inline GPUImage<PixelT> NAME(const GPUImage<PixelT>& image1, float F1, float F2) {  \
    return NAME((GPUImageBase&) image1, F1, F2);                                     \
  }

#define GENERIC_FRAGMENT_SHADER_FUNCTION_1i3f(NAME, F1, F2, F3, FILENAME)      \
  /* \cond INTERNAL */ \
  inline GPUImageBase NAME(const GPUImageBase& image1, float F1, float F2, float F3) {       \
    static GenericFragmentShader_1i3f shader(FILENAME);      \
    return shader(image1, F1, F2, F3);                                    \
  }    \
  /* \endcond  */ \
  template<class PixelT>                                     \
  inline GPUImage<PixelT> NAME(const GPUImage<PixelT>& image1, float F1, float F2, float F3) {  \
    return NAME((GPUImageBase&) image1, F1, F2, F3);                                     \
  }

#define GENERIC_FRAGMENT_SHADER_FUNCTION_1i4f(NAME, F1, F2, F3, F4, FILENAME)      \
  /* \cond INTERNAL */ \
  inline GPUImageBase NAME(const GPUImageBase& image1, float F1, float F2, float F3, float F4) {       \
    static GenericFragmentShader_1i4f shader(FILENAME);      \
    return shader(image1, F1, F2, F3, F4);                                    \
  }                \
  /* \endcond */ \
  template<class PixelT>                                     \
  inline GPUImage<PixelT> NAME(const GPUImage<PixelT>& image1, float F1, float F2, float F3, float F4) {  \
    return NAME((GPUImageBase&)image1, F1, F2, F3, F4);                                     \
  }

  } // namespace vw
} // namespace GPU

#endif
