
#ifndef Filter_H
#define Filter_H

#include "OpenGLManager.h"
#include "GPUImage.h"
#include "GPUProgram.h"

#include <vw/Filter.h>

#include <vw/ImageView.h>
#include <vw/PixelTypes.h>
#include <vw/Convolution.h>
#include <vw/Fourier.h>
#include <vw/ImageMath.h>
#include <vw/ImageOperators.h>

namespace vw { namespace GPU {


// Private generic functions



  // Public templetized functions

  GPUImageBase convolution_filter(const GPUImageBase& image, 
				 const GPUImageBase& kernel);

  template <class PixelT>
  GPUImage<PixelT> 
  inline convolution_filter(const GPUImage<PixelT>& image, const GPUImageBase& kernel) {
    return convolution_filter((GPUImageBase&) image, kernel);
  }

  GPUImageBase seperable_convolution_filter(const GPUImageBase& src, 
					   const GPUImageBase& hKernel, 
					   const GPUImageBase& vKernel);


  template <class PixelT>
  GPUImage<PixelT> 
  inline seperable_convolution_filter(const GPUImage<PixelT>& image, const GPUImageBase& hKernel, const GPUImageBase& vKernel) {
    return  seperable_convolution_filter((GPUImageBase&) image, hKernel, vKernel);
  }

  template <class PixelT>
  GPUImage<PixelT> 
  inline gaussian_filter(const GPUImage<PixelT>& image, float x_sigma, float y_sigma, int x_dim = 0, int y_dim = 0) {
    std::vector<float> x_kernel, y_kernel;
    generate_gaussian_kernel(x_kernel, x_sigma, x_dim);
    generate_gaussian_kernel(y_kernel, y_sigma, y_dim);
    GPUImageBase x_kernel_tex(x_kernel.size(), 1, TEX_R, TEX_FLOAT32, TEX_R, TEX_FLOAT32, &(x_kernel[0]));
    GPUImageBase y_kernel_tex(y_kernel.size(), 1, TEX_R, TEX_FLOAT32, TEX_R, TEX_FLOAT32, &(y_kernel[0]));
    return seperable_convolution_filter((GPUImageBase&) image, x_kernel_tex, y_kernel_tex);
  }

  template <class PixelT>
  GPUImage<PixelT> 
  inline derivative_filter(const GPUImage<PixelT>& image, int x_deriv, int y_deriv,  int x_dim = 0, int y_dim = 0) {
    std::vector<float> x_kernel, y_kernel;
    generate_derivative_kernel(x_kernel, x_deriv, x_dim);
    generate_derivative_kernel(y_kernel, y_deriv, y_dim);
    GPUImageBase x_kernel_tex(x_kernel.size(), 1, TEX_R, TEX_FLOAT32, TEX_R, TEX_FLOAT32, &(x_kernel[0]));
    GPUImageBase y_kernel_tex(y_kernel.size(), 1, TEX_R, TEX_FLOAT32, TEX_R, TEX_FLOAT32, &(y_kernel[0]));
    return seperable_convolution_filter((GPUImageBase&) image, x_kernel_tex, y_kernel_tex);
  }

  template <class PixelT>
  GPUImage<PixelT> 
  inline laplacian_filter(const GPUImage<PixelT>& image) {
    ImageView<float> kernel(3, 3);
    kernel(0,0)=0; kernel(1,0)=1;  kernel(2,0)=0;
    kernel(0,1)=1; kernel(1,1)=-4; kernel(2,1)=1;
    kernel(0,2)=0; kernel(1,2)=1;  kernel(2,2)=0;
    GPUImageBase kernel_tex(3, 3, TEX_R, TEX_FLOAT32, TEX_R, TEX_FLOAT32, &(kernel(0, 0)));
    return convolution_filter((GPUImageBase&) image, kernel_tex);
  }

  template <class PixelT>
  GPUImage<PixelT> 
    inline gaussian_derivative_filter(const GPUImage<PixelT>& image, float x_sigma, int x_deriv, float y_sigma, int y_deriv, float angle, int size) {
    ImageView<float> kernel;
    generate_gaussian_derivative_kernel( kernel, x_sigma, x_deriv, y_sigma, y_deriv, angle, size );
    GPUImageBase kernel_tex(kernel.cols(), kernel.rows(), TEX_R, TEX_FLOAT32, TEX_R, TEX_FLOAT32, &(kernel(0, 0)));
    return convolution_filter((GPUImageBase&) image, kernel_tex);
  }

} } // namespaces GPU, vw

#endif

