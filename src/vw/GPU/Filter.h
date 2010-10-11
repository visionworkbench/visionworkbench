// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef Filter_H
#define Filter_H

#include <vector>

#include <vw/GPU/Setup.h>
#include <vw/GPU/GPUImage.h>
#include <vw/GPU/GPUProgram.h>

#include <vw/Image.h>
#include <vw/Math.h>


namespace vw { namespace GPU {

  // General 2D convolution filter functions

  /// This function computes the convolution of an image with a kernel stored in another image.
  /// It assumes the origin of the kernel is at the point <CODE>(cx,cy)</CODE> and uses the given
  /// edge extension mode to extend the source image as needed.
  /// \see FFTConvolutionFilter

  GPUImageBase convolution_filter(const GPUImageBase& image, const GPUImageBase& kernel);

  template <class PixelT>
  inline GPUImage<PixelT>
  convolution_filter(const GPUImage<PixelT>& image, const GPUImageBase& kernel) {
    return convolution_filter((GPUImageBase&) image, kernel);
  }



  // Separable convolution filter functions


  GPUImageBase seperable_convolution_filter(const GPUImageBase& src,
                       const GPUImageBase& hKernel,
                       const GPUImageBase& vKernel);


  template <class PixelT>
  inline GPUImage<PixelT>
  seperable_convolution_filter(const GPUImage<PixelT>& image,
                                        const GPUImageBase& hKernel,
                                        const GPUImageBase& vKernel) {
    return  seperable_convolution_filter((GPUImageBase&) image, hKernel, vKernel);
  }


  template <class PixelT>
  inline GPUImage<PixelT>
  seperable_convolution_filter(const GPUImage<PixelT>& image,
                                        const std::vector<float>& hKernalVector,
                                        const std::vector<float>& vKernalVector) {
    GPUImage<float> hKernel(hKernalVector.size(), 1);
    hKernel.write(GPU_RED, GPU_FLOAT32, &(hKernalVector[0]));
    GPUImage<float> vKernel(vKernalVector.size(), 1);
    vKernel.write(GPU_RED, GPU_FLOAT32, &(vKernalVector[0]));
    return seperable_convolution_filter((GPUImageBase&) image, hKernel, vKernel);
  }

  // Gaussian convolution functions

  /// This function applies a Gaussian smoothing filter to an image.
  /// It uses an axis-aligned Guassian kernel with standard deviations
  /// of x_sigma and y_sigma, kernel dimensions of x_dim and y_dim,
  /// and the given edge extension mode to extend the source image as
  /// needed.  Specifying a zero falue of x_dim or y_dim causes the
  /// corresponding dimension to be chose automatically as appropriate
  /// for the requested standard deviation.

  template <class PixelT>
  inline GPUImage<PixelT>
  gaussian_filter(const GPUImage<PixelT>& image, float x_sigma, float y_sigma, int x_dim = 0, int y_dim = 0) {
    std::vector<float> x_kernel, y_kernel;
    generate_gaussian_kernel(x_kernel, x_sigma, x_dim);
    generate_gaussian_kernel(y_kernel, y_sigma, y_dim);
    GPUImageBase x_kernel_tex(x_kernel.size(), 1, GPU_RED, GPU_FLOAT32, GPU_RED, GPU_FLOAT32, &(x_kernel[0]));
    GPUImageBase y_kernel_tex(y_kernel.size(), 1, GPU_RED, GPU_FLOAT32, GPU_RED, GPU_FLOAT32, &(y_kernel[0]));
    return seperable_convolution_filter((GPUImageBase&) image, x_kernel_tex, y_kernel_tex);
  }

  // Image differentiation functions

  /// Applies a differentiation filter to an image.  This function
  /// performs differentiation of order x_deriv in the x direction and
  /// order y_deriv in the y direction.  The dimensions of the
  /// differentiation kernel are x_dim and y_dim, and must be odd
  /// numbers large enough to accomodate a kernel of the requested
  /// order.  You may also specify zero for either kernel dimension,
  /// in which case it is chosen automatically to be the smallest
  /// suitable dimension.  The source image is edge-exetnded using the
  /// given edge extension mode as needed.

  template <class PixelT>
  inline GPUImage<PixelT>
  derivative_filter(const GPUImage<PixelT>& image, int x_deriv, int y_deriv,  int x_dim = 0, int y_dim = 0) {
    std::vector<float> x_kernel, y_kernel;
    generate_derivative_kernel(x_kernel, x_deriv, x_dim);
    generate_derivative_kernel(y_kernel, y_deriv, y_dim);
    GPUImageBase x_kernel_tex(x_kernel.size(), 1, GPU_RED, GPU_FLOAT32, GPU_RED, GPU_FLOAT32, &(x_kernel[0]));
    GPUImageBase y_kernel_tex(y_kernel.size(), 1, GPU_RED, GPU_FLOAT32, GPU_RED, GPU_FLOAT32, &(y_kernel[0]));
    return seperable_convolution_filter((GPUImageBase&) image, x_kernel_tex, y_kernel_tex);
  }

  template <class PixelT>
  GPUImage<PixelT>
  inline laplacian_filter(const GPUImage<PixelT>& image) {
    ImageView<float> kernel(3, 3);
    kernel(0,0)=0; kernel(1,0)=1;  kernel(2,0)=0;
    kernel(0,1)=1; kernel(1,1)=-4; kernel(2,1)=1;
    kernel(0,2)=0; kernel(1,2)=1;  kernel(2,2)=0;
    GPUImageBase kernel_tex(3, 3, GPU_RED, GPU_FLOAT32, GPU_RED, GPU_FLOAT32, &(kernel(0, 0)));
    return convolution_filter((GPUImageBase&) image, kernel_tex);
  }

  template <class PixelT>
  GPUImage<PixelT>
    inline gaussian_derivative_filter(const GPUImage<PixelT>& image, float x_sigma, int x_deriv, float y_sigma, int y_deriv, float angle, int size) {
     ImageView<float> kernel;
    generate_gaussian_derivative_kernel( kernel, x_sigma, x_deriv, y_sigma, y_deriv, angle, size );
    GPUImageBase kernel_tex(kernel.cols(), kernel.rows(), GPU_RED, GPU_FLOAT32, GPU_RED, GPU_FLOAT32, &(kernel(0, 0)));
    return convolution_filter((GPUImageBase&) image, kernel_tex);
  }

} } // namespaces GPU, vw

#endif

