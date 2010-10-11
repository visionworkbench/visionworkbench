// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef ImageMath_H
#define ImageMath_H

#include <vw/GPU/GPUImage.h>
#include <vw/GPU/GPUProgram.h>
#include <vw/GPU/GenericShaders.h>




namespace vw {
  namespace GPU {

 // Math Operators

 // operator+
    GENERIC_FRAGMENT_SHADER_FUNCTION_2i0f(operator+, "ImageMath/add-II-2i0f")
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i1f_IF(operator+, scalar, "ImageMath/add-IF-1i1f")
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i1f_FI(operator+, scalar, "ImageMath/add-IF-1i1f")

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
    GENERIC_FRAGMENT_SHADER_FUNCTION_2i0f(operator-, "ImageMath/subtract-II-2i0f")
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i1f_IF(operator-, scalar, "ImageMath/subtract-IF-1i1f")
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i1f_FI(operator-, scalar, "ImageMath/subtract-FI-1i1f")

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
    GENERIC_FRAGMENT_SHADER_FUNCTION_2i0f(operator*, "ImageMath/multiply-II-2i0f")
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i1f_IF(operator*, scalar, "ImageMath/multiply-IF-1i1f")
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i1f_FI(operator*, scalar, "ImageMath/multiply-IF-1i1f")

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
    GENERIC_FRAGMENT_SHADER_FUNCTION_2i0f(operator/, "ImageMath/divide-II-2i0f")
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i1f_IF(operator/, scalar, "ImageMath/divide-IF-1i1f")
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i1f_FI(operator/, scalar, "ImageMath/divide-FI-1i1f")

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


   // Math Functions

  /// Computes the absolute value of each channel of each pixel in an image.
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i0f(abs, "ImageMath/abs-1i0f")
  /// Computes the arccosine, \f$\cos^{-1} x\f$, of each channel of each pixel in an image.
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i0f(acos, "ImageMath/acos-1i0f")
  /// Computes the arcsine, \f$\sin^{-1} x\f$, of each channel of each pixel in an image.
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i0f(asin, "ImageMath/asin-1i0f")
  /// Computes the arctangent, \f$\tan^{-1} x\f$, of each channel of each pixel in an image.
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i0f(atan, "ImageMath/atan-1i0f")
  /// Computes the smallest integer not less than each channel of each pixel in an image.
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i0f(ceil, "ImageMath/ceil-1i0f")
  /// Computes the cosine, \f$\cos x\f$, of each channel of each pixel in an image.
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i0f(cos, "ImageMath/cos-1i0f")
  /// Computes the base-2 exponential, \f$2^x\f$, of each channel of each pixel in an image.
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i0f(exp2, "ImageMath/exp2-1i0f")
  /// Computes the largest integer not greater than each channel of each pixel in an image.
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i0f(floor, "ImageMath/floor-1i0f")
  /// Computes the nearest integer to each channel of each pixel in an image, rounding half-way cases away from zero.
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i0f(sin, "ImageMath/sin-1i0f")
  /// Computes the square root of each channel of each pixel in an image.
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i0f(sqrt, "ImageMath/sqrt-1i0f")
  /// Computes the tangent, \f$\tan x\f$, of each channel of each pixel in an image.
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i0f(tan, "ImageMath/tan-1i0f")
  /// Computes the cube root of each channel of each pixel in an image.
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i0f(cbrt, "ImageMath/cbrt-1i0f")
  /// Computes the base-e exponential, \f$e^x\f$, of each channel of each pixel in an image.
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i0f(exp, "ImageMath/exp-1i0f")
  /// Computes the base-e exponential minus one, \f$e^x-1\f$, of each channel of each pixel in an image.
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i0f(expm1, "ImageMath/expm1-1i0f")
  /// Computes the natural logarithm, \f$\ln x\f$, of each channel of each pixel in an image.
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i0f(log, "ImageMath/log-1i0f")
  /// Computes the natural logarithm of one plus each channel of each pixel in an image, \f$\ln(1+x)\f$.
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i0f(log1p, "ImageMath/log1p-1i0f")
  /// Computes the base-10 logarithm, \f$\log_{10} x\f$, of each channel of each pixel in an image.
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i0f(log10, "ImageMath/log10-1i0f")
  /// Computes the largest integer not greater than each channel of each pixel in an image.
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i0f(trunc, "ImageMath/trunc-1i0f")
  /// Computes the nearest integer to each channel of each pixel in an image, rounding half-way cases away from zero.
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i0f(round, "ImageMath/round-1i0f")
  /// Computes the power function, \f$x^y\f$, of each channel of each pixel from two images.
    GENERIC_FRAGMENT_SHADER_FUNCTION_2i0f(pow, "ImageMath/pow-II-2i0f")
  /// Computes the given power of each channel of each pixel of an image.
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i1f_IF(pow, scalar, "ImageMath/pow-IF-1i1f")
  /// Raises a scalar to the power of each channel of each pixel of an image.
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i1f_FI(pow, scalar, "ImageMath/pow-FI-1i1f")
  /// Computes the base-2 exponential, \f$2^x\f$, of each channel of each pixel in an image.
    GENERIC_FRAGMENT_SHADER_FUNCTION_2i0f(hypot, "ImageMath/hypot-II-2i0f")
  /// Computes the hypotenuse, \f$\sqrt{x^2+y^2}\f$, of each channel of each pixel of an image with a scalar.
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i1f_IF(hypot, scalar, "ImageMath/hypot-IF-1i1f")
  /// Computes the hypotenuse, \f$\sqrt{x^2+y^2}\f$, of a scalar with each channel of each pixel of an image.
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i1f_FI(hypot, scalar, "ImageMath/hypot-FI-1i1f")

  ///
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i0f(cosh, "ImageMath/cosh-1i0f")
  ///
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i0f(sinh, "ImageMath/sinh-1i0f")
  ///
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i0f(tanh, "ImageMath/tanh-1i0f")
  ///
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i0f(acosh, "ImageMath/acosh-1i0f")
  ///
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i0f(asinh, "ImageMath/asinh-1i0f")
  ///
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i0f(atanh, "ImageMath/atanh-1i0f")
  ///
    GENERIC_FRAGMENT_SHADER_FUNCTION_2i0f(atan2, "ImageMath/atan2-II-2i0f")
  ///
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i1f_IF(atan2, scalar, "ImageMath/atan2-IF-1i1f")
  ///
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i1f_FI(atan2, scalar, "ImageMath/atan2-FI-1i1f")
  ///
    GENERIC_FRAGMENT_SHADER_FUNCTION_2i0f(fdim, "ImageMath/fdim-II-2i0f")
  ///
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i1f_IF(fdim, scalar, "ImageMath/fdim-IF-1i1f")
  ///
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i1f_FI(fdim, scalar, "ImageMath/fdim-FI-1i1f")
  ///
    GENERIC_FRAGMENT_SHADER_FUNCTION_2i0f(copysign, "ImageMath/copysign-II-2i0f")
  ///
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i1f_IF(copysign, scalar, "ImageMath/copysign-IF-1i1f")
  ///
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i1f_FI(copysign, scalar, "ImageMath/copysign-FI-1i1f")
 ///
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i0f(log2, "ImageMath/log2-1i0f")


  } // namespace vw
} // namespace GPU

#endif



