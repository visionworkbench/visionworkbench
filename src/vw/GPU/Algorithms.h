// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef ImageAlgorithms_H
#define ImageAlgorithms_H

#include <vw/GPU/GPUImage.h>
#include <vw/GPU/GenericShaders.h>
#include <vw/GPU/GPUProgram.h>
#include <vw/GPU/Statistics.h>

namespace vw {
  namespace GPU {

  // *******************************************************************
  // fill()
  // *******************************************************************

  /// \cond INTERNAL
    void fill_impl(GPUImageBase& image, float red, float green, float blue, float alpha, const Rectangle2D<int>& bounds);
  /// \endcond

  /// Fill an image with a constant pixel value
    inline void fill(GPUImageBase& image, float red, float green, float blue, float alpha) {
        Rectangle2D<int> bounds(0, 0, image.width(), image.height());
        fill_impl(image, red, green, blue, alpha, bounds);
    }

  /// Fill the specified rectangular region of an image with a constant pixel value,
    inline void fill(GPUImageBase& image, float red, float green, float blue, float alpha, const Rectangle2D<int>& bounds) {
        fill_impl(image, red, green, blue, alpha, bounds);
    }

  /// Fill an image with a constant pixel value
    template <class ChannelT>
    inline void fill(GPUImage<PixelGray<ChannelT> >& image, const PixelGray<ChannelT>& color) {
      fill(image, color.v(), 0.0, 0.0, 0.0);
    }

  /// Fill an image with a constant pixel value
    template <class ChannelT>
    inline void fill(GPUImage<PixelRGB<ChannelT> >& image, const PixelRGB<ChannelT>& color) {
      fill(image, color.r(), color.g(), color.b(), 0.0);
    }

  /// Fill an image with a constant pixel value
    template <class ChannelT>
    inline void fill(GPUImage<PixelRGBA<ChannelT> >& image, const PixelRGBA<ChannelT>& color) {
      fill(image, color.r(), color.g(), color.b(), color.a());
    }

  /// Fill an image with a constant pixel value
    template <class ChannelT>
    inline void fill(GPUImage<PixelGray<ChannelT> >& image, const PixelGray<ChannelT>& color, const Rectangle2D<int>& bounds) {
      fill(image, color.v(), 0.0, 0.0, 0.0);
    }

  /// Fill an image with a constant pixel value
    template <class ChannelT>
    inline void fill(GPUImage<PixelRGB<ChannelT> >& image, const PixelRGB<ChannelT>& color, const Rectangle2D<int>& bounds) {
      fill(image, color.r(), color.g(), color.b(), 0.0);
    }

  /// Fill an image with a constant pixel value
    template <class ChannelT>
    inline void fill(GPUImage<PixelRGBA<ChannelT> >& image, const PixelRGBA<ChannelT>& color, const Rectangle2D<int>& bounds) {
      fill(image, color.r(), color.g(), color.b(), color.a());
    }

  // *******************************************************************
  // copy()
  // *******************************************************************

    /// Create a deep copy of an image.
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i0f(copy, "Algorithms/copy-1i0f")

  // *******************************************************************
  // clamp()
  // *******************************************************************

  /// Clamp the values in an image to fall within the range [low,high].
    GENERIC_FRAGMENT_SHADER_FUNCTION_1i2f(clamp, low, high, "Algorithms/clamp-1i2f")

  /// \cond internal
    inline GPUImageBase clamp(const GPUImageBase& image, float high) {
      static GenericFragmentShader_1i2f shader("Algorithms/clamp-1i2f");
      return shader(image, 0.0, high);
    }
  /// \endcond

  /// Clamp the values in an image to fall within the range [0,high].
  /// The low end of the range is actually determined by the
  /// ChannelRange type trait but is generally zero.
    template <class PixelT>
    inline GPUImage<PixelT> clamp(const GPUImage<PixelT>& image, float high) {
      return clamp((GPUImageBase&) image, high);
    }

  /// \cond internal
    inline GPUImageBase clamp(const GPUImageBase& image) {
      static GenericFragmentShader_1i2f shader("Algorithms/clamp-1i2f");
      return shader(image, 0.0, 1.0);
    }
  /// \endcond

  /// Clamp the values in an image to fall within the range [min,max],
  /// where min and max are determined by the ChannelRange type trait
  /// and are generally equal to 0.0 and 1.0 for floating point types
  /// and 0 and the largest positve value for integral types.
    template <class PixelT>
    inline GPUImage<PixelT> clamp(const GPUImage<PixelT>& image) {
      return clamp((GPUImageBase&) image);
    }

  // *******************************************************************
  // normalize()
  // *******************************************************************

  /// \cond INTERNAL
    inline GPUImageBase normalize(const GPUImageBase& image, float low, float high) {
      static GenericFragmentShader_1i3f shader("Algorithms/normalize-1i3f");
      float old_min = 0.0;
      float old_max = 1.0;
      min_max_channel_values(image, old_min, old_max);
      return shader(image, old_min, low, (high - low) / (old_max - old_min));
    }
  /// \endcond

  /// Renormalize the values in an image to fall within the range [low,high).
    template <class PixelT>
    inline GPUImage<PixelT> normalize(const GPUImage<PixelT>& image, float low, float high) {
      return normalize((GPUImageBase&) image, low, high);
    }

  /// \cond INTERNAL
    inline GPUImageBase normalize(const GPUImageBase& image, float high) {
        normalize(image, 0, high);
    }
  /// \endcond

  /// Renormalize the values in an image to fall within the range [0,high).
    template <class PixelT>
    inline GPUImage<PixelT> normalize(const GPUImage<PixelT>& image, float high) {
      return normalize((GPUImageBase&) image, 0, high);
    }

  /// \cond INTERNAL
    inline GPUImageBase normalize(const GPUImageBase& image) {
        normalize(image, 0, 1);
    }
  /// \endcond

  /// Renormalize the values in an image to fall within the range [0, 1).
    template <class PixelT>
    inline GPUImage<PixelT> normalize(const GPUImage<PixelT>& image) {
      return normalize((GPUImageBase&) image, 0, 1);
    }



  // *******************************************************************
  // threshold()
  // *******************************************************************

  /// Threshold the values in an image, generating a two-valued output
  /// image with values low and high.
    GPUImageBase threshold(const GPUImageBase& image, float thresh, float low, float high);

  /// \cond internal
    inline GPUImageBase threshold(const GPUImageBase& image, float thresh, float high) {
        threshold(image, thresh, 0.0, high);
    }
  /// \endcond

  /// Threshold the values in an image, generating a two-valued output
  /// image with values 0 and high.
    template <class PixelT>
    inline GPUImage<PixelT> threshold(const GPUImage<PixelT>& image, float thresh, float high) {
      return threshold((GPUImageBase&) image, thresh, 0.0, high);
    }

  /// \cond internal
    inline GPUImageBase threshold(const GPUImageBase& image, float thresh) {
        threshold(image, thresh, 0.0, 1.0);
    }
  /// \endcond

  /// Threshold the values in an image, generating a two-valued output
  /// where the values are either 0 or 1.
    template <class PixelT>
    inline GPUImage<PixelT> threshold(const GPUImage<PixelT>& image, float thresh) {
      return threshold((GPUImageBase&) image, thresh, 0.0, 1.0);
    }

  /// \cond internal
    inline GPUImageBase threshold(const GPUImageBase& image) {
        threshold(image, 0.0, 0.0, 1.0);
    }
  /// \endcond

  /// Threshold the values in an image, generating a two-valued output
  /// where the values are either 0 or 1.
    template <class PixelT>
    inline GPUImage<PixelT> threshold(const GPUImage<PixelT>& image) {
      return threshold((GPUImageBase&) image, 0.0, 0.0, 1.0);
    }

  } // namespace vw
} // namespace GPU

#endif

