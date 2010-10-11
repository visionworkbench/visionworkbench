// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef Manipulation_H
#define Manipulation_H

#include <vw/Image.h>
#include <vw/Math.h>

#include <vw/GPU/Setup.h>
#include <vw/GPU/GPUImage.h>
#include <vw/GPU/GPUProgram.h>
#include <vw/GPU/Transform.h>
#include <vw/GPU/Algorithms.h>


namespace vw {
  namespace GPU {

    // crop

    GPUImageBase crop(const GPUImageBase& image,
                        unsigned int upper_left_i,
                        unsigned int upper_left_j,
                        unsigned int width,
                      unsigned int height);


    template <class PixelT>
    inline GPUImage<PixelT> crop(const GPUImage<PixelT>& src,
                        unsigned int upper_left_i,
                        unsigned int upper_left_j,
                        unsigned int width,
                        unsigned int height)
      {
        return crop((GPUImageBase&) src, upper_left_i, upper_left_j, width, height);
      }

    //  flip_horizontal

    GPUImageBase flip_horizontal(const GPUImageBase& image);

    template <class PixelT>
    inline GPUImage<PixelT> flip_horizontal(const GPUImage<PixelT>& image) {
      return flip_horizontal((GPUImageBase&) image);
    }


    //  flip_vertical

    GPUImageBase flip_vertical(const GPUImageBase& image);

    template <class PixelT>
    inline GPUImage<PixelT> flip_vertical(const GPUImage<PixelT>& image) {
      return flip_vertical((GPUImageBase&) image);
    }

    //  rotate_180

    GPUImageBase rotate_180(const GPUImageBase& image);

    template <class PixelT>
    inline GPUImage<PixelT> rotate_180(const GPUImage<PixelT>& image) {
      return rotate_180((GPUImageBase&) image);
    }

    //  rotate_90_cw

    GPUImageBase rotate_90_cw(const GPUImageBase& image);

    template <class PixelT>
    inline GPUImage<PixelT> rotate_90_cw(const GPUImage<PixelT>& image) {
      return rotate_90_cw((GPUImageBase&) image);
    }

    //  rotate_90_ccw

    GPUImageBase rotate_90_ccw(const GPUImageBase& image);

    template <class PixelT>
    inline GPUImage<PixelT> rotate_90_ccw(const GPUImage<PixelT>& image) {
      return rotate_90_ccw((GPUImageBase&) image);
    }

    // subsample

    inline GPUImageBase subsample(const GPUImageBase& image, unsigned xfactor, unsigned yfactor) {
      float x_resample_ratio = 1.0 / xfactor;
      float y_resample_ratio = 1.0 / yfactor;
      // return resample(image, x_resample_ratio, y_resample_ratio, (int) ceilf(image.width() * x_resample_ratio),
      //              (int) ceilf(image.height() * y_resample_ratio), interpolation::NearestPixel(), ConstantEdgeExtend());
    }

    inline GPUImageBase subsample(const GPUImageBase& image, unsigned subsampling_factor) {
      subsample(image, subsampling_factor, subsampling_factor);
    }

    template <class PixelT>
    inline GPUImage<PixelT> subsample(const GPUImage<PixelT>& image, unsigned subsampling_factor) {
      return subsample((GPUImageBase&) image, subsampling_factor, subsampling_factor);
    }

    template <class PixelT>
    inline GPUImage<PixelT> subsample(const GPUImage<PixelT>& image, unsigned xfactor, unsigned yfactor) {
      return subsample((GPUImageBase&) image, xfactor, yfactor);
    }

    // transpose

    GPUImageBase transpose(const GPUImageBase& image);

    template <class PixelT>
    inline GPUImage<PixelT> transpose(const GPUImage<PixelT>& image) {
      return transpose((GPUImageBase&) image);
    }

    // inset

    void inset(GPUImageBase& image1, const GPUImageBase& image2, int x, int y, int width, int height);

    inline void inset(GPUImageBase& image1, const GPUImageBase& image2, int x = 0, int y = 0) {
      inset(image1, image2, x, y, image2.width(), image2.height());
    }

    template <class PixelT>
    inline void inset(GPUImage<PixelT>& image1, const GPUImage<PixelT>& image2, int x, int y, int width, int height) {
      return inset((GPUImageBase&) image1, (GPUImageBase&) image2, x, y, width, height);
    }

    template <class PixelT>
    inline void inset(GPUImage<PixelT>& image1, const GPUImage<PixelT>& image2, int x = 0, int y = 0) {
      return inset((GPUImageBase&) image1, (GPUImageBase&) image2, x, y, image2.width(), image2.height());
    }


  // Image Packing - Gray -> RGBA

    template <class ChannelT>
    GPUImage<PixelRGBA<ChannelT> > pack_gray_into_rgba(const GPUImage<PixelGray<ChannelT> >& image, int overlap) {
      // Static
      // GLState - Setup
          ((GPUImageBase&) image).rasterize_homography();
      ShaderInvocation_SetupGLState(image.width(), image.height());
      // Program - Install
      GPUProgram* program = create_gpu_program("Manipulation/PackGrayIntoRGBA");
      program->install();
      // OUTPUT
      int outputWidth = overlap + (int) ceilf(image.width() / 2.0);
      int outputHeight = overlap + (int) ceilf(image.height() / 2.0);
      GPUImage<PixelRGBA<ChannelT> > output(outputWidth, outputHeight);
      glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, output.target(), output.name(), 0);
      // INPUT
      program->set_input_image("image", 0, image);
      program->set_input_float("xCellOffset", outputWidth - 2*overlap);
      program->set_input_float("yCellOffset", outputHeight - 2*overlap);
      // DRAW
      ShaderInvocation_DrawRectOneTexture(image);
      // CleanUp State
      program->uninstall();
      ShaderInvocation_CleanupGLState();

      return output;

    }

    template <class ChannelT>
    GPUImage<PixelGray<ChannelT> > unpack_gray_from_rgba(const GPUImage<PixelRGBA<ChannelT> >& image, int overlap) {
      // GLState - Setup
          ((GPUImageBase&) image).rasterize_homography();
      int inputWidth = image.width();
      int inputHeight = image.height();
      int outputWidth =  2 * (inputWidth - overlap);
      int outputHeight =  2 * (inputHeight - overlap);
      ShaderInvocation_SetupGLState(outputWidth, outputHeight);
      // OUTPUT
      GPUImage<PixelRGBA<ChannelT> > output(outputWidth, outputHeight);
      glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, output.target(), output.name(), 0);
      float tImage_x1 = image.x(0);
      float tImage_x2 = image.x(inputWidth);
      float tImage_y1 = image.y(0);
      float tImage_y2 = image.y(inputHeight);
      float t_x1, t_x2, t_y1, t_y2;
      float v_x1, v_x2, v_y1, v_y2;
      // 1. R-Channel
      GPUProgram* program = create_gpu_program("Manipulation/SelectChannel-R");
      program->install();
      program->set_input_image("image1", 0, image);
      t_x1 = tImage_x1;
      t_x2 = tImage_x2 - overlap;
      t_y1 = tImage_y1;
      t_y2 = tImage_y2 - overlap;
      v_x1 = 0;
      v_x2 = outputWidth / 2;
      v_y1 = 0;
      v_y2 = outputHeight / 2;
      glBegin(GL_QUADS);
      glTexCoord2f(t_x1, t_y1);           glVertex2f(v_x1, v_y1);
      glTexCoord2f(t_x2, t_y1);           glVertex2f(v_x2, v_y1);
      glTexCoord2f(t_x2, t_y2);           glVertex2f(v_x2, v_y2);
      glTexCoord2f(t_x1, t_y2);           glVertex2f(v_x1, v_y2);
      glEnd();
      // 2. G-Channel
      program = create_gpu_program("Manipulation/SelectChannel-G");
      program->install();
      program->set_input_image("image1", 0, image);
      t_x1 = tImage_x1 + overlap;
      t_x2 = tImage_x2;
      t_y1 = tImage_y1;
      t_y2 = tImage_y2 - overlap;
      v_x1 = outputWidth / 2;
      v_x2 = outputWidth;
      v_y1 = 0;
      v_y2 = outputHeight / 2;
      glBegin(GL_QUADS);
      glTexCoord2f(t_x1, t_y1);           glVertex2f(v_x1, v_y1);
      glTexCoord2f(t_x2, t_y1);           glVertex2f(v_x2, v_y1);
      glTexCoord2f(t_x2, t_y2);           glVertex2f(v_x2, v_y2);
      glTexCoord2f(t_x1, t_y2);           glVertex2f(v_x1, v_y2);
      glEnd();
      // 3. B-Channel
      program = create_gpu_program("Manipulation/SelectChannel-B");
      program->install();
      program->set_input_image("image1", 0, image);
      t_x1 = tImage_x1;
      t_x2 = tImage_x2 - overlap;
      t_y1 = tImage_y1 + overlap;
      t_y2 = tImage_y2;
      v_x1 = 0;
      v_x2 = outputWidth / 2;
      v_y1 = outputHeight / 2;
      v_y2 = outputHeight;
      glBegin(GL_QUADS);
      glTexCoord2f(t_x1, t_y1);           glVertex2f(v_x1, v_y1);
      glTexCoord2f(t_x2, t_y1);           glVertex2f(v_x2, v_y1);
      glTexCoord2f(t_x2, t_y2);           glVertex2f(v_x2, v_y2);
      glTexCoord2f(t_x1, t_y2);           glVertex2f(v_x1, v_y2);
      glEnd();
      // 4. A-Channel
      program = create_gpu_program("Manipulation/SelectChannel-A");
      program->install();
      program->set_input_image("image1", 0, image);
      t_x1 = tImage_x1 + overlap;
      t_x2 = tImage_x2;
      t_y1 = tImage_y1 + overlap;
      t_y2 = tImage_y2;
      v_x1 = outputWidth / 2;
      v_x2 = outputWidth;
      v_y1 = outputHeight / 2;
      v_y2 = outputHeight;
      glBegin(GL_QUADS);
      glTexCoord2f(t_x1, t_y1);           glVertex2f(v_x1, v_y1);
      glTexCoord2f(t_x2, t_y1);           glVertex2f(v_x2, v_y1);
      glTexCoord2f(t_x2, t_y2);           glVertex2f(v_x2, v_y2);
      glTexCoord2f(t_x1, t_y2);           glVertex2f(v_x1, v_y2);
      glEnd();
      // GLState CleanUp
      program->uninstall();
      ShaderInvocation_CleanupGLState();

      return output;

    }

  } // namespace vw
} // namespace GPU

#endif
