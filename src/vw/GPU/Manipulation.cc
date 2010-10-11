// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/GPU/Manipulation.h>

namespace vw {
  namespace GPU {

    GPUImageBase crop(const GPUImageBase& image,
                        unsigned int upper_left_i,
                        unsigned int upper_left_j,
                        unsigned int width,
                        unsigned int height)
      {
// GLState - Setup
        ((GPUImageBase&) image).rasterize_homography();
        ShaderInvocation_SetupGLState(image.width(), image.height());
// Program - Install
        GPUProgram* program = create_gpu_program("Interp/identity");
        program->install();
// OUTPUT
        GPUImageBase temp(image.width(), image.height(), image.format(), image.type());
        glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, temp.target(), temp.name(), 0);
// INPUT
        program->set_input_image("image", image);
// DRAW
        float t_x1 = image.x(upper_left_i);
        float t_x2 = image.x(upper_left_i + width);
        float t_y1 = image.y(upper_left_j);
        float t_y2 = image.y(upper_left_j + height);

        glBegin(GL_QUADS);
        glTexCoord2f(t_x1, t_y1);         glVertex2f(0, 0);
        glTexCoord2f(t_x2, t_y1);         glVertex2f(width, 0);
        glTexCoord2f(t_x2, t_y2);         glVertex2f(width, height);
        glTexCoord2f(t_x1, t_y2);         glVertex2f(0, height);
        glEnd();
// CleanUp State
        program->uninstall();
        ShaderInvocation_CleanupGLState();
        return temp;
      }


    GPUImageBase flip_horizontal(const GPUImageBase& image) {
// GLState - Setup
        ((GPUImageBase&) image).rasterize_homography();
        ShaderInvocation_SetupGLState(image.width(), image.height());
// Program - Install
        GPUProgram* program = create_gpu_program("Interp/identity");
        program->install();
// OUTPUT
        GPUImageBase temp(image.width(), image.height(), image.format(), image.type());
        glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, temp.target(), temp.name(), 0);
// INPUT
        program->set_input_image("image", image);
// DRAW
        int width = image.width();
        int height = image.height();
        int x = 0;
        int y = 0;
        float t0_x1 = image.x(width);
        float t0_x2 = image.x(0);
        float t0_y1 = image.y(0);
        float t0_y2 = image.y(height);

        glBegin(GL_QUADS);
        glTexCoord2f(t0_x1, t0_y1);       glVertex2f(x, y);
        glTexCoord2f(t0_x2, t0_y1);       glVertex2f(x + width, y);
        glTexCoord2f(t0_x2, t0_y2);       glVertex2f(x + width, y + height);
        glTexCoord2f(t0_x1, t0_y2);       glVertex2f(x, y + height);
        glEnd();
// CleanUp State
        program->uninstall();
        ShaderInvocation_CleanupGLState();

        return temp;
    }



    GPUImageBase flip_vertical(const GPUImageBase& image) {
// GLState - Setup
        ((GPUImageBase&) image).rasterize_homography();
        ShaderInvocation_SetupGLState(image.width(), image.height());
// Program - Install
        GPUProgram* program = create_gpu_program("Interp/identity");
        program->install();
// OUTPUT
        GPUImageBase temp(image.width(), image.height(), image.format(), image.type());
        glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, temp.target(), temp.name(), 0);
// INPUT
        glColor3f(1.0, 0.0, 0.0);

        program->set_input_image("image", image);
// DRAW
        int width = image.width();
        int height = image.height();
        int x = 0;
        int y = 0;
        float t0_x1 = image.x(0);
        float t0_x2 = image.x(width);
        float t0_y1 = image.y(height);
        float t0_y2 = image.y(0);

        glBegin(GL_QUADS);
        glTexCoord2f(t0_x1, t0_y1);       glVertex2f(x, y);
        glTexCoord2f(t0_x2, t0_y1);       glVertex2f(x + width, y);
        glTexCoord2f(t0_x2, t0_y2);       glVertex2f(x + width, y + height);
        glTexCoord2f(t0_x1, t0_y2);       glVertex2f(x, y + height);
        glEnd();
// CleanUp State
        program->uninstall();
        ShaderInvocation_CleanupGLState();

        return temp;
    }

 GPUImageBase rotate_180(const GPUImageBase& image)
  {
      float theta = M_PI;
      int width = image.width();
      int height = image.height();

      Matrix<float> h_translate(3, 3);
      h_translate.set_identity();
      h_translate(0,2) = width/2.0;
      h_translate(1,2) = height/2.0;

      Matrix<float> h_rotate(3, 3);
      h_rotate.set_identity();
      h_rotate(0,0) = cos(theta);
      h_rotate(0,1) = sin(theta);
      h_rotate(1,0) = -sin(theta);
      h_rotate(1,1) = cos(theta);

      Matrix<float> h_translate_back(3, 3);
      h_translate_back.set_identity();
      h_translate_back(0,2) = -width/2.0;
      h_translate_back(1,2) = -height/2.0;

      Matrix<float> h(3, 3);
      h = h_translate * h_rotate * h_translate_back;

      return fixed_homography_transform(image, h);
  }

 GPUImageBase rotate_90_cw(const GPUImageBase& image)
  {
      float theta = 0.5*M_PI;
      int width = image.width();
      int height = image.height();

      Matrix<float> h_translate(3, 3);
      h_translate.set_identity();
      h_translate(0,2) = width/2.0;
      h_translate(1,2) = height/2.0;

      Matrix<float> h_rotate(3, 3);
      h_rotate.set_identity();
      h_rotate(0,0) = cos(theta);
      h_rotate(0,1) = sin(theta);
      h_rotate(1,0) = -sin(theta);
      h_rotate(1,1) = cos(theta);

      Matrix<float> h_translate_back(3, 3);
      h_translate_back.set_identity();
      h_translate_back(0,2) = -width/2.0;
      h_translate_back(1,2) = -height/2.0;

      Matrix<float> h(3, 3);
      h = h_translate * h_rotate * h_translate_back;

      Vector2 pt_min = Vector2((width - height) / 2.0, (height - width) / 2.0);
      Vector2 pt_max = Vector2(width - (width - height) / 2.0, height -(height - width) / 2.0);
      // return fixed_homography_transform(image, h, pt_min, pt_max);
  }

GPUImageBase rotate_90_ccw(const GPUImageBase& image)
  {
      float theta = -0.5*M_PI;
      int width = image.width();
      int height = image.height();

      Matrix<float> h_translate(3, 3);
      h_translate.set_identity();
      h_translate(0,2) = width/2.0;
      h_translate(1,2) = height/2.0;

      Matrix<float> h_rotate(3, 3);
      h_rotate.set_identity();
      h_rotate(0,0) = cos(theta);
      h_rotate(0,1) = sin(theta);
      h_rotate(1,0) = -sin(theta);
      h_rotate(1,1) = cos(theta);

      Matrix<float> h_translate_back(3, 3);
      h_translate_back.set_identity();
      h_translate_back(0,2) = -width/2.0;
      h_translate_back(1,2) = -height/2.0;

      Matrix<float> h(3, 3);
      h = h_translate * h_rotate * h_translate_back;

      Vector2 pt_min = Vector2((width - height) / 2.0, (height - width) / 2.0);
      Vector2 pt_max = Vector2(width - (width - height) / 2.0, height -(height - width) / 2.0);
      //  return fixed_homography_transform(image, h, pt_min, pt_max);
  }

GPUImageBase transpose(const GPUImageBase& image) {
// GLState - Setup
        ((GPUImageBase&) image).rasterize_homography();
        int out_width = image.height();
        int out_height = image.width();
        ShaderInvocation_SetupGLState(out_width, out_height);
// Program - Install
        GPUProgram* program = create_gpu_program("Interp/identity");
        program->install();
// OUTPUT
        GPUImageBase temp(out_width, out_height, image.format(), image.type());
        glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, temp.target(), temp.name(), 0);
// INPUT
        program->set_input_image("image1", image);
// DRAW
        ShaderInvocation_DrawRectOneTexture(image, Rectangle2D<int>(0, 0, out_width, out_height),
                                            Rectangle2D<int>(0, 0, out_width, out_height));
// CleanUp State
        program->uninstall();
        ShaderInvocation_CleanupGLState();

        return temp;
    }

    void inset(GPUImageBase& image1, const GPUImageBase& image2, int x, int y, int width, int height) {
// Copy tex if necessary
      GPUImageBase inputTex;
      if(image1.same_real_object(image2))
        inputTex = copy(image2);
      else
        inputTex = image2;
// Setup
      ShaderInvocation_SetupGLState(image1);
// Program - Install
      GPUProgram* program = create_gpu_program("Interp/identity");
      program->install();
// OUTPUT
      glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, image1.target(), image1.name(), 0);
// INPUT
      program->set_input_image("image1", inputTex);
// DRAW
      ShaderInvocation_DrawRectOneTexture(inputTex, Rectangle2D<int>(x, y, width, height),
                                          Rectangle2D<int>(0, 0, width, height));
// CleanUp State
      program->uninstall();
      ShaderInvocation_CleanupGLState();
    }


  } // namespace vw
} // namespace GPU
