// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/GPU/Algorithms.h>

namespace vw { namespace GPU {

  void fill_impl(GPUImageBase& image, float red, float green, float blue, float alpha, const Rectangle2D<int>& bounds)
  {
    // GLState - Setup
    ((GPUImageBase&) image).rasterize_homography();
    ShaderInvocation_SetupGLState(image);
    // Program - Install
    GPUProgram* program = create_gpu_program("Algorithms/fill-1i4f");
    program->install();
    glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, image.target(), image.name(), 0);
    // INPUTS
    bool is_uint8 = image.type() == GPU_UINT8;
    program->set_input_float("f1", red, is_uint8);
    program->set_input_float("f2", green, is_uint8);
    program->set_input_float("f3", blue, is_uint8);
    program->set_input_float("f4", alpha, is_uint8);
    // DRAW
    int x1 = (int) image.x(bounds.x);
    int x2 = (int) x1 + bounds.width;
    int y1 = (int) image.y(bounds.y);
    int y2 = (int) (y1 + bounds.height);
    glBegin(GL_QUADS);
    glVertex2f(x1, y1);
    glVertex2f(x2, y1);
    glVertex2f(x2, y2);
    glVertex2f(x1, y2);
    glEnd();
    // CleanUp State
    program->uninstall();
    ShaderInvocation_CleanupGLState();
  }


  GPUImageBase threshold(const GPUImageBase& image, float thresh, float low, float high) {
    ShaderInvocation_SetupGLState(image);
    // Program - Install
    GPUProgram* program;
    switch(image.num_channels()) {
    case 4:
      program = create_gpu_program("Algorithms/threshold-1i3f-rgba"); break;
    case 3:
      program = create_gpu_program("Algorithms/threshold-1i3f-rgb"); break;
    case 1:
      program = create_gpu_program("Algorithms/threshold-1i3f-r"); break;
    }
    program->install();
    // OUTPUT
    GPUImageBase temp;
    temp.copy_attributes(image);
    ShaderInvocation_SetOutputImage(temp);
    // INPUT
    bool is_int8 = image.type() == GPU_UINT8;
    program->set_input_image("i1", image);
    program->set_input_float("f1", thresh, is_int8);
    program->set_input_float("f2", low, is_int8);
    program->set_input_float("f3", high, is_int8);
    // DRAW
    ShaderInvocation_DrawRectOneTexture(image);
    // CleanUp State
    program->uninstall();
    ShaderInvocation_CleanupGLState();

    return temp;
  }

} // namespace GPU
} // namespace vw

