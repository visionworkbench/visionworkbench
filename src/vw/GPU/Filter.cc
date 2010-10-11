// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/GPU/Filter.h>
#include <vector>

using std::vector;

namespace vw { namespace GPU {


//########################################################################
//#    convolution_filter
//########################################################################


GPUImageBase convolution_filter(const GPUImageBase& image,
                               const GPUImageBase& kernel)
{
// Static
  static vector<int> fAttributes(2);
// GLState - Setup
  ((GPUImageBase&) image).rasterize_homography();
  ShaderInvocation_SetupGLState(image.width(), image.height());
// Program - Install
  fAttributes[0] = kernel.width();
  fAttributes[1] = kernel.height();
  GPUProgram* program = create_gpu_program("Filter/convolution", fAttributes);
  program->install();
  // OUTPUT
  GPUImageBase temp(image.width(), image.height(), image.format(), image.type());
  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, temp.target(), temp.name(), 0);
  // INPUT
  program->set_input_image("image", image);
  program->set_input_image("kernel", kernel);
  program->set_input_float("hHalfSize", (kernel.width() / 2));
  program->set_input_float("vHalfSize", (kernel.height() / 2));
  // DRAW
  ShaderInvocation_DrawRectOneTexture(image);
  // CleanUp State
  program->uninstall();
  ShaderInvocation_CleanupGLState();

  return temp;
}


//########################################################################
//#    seperable_convolution_filter
//########################################################################

GPUImageBase seperable_convolution_filter(const GPUImageBase& image,
                                               const GPUImageBase& hKernel,
                                               const GPUImageBase& vKernel)
{
// Static
  static vector<int> fAttributes(1);
// GLState - Setup
  ((GPUImageBase&) image).rasterize_homography();
  GPUImageBase temp1(image.width(), image.height(), image.format(), image.type());
  GPUImageBase temp2(image.width(), image.height(), image.format(), image.type());
  GPUImageBase* input = (GPUImageBase*) &image;
  GPUImageBase* output = &temp1;
// *** STAGE 1 - Rows ***
  int h_kernel_size = hKernel.width();
  if(h_kernel_size) {
      ShaderInvocation_SetupGLState(input->width(), input->height());
          // Program - Install
          fAttributes[0] = h_kernel_size;
          GPUProgram* program = create_gpu_program("Filter/convolution-rows", fAttributes);
          program->install();
          // OUTPUT
          glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, output->target(), output->name(), 0);
          // INPUT
          program->set_input_image("image", *input);
          program->set_input_image("kernel", hKernel);
          program->set_input_float("halfSize", hKernel.width() / 2);
          // DRAW
          ShaderInvocation_DrawRectOneTexture(*input);
          program->uninstall();
          // SWAP Textures
          input = output;
          output = &temp2;
  }

// *** STAGE 2 - Columns ***
  int v_kernel_size = vKernel.width();
  if(v_kernel_size) {
      ShaderInvocation_SetupGLState(input->width(), input->height());
          // Program - Install
          fAttributes[0] = h_kernel_size;
          GPUProgram* program = create_gpu_program("Filter/convolution-rows", fAttributes);
          program->install();
          // OUTPUT
          glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, output->target(), output->name(), 0);
          // INPUT
          program->set_input_image("image", *input);
          program->set_input_image("kernel", hKernel);
          program->set_input_float("halfSize", hKernel.width() / 2);
          // DRAW
          ShaderInvocation_DrawRectOneTexture(*input);
          program->uninstall();
          // SWAP Textures
          GPUImageBase* old_input = input;
          input = output;
          output = old_input;
  }
  ShaderInvocation_CleanupGLState();
  return *input;
}

} } // namespaces GPU, vw
