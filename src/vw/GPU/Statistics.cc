// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/GPU/Statistics.h>
#include <vw/GPU/Algorithms.h>
#include <vw/GPU/Manipulation.h>
#include <vw/GPU/ImageMath.h>
#include <cmath>

//#include <vw/Math.h>
//#include <vw/Geometry/VectorFixed.h>
//#include <vw/Geometry/PointND.h>

namespace vw { namespace GPU {

#define MIN(n1, n2) n1 < n2 ? n1 : n2
#define MAX(n1, n2) n1 > n2 ? n1 : n2


// *******************************************************************
//    min_channel_value()
// *******************************************************************

float min_channel_value(const GPUImageBase& image)
{
// Output Image Size
  int POT_Level = MAX((int) ceilf(log2f(image.width())), (int) ceilf(log2f(image.height())));
  int POT_Dimension = (int) powf(2.0, POT_Level);
// Reduce to 1 Channell; pad to POT dimensions; fill padded area with high value
  GPUImageBase tex_StartLevel(MAX(1, POT_Dimension), MAX(1, POT_Dimension), GPU_RED, image.type());
  GPUProgram* program_Level1;
  if(image.num_channels() == 4)
    program_Level1 = create_gpu_program("Statistics/min-channels-rgba");
  else if(image.num_channels() == 3)
    program_Level1 = create_gpu_program("Statistics/min-channels-rgb");
  if(image.num_channels() == 1)
    program_Level1 = create_gpu_program("Statistics/min-channels-r");
  program_Level1->install();

  ((GPUImageBase&) image).rasterize_homography();
  ShaderInvocation_SetupGLState(tex_StartLevel);
  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, tex_StartLevel.target(), tex_StartLevel.name(), 0);
  program_Level1->set_input_image("image", image);
  glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  int draw_size = image.width();
  ShaderInvocation_DrawRectOneTexture(image, Rectangle2D<int>(0, 0, image.width(), image.height()),
                              Rectangle2D<int>(0, 0, image.width(), image.height()));
  program_Level1->uninstall();
  // Fill edge extend area
  float high_value = 10000000.0;
  int xEdgeSize = tex_StartLevel.width() - image.width();
  int yEdgeSize = tex_StartLevel.height() - image.height();
  GPUImageBase tex_Cropped;
  if(xEdgeSize) {
    fill(tex_StartLevel, high_value, 0, 0, 0, Rectangle2D<int>(image.width(), 0, xEdgeSize, tex_StartLevel.height()));
  }
  if(yEdgeSize) {
    fill(tex_StartLevel, high_value, 0, 0, 0, Rectangle2D<int>(0, image.height(), image.width(), yEdgeSize));
  }
// ITERATIONS: 1 to n
  GPUImageBase tex_PreviousLevel = tex_StartLevel;
  float output;
  GPUProgram* program = create_gpu_program("Statistics/min-quad");

  for(int iLevel = 1; iLevel <= POT_Level; iLevel++) {
    program->install();
    int reductionRatio = (int) powf(2.0, iLevel);
    GPUImageBase tex_CurrentLevel(MAX(1, POT_Dimension / reductionRatio), MAX(1, POT_Dimension / reductionRatio), GPU_RED, image.type());
    ShaderInvocation_SetupGLState(tex_CurrentLevel);
    glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, tex_CurrentLevel.target(), tex_CurrentLevel.name(), 0);
    program->set_input_image("image", tex_PreviousLevel);
    glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_T, GL_CLAMP);
    ShaderInvocation_DrawRectOneTexture(tex_PreviousLevel, Vector2(0, 0), Vector2(tex_CurrentLevel.width(), tex_CurrentLevel.height()),
                                                                        Vector2(-0.5, -0.5), Vector2((tex_CurrentLevel.width() * 2) - 0.5, (tex_CurrentLevel.height() * 2) - 0.5));
    tex_PreviousLevel = tex_CurrentLevel;
  }
  program->uninstall();
  // CleanUp State
  ShaderInvocation_CleanupGLState();
  tex_PreviousLevel.read(GPU_RED, GPU_FLOAT32, &output);
  return output;
}

// *******************************************************************
//    max_channel_value()
// *******************************************************************

float max_channel_value(const GPUImageBase& image)
{
// Output Image Size
  int POT_Level = MAX((int) ceilf(log2f(image.width())), (int) ceilf(log2f(image.height())));
  int POT_Dimension = (int) powf(2.0, POT_Level);
// GLState - Setup
  ((GPUImageBase&) image).rasterize_homography();
  ShaderInvocation_SetupGLState(POT_Dimension, POT_Dimension);
// Reduce to 1 Channell; pad to POT dimensions; fill padded area with high value
  GPUProgram* program_Level1;
  if(image.num_channels() == 4)
    program_Level1 = create_gpu_program("Statistics/max-channels-rgba");
  else if(image.num_channels() == 3)
    program_Level1 = create_gpu_program("Statistics/max-channels-rgb");
  if(image.num_channels() == 1)
    program_Level1 = create_gpu_program("Statistics/max-channels-r");
  program_Level1->install();

  GPUImageBase tex_StartLevel(MAX(1, POT_Dimension), MAX(1, POT_Dimension), GPU_RED, image.type());
  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, tex_StartLevel.target(), tex_StartLevel.name(), 0);
  program_Level1->set_input_image("image", image);
  ShaderInvocation_DrawRectOneTexture(image, Rectangle2D<int>(0, 0, tex_StartLevel.width(), tex_StartLevel.height()),
                                      Rectangle2D<int>(0, 0, tex_StartLevel.width(), tex_StartLevel.height()));
  program_Level1->uninstall();
  // Fill edge extend area
  float pad_value = -10000000.0;
  int xEdgeSize = tex_StartLevel.width() - image.width();
  int yEdgeSize = tex_StartLevel.height() - image.height();
  GPUImageBase tex_Cropped;
  if(xEdgeSize) {
    fill(tex_StartLevel, pad_value, 0, 0, 0, Rectangle2D<int>(image.width(), 0, xEdgeSize, tex_StartLevel.height()));
  }
  if(yEdgeSize) {
    fill(tex_StartLevel, pad_value, 0, 0, 0, Rectangle2D<int>(0, image.height(), image.width(), yEdgeSize));
  }
// ITERATIONS: 1 to n
  GPUImageBase tex_PreviousLevel = tex_StartLevel;
  GPUProgram* program = create_gpu_program("Statistics/max-quad");

  for(int iLevel = 1; iLevel <= POT_Level; iLevel++) {
    program->install();
    int reductionRatio = (int) powf(2.0, iLevel);
    GPUImageBase tex_CurrentLevel(MAX(1, POT_Dimension / reductionRatio), MAX(1, POT_Dimension / reductionRatio), GPU_RED, image.type());
    ShaderInvocation_SetupGLState(tex_CurrentLevel);
    glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, tex_CurrentLevel.target(), tex_CurrentLevel.name(), 0);
    program->set_input_image("image", tex_PreviousLevel);
    glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_T, GL_CLAMP);
    ShaderInvocation_DrawRectOneTexture(tex_PreviousLevel, Vector2(0, 0), Vector2(tex_CurrentLevel.width(), tex_CurrentLevel.height()),
                                                                        Vector2(-0.5, -0.5), Vector2((tex_CurrentLevel.width() * 2) - 0.5, (tex_CurrentLevel.height() * 2) - 0.5));
    tex_PreviousLevel = tex_CurrentLevel;
  }
  program->uninstall();
  // CleanUp State
  ShaderInvocation_CleanupGLState();
  float output;
  tex_PreviousLevel.read(GPU_RED, GPU_FLOAT32, &output);
  return output;
}

// *******************************************************************
//    sum_channel_value()
// *******************************************************************

float sum_channel_value(const GPUImageBase& image)
{
// Output Image Size
  int POT_Level = MAX((int) ceilf(log2f(image.width())), (int) ceilf(log2f(image.height())));
  int POT_Dimension = (int) powf(2.0, POT_Level);
// GLState - Setup
        ((GPUImageBase&) image).rasterize_homography();
  ShaderInvocation_SetupGLState(POT_Dimension, POT_Dimension);
// Reduce to 1 Channell; pad to POT dimensions; fill padded area with high value
  GPUProgram* program_Level1;
  if(image.num_channels() == 4)
    program_Level1 = create_gpu_program("Statistics/sum-channels-rgba");
  else if(image.num_channels() == 3)
    program_Level1 = create_gpu_program("Statistics/sum-channels-rgb");
  if(image.num_channels() == 1)
    program_Level1 = create_gpu_program("Statistics/sum-channels-r");
  program_Level1->install();

  GPUImageBase tex_StartLevel(MAX(1, POT_Dimension), MAX(1, POT_Dimension), GPU_RED, image.type());
  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, tex_StartLevel.target(), tex_StartLevel.name(), 0);
  program_Level1->set_input_image("image", image);
  ShaderInvocation_DrawRectOneTexture(image, Rectangle2D<int>(0, 0, tex_StartLevel.width(), tex_StartLevel.height()),
                                      Rectangle2D<int>(0, 0, tex_StartLevel.width(), tex_StartLevel.height()));
  program_Level1->uninstall();
  // Fill edge extend area
  float pad_value = 0.0;
  int xEdgeSize = tex_StartLevel.width() - image.width();
  int yEdgeSize = tex_StartLevel.height() - image.height();
  GPUImageBase tex_Cropped;
  if(xEdgeSize) {
    fill(tex_StartLevel, pad_value, 0, 0, 0, Rectangle2D<int>(image.width(), 0, xEdgeSize, tex_StartLevel.height()));
  }
  if(yEdgeSize) {
    fill(tex_StartLevel, pad_value, 0, 0, 0, Rectangle2D<int>(0, image.height(), image.width(), yEdgeSize));
  }
// ITERATIONS: 1 to n
  GPUImageBase tex_PreviousLevel = tex_StartLevel;
  GPUProgram* program = create_gpu_program("Statistics/sum-quad");

  for(int iLevel = 1; iLevel <= POT_Level; iLevel++) {
    program->install();
    int reductionRatio = (int) powf(2.0, iLevel);
    GPUImageBase tex_CurrentLevel(MAX(1, POT_Dimension / reductionRatio), MAX(1, POT_Dimension / reductionRatio), GPU_RED, image.type());
    ShaderInvocation_SetupGLState(tex_CurrentLevel);
    glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, tex_CurrentLevel.target(), tex_CurrentLevel.name(), 0);
    program->set_input_image("image", tex_PreviousLevel);
    glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_T, GL_CLAMP);
    ShaderInvocation_DrawRectOneTexture(tex_PreviousLevel, Vector2(0, 0), Vector2(tex_CurrentLevel.width(), tex_CurrentLevel.height()),
                                                                        Vector2(-0.5, -0.5), Vector2((tex_CurrentLevel.width() * 2) - 0.5, (tex_CurrentLevel.height() * 2) - 0.5));
    tex_PreviousLevel = tex_CurrentLevel;
  }
  program->uninstall();
  // CleanUp State
  ShaderInvocation_CleanupGLState();
  PixelGray<float> output_pixel = ((GPUImage<PixelGray<float> >&) tex_PreviousLevel)(0,0);
  return output_pixel.v();
}

// *******************************************************************
//    stddev_channel_value()
// *******************************************************************

float stddev_channel_value(const GPUImageBase& image)
{
  float mean = mean_channel_value(image);
  float sum = sum_channel_value(pow(image - mean, 2.0));
  return sqrtf(sum) / (image.width() * image.height() * image.num_channels() - 1);
}


  } // namespace GPU
} // namespace vw
