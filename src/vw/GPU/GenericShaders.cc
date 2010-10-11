// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/GPU/GenericShaders.h>
using std::vector;

namespace vw { namespace GPU {


  // ##################### GenericFragmentShader_1i0f Members #####################

  GenericFragmentShader_1i0f::GenericFragmentShader_1i0f(const char* fragmentShaderBaseName)
  { path = fragmentShaderBaseName; }


  GPUImageBase GenericFragmentShader_1i0f::operator()(const GPUImageBase& image1) {
    // Setup
    image1.rasterize_homography();
    ShaderInvocation_SetupGLState(image1.width(), image1.height());
    // Program
    GPUProgram* program = create_gpu_program(path);
    program->install();
    // Output
    GPUImageBase temp;
    temp.copy_attributes(image1);
    ShaderInvocation_SetOutputImage(temp);
    // Input
    program->set_input_image("i1", image1);
    // Drawing
    ShaderInvocation_DrawRectOneTexture(image1);
    // Cleanup
    program->uninstall();
    ShaderInvocation_CleanupGLState();

    return temp;
  }


  // ##################### GenericFragmentShader_2i0f Members #####################

  GenericFragmentShader_2i0f::GenericFragmentShader_2i0f(const char* fragmentShaderBaseName)
  { path = fragmentShaderBaseName; }

  GPUImageBase GenericFragmentShader_2i0f::operator()(const GPUImageBase& image1, const GPUImageBase& image2) {
    // GLState - Setup
    ((GPUImageBase&) image1).rasterize_homography();
    ((GPUImageBase&) image2).rasterize_homography();
    ShaderInvocation_SetupGLState(image1.width(), image1.height());
    // Program - Install
    GPUProgram* program = create_gpu_program(path);
    program->install();
    // OUTPUT
    GPUImageBase temp(image1.width(), image1.height(), image1.format(), image1.type());
    glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, temp.target(), temp.name(), 0);
    // INPUT
    program->set_input_image("i1", image1);
    program->set_input_image("i2", image2);
    // DRAW
    ShaderInvocation_DrawRectOneTexture(image1);
    // CleanUp State
    program->uninstall();
    ShaderInvocation_CleanupGLState();

    return temp;
  }


  // ##################### GenericFragmentShader_1i1f Members #####################

  GenericFragmentShader_1i1f::GenericFragmentShader_1i1f(const char* fragmentShaderBaseName)
  {
    path = fragmentShaderBaseName;
  }

  GPUImageBase GenericFragmentShader_1i1f::operator()(const GPUImageBase& image1, float float1) {
    // GLState - Setup
    ((GPUImageBase&) image1).rasterize_homography();
    ShaderInvocation_SetupGLState(image1.width(), image1.height());
    // Program - Install
    GPUProgram* program = create_gpu_program(path);
    program->install();
    // OUTPUT
    GPUImageBase temp(image1.width(), image1.height(), image1.format(), image1.type());
    glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, temp.target(), temp.name(), 0);
    // INPUT
    bool is_int8 = image1.type() == GPU_UINT8;
    program->set_input_image("i1", image1);
    program->set_input_float("f1", float1, is_int8);
    // DRAW
    ShaderInvocation_DrawRectOneTexture(image1);
    // CleanUp State
    program->uninstall();
    ShaderInvocation_CleanupGLState();

    return temp;
  }

  // ##################### GenericFragmentShader_1i2f Members #####################

  GenericFragmentShader_1i2f::GenericFragmentShader_1i2f(const char* fragmentShaderBaseName)
  { path = fragmentShaderBaseName; }


  GPUImageBase GenericFragmentShader_1i2f::operator()(const GPUImageBase& image1, float float1, float float2) {
    // Static
    static vector<int> fAttributes(1);
    static vector<int> emptyVector;
    // GLState - Setup
    ((GPUImageBase&) image1).rasterize_homography();
    ShaderInvocation_SetupGLState(image1.width(), image1.height());
    // Program - Install
    GPUProgram* program = create_gpu_program(path);
    program->install();
    // OUTPUT
    GPUImageBase temp(image1.width(), image1.height(), image1.format(), image1.type());
    glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, temp.target(), temp.name(), 0);
    // INPUT
    bool is_int8 = image1.type() == GPU_UINT8;
    program->set_input_image("i1", image1);
    program->set_input_float("f1", float1, is_int8);
    program->set_input_float("f2", float2, is_int8);
    // DRAW
    ShaderInvocation_DrawRectOneTexture(image1);
    // CleanUp State
    program->uninstall();
    ShaderInvocation_CleanupGLState();

    return temp;
  }

  // ##################### GenericFragmentShader_1i3f Members #####################

  GenericFragmentShader_1i3f::GenericFragmentShader_1i3f(const char* fragmentShaderBaseName)
  { path = fragmentShaderBaseName; }


  GPUImageBase GenericFragmentShader_1i3f::operator()(const GPUImageBase& image1, float float1, float float2, float float3) {
    // GLState - Setup
    ((GPUImageBase&) image1).rasterize_homography();
    ShaderInvocation_SetupGLState(image1.width(), image1.height());
    // Program - Install
    GPUProgram* program = create_gpu_program(path);
    program->install();
    // OUTPUT
    GPUImageBase temp(image1.width(), image1.height(), image1.format(), image1.type());
    glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, temp.target(), temp.name(), 0);
    // INPUT
    bool is_int8 = image1.type() == GPU_UINT8;
    program->set_input_image("i1", image1);
    program->set_input_float("f1", float1, is_int8);
    program->set_input_float("f2", float2, is_int8);
    program->set_input_float("f3", float3, is_int8);
    // DRAW
    ShaderInvocation_DrawRectOneTexture(image1);
    // CleanUp State
    program->uninstall();
    ShaderInvocation_CleanupGLState();

    return temp;
  }

  // ##################### GenericFragmentShader_1i4f Members #####################

  GenericFragmentShader_1i4f::GenericFragmentShader_1i4f(const char* fragmentShaderBaseName)
  { path = fragmentShaderBaseName; }

  GPUImageBase GenericFragmentShader_1i4f::operator()(const GPUImageBase& image1, float float1, float float2, float float3, float float4) {
    // GLState - Setup
    ((GPUImageBase&) image1).rasterize_homography();
    ShaderInvocation_SetupGLState(image1.width(), image1.height());
    // Program - Install
    GPUProgram* program = create_gpu_program(path);
    program->install();
    // OUTPUT
    GPUImageBase temp(image1.width(), image1.height(), image1.format(), image1.type());
    glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, temp.target(), temp.name(), 0);
    // INPUT
    bool is_int8 = image1.type() == GPU_UINT8;
    program->set_input_image("i1", image1);
    program->set_input_float("f1", float1, is_int8);
    program->set_input_float("f2", float2, is_int8);
    program->set_input_float("f3", float3, is_int8);
    program->set_input_float("f4", float4, is_int8);
    // DRAW
    ShaderInvocation_DrawRectOneTexture(image1);
    // CleanUp State
    program->uninstall();
    ShaderInvocation_CleanupGLState();

    return temp;
  }
} // namespace vw
} // namespace GPU


