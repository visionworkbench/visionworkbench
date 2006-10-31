
#include "vw/GPU/GenericShaders.h"

namespace vw { namespace GPU {


  // ##################### GenericFragmentShader_1i0f Members #####################

     GenericFragmentShader_1i0f::GenericFragmentShader_1i0f(const char* fragmentShaderBaseName) 
     { programSet.set_base_paths("", fragmentShaderBaseName); }


  GPUImageBase GenericFragmentShader_1i0f::operator()(const GPUImageBase& image1) {
// Static
	static vector<int> fAttributes(1);
	static vector<int> emptyVector;
// GLState - Setup
    ((GPUImageBase&) image1).rasterize_homography();
	ShaderInvocation_SetupGLState(image1.width(), image1.height());
// Program - Install
	fAttributes[0] = 4;
	GPUProgram* program = programSet.get_program(emptyVector, fAttributes, false);
	program->install();
// OUTPUT
	GPUImageBase temp(image1.width(), image1.height(), image1.format(), image1.type());
	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, temp.target(), temp.name(), 0);	
// INPUT
	program->set_uniform_texture("i1", 0, image1);
// DRAW
	ShaderInvocation_DrawRectOneTexture(image1);
// CleanUp State
	program->uninstall();
	ShaderInvocation_CleanupGLState();

	return temp;
  }


  // ##################### GenericFragmentShader_2i0f Members #####################

     GenericFragmentShader_2i0f::GenericFragmentShader_2i0f(const char* fragmentShaderBaseName) 
     { programSet.set_base_paths("", fragmentShaderBaseName); }

  GPUImageBase GenericFragmentShader_2i0f::operator()(const GPUImageBase& image1, const GPUImageBase& image2) {
// Static
	static vector<int> fAttributes(1);
	static vector<int> emptyVector;
// GLState - Setup
	((GPUImageBase&) image1).rasterize_homography();	
	((GPUImageBase&) image2).rasterize_homography();	
	ShaderInvocation_SetupGLState(image1.width(), image1.height());
// Program - Install
	fAttributes[0] = 4;
	GPUProgram* program = programSet.get_program(emptyVector, fAttributes, false);
	program->install();
// OUTPUT
	GPUImageBase temp(image1.width(), image1.height(), image1.format(), image1.type());
	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, temp.target(), temp.name(), 0);	
// INPUT
	program->set_uniform_texture("i1", 0, image1);
	program->set_uniform_texture("i2", 0, image2);
// DRAW
	ShaderInvocation_DrawRectOneTexture(image1);
// CleanUp State
	program->uninstall();
	ShaderInvocation_CleanupGLState();

	return temp;
  }


  // ##################### GenericFragmentShader_1i1f Members #####################

     GenericFragmentShader_1i1f::GenericFragmentShader_1i1f(const char* fragmentShaderBaseName) 
     { programSet.set_base_paths("", fragmentShaderBaseName); }

  GPUImageBase GenericFragmentShader_1i1f::operator()(const GPUImageBase& image1, float float1) {
// Static
	static vector<int> fAttributes(1);
	static vector<int> emptyVector;
// GLState - Setup
	((GPUImageBase&) image1).rasterize_homography();	
	ShaderInvocation_SetupGLState(image1.width(), image1.height());
// Program - Install
	fAttributes[0] = 4;
	GPUProgram* program = programSet.get_program(emptyVector, fAttributes, false);
	program->install();
// OUTPUT
	GPUImageBase temp(image1.width(), image1.height(), image1.format(), image1.type());
	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, temp.target(), temp.name(), 0);	
// INPUT
	bool is_int8 = image1.type() == TEX_UINT8;
	program->set_uniform_texture("i1", 0, image1);
	program->set_uniform_float("f1", float1, is_int8);	
// DRAW
	ShaderInvocation_DrawRectOneTexture(image1);
// CleanUp State
	program->uninstall();
	ShaderInvocation_CleanupGLState();

	return temp;
  }

  // ##################### GenericFragmentShader_1i2f Members #####################

     GenericFragmentShader_1i2f::GenericFragmentShader_1i2f(const char* fragmentShaderBaseName) 
     { programSet.set_base_paths("", fragmentShaderBaseName); }

  GPUImageBase GenericFragmentShader_1i2f::operator()(const GPUImageBase& image1, float float1, float float2) {
// Static
	static vector<int> fAttributes(1);
	static vector<int> emptyVector;
// GLState - Setup
	((GPUImageBase&) image1).rasterize_homography();	
	ShaderInvocation_SetupGLState(image1.width(), image1.height());
// Program - Install
	fAttributes[0] = 4;
	GPUProgram* program = programSet.get_program(emptyVector, fAttributes, false);
	program->install();
// OUTPUT
	GPUImageBase temp(image1.width(), image1.height(), image1.format(), image1.type());
	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, temp.target(), temp.name(), 0);	
// INPUT
	bool is_int8 = image1.type() == TEX_UINT8;
	program->set_uniform_texture("i1", 0, image1);
	program->set_uniform_float("f1", float1, is_int8);
	program->set_uniform_float("f2", float2, is_int8);		
// DRAW
	ShaderInvocation_DrawRectOneTexture(image1);
// CleanUp State
	program->uninstall();
	ShaderInvocation_CleanupGLState();

	return temp;
  }

  // ##################### GenericFragmentShader_1i3f Members #####################

     GenericFragmentShader_1i3f::GenericFragmentShader_1i3f(const char* fragmentShaderBaseName) 
     { programSet.set_base_paths("", fragmentShaderBaseName); }


  GPUImageBase GenericFragmentShader_1i3f::operator()(const GPUImageBase& image1, float float1, float float2, float float3) {
// Static
	static vector<int> fAttributes(1);
	static vector<int> emptyVector;
// GLState - Setup
	((GPUImageBase&) image1).rasterize_homography();	
	ShaderInvocation_SetupGLState(image1.width(), image1.height());
// Program - Install
	fAttributes[0] = 4;
	GPUProgram* program = programSet.get_program(emptyVector, fAttributes, false);
	program->install();
// OUTPUT
	GPUImageBase temp(image1.width(), image1.height(), image1.format(), image1.type());
	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, temp.target(), temp.name(), 0);	
// INPUT
	bool is_int8 = image1.type() == TEX_UINT8;
	program->set_uniform_texture("i1", 0, image1);
	program->set_uniform_float("f1", float1, is_int8);
	program->set_uniform_float("f2", float2, is_int8);	
	program->set_uniform_float("f3", float3, is_int8);		
// DRAW
	ShaderInvocation_DrawRectOneTexture(image1);
// CleanUp State
	program->uninstall();
	ShaderInvocation_CleanupGLState();

	return temp;
  }

  // ##################### GenericFragmentShader_1i4f Members #####################

     GenericFragmentShader_1i4f::GenericFragmentShader_1i4f(const char* fragmentShaderBaseName) 
     { programSet.set_base_paths("", fragmentShaderBaseName); }

  GPUImageBase GenericFragmentShader_1i4f::operator()(const GPUImageBase& image1, float float1, float float2, float float3, float float4) {
// Static
	static vector<int> fAttributes(1);
	static vector<int> emptyVector;
// GLState - Setup
	((GPUImageBase&) image1).rasterize_homography();	
	ShaderInvocation_SetupGLState(image1.width(), image1.height());
// Program - Install
	fAttributes[0] = 4;
	GPUProgram* program = programSet.get_program(emptyVector, fAttributes, false);
	program->install();
// OUTPUT
	GPUImageBase temp(image1.width(), image1.height(), image1.format(), image1.type());
	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, temp.target(), temp.name(), 0);	
// INPUT
	bool is_int8 = image1.type() == TEX_UINT8;
	program->set_uniform_texture("i1", 0, image1);
	program->set_uniform_float("f1", float1, is_int8);
	program->set_uniform_float("f2", float2, is_int8);	
	program->set_uniform_float("f3", float3, is_int8);	
	program->set_uniform_float("f4", float4, is_int8);	
// DRAW
	ShaderInvocation_DrawRectOneTexture(image1);
// CleanUp State
	program->uninstall();
	ShaderInvocation_CleanupGLState();

	return temp;
  }
  } // namespace vw
} // namespace GPU


