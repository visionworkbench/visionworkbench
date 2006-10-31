 
#include <vw/Exception.h>

#include <vw/GPU/TexObj.h>
#include <vw/GPU/OpenGLManager.h>
#include <vw/GPU/TexAlloc.h>
#include <vw/GPU/TexBlock.h>


using namespace std;

namespace vw { namespace GPU {

#define USE_PBO false

//#################################################################################################################
//                                        TexObj: Instance Functions
//#################################################################################################################

void 
TexObj::retain(TexBlock* texBlock) {
  texBlockList.push_back(texBlock);
}


void 
TexObj::release(TexBlock* texBlock) {
  texBlockList.remove(texBlock);
  if(!texBlockList.size()) {
    TexAlloc::release(this);
  }
}


TexObj::TexObj(int w, int h, Tex_Format internalFormat, Tex_Type internalType) {
  GLuint format_gl;
  if(internalFormat == TEX_DEPTH) {
    if(internalType == TEX_FLOAT32)
      format_gl == GL_DEPTH_COMPONENT24;
    else if(internalType == TEX_FLOAT16)
      format_gl == GL_DEPTH_COMPONENT16;
    else
      format_gl == GL_DEPTH_COMPONENT;
    
    glGenRenderbuffersEXT(1, &texName);
    glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, format_gl, w, h);
    glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, GL_RENDERBUFFER_EXT, texName);
    return;
  }

  if(internalType == TEX_FLOAT32) {
    if(internalFormat == TEX_DEPTH)
      format_gl = GL_DEPTH_COMPONENT24;
    if(internalFormat == TEX_RGBA)
      format_gl = GL_FLOAT_RGBA32_NV;
    else if(internalFormat == TEX_RGB)
      format_gl = GL_FLOAT_RGB32_NV;
    else 
      format_gl = GL_FLOAT_R32_NV;
  }
  else if(internalType == TEX_FLOAT16) {
    if(internalFormat == TEX_DEPTH)
      format_gl = GL_DEPTH_COMPONENT16;
    if(internalFormat == TEX_RGBA)
      format_gl = GL_FLOAT_RGBA16_NV;
    else if(internalFormat == TEX_RGB)
      format_gl = GL_FLOAT_RGB16_NV;
    else 
      format_gl = GL_FLOAT_R16_NV;
  }
  else {
    if(internalFormat == TEX_DEPTH)
      format_gl = GL_DEPTH_COMPONENT;
    if(internalFormat == TEX_RGBA)
      format_gl = GL_RGBA;
    else if(internalFormat == TEX_RGB)
      format_gl = GL_RGB;
    else 
      format_gl = GL_LUMINANCE;
  }  
//
  _width = w;
  _height = h;
  _format = internalFormat;
  _type = internalType;

  glGenTextures(1, &texName);
  glBindTexture(GL_TEXTURE_RECTANGLE_ARB, texName);
  glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER_ARB);
  glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER_ARB);
  GLfloat black_color[4] = { 0, 0, 0, 1 };
  glTexParameterfv(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_BORDER_COLOR, black_color);

  glGetError();
  glTexImage2D(GL_TEXTURE_RECTANGLE_ARB, 0, format_gl, w, h, 0, GL_LUMINANCE, GL_FLOAT, NULL); 
  if(glGetError()) {
    _width = 0;
    _height = 0;
    throw(Exception("[vw::GPU::TexObj()] Error allocating texture."));
  }
}

TexObj::~TexObj() {
  if(texName) {
    glDeleteTextures(1, &texName);
  }
}

float TexObj::x(float x) {
	return x;
}

float TexObj::y(float y) {
	return y;
}


void 
TexObj::write(int x, int y, int w, int h, Tex_Format inputFormat, Tex_Type inputType, void* data) {
// FORMAT
  GLuint inputFormat_gl;
  if(inputFormat == TEX_R && _format == TEX_R) 
    inputFormat_gl = GL_LUMINANCE;
  else
    inputFormat_gl = inputFormat;  
// TYPE
  GLuint type_gl;
  if(inputType == TEX_FLOAT32)
    type_gl = GL_FLOAT;
  else if(inputType == TEX_FLOAT16) {
    type_gl = GL_FLOAT;
  }
  else
    type_gl = GL_UNSIGNED_BYTE;  
// WRITE
  if(USE_PBO) { // This code is supposed to speed up the CPU->GPU transfer, but it never quite seemed to work.
		int nComps;
		if(inputFormat == TEX_RGBA)	 
			nComps = 4;
		else if(inputFormat == TEX_RGB)	 
			nComps = 3;
		else
			nComps = 1;
		int texSize = w * h * inputType * nComps;
		GLuint pboBuffer;
		glGenBuffers(1, &pboBuffer);
		glBindBuffer(GL_PIXEL_UNPACK_BUFFER_ARB, pboBuffer);
		ClearGLError();

		glBufferData(GL_PIXEL_UNPACK_BUFFER_ARB, texSize, NULL, GL_STREAM_DRAW);
		CheckGLError("glBufferData");
		void* mappedTexPtr = glMapBuffer(GL_PIXEL_UNPACK_BUFFER_ARB, GL_WRITE_ONLY);
		CheckGLError("glMapBuffer");
		memcpy(mappedTexPtr, data, texSize);
		glUnmapBuffer(GL_PIXEL_UNPACK_BUFFER_ARB);
		CheckGLError("glUnmapBuffer");
		
		glBindTexture(GL_TEXTURE_RECTANGLE_ARB, texName);
		ClearGLError();
		glTexSubImage2D(GL_TEXTURE_RECTANGLE_ARB, 0, x, y, w, h, inputFormat_gl, type_gl, 0);
		CheckGLError("TexObj::write(PBO)");
		glBindBuffer(GL_PIXEL_UNPACK_BUFFER_ARB, 0);
		glDeleteBuffers(1, &pboBuffer);
	}
	{
	        glBindTexture(GL_TEXTURE_RECTANGLE_ARB, texName);
		glGetError();
		glTexSubImage2D(GL_TEXTURE_RECTANGLE_ARB, 0, x, y, w, h, inputFormat_gl, type_gl, data);
		if(glGetError())
		  throw(Exception("[vw::GPU::TexObj::write()] gl error from glTexSubImage2D."));
	}
}


void 
TexObj::read(int x, int y, int w, int h, Tex_Format outputFormat, Tex_Type outputType, void* data) {
	
// Format
  GLuint outFormat_gl;
  if(outputFormat == TEX_R &&  _format == TEX_R) 
    outFormat_gl = GL_LUMINANCE;
  else
    outFormat_gl = outputFormat;
// Type
  GLuint outType_gl;
  if(outputType == TEX_FLOAT32)
    outType_gl = GL_FLOAT;
  else if(outputType == TEX_FLOAT16)
    outType_gl = GL_FLOAT;
  else
    outType_gl = GL_UNSIGNED_BYTE;  

  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, g_framebuffer);
  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_RECTANGLE_ARB, texName, 0);
  glDrawBuffer(GL_COLOR_ATTACHMENT0_EXT);
  glReadBuffer(GL_COLOR_ATTACHMENT0_EXT);
  CheckFramebuffer(true);

  glGetError();
  glReadPixels(x, y, w, h, outFormat_gl, outType_gl, data);
  if(glGetError())
    throw(Exception("[vw::GPU::TexObj::read()] gl error from glReadPixels."));

  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
}

GLuint 
TexObj::target() {
  return GL_TEXTURE_RECTANGLE_ARB;
}

void
TexObj::bind() {
  glBindTexture(GL_TEXTURE_RECTANGLE_ARB, texName);	
}



} } // namespaces GPU, vw
