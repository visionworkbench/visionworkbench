// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Core.h>

#include <vw/GPU/TexObj.h>
#include <vw/GPU/Setup.h>
#include <vw/GPU/TexAlloc.h>

namespace vw {
namespace GPU {

#define USE_PBO false

  //###################################################################
  //                     TexObj: Instance Functions
  //###################################################################

  TexObj::TexObj(int w, int h, Tex_Format internalFormat, Tex_Type internalType) {
    m_refCount = 0;
    GLuint format_gl;

    if(internalFormat == GPU_DEPTH) {
      if(internalType == GPU_FLOAT32)
        format_gl == GL_DEPTH_COMPONENT24;
      else if(internalType == GPU_FLOAT16)
        format_gl == GL_DEPTH_COMPONENT16;
      else
        format_gl == GL_DEPTH_COMPONENT;

      glGenRenderbuffersEXT(1, &texName);
      glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, format_gl, w, h);
      glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, GL_RENDERBUFFER_EXT, texName);
      return;
    }

    if(internalType == GPU_FLOAT32) {
      if(internalFormat == GPU_DEPTH)
        format_gl = GL_DEPTH_COMPONENT24;
      if(internalFormat == GPU_RGBA)
        format_gl = GL_FLOAT_RGBA32_NV;
      else if(internalFormat == GPU_RGB)
        format_gl = GL_FLOAT_RGB32_NV;
      else
        format_gl = GL_FLOAT_R32_NV;
    }
    else if(internalType == GPU_FLOAT16) {
      if(internalFormat == GPU_DEPTH)
        format_gl = GL_DEPTH_COMPONENT16;
      if(internalFormat == GPU_RGBA)
        format_gl = GL_FLOAT_RGBA16_NV;
      else if(internalFormat == GPU_RGB)
        format_gl = GL_FLOAT_RGB16_NV;
      else
        format_gl = GL_FLOAT_R16_NV;
    }

    else {
      if(internalFormat == GPU_DEPTH)
        format_gl = GL_DEPTH_COMPONENT;
      if(internalFormat == GPU_RGBA)
        format_gl = GL_RGBA;
      else if(internalFormat == GPU_RGB)
        format_gl = GL_RGB;
      else
        format_gl = GL_LUMINANCE;
    }

    m_width = w;
    m_height = h;
    m_format = internalFormat;
    m_type = internalType;

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
      m_width = 0;
      m_height = 0;
      throw(Exception("[vw::GPU::TexObj()] Error allocating texture."));
    }
  }

  TexObj::~TexObj() {
    if(texName) {
      glDeleteTextures(1, &texName);
    }
  }

  void TexObj::retain() {
    m_refCount++;
  }

  void TexObj::release() {
    m_refCount--;
    if(!m_refCount) TexAlloc::release(this);
  }

  void TexObj::write(int x, int y, int w, int h,
                     Tex_Format inputFormat,
                     Tex_Type inputType, void* data) {
    // FORMAT
    GLuint inputFormat_gl;
    if(inputFormat == GPU_RED && m_format == GPU_RED)
      inputFormat_gl = GL_LUMINANCE;
    else
      inputFormat_gl = inputFormat;
    // TYPE
    GLuint type_gl;
    if(inputType == GPU_FLOAT32)
      type_gl = GL_FLOAT;
    else if(inputType == GPU_FLOAT16) {
      type_gl = GL_FLOAT;
    }
    else
      type_gl = GL_UNSIGNED_BYTE;
    // WRITE
    if(USE_PBO) { // This code is supposed to speed up the CPU->GPU transfer, but it never quite seemed to work.
      int nComps;
      if(inputFormat == GPU_RGBA)
        nComps = 4;
      else if(inputFormat == GPU_RGB)
        nComps = 3;
      else
        nComps = 1;
      int texSize = w * h * inputType * nComps;
      GLuint pboBuffer;
      glGenBuffers(1, &pboBuffer);
      glBindBuffer(GL_PIXEL_UNPACK_BUFFER_ARB, pboBuffer);

      glBufferData(GL_PIXEL_UNPACK_BUFFER_ARB, texSize, NULL, GL_STREAM_DRAW);
      void* mappedTexPtr = glMapBuffer(GL_PIXEL_UNPACK_BUFFER_ARB, GL_WRITE_ONLY);
      memcpy(mappedTexPtr, data, texSize);
      glUnmapBuffer(GL_PIXEL_UNPACK_BUFFER_ARB);

      glBindTexture(GL_TEXTURE_RECTANGLE_ARB, texName);
      glTexSubImage2D(GL_TEXTURE_RECTANGLE_ARB, 0, x, y, w, h, inputFormat_gl, type_gl, 0);
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

  void TexObj::read(int x, int y, int w, int h,
                    Tex_Format outputFormat,
                    Tex_Type outputType, void* data) {

    // Format
    GLuint outFormat_gl;
    if(outputFormat == GPU_RED &&  m_format == GPU_RED)
      outFormat_gl = GL_LUMINANCE;
    else
      outFormat_gl = outputFormat;
    // Type
    GLuint outType_gl;
    if(outputType == GPU_FLOAT32)
      outType_gl = GL_FLOAT;
    else if(outputType == GPU_FLOAT16)
      outType_gl = GL_FLOAT;
    else
      outType_gl = GL_UNSIGNED_BYTE;

    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, g_framebuffer);
    glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_RECTANGLE_ARB, texName, 0);
    glDrawBuffer(GL_COLOR_ATTACHMENT0_EXT);
    glReadBuffer(GL_COLOR_ATTACHMENT0_EXT);
    if(!CheckFramebuffer(true))
      throw(Exception("[vw::GPU::TexObj::read()] Framebuffer error."));
    glGetError();
    glReadPixels(x, y, w, h, outFormat_gl, outType_gl, data);
    if(glGetError())
      throw(Exception("[vw::GPU::TexObj::read()] GL error from glReadPixels."));

    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
  }

}} // namespaces GPU, vw
