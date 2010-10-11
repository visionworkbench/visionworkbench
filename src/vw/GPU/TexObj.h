// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef TexObj_H
#define TexObj_H

#include <list>
#include <vw/GPU/Setup.h>

namespace vw {
namespace GPU {

#define GPU_DEBUG_PRINT 0
class TexBlock;
class TexAlloc;

//#################################################################################################################
//                                        Enums: Tex_Format and Tex_Type
//#################################################################################################################

enum Tex_Format {
  GPU_RED = GL_RED,
  GPU_GREEN = GL_GREEN,
  GPU_BLUE = GL_BLUE,
  GPU_RGB = GL_RGB,
  GPU_RGBA = GL_RGBA,
  GPU_DEPTH
};

enum Tex_Type {
  GPU_UINT8 = 1,
  GPU_FLOAT16 = 2,
  GPU_FLOAT32 = 4
};

  inline const char* TexFormatToString(Tex_Format format) {
    if(format == GPU_RED)
     return "GPU_RED";
   else if(format == GPU_GREEN)
     return "GPU_GREEN";
   else if(format == GPU_RED)
     return "GPU_RED";
   else if(format == GPU_BLUE)
     return "GPU_BLUE";
   else if(format == GPU_RGB)
     return "GPU_RGB";
   else
     return "GPU_RGBA";
 }

 inline const char* TexTypeToString(Tex_Type type) {
   if(type == GPU_UINT8)
     return "GPU_UINT8";
   else if(type == GPU_FLOAT16)
     return "GPU_FLOAT16";
   else
     return "GPU_FLOAT32";
 }

//#################################################################################################################
//                                               Class: TexObj
//#################################################################################################################

class TexObj {
  friend class TexAlloc;
protected:
  GLuint texName;
  Tex_Type m_type;
  Tex_Format m_format;
  int m_width;
  int m_height;
  int m_refCount;
public:
  TexObj(int w, int h, Tex_Format internalFormat, Tex_Type internalType);

  ~TexObj();

  void retain();

  void release();
// INLINE
  GLuint name() { return texName; }

  Tex_Format format() { return m_format; }

  Tex_Type type() { return m_type; }

  int MemorySize() {
    int nComponents;
    if(m_format == GPU_RGB)
      nComponents = 3;
    else if(m_format == GPU_RGBA)
      nComponents = 4;
    else
      nComponents = 1;
    return m_width * m_height * m_type * nComponents;
  }

  int Components() {
    if(m_format == GPU_RGB)
      return 3;
    else if(m_format == GPU_RGBA)
      return 4;
    else
      return 1;
  }

  int width() { return m_width; }

  int height() { return m_height; }

  int real_width() { return m_width; }

  int real_height() { return m_height; }

  float x(float x) { return x; }

  float y(float y) { return y; }

  GLuint target() { return GL_TEXTURE_RECTANGLE_ARB; }

  void bind() { glBindTexture(GL_TEXTURE_RECTANGLE_ARB, texName); }

// MEMBER
  void write(int x, int y, int w, int h, Tex_Format inputFormat, Tex_Type inputType, void* data);

  void read(int x, int y, int w, int h, Tex_Format outputFormat, Tex_Type outputType, void* data);

};

} } // namespaces GPU, vw


#endif
