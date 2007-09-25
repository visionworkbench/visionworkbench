
#ifndef TexObj_H
#define TexObj_H

#include <list>
#include <vw/GPU/Setup.h>

namespace vw { namespace GPU {

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
  Tex_Type _type;
  Tex_Format _format;
  int _width;
  int _height;
  int _refCount;
public:
  TexObj(int w, int h, Tex_Format internalFormat, Tex_Type internalType);

  ~TexObj();

  void retain();

  void release();
// INLINE
  GLuint name() { return texName; }

  Tex_Format format() { return _format; }

  Tex_Type type() { return _type; }

  int MemorySize() { 
    int nComponents; 
    if(_format == GPU_RGB) 
      nComponents = 3;
    else if(_format == GPU_RGBA) 
      nComponents = 4;
    else
      nComponents = 1;
    return _width * _height * _type * nComponents;
  }
  
  int Components() {
    if(_format == GPU_RGB) 
		return 3;
	else if(_format == GPU_RGBA) 
		return 4;
    else
		return 1;
  }

  int width() { return _width; }

  int height() { return _height; }

  int real_width() { return _width; }

  int real_height() { return _height; }

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
