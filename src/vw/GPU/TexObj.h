
#ifndef TexObj_H
#define TexObj_H

#include <list>
#include <vw/GPU/OpenGLManager.h>

namespace vw { namespace GPU {

#define TEX_DEBUG_PRINT 0
class TexBlock;
class TexAlloc;

//#################################################################################################################
//                                        Enums: Tex_Format and Tex_Type
//#################################################################################################################

enum Tex_Format {
  TEX_R = GL_RED,
  TEX_G = GL_GREEN,
  TEX_B = GL_BLUE,
  TEX_RGB = GL_RGB,
  TEX_RGBA = GL_RGBA,
  TEX_DEPTH
};

enum Tex_Type {
  TEX_UINT8 = 1,
  TEX_FLOAT16 = 2,
  TEX_FLOAT32 = 4
};

 inline const char* TexFormatToString(Tex_Format format) {
   if(format == TEX_R)
     return "TEX_R";
   else if(format == TEX_G)
     return "TEX_G";
   else if(format == TEX_R)
     return "TEX_R";
   else if(format == TEX_B)
     return "TEX_B";
   else if(format == TEX_RGB)
     return "TEX_RGB";
   else
     return "TEX_RGBA";
 }

 inline const char* TexTypeToString(Tex_Type type) {
   if(type == TEX_UINT8)
     return "TEX_UINT8";
   else if(type == TEX_FLOAT16)
     return "TEX_FLOAT16";
   else 
     return "TEX_FLOAT32";
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
  list<TexBlock*> texBlockList;
public:
  TexObj(int w, int h, Tex_Format internalFormat, Tex_Type internalType);

  ~TexObj();

  void retain(TexBlock* texBlock);

  void release(TexBlock* texBlock);
// INLINE
  GLuint name() { return texName; }

  Tex_Format format() { return _format; }

  Tex_Type type() { return _type; }

  int MemorySize() { 
    int nComponents; 
    if(_format == TEX_RGB) 
      nComponents = 3;
    else if(_format == TEX_RGBA) 
      nComponents = 4;
    else
      nComponents = 1;
    return _width * _height * _type * nComponents;
  }
  
  int Components() {
    if(_format == TEX_RGB) 
		return 3;
	else if(_format == TEX_RGBA) 
		return 4;
    else
		return 1;
  }

  int width() { return _width; }

  int height() { return _height; }

// VIRTUAL
  float x(float x);
  float y(float y);
  void write(int x, int y, int w, int h, Tex_Format inputFormat, Tex_Type inputType, void* data);
  void read(int x, int y, int w, int h, Tex_Format outputFormat, Tex_Type outputType, void* data);
  GLuint target();
  void bind();

};

} } // namespaces GPU, vw


#endif
