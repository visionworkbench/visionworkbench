
#ifndef TexBlock_H
#define TexBlock_H

#include <vw/GPU/TexObj.h>

namespace vw { namespace GPU {

class TexBlock {
 private:
 // Instance Variables - Private
  TexObj* _texObj;
  int _width;
  int _height;
  int _xOffset;
  int _yOffset;
  int _refCount;
  Tex_Format _format;
  Tex_Type _type;
 public:  
 // Instance Functions
  TexBlock(TexObj* obj, int w, int h, int xOff, int yOff, Tex_Format inFormat, Tex_Type inType);
  ~TexBlock();
 // Instance Functions - Inline
  void retain() { _refCount++; }

  void release() { 
    _refCount--; 
    if(!_refCount) delete this; 
  }

  void write(int x, int y, int w, int h, Tex_Format inputFormat, Tex_Type inputType, void* data) 
    { _texObj->write(x + _xOffset, y + _yOffset, w, h, inputFormat, inputType, data); }

  void read(int x, int y, int w, int h, Tex_Format outputFormat, Tex_Type outputType, void* data) 
    { _texObj->read(x + _xOffset, y + _yOffset, w, h, outputFormat, outputType, data); }
 
  GLuint target() { return _texObj->target(); }

  GLuint name() { return _texObj->name(); }

  void bind() { _texObj->bind(); }

  float x(float x) { return _texObj->x(x + _xOffset); }

  float y(float y) { return _texObj->y(y + _yOffset); }

  int width() { return _width; }

  int height() { return _height; }

  int real_width() const { return _texObj->width(); }

  int real_height() const { return _texObj->height(); }

  bool same_real_object(const TexBlock& other) { return _texObj == other._texObj; }

  Tex_Format format() { return _format; } 

  Tex_Type type() { return _type; }

};

} } // namespaces GPU, vw



#endif
