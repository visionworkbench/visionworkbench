
#include "TexBlock.h"
#include "TexObj.h"


namespace vw { namespace GPU {

//#################################################################################################################
//                                        TexBlock: Instance Functions
//#################################################################################################################

TexBlock::TexBlock(TexObj* obj, int w, int h, int xOff, int yOff, Tex_Format inFormat, Tex_Type inType) {
  _refCount = 0;
  _width = w;
  _height = h;
  _xOffset = xOff; 
  _yOffset = yOff;
  _format = inFormat;
  _type = inType;
  _texObj = obj;
  _texObj->retain(this);
}
   
TexBlock::~TexBlock() { 
  if(TEX_DEBUG_PRINT) 
    printf("   *** ~TexBlock %p\n", this); 
  _texObj->release(this); 
}   
 
} } // namespaces GPU, vw
