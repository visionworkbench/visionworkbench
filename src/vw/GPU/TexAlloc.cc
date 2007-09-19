
#include <vw/GPU/TexAlloc.h>
#include <vw/GPU/TexBlock.h>

using namespace std;

namespace vw { namespace GPU {

//#################################################################################################################
//                                        TexAlloc: Class Variables
//#################################################################################################################

bool TexAlloc::isInit = false;
list<TexObj*> TexAlloc::texRecycleList;
int TexAlloc::allocatedCount = 0;
int TexAlloc::allocatedSize = 0;
bool TexAlloc::recylingEnabled = false;
map<pair<Tex_Format, Tex_Type>, pair<Tex_Format, Tex_Type> > TexAlloc::textureSubstitutesMap;

//#################################################################################################################
//                                        TexAlloc: Class Functions
//#################################################################################################################

static char buffer[256];

TexBlock* 
TexAlloc::alloc(int w, int h, Tex_Format format, Tex_Type type) {
// check for init
  if(!isInit)
    _init();
// Attribute Substition
  Tex_Format realFormat = format;
  Tex_Type realType = type;
  get_texture_substitution(format, type, realFormat, realType);
  bool recycledFound = false;
// Allocate
  TexObj* texObj;
  if(recylingEnabled) {
      TexObj* bestTex = NULL;
      TexObj* cTex;
      int cWidth;
      int cHeight;
      if(TEX_DEBUG_PRINT)	
	printf("   *** TexAlloc: %i items in recycle list.\n", texRecycleList.size());
      std::list<TexObj*>::iterator iter;
      int c=0;
      for(iter = texRecycleList.begin(); iter != texRecycleList.end(); iter++) {
	cTex = *iter;
	if(cTex->format() == realFormat && cTex->type() == realType && cTex->width() == w && cTex->height() == h)
	  bestTex = cTex;	  
      }
      if(bestTex) {
	if(TEX_DEBUG_PRINT)
	  printf("### TexAlloc: RECYCLE [%i, %i]  %p\n", w, h, bestTex);
	texRecycleList.remove(bestTex);
	texObj = bestTex;
	recycledFound = true;
	sprintf(buffer, "   ### TexAlloc: [%i x %i]. Format = %i.   Found in recycling (%i x %i)\n", w, h, (int) format, bestTex->width(), bestTex->height());
	gpu_log(buffer);
      }	  
    }
// if recycled not found, create new TexObj
  if(!recycledFound) {
    texObj = new TexObj(w, h, realFormat, realType);
    allocatedCount++;
    allocatedSize += texObj->MemorySize();
    if(TEX_DEBUG_PRINT)
      printf("### TexAlloc: NEW [%i, %i]  %p\n", w, h, texObj);
	if(format != realFormat || type != realType)
		sprintf(buffer, "   ### TexAlloc: [%i x %i] Texture Format/Type Substitution.  Created New\n", w, h, (int) format);
	else
		sprintf(buffer, "   ### TexAlloc: [%i x %i] Format = %i.  Created New - Total MB: %f\n", w, h, (int) format, allocatedSize / 1000000.0);
    gpu_log(buffer);
	
    if(!texObj->width()) {
      w = 0; h = 0;
    }
  }
// Create texBlock and return it
  TexBlock* texBlock = new TexBlock(texObj, w, h, 0, 0, format, type);
  return texBlock;		  
}



void 
TexAlloc::release(TexObj* texObj) {
  if(TEX_DEBUG_PRINT)
    printf("   *** TexAlloc::release %p\n", texObj);
  if(recylingEnabled)
    texRecycleList.push_back(texObj);
  else {
    allocatedCount--;
    allocatedSize -= texObj->MemorySize();
    sprintf(buffer, "   ### Tex DELETED: [%i x %i]. Format = %i. - Total MB: %f\n", texObj->width(), texObj->height(), (int) texObj->format(), allocatedSize / 1000000.0);
    gpu_log(buffer);
    delete texObj;
  }
}



void 
TexAlloc::clear_recycled() {
  std::list<TexObj*>::iterator iter;
  for(iter = texRecycleList.begin(); iter != texRecycleList.end(); iter++) {
    TexObj* texObj = *iter;
    allocatedCount--;
    allocatedSize -= texObj->MemorySize();
    sprintf(buffer, "   ### Tex DELETED: [%i x %i]. Format = %i. - Total MB: %f\n", 
	    texObj->width(), texObj->height(), (int) texObj->format(), allocatedSize / 1000000.0);
    gpu_log(buffer);
    delete *iter;
  }
  texRecycleList.clear();
}


void TexAlloc::generate_texture_substitutions(bool verbose) {
  textureSubstitutesMap.clear();
  char buffer1[256];
  char buffer2[256];
  bool chart[3][3];
  static Tex_Format textureFormats[] = { TEX_R, TEX_RGB, TEX_RGBA };
  static Tex_Type textureTypes[] = { TEX_UINT8, TEX_FLOAT16, TEX_FLOAT32 };
  for(int i_format=0; i_format < 3;  i_format++) {
    Tex_Format format = textureFormats[i_format];
    if(format == TEX_R) sprintf(buffer1, "TEX_R");
    else if(format == TEX_RGB) sprintf(buffer1, "TEX_RGB");
    else if(format == TEX_RGBA) sprintf(buffer1, "TEX_RGBA");
    for(int i_type=0; i_type < 3; i_type++) {
      bool failed = false;
      Tex_Type type = textureTypes[i_type];
      if(type == TEX_UINT8) sprintf(buffer2, "TEX_UINT8");
      else if(type == TEX_FLOAT16) sprintf(buffer2, "TEX_FLOAT16");
      else if(type == TEX_FLOAT32) sprintf(buffer2, "TEX_FLOAT32");
      if(verbose)
	printf("Checking Texture Support for:  format = %s, type = %s\n", buffer1, buffer2);
      auto_ptr<TexObj> texRef(new TexObj(100, 100, format, type));
      if(!texRef->width()) {
	failed = true;
	if(verbose)
	  printf("  *** Failed on create.\n");
      }
      if(!failed) {
	glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, g_framebuffer);
	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, texRef->target(), texRef->name(), 0);
	glDrawBuffer(GL_COLOR_ATTACHMENT0_EXT);
	glReadBuffer(GL_COLOR_ATTACHMENT0_EXT);
	if(!CheckFramebuffer(false)) {  
	  failed = true;
	  if(verbose)
	    printf("  *** Failed on bind to framebuffer.\n");
        }	
    }
    chart[i_format][i_type] = !failed;
  }
  }
  for(int i_format=0; i_format < 3;  i_format++) {
    for(int i_type=0; i_type < 3; i_type++) {
      bool found = false;
      Tex_Format inFormat = textureFormats[i_format];
      Tex_Type inType = textureTypes[i_type];
      Tex_Format outFormat;
      Tex_Type outType;
      for(int j_format=i_format; j_format < 3 && !found;  j_format++) {
	for(int j_type=i_type; j_type < 3; j_type++) {
	  if(chart[j_format][j_type]) {
	    found = true;
	    outFormat = textureFormats[j_format];
	    outType =  textureTypes[j_type];
	    textureSubstitutesMap[pair<Tex_Format, Tex_Type>(inFormat, inType)] = pair<Tex_Format, Tex_Type>(outFormat, outType);
	    if(verbose) {
	      printf("[%s, %s] implemented as: ", TexFormatToString(inFormat), TexTypeToString(inType));
	      printf("[%s, %s]\n", TexFormatToString(outFormat), TexTypeToString(outType));
	    }
	    break;
	  }
	}
      }
      if(verbose && !found) {
	Tex_Format format = inFormat;
	if(format == TEX_R) sprintf(buffer1, "TEX_R");
	else if(format == TEX_RGB) sprintf(buffer1, "TEX_RGB");
	else if(format == TEX_RGBA) sprintf(buffer1, "TEX_RGBA");
	Tex_Type type = inType;
	if(type == TEX_UINT8) sprintf(buffer2, "TEX_UINT8");
	else if(type == TEX_FLOAT16) sprintf(buffer2, "TEX_FLOAT16");
	else if(type == TEX_FLOAT32) sprintf(buffer2, "TEX_FLOAT32");
	printf("[%s, %s]  NOT IMPLEMENTED!\n", buffer1, buffer2);
      }
    }
  }
}

bool TexAlloc::get_texture_substitution(Tex_Format inFormat, Tex_Type inType, Tex_Format& outFormat, Tex_Type& outType) {
  map<pair<Tex_Format, Tex_Type>, pair<Tex_Format, Tex_Type> >::iterator iter;
  iter = textureSubstitutesMap.find(pair<Tex_Format, Tex_Type>(inFormat, inType));
  if(iter != textureSubstitutesMap.end()) {
    outFormat = (*iter).second.first;
    outType = (*iter).second.second;
    return true;
  }
  else
    return false;
}

//#################################################################################################################
//                                        TexAlloc: Class Functions - Private
//#################################################################################################################

void 
TexAlloc::_init() {
  // Replace with GL Extension query functions
  generate_texture_substitutions();
  isInit = true;
}

} } // namespaces GPU, vw
