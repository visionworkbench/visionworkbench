
#ifndef TexAlloc_H
#define TexAlloc_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <list>
#include <map>

#include <vw/GPU/TexObj.h>

namespace vw { namespace GPU {

class TexAlloc {
  // Class Variables
  static bool isInit;
  static std::list<TexObj*> texRecycleList;
  static int allocatedCount;
  static int allocatedSize;
  static bool recylingEnabled;
  static std::map<std::pair<Tex_Format, Tex_Type>, std::pair<Tex_Format, Tex_Type> > textureSubstitutesMap;
 public:
  // Class Functions
  static TexBlock* alloc(int w, int h, Tex_Format format, Tex_Type type);
  static void release(TexObj* texObj);
  static void clear_recycled();
  static void generate_texture_substitutions(bool verbose = false);
  static bool get_texture_substitution(Tex_Format inFormat, Tex_Type inType, Tex_Format& outFormat, Tex_Type& outType);\
  // Class Functions - Inline
  static void set_recycling(bool value) { recylingEnabled = value; if(!value) clear_recycled(); }
  static int get_allocated_count() { return allocatedCount; }
  static int get_allocated_size() { return allocatedSize; }
private:
  // Class Functions - Private
  static void _init();
};

} } // namespaces GPU, vw


#endif
