// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_GPU_TEXALLOC_H__
#define __VW_GPU_TEXALLOC_H__

#include <list>
#include <map>

#include <vw/GPU/TexObj.h>

namespace vw {
namespace GPU {

  class TexAlloc {

    // Class Variables
    static bool isInit;
    static std::list<TexObj*> texRecycleList;
    static int allocatedCount;
    static int allocatedSize;
    static bool recylingEnabled;
    static std::map<std::pair<Tex_Format, Tex_Type>, std::pair<Tex_Format, Tex_Type> > textureSubstitutesMap;

    // Class Functions - Private
    static void initialize_texalloc();

  public:

    // Class Functions
    static TexObj* alloc(int w, int h, Tex_Format format, Tex_Type type);
    static void release(TexObj* texObj);
    static void clear_recycled();
    static void generate_texture_substitutions(bool verbose = false);
    static bool get_texture_substitution(Tex_Format inFormat, Tex_Type inType, Tex_Format& outFormat, Tex_Type& outType); \

    // Class Functions - Inline
    static void set_recycling(bool value) { recylingEnabled = value; if(!value) clear_recycled(); }
    static int get_allocated_count() { return allocatedCount; }
    static int get_allocated_size() { return allocatedSize; }

  private:
  };

}} // namespaces vw::GPU


#endif // __VW_GPU_TEXALLOC_H__
