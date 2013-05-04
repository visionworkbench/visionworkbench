// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NASA Vision Workbench is licensed under the Apache License,
//  Version 2.0 (the "License"); you may not use this file except in
//  compliance with the License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
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
