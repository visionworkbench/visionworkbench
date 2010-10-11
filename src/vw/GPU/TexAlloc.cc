// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/GPU/TexAlloc.h>

namespace vw { namespace GPU {

  //############################################################
  //                  TexAlloc: Class Variables
  //############################################################

  bool TexAlloc::isInit = false;
  list<TexObj*> TexAlloc::texRecycleList;
  int TexAlloc::allocatedCount = 0;
  int TexAlloc::allocatedSize = 0;
  bool TexAlloc::recylingEnabled = false;
  std::map<pair<Tex_Format, Tex_Type>, pair<Tex_Format, Tex_Type> > TexAlloc::textureSubstitutesMap;

  //#############################################################
  //                  TexAlloc: Class Functions
  //#############################################################

  static char buffer[512];

  TexObj*
  TexAlloc::alloc(int w, int h, Tex_Format format, Tex_Type type) {
    // check for init
    if(!isInit)
      initialize_texalloc();

    // Attribute Substition
    Tex_Format realFormat = format;
    Tex_Type realType = type;
    get_texture_substitution(format, type, realFormat, realType);
    // Try to find it in recycling if enabled...
    TexObj* texObj = NULL;
    if(recylingEnabled) {
      TexObj* bestTex = NULL;
      TexObj* cTex;
      std::list<TexObj*>::iterator iter;
      for(iter = texRecycleList.begin(); iter != texRecycleList.end(); iter++) {
        cTex = *iter;
        if(cTex->format() == realFormat && cTex->type() == realType && cTex->width() == w && cTex->height() == h)
          bestTex = cTex;
      }
      if(bestTex) {
        texRecycleList.remove(bestTex);
        texObj = bestTex;
        if(gpu_log_enabled()) {
          sprintf(buffer, "+++ Creating Texture: { %s, %s, (%i x %i) } - RECYCLED from (%i x %i)    (Total Allocated: %.2fMB)\n",
                  TexFormatToString(format), TexTypeToString(type), w, h, bestTex->width(), bestTex->height(),
                  allocatedSize/1000000.0);
          gpu_log(buffer);
        }
      }
    }
    // if necessary, create new TexObj
    if(!texObj) {
      texObj = new TexObj(w, h, realFormat, realType);
      allocatedCount++;
      allocatedSize += texObj->MemorySize();

      if(gpu_log_enabled()) {
        if(format != realFormat || type != realType)
          sprintf(buffer, "+++ Creating Texture: { %s, %s, (%i x %i) } - NEW, Substituting {%s, %s }    (Total Allocated: %.2fMB)\n",
                  TexFormatToString(format), TexTypeToString(type), w, h, TexFormatToString(realFormat), TexTypeToString(realType), allocatedSize/1000000.0);
        else
          sprintf(buffer, "+++ Creating Texture: { %s, %s, (%i x %i) } - NEW    (Total of %.3fMB)\n",
                  TexFormatToString(format), TexTypeToString(type), w, h, allocatedSize/1000000.0);
        gpu_log(buffer);
      }
    }
    // return
    return texObj;
  }


  void
  TexAlloc::release(TexObj* texObj) {
    if(recylingEnabled) {
      texRecycleList.push_back(texObj);
      if(gpu_log_enabled()) {
        sprintf(buffer, "--- Recycling Texture: { %s, %s, (%i x %i) }    (Total Allocated %.2fMB)\n",
                TexFormatToString(texObj->format()), TexTypeToString(texObj->type()), texObj->width(), texObj->height(), allocatedSize/1000000.0);
        gpu_log(buffer);
      }
    }
    else {
      allocatedCount--;
      allocatedSize -= texObj->MemorySize();
      if(gpu_log_enabled()) {
        sprintf(buffer, "--- Deleting Texture: { %s, %s, (%i x %i) }    (Total Allocated %.2fMB)\n",
                TexFormatToString(texObj->format()), TexTypeToString(texObj->type()), texObj->width(), texObj->height(), allocatedSize/1000000.0);
        gpu_log(buffer);
      }
      delete texObj;
    }
  }

  void TexAlloc::clear_recycled() {
    std::list<TexObj*>::iterator iter;
    for(iter = texRecycleList.begin(); iter != texRecycleList.end(); iter++) {
      TexObj* texObj = *iter;
      allocatedCount--;
      allocatedSize -= texObj->MemorySize();
      if(gpu_log_enabled()) {
        sprintf(buffer, "--- Deleting Texture: { %s, %s, (%i x %i) }    (Total Allocated %.2fMB)\n",
                TexFormatToString(texObj->format()), TexTypeToString(texObj->type()), texObj->width(), texObj->height(), allocatedSize/1000000.0);
        gpu_log(buffer);
      }
      delete *iter;
    }
    texRecycleList.clear();
  }

  void TexAlloc::generate_texture_substitutions(bool verbose) {
    textureSubstitutesMap.clear();
    char buffer1[256];
    char buffer2[256];
    bool chart[3][3];
    static Tex_Format textureFormats[] = { GPU_RED, GPU_RGB, GPU_RGBA };
    static Tex_Type textureTypes[] = { GPU_UINT8, GPU_FLOAT16, GPU_FLOAT32 };
    for(int i_format=0; i_format < 3;  i_format++) {
      Tex_Format format = textureFormats[i_format];
      if(format == GPU_RED) sprintf(buffer1, "GPU_RED");
      else if(format == GPU_RGB) sprintf(buffer1, "GPU_RGB");
      else if(format == GPU_RGBA) sprintf(buffer1, "GPU_RGBA");

      for(int i_type=0; i_type < 3; i_type++) {
        bool failed = false;
        Tex_Type type = textureTypes[i_type];
        if(type == GPU_UINT8) sprintf(buffer2, "GPU_UINT8");
        else if(type == GPU_FLOAT16) sprintf(buffer2, "GPU_FLOAT16");
        else if(type == GPU_FLOAT32) sprintf(buffer2, "GPU_FLOAT32");
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
          if(format == GPU_RED) sprintf(buffer1, "GPU_RED");
          else if(format == GPU_RGB) sprintf(buffer1, "GPU_RGB");
          else if(format == GPU_RGBA) sprintf(buffer1, "GPU_RGBA");
          Tex_Type type = inType;
          if(type == GPU_UINT8) sprintf(buffer2, "GPU_UINT8");
          else if(type == GPU_FLOAT16) sprintf(buffer2, "GPU_FLOAT16");
          else if(type == GPU_FLOAT32) sprintf(buffer2, "GPU_FLOAT32");
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

  //###############################################################
  //               TexAlloc: Class Functions - Private
  //###############################################################

  void
  TexAlloc::initialize_texalloc() {
    generate_texture_substitutions();
    isInit = true;
  }

} } // namespaces GPU, vw
