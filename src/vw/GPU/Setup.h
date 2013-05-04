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


#ifndef Setup_H
#define Setup_H

#include <fstream>

#include <vw/GPU/Utilities.h>
#include <vw/GPU/Shaders.h>

#ifdef HAVE_CONFIG_H
#include <vw/config.h>
#endif


namespace vw { namespace GPU {

extern GLuint g_framebuffer;

extern bool loggingEnabled;

// Initialization

 void gpu_init(const std::string &log_path = std::string(), bool dummyWindow = true);

 void gpu_cleanup();

// Settings


  enum ShaderLanguageChoiceEnum {
    SHADER_LANGUAGE_CHOICE_CG_GLSL,
    SHADER_LANGUAGE_CHOICE_GLSL_CG,
    SHADER_LANGUAGE_CHOICE_GLSL,
    SHADER_LANGUAGE_CHOICE_CG
  };

  inline const char* get_string_for_shader_language_choice_enum(ShaderLanguageChoiceEnum choice) {
    switch(choice) {
    case SHADER_LANGUAGE_CHOICE_CG_GLSL:
      return "SHADER_LANGUAGE_CHOICE_CG_GLSL";
    case SHADER_LANGUAGE_CHOICE_GLSL_CG:
      return "SHADER_LANGUAGE_CHOICE_GLSL_CG";
    case SHADER_LANGUAGE_CHOICE_GLSL:
      return "SHADER_LANGUAGE_CHOICE_GLSL";
    case SHADER_LANGUAGE_CHOICE_CG:
      return "SHADER_LANGUAGE_CHOICE_CG";
    }
  }

  extern ShaderLanguageChoiceEnum shaderLanguageChoice;


 void set_shader_language_choice(ShaderLanguageChoiceEnum choice);

 void set_gpu_memory_recycling(bool value);

 void set_shader_base_path(const std::string& path);

 void set_shader_assembly_cache_path(const std::string& path);

// Logging

 void gpu_log(const char* str);

 static bool gpu_log_enabled() { return loggingEnabled; }



} } // namespace GPU, namespace vw

#endif

