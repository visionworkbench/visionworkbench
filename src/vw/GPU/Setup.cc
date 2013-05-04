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


#include <vw/GPU/Setup.h>
#include <vw/GPU/GPUProgram.h>
#include <vw/GPU/TexAlloc.h>

using std::string;

namespace vw { namespace GPU {


  // Globals

  GLuint g_framebuffer;

  std::ofstream gpuLogFile;

  bool loggingEnabled;

  ShaderLanguageChoiceEnum shaderLanguageChoice = SHADER_LANGUAGE_CHOICE_GLSL;

  // gpu_init

  void glut_display_dummy() { }

  void gpu_init(const string &log_path, bool dummyWindow) {
    // Check if already init
    static bool isInit = false;
    if(isInit)
      return;
    else
      isInit = true;
    // Setup logging
    if(!log_path.empty()) {
      loggingEnabled = true;
      gpuLogFile.open(log_path.c_str());
      if(!gpuLogFile) {
        loggingEnabled = false;
        throw(Exception("[gpu_init] Error creating log file."));
      }
    }
    else {
      loggingEnabled = false;
    }
    // Init
    init_standard_shaders();

    gpu_log("Trying glutInit...");
    int argc = 1;
    char* argv[1] = { "" };

    glutInit(&argc, argv);
    gpu_log("Success\n");

    if(dummyWindow) {
      gpu_log("Trying to create a GLUT window...");
      glutInitWindowSize(10, 10);
      glutInitDisplayMode(GLUT_RGBA);
      glutCreateWindow("GLUT Dummy Window");
      glutDisplayFunc(glut_display_dummy);
      gpu_log("Success\n");

#ifndef __APPLE__
      gpu_log("Trying glewInit...");
      glewInit();
      gpu_log("Success\n");
#endif
    }

    gpu_log("Trying to setup GL state...");
    glEnable(GL_TEXTURE_2D);
    glEnable(GL_TEXTURE_RECTANGLE_ARB);
    glDisable(GL_DITHER);
    glDisable(GL_BLEND);
    gpu_log("Success\n");

    gpu_log("Trying to create a GL framebuffer...");
    glGenFramebuffersEXT(1, &g_framebuffer);
    gpu_log("Success\n");

  }

  // gpu_cleanup

  void gpu_cleanup() {
    clear_gpu_program_cache();
    set_gpu_memory_recycling(false);
  }

  // Settings

  void set_shader_language_choice(ShaderLanguageChoiceEnum choice) {
    shaderLanguageChoice = choice;
  }


  void set_shader_base_path(const string& path) {
    shader_base_path = path;
  }

  void set_shader_assembly_cache_path(const string& path) {
    shader_assembly_cache_enabled = true;
    shader_assembly_cache_path = path;
  }

  void set_gpu_memory_recycling(bool value) {
    TexAlloc::set_recycling(value);
  }

  void gpu_log(const char* string) {
    if(loggingEnabled) {
      gpuLogFile << string;
      gpuLogFile.flush();
    }
  }

} } // namespace GPU, namespace vw
