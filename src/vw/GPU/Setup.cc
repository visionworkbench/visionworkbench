
#include <vw/GPU/Setup.h>
#include <vw/GPU/GPUProgram.h>
#include <vw/GPU/TexAlloc.h>

namespace vw { namespace GPU {


  // Globals 

  GLuint g_framebuffer;

  ofstream gpuLogFile;

  static bool loggingEnabled;

  ShaderLanguageChoiceEnum shaderLanguageChoice = SHADER_LANGUAGE_CHOICE_CG_GLSL;

  // gpu_init

  void glut_display_dummy() { }

  void gpu_init(bool dummyWindow, bool isLoggingEnabled) {
    init_standard_shaders();
    static bool isInit = false;
    if(!isInit) {
      isInit = true;
      loggingEnabled = isLoggingEnabled;
      // GLUT Init
      int argc = 1;
      char* argv[1];
      argv[0] = "";
      glutInit(&argc, argv);
      printf("\nGlut init succeeded...\n");
      if(dummyWindow) {
	glutInitWindowSize(100, 100);
	glutInitDisplayMode(GLUT_RGBA);
	glutCreateWindow("GLUT Dummy Window");
	glutDisplayFunc(glut_display_dummy);
#ifndef __APPLE__
	glewInit();
#endif
	printf("Glut context creation succeeded...\n");
      }
      glEnable(GL_TEXTURE_2D);
      glEnable(GL_TEXTURE_RECTANGLE_ARB);
      printf("GL Texture Setup Succeeded\n");

      glDisable(GL_DITHER);
      glDisable(GL_BLEND);
      printf("GL Setup Succeeded\n");

      glGenFramebuffersEXT(1, &g_framebuffer);
      printf("GL Framebuffer creation Succeeded\n");

      if(loggingEnabled)
	gpuLogFile.open("GPU.log");
    }
    printf("GPU Init succeeded...\n");

  }

  // gpu_cleanup

  void gpu_cleanup() {
    clear_gpu_program_cache();
  }


  // Setup

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
