// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef GPUProgram_H
#define GPUProgram_H

#include <map>
#include <cstdlib>
#include <cstdio>

#include <vw/config.h>
#include <vw/GPU/Setup.h>
#include <vw/GPU/GPUImage.h>


#if defined(VW_HAVE_PKG_CG) && VW_HAVE_PKG_CG==1
#include "Cg/cg.h"
#include "Cg/cgGL.h"
#endif

namespace vw { namespace GPU {


  extern std::string shader_base_path;

  extern bool shader_assembly_cache_enabled;

  extern std::string shader_assembly_cache_path;


  //########################################################################
  //#
  //########################################################################

  enum ShaderCompilationStatusEnum  {
    SHADER_COMPILATION_STATUS_SUCCESS = 0,
    SHADER_COMPILATION_STATUS_FILE_ERROR,
    SHADER_COMPILATION_STATUS_COMPILE_ERROR,
    SHADER_COMPILATION_STATUS_SUCCESS_FILE,
    SHADER_COMPILATION_STATUS_SUCCESS_CACHE,
    SHADER_COMPILATION_STATUS_ERROR_FILE,
    SHADER_COMPILATION_STATUS_ERROR_COMPILE,
    SHADER_COMPILATION_STATUS_ERROR_LINK
  };

  extern ShaderCompilationStatusEnum shaderCompilationStatus;

  inline ShaderCompilationStatusEnum get_shader_compilation_status() {
    return shaderCompilationStatus;
  }

  //########################################################################
  //#    Class: GPUProgram  (abstract)
  //########################################################################

  class GPUProgram {
  public:
    // Virtual
    virtual ~GPUProgram() { }
    virtual void install() = 0;
    virtual void uninstall() = 0;
    virtual void set_input_image(const char* name, const GPUImageBase& texture) = 0;
    virtual void set_input_float(const char* name, float value, bool is_uint8 = false) = 0;
  };

  //################################################################################################################
  //#    GLSL:    GPUVertexShader_GLSL, GPUFragmentShader_GLSL,  GPUProgram_GLSL (subclass), GPUProgramSet_GLSL (subclass)
  //################################################################################################################


  class GPUVertexShader_GLSL {
    GLhandleARB shader;
  public:
    bool compile(const std::string& vertexString);
    // Inline
    GPUVertexShader_GLSL() { shader = 0; }
    ~GPUVertexShader_GLSL() { if(shader) glDeleteObjectARB(shader); }

    bool is_compiled() { return shader != 0; }
    GLhandleARB get_shader() { return shader; }
  };


  class GPUFragmentShader_GLSL {
    GLhandleARB shader;
  public:
    bool compile(const std::string& fragmentString);
    // Inline
    GPUFragmentShader_GLSL() { shader = 0; }
    ~GPUFragmentShader_GLSL() { if(shader) glDeleteObjectARB(shader); }

    bool is_compiled() { return shader != 0; }
    GLhandleARB get_shader() { return shader; }

  };


  class GPUProgram_GLSL : public GPUProgram {
    GLhandleARB program;
    int bound_texture_count;
  public:
    bool Link(GPUVertexShader_GLSL& vertex, GPUFragmentShader_GLSL& fragment);
    // INLINE
    GPUProgram_GLSL() { program = 0; bound_texture_count = 0; }
    ~GPUProgram_GLSL() { if(program) glDeleteObjectARB(program); }
    GLhandleARB get_program() { return program; }
    // Over-Riden
    void install() {
      glUseProgramObjectARB(program);
      bound_texture_count = 0;
    }
    void uninstall() {
      glUseProgramObjectARB(0);
    }
    void set_input_image(const char* name, const GPUImageBase& texture) {
      glUniform1i(glGetUniformLocationARB(program, name), bound_texture_count);
      glGetError();
      glActiveTexture(GL_TEXTURE0 + bound_texture_count);
      if(glGetError() != GL_NO_ERROR)
        throw(Exception("[GPUProgram::set_input_image] Error: Too many textures bound."));
      texture.bind();
      bound_texture_count++;
    }
    void set_input_float(const char* name, float value, bool is_uint8) {
      if(is_uint8) value /= 255.0;
      glUniform1f(glGetUniformLocationARB(program, name), value);
    }
  };


  //#############################################################################################
  //#    CG Classes:  GPUShader_CG, GPUProgram_CG (subclass), GPUProgramSet_CG (subclass)
  //#############################################################################################

#if defined(VW_HAVE_PKG_CG) && VW_HAVE_PKG_CG==1

  class GPUShader_CG {
    // Variables - Member
    CGprogram program;
    CGprofile profile;
  public:
    // INLINE
    bool is_compiled() { return (program && cgIsProgramCompiled(program)); }
    CGprogram get_cg_program() { return program; }
    // MEMBER
    GPUShader_CG();
    ~GPUShader_CG();
    bool compile_source_with_string(const char *sourceString, CGprofile profile, const char *entry, const char **args);
    bool compile_source_with_file(const char *sourceFile, CGprofile profile, const char *entry, const char **args);

    bool save_compiled_file(const char *file);
    bool load_compiled_file(const char *file);

    void install();
    void uninstall();
  };


  class GPUProgram_CG : public GPUProgram {
    GPUShader_CG* vertexShader;
    GPUShader_CG* fragmentShader;
    int bound_texture_count;
  public:
    // INLINE
    GPUProgram_CG(GPUShader_CG* inVertexShader, GPUShader_CG* inFragmentShader)
      : vertexShader(inVertexShader), fragmentShader(inFragmentShader) { }
    ~GPUProgram_CG() {
      if(vertexShader) delete vertexShader;
      if(fragmentShader) delete fragmentShader;
    }
    // Over-Riden
    void install() {
      if(vertexShader)
        vertexShader->install();
      if(fragmentShader)
        fragmentShader->install();
      bound_texture_count = 0;
    }
    void uninstall() {
      if(vertexShader)
        vertexShader->uninstall();
      if(fragmentShader)
        fragmentShader->uninstall();
    }
    void set_input_image(const char* name, const GPUImageBase& texture) {
      if(vertexShader) {
        cgGLSetTextureParameter(cgGetNamedParameter(vertexShader->get_cg_program(), name), texture.name());
        cgGLEnableTextureParameter(cgGetNamedParameter(vertexShader->get_cg_program(), name));
      }
      if(fragmentShader) {
        cgGLSetTextureParameter(cgGetNamedParameter(fragmentShader->get_cg_program(), name), texture.name());
        cgGLEnableTextureParameter(cgGetNamedParameter(fragmentShader->get_cg_program(), name));
      }
      glGetError();
      glActiveTexture(GL_TEXTURE0 + bound_texture_count);
      if(glGetError() != GL_NO_ERROR)
        throw(Exception("[GPUProgram::set_input_image] Error: Too many textures bound."));
      texture.bind();
      bound_texture_count++;
    }
    void set_input_float(const char* name, float value, bool is_uint8) {
      if(is_uint8) value /= 255.0;
      if(vertexShader) cgGLSetParameter1f(cgGetNamedParameter(vertexShader->get_cg_program(), name), value);
      if(fragmentShader) cgGLSetParameter1f(cgGetNamedParameter(fragmentShader->get_cg_program(), name), value);
    }
  };

#endif // ifdef VW_HAVE_PKG_CG___



  //#############################################################################################
  //#    Creation Free Functions
  //#############################################################################################



  GPUProgram* create_gpu_program(const std::string& fragmentPath, const std::vector<int>& fragmentAttributes = std::vector<int>(),
                                 const std::string& vertexPath = "", const std::vector<int>& vertexAttributes = std::vector<int>());

  GPUProgram_GLSL* create_gpu_program_glsl_string(const std::string& fragmentString, const std::vector<int>& fragmentAttributes = std::vector<int>(),
                                                  const std::string& vertexString = "", const std::vector<int>& vertexAttributes = std::vector<int>());

  GPUProgram_GLSL* create_gpu_program_glsl(const std::string& fragmentPath, const std::vector<int>& fragmentAttributes = std::vector<int>(),
                                           const std::string& vertexPath = "", const std::vector<int>& vertexAttributes = std::vector<int>());

#if defined(VW_HAVE_PKG_CG) && VW_HAVE_PKG_CG==1

  GPUProgram_CG* create_gpu_program_cg_string(const std::string& fragmentString, const std::vector<int>& fragmentAttributes = std::vector<int>(),
                                              const std::string& vertexString = "", const std::vector<int>& vertexAttributes = std::vector<int>());

  GPUProgram_CG* create_gpu_program_cg(const std::string& fragmentPath, const std::vector<int>& fragmentAttributes = std::vector<int>(),
                                       const std::string& vertexPath = "", const std::vector<int>& vertexAttributes = std::vector<int>());

#endif

  void clear_gpu_program_cache();

} } // namespaces GPU, vw

#endif

