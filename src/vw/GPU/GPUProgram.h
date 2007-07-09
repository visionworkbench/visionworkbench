

#ifndef GPUProgram_H
#define GPUProgram_H

#include <map>
#include <stdlib.h>
#include <stdio.h>

#include <vw/GPU/Setup.h>
#include <vw/GPU/GPUImage.h>


#ifdef HAVE_PKG_CG
  #include "Cg/cg.h"
  #include "Cg/cgGL.h"
#endif

using namespace std;

namespace vw { namespace GPU {



//########################################################################
//#                  
//########################################################################

  enum ShaderCompilationStatusEnum  {
    SHADER_COMPILATION_STATUS_SUCCESS = 0,
    SHADER_COMPILATION_STATUS_FILE_ERROR,
    SHADER_COMPILATION_STATUS_COMPILE_ERROR
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
  virtual void set_uniform_texture(const char* name, int id, const GPUImageBase& texture) = 0;
  virtual void set_uniform_float(const char* name, float value, bool is_uint8 = false) = 0;
};

//################################################################################################################
//#    GLSL:    GPUVertexShader_GLSL, GPUFragmentShader_GLSL,  GPUProgram_GLSL (subclass), GPUProgramSet_GLSL (subclass)      
//################################################################################################################


class GPUVertexShader_GLSL {
	GLhandleARB shader;
public:
	bool compile(string& vertexString);
// Inline	
	GPUVertexShader_GLSL() { shader = 0; }
	~GPUVertexShader_GLSL() { if(shader) glDeleteObjectARB(shader); }

	bool is_compiled() { return shader != 0; }
	GLhandleARB get_shader() { return shader; }
};


class GPUFragmentShader_GLSL {
	GLhandleARB shader;
public:
	bool compile(string& fragmentString);
// Inline	
	GPUFragmentShader_GLSL() { shader = 0; }
	~GPUFragmentShader_GLSL() { if(shader) glDeleteObjectARB(shader); }

	bool is_compiled() { return shader != 0; }
	GLhandleARB get_shader() { return shader; }
	
};


class GPUProgram_GLSL : public GPUProgram {
	GLhandleARB program;
public:
	bool Link(GPUVertexShader_GLSL& vertex, GPUFragmentShader_GLSL& fragment);
// INLINE
	GPUProgram_GLSL() { program = 0; }
	~GPUProgram_GLSL() { if(program) glDeleteObjectARB(program); }
	GLhandleARB get_program() { return program; }
// Over-Riden
	void install() { 
		glUseProgramObjectARB(program); 
	}
	void uninstall() {
		glUseProgramObjectARB(0); 
	}
	void set_uniform_texture(const char* name, int n, const GPUImageBase& texture) { 
		glUniform1i(glGetUniformLocationARB(program, name), n);
		glActiveTexture(GL_TEXTURE0 + n);
		texture.bind();
	}
	void set_uniform_float(const char* name, float value, bool is_uint8) {
	  if(is_uint8) value /= 255.0;
	  glUniform1f(glGetUniformLocationARB(program, name), value);	
	}
};


class GPUProgramSet_GLSL {
	map<pair<vector<int>, vector<int> >, GPUProgram_GLSL*> programMap;
	map<vector<int>, GPUVertexShader_GLSL*> vertexMap;
	map<vector<int>, GPUFragmentShader_GLSL*> fragmentMap;
	string vertexBasePath;
	string fragmentBasePath;
public:
	GPUProgramSet_GLSL();
	~GPUProgramSet_GLSL();
	
	GPUProgram_GLSL* get_program(const vector<int>& vertexAttributes, 
				    const vector<int>& fragmentAttributes, 
				    bool verbose = false );
	  
// Inline
	void set_base_paths(const string& inVertexBasePath, const string& inFragmentBasePath) {
		vertexBasePath = inVertexBasePath;
		fragmentBasePath = inFragmentBasePath;
	}
	void set_base_paths(const char* inVertexBasePath, const char* inFragmentBasePath) {
		vertexBasePath = inVertexBasePath;
		fragmentBasePath = inFragmentBasePath;
	}
};


//#############################################################################################
//#    CG Classes:  GPUShader_CG, GPUProgram_CG (subclass), GPUProgramSet_CG (subclass)         
//#############################################################################################


#ifdef HAVE_PKG_CG

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
public:
// INLINE
  GPUProgram_CG(GPUShader_CG* inVertexShader, GPUShader_CG* inFragmentShader) 
    : vertexShader(inVertexShader), fragmentShader(inFragmentShader) { }
  ~GPUProgram_CG() { if(vertexShader) delete vertexShader; if(fragmentShader) delete fragmentShader; }
// Over-Riden
  void install() { 
    if(vertexShader)
      vertexShader->install();
    if(fragmentShader)
      fragmentShader->install();
  }
  void uninstall() {
    if(vertexShader)
      vertexShader->uninstall();
    if(fragmentShader)
      fragmentShader->uninstall();
  }
  void set_uniform_texture(const char* name, int n, const GPUImageBase& texture) { 
    if(vertexShader) cgGLSetTextureParameter(cgGetNamedParameter(vertexShader->get_cg_program(), name), texture.name());
    if(fragmentShader) {
      cgGLSetTextureParameter(cgGetNamedParameter(fragmentShader->get_cg_program(), name), texture.name());
      cgGLEnableTextureParameter(cgGetNamedParameter(fragmentShader->get_cg_program(), name));
    }
    texture.bind();
  }
  void set_uniform_texture(int id, int n, const GPUImageBase& texture) { 
    //cgGLSetTextureParameter(id, inputTex.name());
  }
  void set_uniform_float(const char* name, float value, bool is_uint8) {
    if(is_uint8) value /= 255.0;
    if(vertexShader) cgGLSetParameter1f(cgGetNamedParameter(vertexShader->get_cg_program(), name), value);
    if(fragmentShader) cgGLSetParameter1f(cgGetNamedParameter(fragmentShader->get_cg_program(), name), value);
  }
};


class GPUProgramSet_CG {
	map<pair<vector<int>, vector<int> >, GPUProgram_CG*> programMap;
	string vertexBasePath;
	string fragmentBasePath;
public:
	GPUProgramSet_CG();
	~GPUProgramSet_CG();
	
	GPUProgram_CG* get_program(const vector<int>& vertexAttributes, 
				    const vector<int>& fragmentAttributes, 
				    bool verbose = false );
	  
// Inline
	void set_base_paths(const string& inVertexBasePath, const string& inFragmentBasePath) {
		vertexBasePath = inVertexBasePath;
		fragmentBasePath = inFragmentBasePath;
	}
	void set_base_paths(const char* inVertexBasePath, const char* inFragmentBasePath) {
		vertexBasePath = inVertexBasePath;
		fragmentBasePath = inFragmentBasePath;
	}
};

#endif // ifdef HAVE_PKG_CG___

//#############################################################################################
//#    Class:  GPUProgramSet   
//#############################################################################################

class GPUProgramSet {
  string vertexBasePath;
  string fragmentBasePath;
  static bool useAssemblyCaching;
  GPUProgramSet_GLSL programSet_GLSL;
#ifdef HAVE_PKG_CG
  GPUProgramSet_CG programSet_CG;
#endif  	
 public:
  static void set_use_assembly_caching(bool value) { useAssemblyCaching = value; }

  static bool get_use_assembly_caching() { return useAssemblyCaching; }

  GPUProgram* get_program(const vector<int>& vertexAttributes, 
			 const vector<int>& fragmentAttributes, 
			 bool verbose = false );
	  
// Inline
  GPUProgramSet();

  GPUProgramSet(const char* inVertexBasePath, const char* inFragmentBasePath) {
    set_base_paths(inVertexBasePath, inFragmentBasePath);
  }

  ~GPUProgramSet();
  void set_base_paths(const string& inVertexBasePath, const string& inFragmentBasePath) {
    vertexBasePath = inVertexBasePath;
    fragmentBasePath = inFragmentBasePath;
    programSet_GLSL.set_base_paths(inVertexBasePath, inFragmentBasePath);
#ifdef HAVE_PKG_CG
    programSet_CG.set_base_paths(inVertexBasePath, inFragmentBasePath);
#endif  	
  }

  void set_base_paths(const char* inVertexBasePath, const char* inFragmentBasePath) {
    vertexBasePath = inVertexBasePath;
    fragmentBasePath = inFragmentBasePath;
    programSet_GLSL.set_base_paths(inVertexBasePath, inFragmentBasePath);
#ifdef HAVE_PKG_CG
    programSet_CG.set_base_paths(inVertexBasePath, inFragmentBasePath);
#endif  	
  }

};

} } // namespaces GPU, vw

#endif

