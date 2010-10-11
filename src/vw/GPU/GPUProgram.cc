// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/*
  *  GPUProgram.cpp
  *  VWGPU
  *
  *  Created by Ian Saxton on 5/17/06.
  *
  */

 #include <fstream>
 #include <iostream>

 #include <boost/algorithm/string.hpp>

 #include <vw/Core/Exception.h>
 #include <vw/GPU/GPUProgram.h>
 #include <vw/GPU/Utilities.h>

using std::string;
using std::pair;
using std::vector;
using std::map;


namespace vw { namespace GPU {

  int maxErrorLength = 2048;

  char errorString[2048];



  string shader_base_path = "";

  bool shader_assembly_cache_enabled = false;

  string shader_assembly_cache_path = "";


  ShaderCompilationStatusEnum shaderCompilationStatus = SHADER_COMPILATION_STATUS_SUCCESS;


  //########################################################################
  //#    GLSL - Program
  //########################################################################

  bool GPUProgram_GLSL::Link(GPUVertexShader_GLSL& vertex, GPUFragmentShader_GLSL& fragment) {
    program = glCreateProgramObjectARB();
    if(vertex.is_compiled())
      glAttachObjectARB(program, vertex.get_shader());
    if(fragment.is_compiled())
      glAttachObjectARB(program, fragment.get_shader());
    glLinkProgramARB(program);
    GLsizei errorStringLength;
    GLint isCompiled;
    glGetInfoLogARB(program, maxErrorLength, &errorStringLength, errorString);
    glGetObjectParameterivARB(program, GL_OBJECT_LINK_STATUS_ARB, &isCompiled);
    if(!isCompiled) {
      printf("***PROGRAM LINKER***\n%s\n", errorString);
      program = 0;
      return false;
    }
    return true;
  }

  //########################################################################
  //#    GLSL - Vertex Shader
  //########################################################################


  bool GPUVertexShader_GLSL::compile(const string& vertexString) {
    const char* cStr = vertexString.c_str();
    shader = glCreateShaderObjectARB(GL_VERTEX_SHADER_ARB);
    glShaderSourceARB(shader, 1, (const GLcharARB**) &cStr, NULL);
    glCompileShaderARB(shader);
    GLsizei errorStringLength;
    GLint isCompiled;
    glGetInfoLogARB(shader, maxErrorLength, &errorStringLength, errorString);
    glGetObjectParameterivARB(shader, GL_OBJECT_COMPILE_STATUS_ARB, &isCompiled);
    if(!isCompiled) {
      gpu_log("\n*********GLSL Vertex Shader Compilation Error*********\n");
      gpu_log(errorString);
      gpu_log("********************************************************\n");
      gpu_log(vertexString.c_str());
      gpu_log("******************************************************\n");
      return false;
    }
    return true;
  }

  //########################################################################
  //#    GLSL - Fragment Shader
  //########################################################################



  bool GPUFragmentShader_GLSL::compile(const string& fragmentString) {
    const char* cStr = fragmentString.c_str();
    shader = glCreateShaderObjectARB(GL_FRAGMENT_SHADER_ARB);
    glShaderSourceARB(shader, 1, (const GLcharARB**) &cStr, NULL);

    glCompileShaderARB(shader);
    GLint isCompiled;
    GLsizei errorStringLength;
    glGetInfoLogARB(shader, maxErrorLength, &errorStringLength, errorString);
    glGetObjectParameterivARB(shader, GL_OBJECT_COMPILE_STATUS_ARB, &isCompiled);
    if(!isCompiled) {
      gpu_log("\n*********GLSL Fragment Shader Compilation Error*********\n");
      gpu_log(errorString);
      gpu_log("********************************************************\n");
      gpu_log(fragmentString.c_str());
      gpu_log("\n********************************************************\n");

      return false;
    }
    return true;
  }



  //########################################################################
  //#    CG Specific Code
  //########################################################################

#if defined(VW_HAVE_PKG_CG) && VW_HAVE_PKG_CG==1

  CGerror CGCheckError();

  //########################################################################
  //#    GPUShader_CG
  //########################################################################


  CGcontext cgContext = NULL;
  void InitCGContext();
  void CGErrorCallback();


  GPUShader_CG::GPUShader_CG() {
    program = NULL;
  }

  GPUShader_CG::~GPUShader_CG() {
    if(program)
      cgDestroyProgram(program);
  }

  bool GPUShader_CG::compile_source_with_string(const char *sourceString, CGprofile profile, const char *entry, const char **args) {
    InitCGContext();
    // Delete old program object
    if(program)
      cgDestroyProgram(program);
    // Compile
    GPUShader_CG::profile = profile;
    program = cgCreateProgram(cgContext, CG_SOURCE, sourceString, profile, entry, args);
    // Check for errors
    CGerror error = cgGetError();
    if(!cgIsProgramCompiled(program)) {
      program = NULL;
      return false;
    }
    // Load into cgGL
    cgGLLoadProgram(program);
    return true;
  }


  bool GPUShader_CG::compile_source_with_file(const char *sourceString, CGprofile profile, const char *entry, const char **args) {
    InitCGContext();
    // Delete old program object
    if(program)
      cgDestroyProgram(program);
    // Compile
    program = cgCreateProgramFromFile(cgContext, CG_SOURCE, sourceString, profile, entry, args);
    GPUShader_CG::profile = profile;
    // printf("%s", cgGetLastListing(cgContext));
    // Check for errors
    CGerror error = cgGetError();
    if(!cgIsProgramCompiled(program)) {
      printf("*** ERROR (GPUShader_CG::compile_source_with_string()  Compile Failed. Exiting...\n");
      program = NULL;
      exit(0);
    }
    // Load into cgGL
    cgGLLoadProgram(program);
    return true;
  }

  bool GPUShader_CG::save_compiled_file(const char *file) {
    std::ofstream outFile(file);
    if(!outFile) {
      return false;
    }
    if(!is_compiled()) {
      return false;
    }
    const char *entryString = cgGetProgramString(program, CG_PROGRAM_ENTRY);
    outFile << entryString << " ";
    const char *profileString = cgGetProgramString(program, CG_PROGRAM_PROFILE);
    outFile << profileString << " ";
    const char *objectString = cgGetProgramString(program, CG_COMPILED_PROGRAM);
    outFile << objectString;

    outFile.close();
    return true;
  }

  bool GPUShader_CG::load_compiled_file(const char *file) {
    InitCGContext();
    // Input File
    std::ifstream inFile(file);
    if(!inFile) {
      return false;
    }
    // Read Entry and Profile Strings
    string entryString;
    inFile >> entryString;
    string profileString;
    inFile >> profileString;

    CGprofile profile = cgGetProfile(profileString.c_str());
    if(profile == CG_PROFILE_UNKNOWN) {
      return false;
    }
    // Read Object String
    int objectStringStart = ((int) inFile.tellg()) + 1;
    inFile.seekg(0, std::ios::end);
    int objectStringLength = ((int) inFile.tellg()) - objectStringStart;
    inFile.seekg(objectStringStart, std::ios::beg);
    char *objectString = (char*) malloc(objectStringLength + 1);
    inFile.read(objectString, objectStringLength);
    inFile.close();
    objectString[objectStringLength] = 0;
    // Delete old program object
    if(program)
      cgDestroyProgram(program);
    CGCheckError();
    program = cgCreateProgram(cgContext, CG_OBJECT, objectString, profile, entryString.c_str(), NULL);
    GPUShader_CG::profile = profile;
    // Check for errors
    if(!is_compiled() || CGCheckError()) {
      program = NULL;
      return false;
    }
    CGCheckError();
    cgGLLoadProgram(program);
    if(CGCheckError()) {
      return false;
    }
    return true;
  }




  void GPUShader_CG::install() {
    //cgGLEnableProfile(CG_PROFILE_FP30);
    cgGLEnableProfile(profile);
    cgGLBindProgram(program);
  }


  void GPUShader_CG::uninstall() {
    cgGLDisableProfile(profile);
  }


  //########################################################################
  //#                           CG - Internal Globals                      #
  //########################################################################


  CGerror CGCurrentError = (CGerror) 0;

  void InitCGContext() {
    if(!cgContext) {
      //cgSetErrorCallback(&CGErrorCallback);
      cgContext = cgCreateContext();
    }
  }


  void CGErrorCallback() {
    while(1) {
      CGerror error = cgGetError();
      if(error == 0)
        break;
      else {
        CGCurrentError = error;
        printf("***CG ERROR:   %s \n", cgGetErrorString(error));
      }
    }
  }

  CGerror CGCheckError() {
    CGerror temp = CGCurrentError;
    CGCurrentError = (CGerror) 0;
    return temp;
  }



#endif // ifdef VW_HAVE_PKG_CG___

  //#############################################################################################
  //#    Creation Free Functions - GLSL
  //#############################################################################################

  typedef pair<pair<string, string>, pair<vector<int>, vector<int> > > GPUProgramKey;

  GPUProgram_GLSL* create_gpu_program_glsl_string(const string& fragmentString, const vector<int>& fragmentAttributes,
                                                  const string& vertexString, const vector<int>& vertexAttributes)
  {
    shaderCompilationStatus = SHADER_COMPILATION_STATUS_SUCCESS_FILE;
    static char charBuffer1[32];
    static char charBuffer2[32];
    // Vertex
    GPUVertexShader_GLSL vertexShader;
    if(!vertexString.empty()) {
      string vertReplacedString;
      const string* sourceString = &vertexString;
      if(vertexAttributes.size()) {
        TokenReplacer tr;
        for(unsigned i=0; i < vertexAttributes.size(); i++) {
          sprintf(charBuffer1, "%i", i+1);
          sprintf(charBuffer2, "%i", vertexAttributes[i]);
          tr.AddVariable(charBuffer1, charBuffer2);
        }
        tr.Replace(vertexString, vertReplacedString);
        sourceString = &vertReplacedString;
      }
      if(!vertexShader.compile(*sourceString)) {
        shaderCompilationStatus = SHADER_COMPILATION_STATUS_ERROR_COMPILE;
        throw(Exception("GPUProgram creation failed."));
      }
    }
    // Fragment
    GPUFragmentShader_GLSL fragmentShader;
    if(!fragmentString.empty()) {
      string fragReplacedString;
      const string* sourceString = &fragmentString;
      if(fragmentAttributes.size()) {
        TokenReplacer tr;
        string variable;
        for(unsigned i=0; i < fragmentAttributes.size(); i++) {
          sprintf(charBuffer1, "%i", i+1);
          sprintf(charBuffer2, "%i", fragmentAttributes[i]);
          tr.AddVariable(charBuffer1, charBuffer2);
        }
        tr.Replace(fragmentString, fragReplacedString);
        sourceString = &fragReplacedString;
      }
      if(!fragmentShader.compile(*sourceString)) {
        shaderCompilationStatus = SHADER_COMPILATION_STATUS_ERROR_COMPILE;
        throw(Exception("GPUProgram creation failed."));
      }
    }
    // Program
    GPUProgram_GLSL* program = new GPUProgram_GLSL;
    if(!program->Link(vertexShader, fragmentShader)) {
      shaderCompilationStatus = SHADER_COMPILATION_STATUS_ERROR_LINK;
      delete program;
      throw(Exception("GPUProgram creation failed."));
    }
    return program;
  }




  static std::map<GPUProgramKey, GPUProgram_GLSL*> program_cache_glsl;



  GPUProgram_GLSL* create_gpu_program_glsl(const string& fragmentPath, const vector<int>& fragmentAttributes,
                                           const string& vertexPath, const vector<int>& vertexAttributes)
  {
    // Static
    // Check cache
    GPUProgramKey programKey = GPUProgramKey(pair<string, string>(fragmentPath, vertexPath), pair<vector<int>, vector<int> >(fragmentAttributes, vertexAttributes));
    map<GPUProgramKey, GPUProgram_GLSL*>::iterator iter_prog;
    iter_prog = program_cache_glsl.find(programKey);
    if(iter_prog != program_cache_glsl.end()) {
      shaderCompilationStatus = SHADER_COMPILATION_STATUS_SUCCESS_CACHE;
      return (*iter_prog).second;
    }
    // Vertex File Read
    string vertexString;
    GPUVertexShader_GLSL vertexShader;
    if(!vertexPath.empty()) {
      std::map<std::string, const char*>::iterator iter_map = standard_shaders_map.find((vertexPath).c_str());
      if(iter_map != standard_shaders_map.end()) {
        vertexString = (*iter_map).second;
      }
      else if(!ReadFileAsString(shader_base_path + vertexPath, vertexString)) {
        shaderCompilationStatus = SHADER_COMPILATION_STATUS_ERROR_FILE;
        throw(Exception("GPUProgram creation failed."));
      }
    }
    // Fragment File Read
    string fragmentString;
    if(!fragmentPath.empty()) {
      std::map<std::string, const char*>::iterator iter_map = standard_shaders_map.find((fragmentPath).c_str());
      if(iter_map != standard_shaders_map.end()) {
        fragmentString = (*iter_map).second;
      }
      else if(!ReadFileAsString(shader_base_path + fragmentPath, fragmentString)) {
        shaderCompilationStatus = SHADER_COMPILATION_STATUS_ERROR_FILE;
        throw(Exception("GPUProgram creation failed."));
      }
    }
    // Compile Strings
    GPUProgram_GLSL* program = create_gpu_program_glsl_string(fragmentString, fragmentAttributes, vertexString, vertexAttributes);
    if(program)
      program_cache_glsl[programKey] = program;
    return program;
  }






  //#############################################################################################
  //#    Creation Free Functions - CG
  //#############################################################################################

#if defined(VW_HAVE_PKG_CG) && VW_HAVE_PKG_CG==1

  GPUProgram_CG* create_gpu_program_cg_string(const string& fragmentString, const vector<int>& fragmentAttributes,
                                              const string& vertexString, const vector<int>& vertexAttributes)
  {
    std::auto_ptr<GPUShader_CG> vertexShader(NULL);
    std::auto_ptr<GPUShader_CG> fragmentShader(NULL);
    static char charBuffer1[32];
    static char charBuffer2[32];
    // VERTEX
    if(!vertexString.empty()) {
      vertexShader.reset(new GPUShader_CG());
      string vertReplacedString;
      const string* sourceString = &vertexString;
      if(vertexAttributes.size()) {     // Specialize String
        TokenReplacer tr;
        string variable;
        for(int i=0; i < vertexAttributes.size(); i++) {
          sprintf(charBuffer1, "%i", i+1);
          sprintf(charBuffer2, "%i", vertexAttributes[i]);
          tr.AddVariable(charBuffer1, charBuffer2);
        }
        tr.Replace(vertexString, vertReplacedString);
        sourceString = &vertReplacedString;
      }
      if(!vertexShader->compile_source_with_string(sourceString->c_str(), CG_PROFILE_VP30, "main", NULL)) {
        shaderCompilationStatus = SHADER_COMPILATION_STATUS_ERROR_COMPILE;
        throw(Exception("GPUProgram creation failed."));
      }
    }
    // FRAGMENT
    if(!fragmentString.empty()) {
      fragmentShader.reset(new GPUShader_CG());
      string fragReplacedString;
      const string* sourceString = &fragmentString;
      if(fragmentAttributes.size()) {     // Specialize String
        TokenReplacer tr;
        string variable;
        for(int i=0; i < fragmentAttributes.size(); i++) {
          sprintf(charBuffer1, "%i", i+1);
          sprintf(charBuffer2, "%i", fragmentAttributes[i]);
          tr.AddVariable(charBuffer1, charBuffer2);
        }
        tr.Replace(fragmentString, fragReplacedString);
        sourceString = &fragReplacedString;
      }
      if(!fragmentShader->compile_source_with_string(sourceString->c_str(), CG_PROFILE_FP30, "main", NULL)) {
        shaderCompilationStatus = SHADER_COMPILATION_STATUS_ERROR_COMPILE;
        throw(Exception("GPUProgram creation failed."));
      }
    }
    // PROGRAM - Make new program, put in cache and return it
    GPUProgram_CG* program = new GPUProgram_CG(vertexShader.get(), fragmentShader.get());
    vertexShader.release();
    fragmentShader.release();
    shaderCompilationStatus = SHADER_COMPILATION_STATUS_SUCCESS_FILE;
    return program;
  }


  std::map<GPUProgramKey, GPUProgram_CG*> program_cache_cg;


  GPUProgram_CG* create_gpu_program_cg(const string& fragmentPath, const vector<int>& fragmentAttributes,
                                       const string& vertexPath, const vector<int>& vertexAttributes)
  {

    char charBuffer1[256];
    char charBuffer2[32];
    shaderCompilationStatus = SHADER_COMPILATION_STATUS_SUCCESS_FILE;
    // Check cache for program
    GPUProgramKey programKey = GPUProgramKey(pair<string, string>(fragmentPath, vertexPath), pair<vector<int>, vector<int> >(fragmentAttributes, vertexAttributes));
    map<GPUProgramKey, GPUProgram_CG*>::iterator iter_prog;
    iter_prog = program_cache_cg.find(programKey);
    if(iter_prog != program_cache_cg.end()) {
      shaderCompilationStatus = SHADER_COMPILATION_STATUS_SUCCESS_CACHE;
      return (*iter_prog).second;
    }

    //
    std::auto_ptr<GPUShader_CG> vertexShader(NULL);
    std::auto_ptr<GPUShader_CG> fragmentShader(NULL);
    // Get Source Strings - Try to find in STD virtual directory, then try the real path
    string vertRawString;
    string fragRawString;
    string fragAssemblyFilePath;
    bool fragComplete = false;
    bool fail = false;
    // VERTEX
    if(!vertexPath.empty()) {
      std::map<std::string, const char*>::iterator iter_map = standard_shaders_map.find((vertexPath).c_str());
      if(iter_map != standard_shaders_map.end()) {
        vertRawString = (*iter_map).second;
      }
      else if(!ReadFileAsString(shader_base_path + vertexPath, vertRawString)) {
        shaderCompilationStatus = SHADER_COMPILATION_STATUS_ERROR_FILE;
        fail = true;
      }
    }
    // FRAGMENT
    if(!fail && !fragmentPath.empty()) {
      // If using assembly caching, create assembly path and attempt to load it.
      if(shader_assembly_cache_enabled) {
        string modifiedFragmentPath = fragmentPath;
        boost::algorithm::replace_all(modifiedFragmentPath, "/", "_");
        fragAssemblyFilePath = shader_assembly_cache_path + modifiedFragmentPath;
        for(int i=0; i < fragmentAttributes.size(); i++) {
          sprintf(charBuffer1, "_%i", fragmentAttributes[i]);
          fragAssemblyFilePath += charBuffer1;
        }
        fragAssemblyFilePath += ".cache";
        if(fragmentShader->load_compiled_file(fragAssemblyFilePath.c_str())) {
          fragComplete = true;
        }
      }
      // If necessary, read source string from file
      if(!fragComplete) {
        std::map<std::string, const char*>::iterator iter_map = standard_shaders_map.find((fragmentPath).c_str());
        if(iter_map != standard_shaders_map.end()) {
          fragRawString = (*iter_map).second;
        }
        else if(!ReadFileAsString(shader_base_path + fragmentPath, fragRawString)) {
          shaderCompilationStatus = SHADER_COMPILATION_STATUS_ERROR_FILE;
          fail = true;
        }
      }
    }
    // Specialize Strings and Compile Shaders
    // VERTEX
    if(!fail && !vertexPath.empty()) {
      vertexShader.reset(new GPUShader_CG);
      string vertReplacedString;
      string& sourceString = vertRawString;
      if(vertexAttributes.size()) {     // Specialize String
        TokenReplacer tr;
        string variable;
        for(int i=0; i < vertexAttributes.size(); i++) {
          sprintf(charBuffer1, "%i", i+1);
          sprintf(charBuffer2, "%i", vertexAttributes[i]);
          tr.AddVariable(charBuffer1, charBuffer2);
        }
        tr.Replace(vertRawString, vertReplacedString);
        sourceString = vertReplacedString;
      }
      if(!vertexShader->compile_source_with_string(sourceString.c_str(), CG_PROFILE_VP30, "main", NULL)) {
        shaderCompilationStatus = SHADER_COMPILATION_STATUS_ERROR_COMPILE;
        return NULL;
      }
    }
    // FRAGMENT
    if(!fail && !fragmentPath.empty() && !fragComplete) {
      fragmentShader.reset(new GPUShader_CG);
      string fragReplacedString;
      string& sourceString = fragRawString;
      if(fragmentAttributes.size()) {     // Specialize String
        TokenReplacer tr;
        string variable;
        for(int i=0; i < fragmentAttributes.size(); i++) {
          sprintf(charBuffer1, "%i", i+1);
          sprintf(charBuffer2, "%i", fragmentAttributes[i]);
          tr.AddVariable(charBuffer1, charBuffer2);
        }
        tr.Replace(fragRawString, fragReplacedString);
        sourceString = fragReplacedString;
      }
      if(!fragmentShader->compile_source_with_string(sourceString.c_str(), CG_PROFILE_FP30, "main", NULL)) {
        shaderCompilationStatus = SHADER_COMPILATION_STATUS_ERROR_COMPILE;
        fail = true;
      }
      if(shader_assembly_cache_enabled) {
        if(fragmentShader->save_compiled_file(fragAssemblyFilePath.c_str())) {
        }
      }
    }

    if(fail)
      throw(Exception("GPUProgram creation failed."));
    // PROGRAM - Make new program, put in cache and return it
    GPUProgram_CG* program = new GPUProgram_CG(vertexShader.get(), fragmentShader.get());
    vertexShader.release();
    fragmentShader.release();

    program_cache_cg[programKey] = program;
    return program;
  }


#endif

  //#############################################################################################
  //#    Creation Free Functions - Base
  //#############################################################################################


  GPUProgram* create_gpu_program(const string& fragmentPath, const vector<int>& fragmentAttributes,
                                 const string& vertexPath, const vector<int>& vertexAttributes) {
    // LOGGING
    static char  buffer1[256];
    static char buffer2[256];
    buffer1[0] = 0;
    if(vertexAttributes.size()) {
      sprintf(buffer1, "< ");
      for(size_t i=0; i<vertexAttributes.size(); i++)
        sprintf(buffer1, "%s%i ", buffer1, vertexAttributes[i]);
      sprintf(buffer1, "%s>", buffer1);
    }
    buffer2[0] = 0;
    if(fragmentAttributes.size()) {
      sprintf(buffer2, "< ");
      for(size_t i=0; i<fragmentAttributes.size(); i++)
        sprintf(buffer2, "%s%i ", buffer2, fragmentAttributes[i]);
      sprintf(buffer2, "%s>", buffer2);
    }

    static char buffer3[512];
    GPUProgram* program = NULL;
    sprintf(buffer3, "[create_gpu_program] FRAGMENT: %s%s, VERTEX: %s%s }    ", fragmentPath.c_str(), buffer2, vertexPath.c_str(), buffer1 );
    gpu_log(buffer3);
    // FIND PROGRAM

    gpu_log(get_string_for_shader_language_choice_enum(shaderLanguageChoice));
    gpu_log("  ");
#if defined(VW_HAVE_PKG_CG) && VW_HAVE_PKG_CG==1
    if(shaderLanguageChoice == SHADER_LANGUAGE_CHOICE_CG_GLSL || shaderLanguageChoice == SHADER_LANGUAGE_CHOICE_CG) {
      try {
        program = create_gpu_program_cg(fragmentPath  + ".cg", fragmentAttributes,
                                        vertexPath.empty() ? "" : vertexPath + ".cg", vertexAttributes);
      }
      catch(...) { program = NULL; }
      gpu_log("Trying CG... ");
    }
#endif

    if(!program && shaderLanguageChoice != SHADER_LANGUAGE_CHOICE_CG) {
      try {
        program = create_gpu_program_glsl(fragmentPath + ".glsl", fragmentAttributes,
                                          vertexPath.empty() ? "" : vertexPath + ".glsl", vertexAttributes);
      }
      catch(...) {  program = NULL; }
      gpu_log("Trying GLSL... ");
    }

#if defined(VW_HAVE_PKG_CG) && VW_HAVE_PKG_CG==1

    if(!program && shaderLanguageChoice == SHADER_LANGUAGE_CHOICE_GLSL_CG) {
      try {
        program = create_gpu_program_cg(fragmentPath  + ".cg", fragmentAttributes,
                                        vertexPath.empty() ? "" : vertexPath + ".cg", vertexAttributes);
      }
      catch(...) {  program = NULL; }

      gpu_log("Trying CG... ");
    }
#endif

    if(program) {
      if(get_shader_compilation_status() == SHADER_COMPILATION_STATUS_SUCCESS_FILE)
        gpu_log("SUCCESS (From File)\n");
      if(get_shader_compilation_status() == SHADER_COMPILATION_STATUS_SUCCESS_CACHE)
        gpu_log("SUCCESS (From Cache)\n");
    }
    else {
      if(get_shader_compilation_status() == SHADER_COMPILATION_STATUS_ERROR_FILE)
        gpu_log("*** FAILED (File Error) ***\n");
      if(get_shader_compilation_status() == SHADER_COMPILATION_STATUS_ERROR_COMPILE)
        gpu_log("*** FAILED (Compile Error) ***\n");
      if(get_shader_compilation_status() == SHADER_COMPILATION_STATUS_ERROR_LINK)
        gpu_log("*** FAILED (Link Error) ***\n");

      throw(Exception("GPUProgram creation failed."));
    }

    return program;
  }

  void clear_gpu_program_cache() {
    for(map<GPUProgramKey, GPUProgram_GLSL*>::iterator iter = program_cache_glsl.begin(); iter != program_cache_glsl.end(); iter++)
      delete (*iter).second;
    program_cache_glsl.clear();

#if defined(VW_HAVE_PKG_CG) && VW_HAVE_PKG_CG==1
    for(map<GPUProgramKey, GPUProgram_CG*>::iterator iter = program_cache_cg.begin(); iter != program_cache_cg.end(); iter++)
      delete (*iter).second;
    program_cache_cg.clear();
#endif

  }



} } // namespaces GPU, vw


