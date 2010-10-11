// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef OpenGLManager_H
#define OpenGLManager_H

#include <vw/Math/Vector.h>

#include <cstdlib>
#include <cstdio>
#include <sys/time.h>
#include <cmath>
#include <string>
#include <map>
#include <vector>
#include <list>

#ifdef __APPLE__
  #include "OpenGL/gl.h"
  #include "GLUT/glut.h"
#else
  #include "GL/glew.h"
  #include <GL/glut.h>
#endif


// Float Buffer Enums
#ifdef __APPLE__
  #define GL_FLOAT_R32_NV  GL_LUMINANCE_FLOAT32_ATI
  #define GL_FLOAT_RGB32_NV  GL_RGB_FLOAT32_ATI
  #define GL_FLOAT_RGBA32_NV  GL_RGBA_FLOAT32_ATI

  #define GL_FLOAT_R_NV  GL_FLOAT_R32_NV
  #define GL_FLOAT_RGB_NV  GL_FLOAT_RGB32_NV
  #define GL_FLOAT_RGBA_NV  GL_FLOAT_RGBA32_NV

  #define GL_FLOAT_R16_NV  GL_LUMINANCE_FLOAT16_ATI
  #define GL_FLOAT_RGB16_NV  GL_RGB_FLOAT16_ATI
  #define GL_FLOAT_RGBA16_NV  GL_RGBA_FLOAT16_ATI
#endif


namespace vw { namespace GPU {

// Class: Rectangle2D

template <class T>
  class Rectangle2D {
 public:
  T x, y, width, height;
  Rectangle2D(T _x, T _y, T _width, T _height) : x(_x), y(_y), width(_width), height(_height) { }
 };

 // Class: UtilityTimer

class UtilityTimer  {
protected:
  timeval  startTime;
        timeval  endTime;
  long averageCount;
  double  averageValue;
  bool enabled;

public:
  UtilityTimer() { Clear(); enabled = true; }
  void Clear() { averageCount = 0; averageValue = 0; }
  void Start() { if(enabled) gettimeofday(&startTime, NULL);  }
  void Stop() { if(enabled) gettimeofday(&endTime, NULL); }
  void SetEnabled(bool value) { enabled = value; }

  void StopAndPrint(char* text);
  void AddToAverage();
  void Print(double scalar = 1.0);
  double ElapsedSeconds();

};

// Class: TokenReplacer

class TokenReplacer {
  std::map<std::string, std::string> stringMap;
public:
  void AddVariable(std::string& variableName, std::string& variableValue);
  void AddVariable(const char* variableName, const char*  variableValue);
  void Replace(const std::string& inText, std::string& outText);
// INLINE
  void Clear() {
    stringMap.clear();
  };
};

// File I/O

bool ReadFileAsString(const std::string& path, std::string& outString);

// OpenGL

 bool CheckFramebuffer(bool enablePrint = false);




  // ##################### Shader Invocation - Helper Functions ########################

class GPUImageBase;

void ShaderInvocation_SetupGLState(int width, int height);

void ShaderInvocation_SetupGLState(const GPUImageBase& image);

void ShaderInvocation_SetOutputImage(const GPUImageBase& image);

void ShaderInvocation_DrawRectOneTexture(const GPUImageBase& image,  int x = 0, int y = 0);

void ShaderInvocation_DrawRectOneTexture(const GPUImageBase& image, int x, int y, int width, int height);

void ShaderInvocation_DrawRectOneTexture(const GPUImageBase& image, const Rectangle2D<int>& rect, const Rectangle2D<int>& texCoord);

void ShaderInvocation_DrawRectOneTexture(const GPUImageBase& image, Vector2 bound_left_upper, Vector2 bound_right_lower,
                      Vector2 tex_left_upper, Vector2 tex_right_lower);

void ShaderInvocation_CleanupGLState();




} } // namespaces GPU, vw


#endif
