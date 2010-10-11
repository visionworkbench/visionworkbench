// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/GPU/Setup.h>
#include <vw/GPU/GPUImage.h>
#include <fstream>

using std::string;
using std::map;

namespace vw { namespace GPU {

// Class: UtilityTimer

void UtilityTimer::StopAndPrint(char* text) {
  if(!enabled) return;
  Stop();
  if(text)
    printf("%s    ", text);
  Print();
}


void UtilityTimer::Print(double scalar) {
  if(!enabled) return;

  long long int deltaTimeNanos;
  double endTimeFloat;
  double startTimeFloat;
  if(averageCount==0) {
    endTimeFloat = (endTime.tv_sec + endTime.tv_usec/1000000.0);
    startTimeFloat = (startTime.tv_sec + startTime.tv_usec/1000000.0);
    deltaTimeNanos = (long long int) (1000000000 * scalar * (endTimeFloat - startTimeFloat));
  }
  else {
    printf("Average Time, %li runs: ", averageCount);
    deltaTimeNanos = (long long int) (1000000000.0 * (scalar * averageValue));
  }
  int seconds = (int) (deltaTimeNanos / 1000000000);
  int ms = (int) (((long long)  deltaTimeNanos % 1000000000) / 1000000);
  int us = (int) (((long long) deltaTimeNanos % 1000000) / 1000);
  int ns = (int) ((long long)  deltaTimeNanos % 1000);
  if (scalar != 1.0)
    printf("(Scaled by %.6f) ", scalar);
  printf("%i.%03i,%03i Seconds\n", seconds, ms, us);
}

double UtilityTimer::ElapsedSeconds() {
  double end = (endTime.tv_sec + endTime.tv_usec/1000000.0);
  double start = (startTime.tv_sec + startTime.tv_usec/1000000.0);
  return (end - start);
}


void UtilityTimer::AddToAverage() {
  //  averageValue = ((averageValue*averageCount) + (endTime - startTime)) / (double) (averageCount + 1);
  //  averageCount++;
}

// Class: TokenReplacer

void TokenReplacer::AddVariable(string& variableName, string& variableValue) {
  stringMap[variableName] = variableValue;
}

void TokenReplacer::AddVariable(const char* variableName, const char*  variableValue) {
  stringMap[string(variableName)] = string(variableValue);
}


void TokenReplacer::Replace(const string& inText, string& outText) {
  int inLength = inText.size();
  outText.resize(0);
  outText.reserve((int) (inLength * 1.5));
  string tempToken;
  int i=0;

  while(i < inLength) {
    char cur = inText[i];
    if(cur == '$') {
      tempToken.resize(0);
      tempToken.reserve(32);
      i++;
      while(i < inLength) {
        cur = inText[i];
        if((cur > 47 && cur < 58) || (cur > 64 && cur < 91) || (cur > 96 && cur < 123) || cur == 95) {
          tempToken.push_back(cur);
          i++;
        }
        else {
          break;
        }
      }
      if(tempToken.size()) {
        map<string, string>::iterator iter = stringMap.find(tempToken);
        if(iter != stringMap.end()) {
          outText.append((*iter).second);  // OK to push string?
        }
      }
    }
    else {
      outText.push_back(cur);
      i++;
    }
  }
}


// OpenGL

bool CheckFramebuffer(bool enablePrint) {
  GLuint status = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
  bool isComplete = (status == GL_FRAMEBUFFER_COMPLETE_EXT);
  if(enablePrint && !isComplete) {
    switch(status) {
    case GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT_EXT:
      printf("Framebuffer Incomplete: FRAMEBUFFER_INCOMPLETE_ATTACHMENT_EXT\n"); break;
    case GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT_EXT:
      printf("Framebuffer Incomplete: FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT_EXT\n"); break;
    case GL_FRAMEBUFFER_INCOMPLETE_DIMENSIONS_EXT:
      printf("Framebuffer Incomplete: FRAMEBUFFER_INCOMPLETE_DIMENSIONS_EXT\n"); break;
    case GL_FRAMEBUFFER_INCOMPLETE_FORMATS_EXT:
      printf("Framebuffer Incomplete: FRAMEBUFFER_INCOMPLETE_FORMATS_EXT\n"); break;
    case GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER_EXT:
      printf("Framebuffer Incomplete: FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER_EXT\n"); break;
    case GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER_EXT:
      printf("Framebuffer Incomplete: FRAMEBUFFER_INCOMPLETE_READ_BUFFER_EXT\n"); break;
    case GL_FRAMEBUFFER_UNSUPPORTED_EXT:
      printf("Framebuffer Incomplete: FRAMEBUFFER_UNSUPPORTED_EXT\n"); break;
    }
  }
  return isComplete;
}

// File I/O

bool ReadFileAsString(const string& path, string& outString) {
  std::ifstream file(path.c_str());
  if(!file)
    return false;
  //while(file)
  //  outString << file;

  int length;
  file.seekg(0, std::ios::end);
  length = file.tellg();
  file.seekg(0, std::ios::beg);
  outString.resize(length);
  file.read((char*) outString.c_str(), length);

  return true;
}




  // ##################### Shader Invocatin - Helper Functions ########################

  void ShaderInvocation_SetupGLState(int width, int height) {
  glEnable(GL_TEXTURE_RECTANGLE_ARB);
  glPolygonMode(GL_FRONT,GL_FILL);
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, g_framebuffer);

  glPushMatrix();
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(0.0,width,0.0, height);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glViewport(0,0,width, height);

  }

  void ShaderInvocation_SetupGLState(const GPUImageBase& image) {
    //glEnable(GL_TEXTURE_RECTANGLE_ARB);
    //glPolygonMode(GL_FRONT,GL_FILL);
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, g_framebuffer);

  int width = image.real_width();
  int height = image.real_height();

  glPushMatrix();
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(0.0,width,0.0, height);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glViewport(0,0,width, height);

  }

  void ShaderInvocation_SetOutputImage(const GPUImageBase& image) {
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, g_framebuffer);
  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, image.target(), image.name(), 0);
  }



  void ShaderInvocation_DrawRectOneTexture(const GPUImageBase& image, int x, int y) { // New Name? DrawRectAligned_OneTexture
  float width = image.width();
  float height = image.height();
  float t0_x1 = image.x(0);
  float t0_x2 = image.x(width);
  float t0_y1 = image.y(0);
  float t0_y2 = image.y(height);

  glBegin(GL_QUADS);
  glTexCoord2f(t0_x1, t0_y1);     glVertex2f(x, y);
  glTexCoord2f(t0_x2, t0_y1);     glVertex2f(x + width, y);
  glTexCoord2f(t0_x2, t0_y2);   glVertex2f(x + width, y + height);
  glTexCoord2f(t0_x1, t0_y2);     glVertex2f(x, y + height);
  glEnd();
  }

  void ShaderInvocation_DrawRectOneTexture(const GPUImageBase& image, int x, int y, int width, int height) {
  float t0_x1 = image.x(0);
  float t0_x2 = image.x(width);
  float t0_y1 = image.y(0);
  float t0_y2 = image.y(height);

  glBegin(GL_QUADS);
  glTexCoord2f(t0_x1, t0_y1);     glVertex2f(x, y);
  glTexCoord2f(t0_x2, t0_y1);     glVertex2f(x + width, y);
  glTexCoord2f(t0_x2, t0_y2);   glVertex2f(x + width, y + height);
  glTexCoord2f(t0_x1, t0_y2);     glVertex2f(x, y + height);
  glEnd();
  }

  void ShaderInvocation_DrawRectOneTexture(const GPUImageBase& image, const Rectangle2D<int>& rect, const Rectangle2D<int>& texCoord) {
  float t0_x1 = image.x(texCoord.x);
  float t0_x2 = image.x(texCoord.x + texCoord.width);
  float t0_y1 = image.y(texCoord.y);
  float t0_y2 = image.y(texCoord.y + texCoord.height);

  glBegin(GL_QUADS);
  glTexCoord2f(t0_x1, t0_y1);     glVertex2f(rect.x, rect.y);
  glTexCoord2f(t0_x2, t0_y1);     glVertex2f(rect.x + rect.width, rect.y);
  glTexCoord2f(t0_x2, t0_y2);   glVertex2f(rect.x + rect.width, rect.y + rect.height);
  glTexCoord2f(t0_x1, t0_y2);     glVertex2f(rect.x, rect.y + rect.height);
  glEnd();
  }

void ShaderInvocation_DrawRectOneTexture(const GPUImageBase& image, Vector2 bound_left_upper, Vector2 bound_right_lower,
                      Vector2 tex_left_upper, Vector2 tex_right_lower)
  {
  float t0_x1 = image.x(tex_left_upper.x());
  float t0_x2 = image.x(tex_right_lower.x());
  float t0_y1 = image.y(tex_left_upper.y());
  float t0_y2 = image.y(tex_right_lower.y());

  glBegin(GL_QUADS);
  glTexCoord2f(t0_x1, t0_y1);     glVertex2f(bound_left_upper.x(), bound_left_upper.y());
  glTexCoord2f(t0_x2, t0_y1);     glVertex2f(bound_right_lower.x(), bound_left_upper.y());
  glTexCoord2f(t0_x2, t0_y2);     glVertex2f(bound_right_lower.x(), bound_right_lower.y());
  glTexCoord2f(t0_x1, t0_y2);     glVertex2f(bound_left_upper.x(), bound_right_lower.y());
  glEnd();
  }




  void ShaderInvocation_CleanupGLState() {
  glDisable(GL_TEXTURE_RECTANGLE_ARB);
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
  glPopMatrix();
  }


} } // namespaces GPU, vw
