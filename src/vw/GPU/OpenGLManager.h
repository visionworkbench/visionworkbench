
#ifndef OpenGLManager_H
#define OpenGLManager_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
//#include <vw/Utilities.h>

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


#include <string>
using namespace std;


namespace vw { namespace GPU {


extern GLuint g_framebuffer;

// INIT

void OpenGLManagerInit(bool dummyWindow = true, bool isLoggingEnabled = false);

void OpenGLManagerExit();

// GLUT Tests

void SetGLUTGPUImage(void* texRef);

void UpdateGLUTWindow();

// UTILITY

bool ReadFileAsString(string& path, string& outString);



 bool CheckFramebuffer(bool enablePrint = false);



inline void CheckGLError(const char* text = NULL) {
	GLenum error = glGetError();
	if(error == GL_INVALID_ENUM)
		printf("GL Error: GL_INVALID_ENUM");
	else if(error == GL_INVALID_VALUE)
		printf("GL Error: GL_INVALID_VALUE");
	else if(error == GL_INVALID_OPERATION)
		printf("GL Error: GL_INVALID_OPERATION");
	else if(error != 0)
		printf("GL Error: 0x%x", error);
	if(error) {
	  if(text)
	    printf("   (%s)", text);
	  printf("\n");
	  }
}


inline void ClearGLError() {
	glGetError();
}

 void WriteToGPULog(const char* string);

  // ##################### Shader Invocatin - Helper Functions ########################

class GPUImageBase;

void ShaderInvocation_SetupGLState(int width, int height);

void ShaderInvocation_SetupGLState(const GPUImageBase& image);

void ShaderInvocation_DrawRectOneTexture(const GPUImageBase& image,  int x = 0, int y = 0);

 void ShaderInvocation_DrawRectOneTexture(const GPUImageBase& image, int x, int y, int width, int height);

 void ShaderInvocation_DrawRectOneTexture(const GPUImageBase& image, const Rectangle2D<int>& rect, const Rectangle2D<int>& texCoord);

void ShaderInvocation_CleanupGLState();

} } // namespaces GPU, vw

#endif
