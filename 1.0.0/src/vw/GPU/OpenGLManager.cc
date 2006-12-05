
#include "OpenGLManager.h"
#include <vw/GPU/GPUImage.h>
#include <fstream>
#include <map>


namespace vw { namespace GPU {


GLuint g_framebuffer;

ofstream gpuLogFile;

static bool loggingEnabled;



// GLUT Tests


GPUImageBase* glutGPUImage;

void glut_display_dummy() {

}
/*
void SetDisplayFunc(void (*func)(void)) {
  glutDisplayFunc(func);
}
*/


// INIT


void OpenGLManagerInit(bool dummyWindow, bool isLoggingEnabled) {
    printf("OpenGLManagerInit\n");
  static bool isInit = false;
  if(!isInit) {
    isInit = true;
    loggingEnabled = isLoggingEnabled;
    // GLUT Init
    int argc = 1;
    char* argv[1];
    argv[0] = "";
    glutInit(&argc, argv);
    if(dummyWindow) {
      glutInitWindowSize(600, 400);
      glutInitDisplayMode(GLUT_RGBA);
      
      glutCreateWindow("GLUT Dummy Window");
      glutDisplayFunc(glut_display_dummy);
#ifndef __APPLE__
      glewInit();
#endif
    }

    glEnable(GL_TEXTURE_2D);
    glEnable(GL_TEXTURE_RECTANGLE_ARB);
    glDisable(GL_DITHER);
    glDisable(GL_BLEND);

    glGenFramebuffersEXT(1, &g_framebuffer);

    if(loggingEnabled)
      gpuLogFile.open("GPU.log");

  }
}

void OpenGLManagerExit() {
     glDeleteFramebuffersEXT(1, &g_framebuffer);
}

  
// UTILITIES

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

bool ReadFileAsString(string& path, string& outString) {
	ifstream file(path.c_str());
	if(!file)
		return false;
	//while(file)
	//	outString << file;
	
	int length;
	file.seekg(0, ios::end);
	length = file.tellg();
	file.seekg(0, ios::beg);
	outString.resize(length);
	file.read((char*) outString.c_str(), length);

	return true;
} 

  void WriteToGPULog(const char* string) {
    if(loggingEnabled) {
      gpuLogFile << string;
      gpuLogFile.flush();
    }
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


  void ShaderInvocation_DrawRectOneTexture(const GPUImageBase& image, int x, int y) { // New Name? DrawRectAligned_OneTexture
	float width = image.width();
	float height = image.height();
	float t0_x1 = image.x(0);
	float t0_x2 = image.x(width);
	float t0_y1 = image.y(0);
	float t0_y2 = image.y(height);

	glBegin(GL_QUADS);							  
	glTexCoord2f(t0_x1, t0_y1); 	  glVertex2f(x, y);	
	glTexCoord2f(t0_x2, t0_y1); 	  glVertex2f(x + width, y);	
	glTexCoord2f(t0_x2, t0_y2);	  glVertex2f(x + width, y + height);	
	glTexCoord2f(t0_x1, t0_y2); 	  glVertex2f(x, y + height);
	glEnd();
  }

  void ShaderInvocation_DrawRectOneTexture(const GPUImageBase& image, int x, int y, int width, int height) { 
	float t0_x1 = image.x(0);
	float t0_x2 = image.x(width);
	float t0_y1 = image.y(0);
	float t0_y2 = image.y(height);

	glBegin(GL_QUADS);							  
	glTexCoord2f(t0_x1, t0_y1); 	  glVertex2f(x, y);	
	glTexCoord2f(t0_x2, t0_y1); 	  glVertex2f(x + width, y);	
	glTexCoord2f(t0_x2, t0_y2);	  glVertex2f(x + width, y + height);	
	glTexCoord2f(t0_x1, t0_y2); 	  glVertex2f(x, y + height);
	glEnd();
  }

  void ShaderInvocation_DrawRectOneTexture(const GPUImageBase& image, const Rectangle2D<int>& rect, const Rectangle2D<int>& texCoord) { 
	float t0_x1 = image.x(texCoord.x);
	float t0_x2 = image.x(texCoord.x + texCoord.width);
	float t0_y1 = image.y(texCoord.y);
	float t0_y2 = image.y(texCoord.y + texCoord.height);

	glBegin(GL_QUADS);							  
	glTexCoord2f(t0_x1, t0_y1); 	  glVertex2f(rect.x, rect.y);	
	glTexCoord2f(t0_x2, t0_y1); 	  glVertex2f(rect.x + rect.width, rect.y);	
	glTexCoord2f(t0_x2, t0_y2);	  glVertex2f(rect.x + rect.width, rect.y + rect.height);	
	glTexCoord2f(t0_x1, t0_y2); 	  glVertex2f(rect.x, rect.y + rect.height);
	glEnd();
  }



  void ShaderInvocation_CleanupGLState() {
	glDisable(GL_TEXTURE_RECTANGLE_ARB); 
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);    
	glPopMatrix();
  }


} } // namespaces GPU, vw
