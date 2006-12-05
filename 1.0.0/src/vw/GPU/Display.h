
#ifndef Display_H
#define Display_H


#include <map>

#include <vw/ImageView.h>
#include <vw/GPU/OpenGLManager.h>
#include <vw/GPU/GPUImage.h>

/*
#ifdef __APPLE__
#include <GLUT/glut.h>
#include <OpenGL/glext.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else 
#include <GL/glut.h>
#include <GL/glext.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif
*/


namespace vw { namespace GPU {
  
  class Display {
  public:
    static map<int, Display*> map_WindowToDisplay;
    GPUImageBase texture;
    int width;
    int height;
  public:
    static Display* GetDisplayForWindow(int windowID) {
      map<int, Display*>::iterator iter = map_WindowToDisplay.find(windowID);
      if(iter != map_WindowToDisplay.end())
	return (*iter).second;
      else
	return NULL; 
    }
    static void RunGLUTMainLoop() { glutMainLoop(); }
    
    Display(GPUImageBase& tex, char* name);
    ~Display();
    void Init(GPUImageBase& tex, char* name);
    

  };

  void glut_display();

  void glut_timer(int value);

  void glut_idle();

} // GPU
} // vw

#endif
