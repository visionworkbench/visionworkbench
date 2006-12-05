

#include "vw/GPU/Display.h"

  // GLUT Callbacks

#include <vw/GPU/GPUProgram.h>
#include <vw/FileIO.h>
#include <vw/GPU/ImageOperators.h>
#include <vw/GPU/Timer.h>
#include <vw/GPU/Transform.h>
#include <vw/GPU/Filter.h>
#include <vw/GPU/ImageMath.h>
#include <vw/GPU/SimpleTransforms.h>

namespace vw { namespace GPU {

  /*
  void glut_timer(int value) {
    static int c = 1;
    glutPostRedisplay();
  }

  void glut_idle() {
    static int c = 1;
    //glutPostRedisplay();
  }

  void glut_display() {
    glutTimerFunc(15, glut_timer, 5);
// Get Display Object
    
    static bool init = false;
    static GPUImageBase image_orig;
    static GPUImageBase image;
    static Timer timer;
    if(!init) {
      timer.Start();
      init = true;
      ImageView<PixelRGB<float> > image_color;
      read_image(image_color, "tests/test_images/lighthouse.jpg");
      //write_image("tests/test_images/out_image_flowers.jpg", image_color);
      printf("image %i %i \n", image_color.cols(), image_color.rows());
      //GPUImageBase tex_new = image_color;
      //GPUImageBase tex_new(image_color.cols(), image_color.rows(), TEX_RGB, TEX_FLOAT32, TEX_RGB, TEX_FLOAT32, &(image_color(0, 0))); 
      //image_orig = tex_new; 
      image_orig = image_color;
    }
    
    timer.Stop();
    float elapsedTime = timer.ElapsedSeconds();
    static int counter=0;

    float value = sinf(elapsedTime) + 1.0;
    //image = gaussian_filter(image_orig, value, value, 0, 0);
    image = gaussian_filter(GPU::crop(fixed_rotate(pow(image_orig, value), elapsedTime), 300, 300, 400, 300), 1, 1, 0, 0);

   
//Static
	static vector<int> fAttributes(1);
	static GPUProgramSet_GLSL programSet_Identity;
	static int needsInit = 1;
	if(needsInit--)
	  programSet_Identity.set_base_paths("", "Identity");
	static vector<int> emptyVector;
// GLState
    glClear(GL_DEPTH_BUFFER_BIT|GL_COLOR_BUFFER_BIT);
    glPolygonMode(GL_FRONT,GL_FILL);
    glColor3f(0,1,0);
// Variables
    int width = image.width();
	int height = image.height();

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0,width,height, 0.0);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glViewport(0,0,width, height);

    
    glEnable(GL_TEXTURE_RECTANGLE_ARB);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
// Resize
	glutReshapeWindow(width, height);

// Viewport
// Program
	fAttributes[0] = 4;
	GPUProgram* program_identity = programSet_Identity.get_program(emptyVector, fAttributes, false);
	ClearGLError();
	program_identity->install();
	CheckGLError();
// OUTPUT
	GPUImageBase temp(width, height, image.format(), image.type());
// INPUT\
	glEnable(GL_TEXTURE_RECTANGLE_ARB); 
	ClearGLError();

	program_identity->set_uniform_texture("image", 0, image);
	CheckGLError();

	glActiveTexture(GL_TEXTURE0);
	glEnable(GL_TEXTURE_RECTANGLE_ARB); 

// DRAW	
	
	float t0_x1 = image.x(0);
	float t0_x2 = image.x(width);
	float t0_y1 = image.y(0);
	float t0_y2 = image.y(height);
	
	glBegin(GL_QUADS);							  
	glTexCoord2f(t0_x1, t0_y1);  glVertex2f(0, 0);
	glTexCoord2f(t0_x2, t0_y1);  glVertex2f(width, 0);
       	glTexCoord2f(t0_x2, t0_y2);  glVertex2f(width, height);	
	glTexCoord2f(t0_x1, t0_y2);  glVertex2f(0, height);		
	glEnd();

// Uninstall
    program_identity->uninstall();
// Flush
    //glFlush();
    glutSwapBuffers();
  }

  */
    //
  /*
  void glut_display() {
    //Static
    static vector<int> fAttributes(1);
    static GPUProgramSet_GLSL programSet_Identity;
    static int needsInit = 1;
    if(needsInit--)
      programSet_Identity.set_base_paths("", "Identity");
    static vector<int> emptyVector;
    // Get Display Object
    Display* display = vw::GPU::Display::GetDisplayForWindow(glutGetWindow());
    GPUImageBase& tex = display->texture;
    // Setup
    int width = tex.width();
    int height = tex.height();
    glutReshapeWindow(width, height);

    glClear(GL_DEPTH_BUFFER_BIT|GL_COLOR_BUFFER_BIT);
    glEnable(GL_TEXTURE_RECTANGLE_ARB);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
    // Program
    fAttributes[0] = 4;
    GPUProgram* program_identity = programSet_Identity.get_program(emptyVector, fAttributes, true);
    program_identity->install();
    // Input
    program_identity->set_uniform_texture("image", 0, tex);
    // Draw
    glBegin(GL_QUADS);
    glTexCoord2f(0, 0);   glVertex2f(0, 0);
    glTexCoord2f(width, 0);   glVertex2f(width, 0);
    glTexCoord2f(width, height);   glVertex2f(width, height);
    glTexCoord2f(0, height);   glVertex2f(0, height);
    glEnd();
    // Uninstall
    program_identity->uninstall();
    // Flush
    glutSwapBuffers();
  }
  */  
  void glut_reshape(int w, int h) {
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0,w,0.0, h);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glViewport(0,0,w, h);

  }

  // Display Members

  map<int, Display*> Display::map_WindowToDisplay;
 
  Display::Display(GPUImageBase& tex, char* name) {
    Init(tex, name);
  }
  
  Display::~Display() {

  }

  void Display::Init(GPUImageBase& tex, char* name) {
    /*
    texture = tex;
    glutInitWindowPosition(30, 30);
    glutInitWindowSize(tex.width(), tex.height());
    int windowID = glutCreateWindow(name);
    map_WindowToDisplay[windowID] = this;
    glutDisplayFunc(glut_display);
    glutReshapeFunc(glut_reshape);
    //glutTimerFunc(30, glut_timer, 1);
    */

  }
    

} // GPU
} // vw
