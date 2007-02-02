
#include "Filter.h"
#include "Timer.h"
#include <vector>


namespace vw { namespace GPU {


//########################################################################
//#    convolution_filter                    
//########################################################################


GPUImageBase convolution_filter(const GPUImageBase& image, 
			       const GPUImageBase& kernel)
{
// Static
  static GPUProgramSet programSet;
  static int init=1;
  if(init--)
    programSet.set_base_paths("", "Filter/convolution");
  static vector<int> fAttributes(3);
  static vector<int> emptyVector;

// GLState - Setup
  ((GPUImageBase&) image).rasterize_homography();
  ShaderInvocation_SetupGLState(image.width(), image.height());
// Program - Install
  fAttributes[0] = 4;
  fAttributes[1] = kernel.width();
  fAttributes[2] = kernel.height();
  GPUProgram* program = programSet.get_program(emptyVector, fAttributes, false);
  program->install();
  // OUTPUT
  GPUImageBase temp(image.width(), image.height(), image.format(), image.type());
  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, temp.target(), temp.name(), 0);	
  // INPUT
  program->set_uniform_texture("image", 0, image);
  program->set_uniform_texture("kernel", 1, kernel);
  program->set_uniform_float("hHalfSize", (kernel.width() / 2) - 1);
  program->set_uniform_float("vHalfSize", (kernel.height() / 2) - 1);
  // DRAW
  ShaderInvocation_DrawRectOneTexture(image);
  // CleanUp State
  program->uninstall();
  ShaderInvocation_CleanupGLState();
  
  return temp;
}


//########################################################################
//#    separable_convolution_filter                    
//########################################################################


GPUImageBase separable_convolution_filter(const GPUImageBase& image, 
					       const GPUImageBase& hKernel, 
					       const GPUImageBase& vKernel)
{
// Static
  static vector<int> fAttributes(2);
  static GPUProgramSet programSet_Rows;
  static GPUProgramSet programSet_Columns;
  static int needsInit = 1;
  if(needsInit--) {
    programSet_Rows.set_base_paths("", "Filter/convolution-rows");
    programSet_Columns.set_base_paths("", "Filter/convolution-columns");
  }
  vector<int> emptyVector;
// GLState - Setup
  ((GPUImageBase&) image).rasterize_homography();
  ShaderInvocation_SetupGLState(image.width(), image.height());
// *** STAGE 1 - Rows ***
  // Program - Install
  fAttributes[0] = 4;
  fAttributes[1] = hKernel.width();
  GPUProgram* program = programSet_Rows.get_program(emptyVector, fAttributes, false);
  program->install();
  // OUTPUT
  GPUImageBase temp(image.width(), image.height(), image.format(), image.type());
  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, temp.target(), temp.name(), 0);	
  // INPUT
  program->set_uniform_texture("image", 0, image);
  program->set_uniform_texture("kernel", 1, hKernel);
  program->set_uniform_float("halfSize", (hKernel.width() / 2) - 1);
  // DRAW
  ShaderInvocation_DrawRectOneTexture(image);

// *** STAGE 2 - Columns ***
  // Program - Install
  fAttributes[0] = 4;
  fAttributes[1] = vKernel.width();
  program = programSet_Columns.get_program(emptyVector, fAttributes, false);
  program->install();
  // OUTPUT
  GPUImageBase temp2(image.width(), image.height(), image.format(), image.type());
  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, temp2.target(), temp2.name(), 0);	
  // INPUT
  program->set_uniform_texture("image", 0, temp);
  program->set_uniform_texture("kernel", 1, vKernel);
  program->set_uniform_float("halfSize", (vKernel.width() / 2) - 1);
  // DRAW
  ShaderInvocation_DrawRectOneTexture(temp);

  // CleanUp State
  program->uninstall();
  ShaderInvocation_CleanupGLState();
  
  return temp2;

}


} } // namespaces GPU, vw
