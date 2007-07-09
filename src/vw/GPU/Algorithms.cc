
#include <vw/GPU/Algorithms.h>

namespace vw { namespace GPU {

void fill_impl(GPUImageBase& image, float red, float green, float blue, float alpha, const Rectangle2D<int>& bounds)
{
// Static
  static GPUProgramSet programSet("", "Algorithms/fill-1i4f");
  static vector<int> fAttributes(1);
// GLState - Setup
  ((GPUImageBase&) image).rasterize_homography();	
  ShaderInvocation_SetupGLState(image);
// Program - Install
  fAttributes[0] = 4;
  GPUProgram* program = programSet.get_program(vector<int>(), fAttributes, false);
  program->install();
  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, image.target(), image.name(), 0);
  // INPUTS
  bool is_uint8 = image.type() == TEX_UINT8;
  program->set_uniform_float("f1", red, is_uint8);
  program->set_uniform_float("f2", green, is_uint8);
  program->set_uniform_float("f3", blue, is_uint8);
  program->set_uniform_float("f4", alpha, is_uint8);
  // DRAW
  int x1 = (int) image.x(bounds.x);
  int x2 = (int) x1 + bounds.width;
  int y1 = (int) image.y(bounds.y);
  int y2 = (int) (y1 + bounds.height);
  glBegin(GL_QUADS);
  glVertex2f(x1, y1);
  glVertex2f(x2, y1);
  glVertex2f(x2, y2);
  glVertex2f(x1, y2);
  glEnd();
  // CleanUp State
  program->uninstall();
  ShaderInvocation_CleanupGLState();
}

  } // namespace GPU
} // namespace vw

