
#include <vw/GPU/ImageStatistics.h>
#include <vw/GPU/ImageAlgorithms.h>
#include <vw/GPU/SimpleTransforms.h>

#include <vw/Vector.h>
//#include <vw/Geometry/VectorFixed.h>
//#include <vw/Geometry/PointND.h>

namespace vw { namespace GPU {

#define MIN(n1, n2) n1 < n2 ? n1 : n2
#define MAX(n1, n2) n1 > n2 ? n1 : n2

float min_channel_value(const GPUImageBase& image)			     
{
// Static
  static GPUProgramSet programSet_Level1("", "ImageStatistics/min-channels");
  static GPUProgramSet programSet_OtherLevels("", "ImageStatistics/min-quad");
  static vector<int> fAttributes(1);
// Output Image Size
  int POT_Level = MAX((int) ceilf(log2f(image.width())), (int) ceilf(log2f(image.height())));
  int POT_Dimension = pow(2.0, POT_Level);
// Reduce to 1 Channell; pad to POT dimensions; fill padded area with high value
  GPUImageBase tex_StartLevel(MAX(1, POT_Dimension), MAX(1, POT_Dimension), TEX_R, image.type());
  fill(tex_StartLevel,  0.9, 0, 0, 0); 

  fAttributes[0] = image.num_channels();
  GPUProgram* program_Level1 = programSet_Level1.get_program(vector<int>(), fAttributes, false);
  program_Level1->install();
  ((GPUImageBase&) image).rasterize_homography();		
  ShaderInvocation_SetupGLState(tex_StartLevel);
  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, tex_StartLevel.target(), tex_StartLevel.name(), 0);	
  program_Level1->set_uniform_texture("image", 0, image);
  glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  int draw_size = image.width();
  ShaderInvocation_DrawRectOneTexture(image, Rectangle2D<int>(0, 0, image.width(), image.height()), 
  			      Rectangle2D<int>(0, 0, image.width(), image.height()));
  
  program_Level1->uninstall();

  // Fill edge extend area
  float high_value = 10000000.0;
  int xEdgeSize = tex_StartLevel.width() - image.width();
  int yEdgeSize = tex_StartLevel.height() - image.height();
  GPUImageBase tex_Cropped;
  if(xEdgeSize) {
    tex_Cropped = GPU::crop(tex_StartLevel, image.width(), 0, xEdgeSize, tex_StartLevel.height());
    fill(tex_Cropped, high_value, 0, 0, 0);
  }
  if(yEdgeSize) {
    tex_Cropped = GPU::crop(tex_StartLevel, 0, image.height(), image.width(), yEdgeSize);
    fill(tex_Cropped, high_value, 0, 0, 0);
  }
// ITERATIONS: 1 to n
  GPUImageBase tex_PreviousLevel = tex_StartLevel;

  float output;
  fAttributes[0] = 4;
  GPUProgram* program = programSet_OtherLevels.get_program(vector<int>(), fAttributes, false);
  program->install();

  for(int iLevel = 1; iLevel <= POT_Level; iLevel++) {
    int reductionRatio = (int) powf(2.0, iLevel);
    GPUImageBase tex_CurrentLevel(MAX(1, POT_Dimension / reductionRatio), MAX(1, POT_Dimension / reductionRatio), TEX_R, image.type());
    ShaderInvocation_SetupGLState(tex_CurrentLevel);
    glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, tex_CurrentLevel.target(), tex_CurrentLevel.name(), 0);	
    program->set_uniform_texture("image", 0, tex_PreviousLevel);
    glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_T, GL_CLAMP);
    ShaderInvocation_DrawRectOneTexture(tex_PreviousLevel, Rectangle2D<int>(0, 0, tex_CurrentLevel.width(), tex_CurrentLevel.height()), 
					Rectangle2D<int>(0, 0, tex_CurrentLevel.width() * 2, tex_CurrentLevel.height() * 2));

    tex_PreviousLevel = tex_CurrentLevel;
  }
  program->uninstall();
  // CleanUp State
  ShaderInvocation_CleanupGLState();
  tex_PreviousLevel.read(TEX_R, TEX_FLOAT32, &output);
  return output;
}



float max_channel_value(const GPUImageBase& image)			     
{
// Static
  static GPUProgramSet programSet_Level1("", "ImageStatistics/max-channels");
  static GPUProgramSet programSet_OtherLevels("", "ImageStatistics/max-quad");
  static vector<int> fAttributes(1);
// Output Image Size
  int POT_Level = MAX((int) ceilf(log2f(image.width())), (int) ceilf(log2f(image.height())));
  int POT_Dimension = pow(2.0, POT_Level);
// GLState - Setup
  ((GPUImageBase&) image).rasterize_homography();	
  ShaderInvocation_SetupGLState(POT_Dimension, POT_Dimension);
// Reduce to 1 Channell; pad to POT dimensions; fill padded area with high value
  fAttributes[0] = image.num_channels();
  GPUProgram* program_Level1 = programSet_Level1.get_program(vector<int>(), fAttributes, false);
  program_Level1->install();
  GPUImageBase tex_StartLevel(MAX(1, POT_Dimension), MAX(1, POT_Dimension), TEX_R, image.type());
  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, tex_StartLevel.target(), tex_StartLevel.name(), 0);	
  program_Level1->set_uniform_texture("image", 0, image);
  ShaderInvocation_DrawRectOneTexture(image, Rectangle2D<int>(0, 0, tex_StartLevel.width(), tex_StartLevel.height()), 
				      Rectangle2D<int>(0, 0, tex_StartLevel.width(), tex_StartLevel.height()));
  program_Level1->uninstall();
  // Fill edge extend area
  float pad_value = -10000000.0;
  int xEdgeSize = tex_StartLevel.width() - image.width();
  int yEdgeSize = tex_StartLevel.height() - image.height();
  GPUImageBase tex_Cropped;
  if(xEdgeSize) {
    tex_Cropped = GPU::crop(tex_StartLevel, image.width(), 0, xEdgeSize, tex_StartLevel.height());
    fill(tex_Cropped, pad_value, 0, 0, 0);
  }
  if(yEdgeSize) {
    tex_Cropped = GPU::crop(tex_StartLevel, 0, image.height(), image.width(), yEdgeSize);
    fill(tex_Cropped, pad_value, 0, 0, 0);
  }
// ITERATIONS: 1 to n
  GPUImageBase tex_PreviousLevel = tex_StartLevel;
  fAttributes[0] = 4;
  GPUProgram* program = programSet_OtherLevels.get_program(vector<int>(), fAttributes, false);
  program->install();

  for(int iLevel = 1; iLevel <= POT_Level; iLevel++) {
    int reductionRatio = (int) powf(2.0, iLevel);
    GPUImageBase tex_CurrentLevel(MAX(1, POT_Dimension / reductionRatio), MAX(1, POT_Dimension / reductionRatio), TEX_R, image.type());
    ShaderInvocation_SetupGLState(tex_CurrentLevel);
    glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, tex_CurrentLevel.target(), tex_CurrentLevel.name(), 0);	
    program->set_uniform_texture("image", 0, tex_PreviousLevel);
    glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_T, GL_CLAMP);
    ShaderInvocation_DrawRectOneTexture(image, Rectangle2D<int>(0, 0, tex_CurrentLevel.width(), tex_CurrentLevel.height()), 
					Rectangle2D<int>(0, 0, tex_CurrentLevel.width() * 2, tex_CurrentLevel.height() * 2));

    tex_PreviousLevel = tex_CurrentLevel;
  }
  program->uninstall();
  // CleanUp State
  ShaderInvocation_CleanupGLState();
  float output;
  tex_PreviousLevel.read(TEX_R, TEX_FLOAT32, &output);
  return output;
}



  void min_max_channel_values(const GPUImageBase& image, float& min, float& max) {
    min = min_channel_value(image);
    max = max_channel_value(image);
  }


float mean_channel_value(const GPUImageBase& image)			     
{
// Static
  static GPUProgramSet programSet_Level1("", "ImageStatistics/sum-channels");
  static GPUProgramSet programSet_OtherLevels("", "ImageStatistics/sum-quad");
  static vector<int> fAttributes(1);
// Output Image Size
  int POT_Level = MAX((int) ceilf(log2f(image.width())), (int) ceilf(log2f(image.height())));
  int POT_Dimension = pow(2.0, POT_Level);
// GLState - Setup
	((GPUImageBase&) image).rasterize_homography();	
  ShaderInvocation_SetupGLState(POT_Dimension, POT_Dimension);
// Reduce to 1 Channell; pad to POT dimensions; fill padded area with high value
  fAttributes[0] = image.num_channels();
  GPUProgram* program_Level1 = programSet_Level1.get_program(vector<int>(), fAttributes, false);
  program_Level1->install();
  GPUImageBase tex_StartLevel(MAX(1, POT_Dimension), MAX(1, POT_Dimension), TEX_R, image.type());
  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, tex_StartLevel.target(), tex_StartLevel.name(), 0);	
  program_Level1->set_uniform_texture("image", 0, image);
  ShaderInvocation_DrawRectOneTexture(image, Rectangle2D<int>(0, 0, tex_StartLevel.width(), tex_StartLevel.height()), 
				      Rectangle2D<int>(0, 0, tex_StartLevel.width(), tex_StartLevel.height()));
  program_Level1->uninstall();
  // Fill edge extend area
  float pad_value = 0.0;
  int xEdgeSize = tex_StartLevel.width() - image.width();
  int yEdgeSize = tex_StartLevel.height() - image.height();
  GPUImageBase tex_Cropped;
  if(xEdgeSize) {
    tex_Cropped = GPU::crop(tex_StartLevel, image.width(), 0, xEdgeSize, tex_StartLevel.height());
    fill(tex_Cropped, pad_value, 0, 0, 0);
  }
  if(yEdgeSize) {
    tex_Cropped = GPU::crop(tex_StartLevel, 0, image.height(), image.width(), yEdgeSize);
    fill(tex_Cropped, pad_value, 0, 0, 0);
  }
// ITERATIONS: 1 to n
  GPUImageBase tex_PreviousLevel = tex_StartLevel;
  fAttributes[0] = 4;
  GPUProgram* program = programSet_OtherLevels.get_program(vector<int>(), fAttributes, false);
  program->install();

  for(int iLevel = 1; iLevel <= POT_Level; iLevel++) {
    int reductionRatio = (int) powf(2.0, iLevel);
    GPUImageBase tex_CurrentLevel(MAX(1, POT_Dimension / reductionRatio), MAX(1, POT_Dimension / reductionRatio), TEX_R, image.type());
    ShaderInvocation_SetupGLState(tex_CurrentLevel);
    glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, tex_CurrentLevel.target(), tex_CurrentLevel.name(), 0);	
    program->set_uniform_texture("image", 0, tex_PreviousLevel);
    glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_T, GL_CLAMP);
    ShaderInvocation_DrawRectOneTexture(image, Rectangle2D<int>(0, 0, tex_CurrentLevel.width(), tex_CurrentLevel.height()), 
					Rectangle2D<int>(0, 0, tex_CurrentLevel.width() * 2, tex_CurrentLevel.height() * 2));

    tex_PreviousLevel = tex_CurrentLevel;
  }
  program->uninstall();
  // CleanUp State
  ShaderInvocation_CleanupGLState();
  float output = ((GPUImage<PixelGray<float> >&) tex_PreviousLevel)(0,0).v();
  return output / (float) (image.width() * image.height() * image.num_channels());
}


  } // namespace GPU
} // namespace vw
