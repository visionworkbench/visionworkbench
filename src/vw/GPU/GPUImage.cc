
#include "GPUImage.h"
#include "TexAlloc.h"
#include <vw/GPU/GPUProgram.h>
#include <vw/GPU/Setup.h>
#include <vw/Image/Transform.h>
#include <map>
#include <string>
#include <vw/Math/Vector.h>


using namespace std;

namespace vw { namespace GPU {

    void GPUImageBase::rasterize_homography() const {
      if(!_isHomography) return;
//Static
	static vector<int> fAttributes(1);
	static map<string, GPUProgramSet*> cachedProgramSets;
// GLState
	ShaderInvocation_SetupGLState(width(), height());
// Program
	GPUProgramSet* programSet;
	map<string, GPUProgramSet*>::iterator iter_map = cachedProgramSets.find(_interpolation_string);
	if(iter_map != cachedProgramSets.end()) {
	  programSet = (*iter_map).second;
	}
	else {
	  programSet = new GPUProgramSet("", _interpolation_string.c_str());
	  cachedProgramSets[_interpolation_string] = programSet;
	}
	fAttributes[0] = 4;
	GPUProgram* program = programSet->get_program(vector<int>(), fAttributes, false);
	program->install();
// OUTPUT
	GPUImageBase temp(width(), height(), format(), type());

	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, g_framebuffer);
	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, temp.target(), temp.name(), 0);	
// INPUT
	program->set_uniform_texture("image", 0, *this);
	EdgeExtension_SetupTexture(_edge_extension_type);
// DRAW
	//Matrix<float> h = inverse(_homography);	
	HomographyTransform h_functor(_homography);
      
	Vector2 t_0_0 = h_functor.forward(Vector2(-0.5, -0.5));
	Vector2 t_1_0 = h_functor.forward(Vector2(width()-0.5, -0.5)); 
	Vector2 t_1_1 = h_functor.forward(Vector2(width()-0.5, height()-0.5));
	Vector2 t_0_1 = h_functor.forward(Vector2(-0.5, height()-0.5));
	
	glBegin(GL_QUADS);							  
	glTexCoord2f(t_0_0[0], t_0_0[1]);  glVertex2f(0, 0);
	glTexCoord2f(t_1_0[0], t_1_0[1]);  glVertex2f(width(), 0);
       	glTexCoord2f(t_1_1[0], t_1_1[1]);  glVertex2f(width(), height());	
	glTexCoord2f(t_0_1[0], t_0_1[1]);  glVertex2f(0, height());		
	glEnd();
    // Clean Up
	program->uninstall();
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
	*((GPUImageBase*) this) = temp;
    }
	
	

// GPUImageBase Members

} } // namespaces GPU, vw




 
