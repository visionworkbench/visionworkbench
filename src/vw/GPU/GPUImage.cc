// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/GPU/GPUImage.h>
#include <vw/GPU/TexAlloc.h>
#include <vw/GPU/GPUProgram.h>
#include <vw/GPU/Setup.h>
#include <vw/Image/Transform.h>
#include <map>
#include <string>
#include <vw/Math/Vector.h>

namespace vw {
namespace GPU {

    void GPUImageBase::rasterize_homography() const {
      if(!m_isHomography) return;
// GLState
      ShaderInvocation_SetupGLState(width(), height());
// Program
      GPUProgram* program = create_gpu_program(m_interpolation_string);
      program->install();
// OUTPUT
      GPUImageBase temp(width(), height(), format(), type());
      ShaderInvocation_SetOutputImage(temp);
// INPUT
      program->set_input_image("image", *this);
      EdgeExtension_SetupTexture(m_edge_extension_type);
// DRAW
      HomographyTransform h_functor(m_homography);

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





