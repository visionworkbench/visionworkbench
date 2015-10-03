// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NASA Vision Workbench is licensed under the Apache License,
//  Version 2.0 (the "License"); you may not use this file except in
//  compliance with the License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
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





