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


#ifndef GPU_EdgeExtend_H
#define GPU_EdgeExtend_H


#include <vw/Image/EdgeExtension.h>

namespace vw { namespace GPU {

  enum EdgeExtensionType {
    ZERO_EDGE_EXTENSION,
    CONSTANT_EDGE_EXTENSION
  };

   template <class EdgeT>
     struct TraitsForEdgeT { };

   template <>
     struct TraitsForEdgeT<ZeroEdgeExtension> {
       static const EdgeExtensionType type = ZERO_EDGE_EXTENSION;
     };

   template <>
     struct TraitsForEdgeT<ConstantEdgeExtension> {
       static const EdgeExtensionType type = CONSTANT_EDGE_EXTENSION;
     };


   inline void EdgeExtension_SetupTexture(EdgeExtensionType type) {
     GLuint gl_wrap_type;
     if(type == ZERO_EDGE_EXTENSION) {
       gl_wrap_type = GL_CLAMP_TO_BORDER_ARB;

     }
     else if(type == CONSTANT_EDGE_EXTENSION) {
       gl_wrap_type = GL_CLAMP;
     }
     glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_S, gl_wrap_type);
     glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_T, gl_wrap_type);

   }

   inline void EdgeExtension_RestoreTexture(EdgeExtensionType type) {
     glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER_ARB);
     glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER_ARB);
   }

   typedef ZeroEdgeExtension DefaultEdgeExtension;

} } // namespaces GPU, vw


#endif
