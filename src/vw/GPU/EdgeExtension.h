 
#ifndef GPU_EdgeExtension_H
#define GPU_EdgeExtension_H


#include <vw/Image/EdgeExtension.h>

namespace vw { namespace GPU {

  enum EdgeExtensionType {
    ZERO_EDGE_EXTEND,
    CONSTANT_EDGE_EXTEND
  };

   template <class EdgeT> 
     struct TraitsForEdgeT { };

   template <> 
     struct TraitsForEdgeT<ZeroEdgeExtension> {
       static const EdgeExtensionType type = ZERO_EDGE_EXTEND;
     };

   template <> 
     struct TraitsForEdgeT<ConstantEdgeExtension> {
       static const EdgeExtensionType type = CONSTANT_EDGE_EXTEND;
     };


   inline void EdgeExtension_SetupTexture(EdgeExtensionType type) {
     GLuint gl_wrap_type;
     if(type == ZERO_EDGE_EXTEND) {
       gl_wrap_type = GL_CLAMP_TO_BORDER_ARB;

     }
     else if(type == CONSTANT_EDGE_EXTEND) {
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
