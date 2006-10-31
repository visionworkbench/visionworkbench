 
#ifndef GPU_EdgeExtend_H
#define GPU_EdgeExtend_H


#include <vw/Image/EdgeExtend.h>

namespace vw { namespace GPU {

  enum EdgeExtensionType {
    ZERO_EDGE_EXTEND,
    CONSTANT_EDGE_EXTEND
  };

   template <class EdgeT> 
     struct TraitsForEdgeT { };

   template <> 
     struct TraitsForEdgeT<ZeroEdgeExtend> {
       static const EdgeExtensionType type = ZERO_EDGE_EXTEND;
     };

   template <> 
     struct TraitsForEdgeT<ConstantEdgeExtend> {
       static const EdgeExtensionType type = CONSTANT_EDGE_EXTEND;
     };


   inline void EdgeExtend_SetupTexture(EdgeExtensionType type) {
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

   inline void EdgeExtend_RestoreTexture(EdgeExtensionType type) {
     glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER_ARB);
     glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER_ARB);
   }

   typedef ZeroEdgeExtend DefaultEdgeExtend;

} } // namespaces GPU, vw


#endif
