// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
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
