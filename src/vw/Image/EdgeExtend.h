// __BEGIN_LICENSE__
// 
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
// 
// Copyright 2006 Carnegie Mellon University. All rights reserved.
// 
// This software is distributed under the NASA Open Source Agreement
// (NOSA), version 1.3.  The NOSA has been approved by the Open Source
// Initiative.  See the file COPYING at the top of the distribution
// directory tree for the complete NOSA document.
// 
// THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
// KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
// LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
// SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
// A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
// THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
// DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
// 
// __END_LICENSE__

// This deprecated header file defines a number of deprecated typdefs 
// for the edge extension modes in EdgeExtension.h.
//
#ifndef __VW_IMAGE_EDGEEXTEND_H__
#define __VW_IMAGE_EDGEEXTEND_H__

#include <vw/Image/EdgeExtension.h>

#warning The <vw/Image/EdgeExtend.h> header is deprecated: use <vw/Image/EdgeExtension.h> instead!

namespace vw {
 
  typedef EdgeExtensionBase EdgeExtendBase;
  typedef NoEdgeExtension NoEdgeExtend;
  typedef ZeroEdgeExtension ZeroEdgeExtend;
  typedef ConstantEdgeExtension ConstantEdgeExtend;
  typedef PeriodicEdgeExtension PeriodicEdgeExtend;
  typedef ReflectEdgeExtension ReflectEdgeExtend;

  template <class ImageT, class ExtensionT>
  class EdgeExtendView : public EdgeExtensionView<ImageT,ExtensionT>
  {
    typedef EdgeExtensionView<ImageT,ExtensionT> impl_type;
  public:
    EdgeExtendView( ImageT const& image )
      : impl_type( image ) {}
    
    EdgeExtendView( ImageT const& image, ExtensionT const& extension )
      : impl_type( image, extension ) {}
    
    EdgeExtendView( ImageT const& image, ptrdiff_t xoffset, ptrdiff_t yoffset, int32 cols, int32 rows )
      : impl_type( image, xoffset, yoffset, cols, rows ) {}
    
    EdgeExtendView( ImageT const& image, ptrdiff_t xoffset, ptrdiff_t yoffset, int32 cols, int32 rows, ExtensionT const& extension )
      : impl_type( image, xoffset, yoffset, cols, rows, extension ) {}

    EdgeExtendView( impl_type const& other ) : impl_type( other ) {}
  };
  
} // namespace vw

#endif // __VW_IMAGE_EDGEEXTEND_H__
