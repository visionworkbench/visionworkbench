// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// This deprecated header file defines a number of deprecated typdefs
// for the edge extension modes in EdgeExtension.h.
//
#ifndef __VW_IMAGE_EDGEEXTEND_H__
#define __VW_IMAGE_EDGEEXTEND_H__

#include <vw/Image/EdgeExtension.h>

#warning The <vw/Image/EdgeExtend.h> header is deprecated: use <vw/Image/EdgeExtension.h> instead!

namespace vw {

  typedef EdgeExtensionBase EdgeExtendBase VW_DEPRECATED;
  typedef NoEdgeExtension NoEdgeExtend VW_DEPRECATED;
  typedef ZeroEdgeExtension ZeroEdgeExtend VW_DEPRECATED;
  typedef ConstantEdgeExtension ConstantEdgeExtend VW_DEPRECATED;
  typedef PeriodicEdgeExtension PeriodicEdgeExtend VW_DEPRECATED;
  typedef ReflectEdgeExtension ReflectEdgeExtend VW_DEPRECATED;

  template <class ImageT, class ExtensionT>
  class VW_DEPRECATED EdgeExtendView : public EdgeExtensionView<ImageT,ExtensionT>
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
