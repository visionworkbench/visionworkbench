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

/// \file UtilityViews.h
/// 
/// Utility objects conforming to the view concept.
///
#ifndef __VW_IMAGE_UTILITYVIEWS_H__
#define __VW_IMAGE_UTILITYVIEWS_H__

#include <boost/smart_ptr.hpp>
#include <boost/type_traits.hpp>

#include <vw/Core/Exception.h>

#include <vw/Image/ImageViewBase.h>
#include <vw/Image/ImageResource.h>
#include <vw/Image/PixelAccessors.h>

namespace vw {

  // ConstantView
  //
  // A view that returns a fixed constant value.

  template <class PixelT>
  class ConstantView : public ImageViewBase<ConstantView<PixelT> > {
    PixelT m_value;
    int32 m_cols, m_rows, m_planes;
  public:
    typedef PixelT pixel_type;
    typedef PixelT const& result_type;
    typedef ProceduralPixelAccessor<ConstantView> pixel_accessor;

    ConstantView( PixelT const& value, int32 cols, int32 rows, int32 planes = 1 )
      : m_value(value), m_cols(cols), m_rows(rows), m_planes(planes) {}

    inline int32 cols() const { return m_cols; }
    inline int32 rows() const { return m_rows; }
    inline int32 planes() const { return m_planes; }

    inline pixel_accessor origin() const { return pixel_accessor( *this ); }

    inline result_type operator()( int32 col, int32 row, int32 plane=0 ) const { return m_value; }

    typedef ConstantView prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i /*bbox*/ ) const { return *this; }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i bbox ) const {
      vw::rasterize( prerasterize(bbox), dest, bbox );
    }
  };

  template <class PixelT>
  struct IsMultiplyAccessible<ConstantView<PixelT> > : public true_type {};

  template <class PixelT>
  ConstantView<PixelT> constant_view( PixelT const& value, int32 cols, int32 rows, int32 planes=1 ) {
    return ConstantView<PixelT>( value, cols, rows, planes );
  }
  

  /// PixelIndexView
  ///
  /// A view that returns a Vector2 representing the pixel location.
  class PixelIndexView : public ImageViewBase<PixelIndexView> {
    int32 m_cols, m_rows, m_planes;
  public:
    typedef Vector2 pixel_type;
    typedef Vector2 result_type;
    typedef ProceduralPixelAccessor<PixelIndexView> pixel_accessor;

    /// Initialize from another view
    template <class ImageT>
    PixelIndexView( ImageViewBase<ImageT> const& view ) : 
      m_rows(view.rows()), m_cols(view.cols()), m_planes(view.planes()) {}

    PixelIndexView( int32 cols, int32 rows, int32 planes = 1 )
      : m_cols(cols), m_rows(rows), m_planes(planes) {}

    inline int32 cols() const { return m_cols; }
    inline int32 rows() const { return m_rows; }
    inline int32 planes() const { return m_planes; }

    inline pixel_accessor origin() const { return pixel_accessor(*this); }

    inline result_type operator()( int32 i, int32 j, int32 p=0 ) const { 
      return Vector2(i,j);
    }

    typedef PixelIndexView prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i /*bbox*/ ) const { return *this; }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i bbox ) const {
      vw::rasterize( prerasterize(bbox), dest, bbox );
    }
  };

  template <>
  struct IsMultiplyAccessible<PixelIndexView> : public true_type {};

  inline PixelIndexView pixel_index_view( int32 cols, int32 rows, int32 planes=1 ) {
    return PixelIndexView( cols, rows, planes );
  }

  template <class ImageT>
  inline PixelIndexView pixel_index_view( ImageViewBase<ImageT> const& image ) {
    return PixelIndexView( image );
  }


  /// PixelIndex3View
  ///
  /// This is a procedurally generated utility view whose pixels are
  /// Vector3's that contain the column, row, and plane of the image
  /// at that index (i,j,p).
  class PixelIndex3View : public ImageViewBase<PixelIndex3View> {
    int32 m_rows, m_cols, m_planes;
  public:
    typedef Vector3 pixel_type;
    typedef Vector3 result_type;
    typedef ProceduralPixelAccessor<PixelIndex3View> pixel_accessor;

    /// Initialize from another view
    template <class ImageT>
    PixelIndex3View( ImageViewBase<ImageT> const& view ) : 
      m_rows(view.impl().rows()), m_cols(view.impl().cols()), m_planes(view.impl().planes()) {}
    
    /// Initialize explicitly
    PixelIndex3View( int rows, int cols, int planes = 1 ) : 
      m_rows(rows), m_cols(cols), m_planes(planes) {}

    inline int32 cols() const { return m_cols; }
    inline int32 rows() const { return m_rows; }
    inline int32 planes() const { return m_planes; }
    
    inline pixel_accessor origin() const { return pixel_accessor(*this); }

    inline result_type operator()( int32 i, int32 j, int32 p=0 ) const { 
      return Vector3(i,j,p);
    } 

    typedef PixelIndex3View prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i const& bbox ) const { return *this; }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const {
      vw::rasterize( prerasterize(bbox), dest, bbox );
    }
  };

  template<>
  struct IsMultiplyAccessible<PixelIndex3View> : public true_type {};

  inline PixelIndex3View pixel_index3_view( int32 cols, int32 rows, int32 planes=1 ) {
    return PixelIndex3View( cols, rows, planes );
  }

  template <class ImageT>
  inline PixelIndex3View pixel_index3_view( ImageViewBase<ImageT> const& image ) {
    return PixelIndex3View( image );
  }

} // namespace vw

#endif // __VW_IMAGE_IMAGEVIEW_H__
