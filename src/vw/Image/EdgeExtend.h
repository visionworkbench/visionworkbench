// __BEGIN_LICENSE__
//
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
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

/// \file EdgeExtend.h
///
/// Supports extension of images beyond their edges.
///
/// This file provides wrapper mechanisms to procedurally extend 
/// images beyond their edges, for use in (for example) filtering 
/// or compositing operations.  The two main facilities provided 
/// by this file are a function (\ref vw::edge_extend) that you 
/// can use to manually wrap an image and a set of classes that 
/// specify the various edge extension modes.  At the moment 
/// there are five such classes:
///
///  - \ref vw::ZeroEdgeExtend extends an image with zeros
///  - \ref vw::ConstantEdgeExtend extends an image with constant (i.e. nearest neighbor) values
///  - \ref vw::ReflectEdgeExtend extends an image by reflecting across its edges
///  - \ref vw::PeriodicEdgeExtend extends an image by repeating it periodically
///  - \ref vw::NoEdgeExtend is a special class that can be used to represent no edge extension
///
#ifndef __VW_IMAGE_EDGEEXTEND_H__
#define __VW_IMAGE_EDGEEXTEND_H__

#include <boost/type_traits.hpp>

#include <vw/Image/ImageViewBase.h>
#include <vw/Image/PixelAccessors.h>
#include <vw/Math/BBox.h>

namespace vw {
 
  // *******************************************************************
  // The edge extension types
  // *******************************************************************

  // Here we define the edge extension modes as polymorphic functors.
  // You can extend the list of supported edge extension modes by
  // creating a new functor type derived from EdgeExtendBase and 
  // implementing the function call operator, as shown in the 
  // examples below.

  /// A base class for the edge extension types that provides the 
  /// common return type deduction logic in case users want to use 
  /// these types in a more general manner.
  struct EdgeExtendBase {
    template <class ArgsT> struct result {};
    template <class FuncT, class ViewT, class IT, class JT, class PT>
    struct result<FuncT(ViewT,IT,JT,PT)> {
      typedef typename boost::remove_reference<ViewT>::type::pixel_type type;
    };
  };

  /// A special type providing no edge extension.  This is mainly 
  /// used as a signal to certain types to adopt fundamentally 
  /// different behavior.
  struct NoEdgeExtend : EdgeExtendBase {
    template <class ViewT>
    inline typename ViewT::pixel_type operator()( const ViewT &view, int i, int j, int p ) const { 
      return view(i,j,p);
    }
  };

  /// An edge extention type that extends the image with zeroes in 
  /// all directions.
  struct ZeroEdgeExtend : EdgeExtendBase {
    template <class ViewT>
    inline typename ViewT::pixel_type operator()( const ViewT &view, int i, int j, int p ) const { 
      if( i>=0 && j>=0 && i<view.cols() && j<view.rows() )
        return view(i,j,p);
      else
        return typename ViewT::pixel_type();
    }
  };

  /// An edge extension type that extends the image using constant 
  /// functions in all directions.  In other words, it returns the 
  /// nearest valid pixel.
  struct ConstantEdgeExtend : EdgeExtendBase {
    template <class ViewT>
    inline typename ViewT::pixel_type operator()( const ViewT &view, int i, int j, int p ) const { 
      return view((i<0) ? 0 : (i>int(view.cols()-1)) ? (view.cols()-1) : i,
                  (j<0) ? 0 : (j>=int(view.rows()-1)) ? (view.rows()-1) : j, p);
    }
  };

  /// A periodic edge extension type.
  struct PeriodicEdgeExtend : EdgeExtendBase {
    template <class ViewT>
    inline typename ViewT::pixel_type operator()( const ViewT &view, int i, int j, int p ) const { 
      int d_i=i, d_j=j;
      d_i %= int(view.cols());
      if( d_i < 0 ) d_i += view.cols();
      d_j %= int(view.rows());
      if( d_j < 0 ) d_j += view.rows();
      return view(d_i,d_j,p);
    }
  };

  /// A reflection edge extension type.
  struct ReflectEdgeExtend : EdgeExtendBase {
    template <class ViewT>
    inline typename ViewT::pixel_type operator()( const ViewT &view, int i, int j, int p ) const { 
      int d_i=i, d_j=j;
      if( d_i < 0 ) d_i = -d_i;
      int vcm1 = view.cols() - 1;
      d_i %= 2*vcm1;
      if( d_i > vcm1 ) d_i = 2*vcm1 - d_i;
      if( d_j<0 ) d_j=-d_j;
      int vrm1 = view.rows() - 1;
      d_j %= 2*vrm1;
      if( d_j > vrm1 ) d_j = 2*vrm1 - d_j;
      return view(d_i,d_j,p);
    }
  };


  // *******************************************************************
  // The edge extension view
  // *******************************************************************

  /// A wrapper view supporting procedural extension of an image beyond its 
  /// boundaries using a user-specified edge extension mode.
  template <class ImageT, class EdgeExtendT>
  class EdgeExtendView : public ImageViewBase<EdgeExtendView<ImageT,EdgeExtendT> >
  {
  private:
    ImageT m_image;
    ptrdiff_t m_xoffset, m_yoffset;
    unsigned m_cols, m_rows;
    EdgeExtendT m_extend;
  public:

    typedef typename ImageT::pixel_type pixel_type;
    typedef ProceduralPixelAccessor<EdgeExtendView<ImageT, EdgeExtendT> > pixel_accessor;

    EdgeExtendView( ImageT const& image )
      : m_image(image), m_xoffset(0), m_yoffset(0), m_cols(image.cols()), m_rows(image.rows()), m_extend() {}
    
    EdgeExtendView( ImageT const& image, EdgeExtendT const& extend )
      : m_image(image), m_xoffset(0), m_yoffset(0), m_cols(image.cols()), m_rows(image.rows()), m_extend(extend) {}
    
    EdgeExtendView( ImageT const& image, ptrdiff_t xoffset, ptrdiff_t yoffset, unsigned cols, unsigned rows )
      : m_image(image), m_xoffset(xoffset), m_yoffset(yoffset), m_cols(cols), m_rows(rows), m_extend() {}
    
    EdgeExtendView( ImageT const& image, ptrdiff_t xoffset, ptrdiff_t yoffset, unsigned cols, unsigned rows, EdgeExtendT const& extend )
      : m_image(image), m_xoffset(xoffset), m_yoffset(yoffset), m_cols(cols), m_rows(rows), m_extend(extend) {}

    inline unsigned cols() const { return m_cols; }
    inline unsigned rows() const { return m_rows; }
    inline unsigned planes() const { return m_image.planes(); }

    inline pixel_accessor origin() const { return pixel_accessor(*this,m_xoffset,m_yoffset); }

    inline pixel_type operator()( int i, int j, int p = 0 ) const { return m_extend(m_image,i,j,p); }

    typedef EdgeExtendView<typename ImageT::prerasterize_type,EdgeExtendT> prerasterize_type;
    inline prerasterize_type prerasterize() const { return prerasterize_type( m_image.prerasterize(), m_xoffset, m_yoffset, m_cols, m_rows ); }
    template <class DestT> inline void rasterize( DestT const& dest ) const { vw::rasterize( prerasterize(), dest ); }
  };


  // *******************************************************************
  // General-purpose edge extension functions
  // *******************************************************************

  /// Procedurally edge-extends an image view and alters its bounding box.
  /// This function wraps an image view in a vw::EdgeExtendView, extending 
  /// its definition beyond the edges of the original image according to 
  /// the user-specified edge extension mode.  It also alters the bounding 
  /// box, i.e. the position of the origin and the values returned by 
  /// <B>cols()</B> and <B>rows()</B>.  The new bounding box is expressed 
  /// in terms of an offset and new dimensions.  The sign of the offset is
  /// consistent with \ref vw::crop, so if you wish to expand the image 
  /// you must specify a negative offset!
  template <class ImageT, class EdgeExtendT>
  EdgeExtendView<ImageT,EdgeExtendT> edge_extend( ImageViewBase<ImageT> const& v, ptrdiff_t x_offset, ptrdiff_t y_offset, unsigned cols, unsigned rows, EdgeExtendT const& extend ) {
    return EdgeExtendView<ImageT,EdgeExtendT>( v.impl(), x_offset, y_offset, cols, rows, extend );
  }

  /// This is an overloaded function provided for convenience; see vw::edge_extend.
  /// It takes a bounding box object (BBox<int,2>) as an argument.
  template <class ImageT, class EdgeExtendT>
  EdgeExtendView<ImageT,EdgeExtendT> edge_extend( ImageViewBase<ImageT> const& v, BBox<int,2> const& bbox, EdgeExtendT const& extend ) {
    return EdgeExtendView<ImageT,EdgeExtendT>( v.impl(), bbox.min().x(), bbox.min().y(), bbox.width(), bbox.height(), extend );
  }

  /// This is an overloaded function provided for convenience; see vw::edge_extend.
  /// It uses the default vw::ConstantEdgeExtend mode.
  template <class ImageT>
  EdgeExtendView<ImageT,ConstantEdgeExtend> edge_extend( ImageViewBase<ImageT> const& v, ptrdiff_t x_offset, ptrdiff_t y_offset, unsigned cols, unsigned rows ) {
    return EdgeExtendView<ImageT,ConstantEdgeExtend>( v.impl(), x_offset, y_offset, cols, rows );
  }
  
  /// This is an overloaded function provided for convenience; see vw::edge_extend.
  /// It takes a bounding box object (BBox<int,2>) as an argument and uses the default 
  /// vw::ConstantEdgeExtend mode.
  template <class ImageT>
  EdgeExtendView<ImageT,ConstantEdgeExtend> edge_extend( ImageViewBase<ImageT> const& v, BBox<int,2> const& bbox ) {
    return EdgeExtendView<ImageT,ConstantEdgeExtend>( v.impl(), bbox.min().x(), bbox.min().y(), bbox.width(), bbox.height() );
  }
  
  /// This is an overloaded function provided for convenience; see vw::EdgeExtend.
  /// It does not alter the bounding box of the image, but still allows access 
  /// beyond the bounding box.
  template <class ImageT, class EdgeExtendT>
  EdgeExtendView<ImageT,EdgeExtendT> edge_extend( ImageViewBase<ImageT> const& v, EdgeExtendT const& extend ) {
    return EdgeExtendView<ImageT,EdgeExtendT>( v.impl(), extend );
  }

  /// This is an overloaded function provided for convenience; see vw::edge_extend.
  /// It does not alter the bounding box of the image, but still allows access 
  /// beyond the bounding box using the default vw::ConstantEdgeExtend mode.
  template <class ImageT>
  EdgeExtendView<ImageT,ConstantEdgeExtend> edge_extend( ImageViewBase<ImageT> const& v ) {
    return EdgeExtendView<ImageT,ConstantEdgeExtend>( v.impl() );
  }
  
} // namespace vw

#endif // __VW_IMAGE_EDGEEXTEND_H__
