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


/// \file EdgeExtension.h
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
///  - \ref vw::ZeroEdgeExtension extends an image with zeros
///  - \ref vw::ConstantEdgeExtension extends an image with constant (i.e. nearest neighbor) values
///  - \ref vw::ReflectEdgeExtension extends an image by reflecting across its edges
///  - \ref vw::PeriodicEdgeExtension extends an image by repeating it periodically
///  - \ref vw::NoEdgeExtension is a special class that can be used to represent no edge extension
///
#ifndef __VW_IMAGE_EDGEEXTENSION_H__
#define __VW_IMAGE_EDGEEXTENSION_H__

#include <boost/type_traits.hpp>

#include <vw/Image/ImageViewBase.h>
#include <vw/Image/PixelAccessors.h>
#include <vw/Image/SparseImageCheck.h>
#include <vw/Core/Log.h>

namespace vw {

  // *******************************************************************
  // The edge extension types
  // *******************************************************************

  // Here we define the edge extension modes as polymorphic functors.
  // You can extend the list of supported edge extension modes by
  // creating a new functor type derived from EdgeExtensionBase and
  // implementing the function call operator, as shown in the
  // examples below.

  /// A base class for the edge extension types that provides the
  /// common return type deduction logic in case users want to use
  /// these types in a more general manner.
  struct EdgeExtensionBase {
    template <class ArgsT> struct result {};
    template <class FuncT, class ViewT, class IT, class JT, class PT>
    struct result<FuncT(ViewT,IT,JT,PT)> {
      typedef typename boost::remove_reference<ViewT>::type::pixel_type type;
    };
  };

  /// A special type providing no edge extension.  This is mainly
  /// used as a signal to certain types to adopt fundamentally
  /// different behavior.
  struct NoEdgeExtension;

  /// An edge extention type that extends the image with zeroes in all directions.
  struct ZeroEdgeExtension;

  /// An edge extention type that extends the image with a
  /// user-supplied value in all directions.
  template <class PixelT>
  struct ValueEdgeExtension;

  /// An edge extension type that extends the image using constant
  /// functions in all directions.  In other words, it returns the
  /// nearest valid pixel.
  struct ConstantEdgeExtension;

  /// A periodic edge extension type.
  struct PeriodicEdgeExtension;

  /// A cylindrical edge extension type: periodic in the x axis, constant in the y axis.
  struct CylindricalEdgeExtension;

  /// A reflection edge extension type.
  struct ReflectEdgeExtension;

  /// A linear extrapolation edge extension type.
  struct LinearEdgeExtension;

  // *******************************************************************
  // The edge extension view
  // *******************************************************************

  /// A wrapper view supporting procedural extension of an image beyond its
  /// boundaries using a user-specified edge extension mode.
  template <class ImageT, class ExtensionT>
  class EdgeExtensionView : public ImageViewBase<EdgeExtensionView<ImageT,ExtensionT> >
  {
  private:
    ImageT     m_image;
    int32      m_xoffset, m_yoffset;
    int32      m_cols, m_rows;
    ExtensionT m_extension_func;
  public:

    typedef typename ImageT::pixel_type pixel_type;
    typedef pixel_type result_type;
    typedef ProceduralPixelAccessor<EdgeExtensionView<ImageT, ExtensionT> > pixel_accessor;

    EdgeExtensionView( ImageT const& image )
      : m_image(image), m_xoffset(0), m_yoffset(0), 
        m_cols(image.cols()), m_rows(image.rows()), m_extension_func() {}

    EdgeExtensionView( ImageT const& image, ExtensionT const& extension )
      : m_image(image), m_xoffset(0), m_yoffset(0), m_cols(image.cols()), 
        m_rows(image.rows()), m_extension_func(extension) {}

    EdgeExtensionView( ImageT const& image, int32 xoffset, int32 yoffset, int32 cols, int32 rows )
      : m_image(image), m_xoffset(xoffset), m_yoffset(yoffset), 
        m_cols(cols), m_rows(rows), m_extension_func() {}

    EdgeExtensionView( ImageT const& image, int32 xoffset, int32 yoffset, int32 cols, int32 rows, ExtensionT const& extension )
      : m_image(image), m_xoffset(xoffset), m_yoffset(yoffset), 
        m_cols(cols), m_rows(rows), m_extension_func(extension) {}

    inline int32 cols  () const { return m_cols;           }
    inline int32 rows  () const { return m_rows;           }
    inline int32 planes() const { return m_image.planes(); }

    inline pixel_accessor origin() const { return pixel_accessor(*this,0,0); }
    inline result_type operator()( int32 i, int32 j, int32 p = 0 ) const { 
      return m_extension_func(m_image,i+m_xoffset,j+m_yoffset,p); 
    }

    ImageT     const& child() const { return m_image;          }
    ExtensionT const& func () const { return m_extension_func; }
    BBox2i source_bbox( BBox2i const& bbox ) const {
      return m_extension_func.source_bbox( m_image, bbox + Vector2i( m_xoffset, m_yoffset ) );
    }

    typedef EdgeExtensionView<typename ImageT::prerasterize_type,ExtensionT> prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i const& bbox ) const {
      BBox2i src_bbox = source_bbox( bbox );
      // Make degenerate bboxes sane
      if( src_bbox.empty() ) src_bbox = BBox2i(0,0,0,0);
      VW_OUT(VerboseDebugMessage, "image") << "EdgeExtensionView: prerasterizing child view with bbox " << src_bbox << ".\n";
      return prerasterize_type(m_image.prerasterize(src_bbox), m_xoffset, m_yoffset, m_cols, m_rows, m_extension_func );
    }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const { vw::rasterize( prerasterize(bbox), dest, bbox ); }
  };

  template <class ImageT, class ExtensionT>
  class SparseImageCheck<EdgeExtensionView<ImageT, ExtensionT> > {
    EdgeExtensionView<ImageT, ExtensionT> const& m_view;
  public:
    SparseImageCheck(EdgeExtensionView<ImageT, ExtensionT> const& view)
      : m_view(view) {}
    bool operator()( BBox2i const& bbox ) const {
      return SparseImageCheck<ImageT>(m_view.child())( m_view.source_bbox( bbox ) );
    }
  };


  // *******************************************************************
  // General-purpose edge extension functions
  // *******************************************************************

  /// Procedurally edge-extends an image view and alters its bounding box.
  /// This function wraps an image view in a vw::EdgeExtensionView, extending
  /// its definition beyond the edges of the original image according to
  /// the user-specified edge extension mode.  It also alters the bounding
  /// box, i.e. the position of the origin and the values returned by
  /// <B>cols()</B> and <B>rows()</B>.  The new bounding box is expressed
  /// in terms of an offset and new dimensions.  The sign of the offset is
  /// consistent with \ref vw::crop, so if you wish to expand the image
  /// you must specify a negative offset!
  template <class ImageT, class ExtensionT>
  EdgeExtensionView<ImageT,ExtensionT> edge_extend( ImageViewBase<ImageT> const& v, int32 x_offset, int32 y_offset, int32 cols, int32 rows, ExtensionT const& extension ) {
    return EdgeExtensionView<ImageT,ExtensionT>( v.impl(), x_offset, y_offset, cols, rows, extension );
  }

  /// This is an overloaded function provided for convenience; see vw::edge_extend.
  /// It takes a bounding box object (BBox2i) as an argument.
  template <class ImageT, class ExtensionT>
  EdgeExtensionView<ImageT,ExtensionT> edge_extend( ImageViewBase<ImageT> const& v, BBox2i const& bbox, ExtensionT const& extension ) {
    return EdgeExtensionView<ImageT,ExtensionT>( v.impl(), bbox.min().x(), bbox.min().y(), bbox.width(), bbox.height(), extension );
  }

  /// This is an overloaded function provided for convenience; see vw::edge_extend.
  /// It uses the default vw::ConstantEdgeExtension mode.
  template <class ImageT>
  EdgeExtensionView<ImageT,ConstantEdgeExtension> edge_extend( ImageViewBase<ImageT> const& v, int32 x_offset, int32 y_offset, int32 cols, int32 rows ) {
    return EdgeExtensionView<ImageT,ConstantEdgeExtension>( v.impl(), x_offset, y_offset, cols, rows );
  }

  /// This is an overloaded function provided for convenience; see vw::edge_extend.
  /// It takes a bounding box object (BBox2i) as an argument and uses the default
  /// vw::ConstantEdgeExtension mode.
  template <class ImageT>
  EdgeExtensionView<ImageT,ConstantEdgeExtension> edge_extend( ImageViewBase<ImageT> const& v, BBox2i const& bbox ) {
    return EdgeExtensionView<ImageT,ConstantEdgeExtension>( v.impl(), bbox.min().x(), bbox.min().y(), bbox.width(), bbox.height() );
  }

  /// This is an overloaded function provided for convenience; see vw::edge_extend.
  /// It does not alter the bounding box of the image, but still allows access
  /// beyond the bounding box.
  template <class ImageT, class ExtensionT>
  EdgeExtensionView<ImageT,ExtensionT> edge_extend( ImageViewBase<ImageT> const& v, ExtensionT const& extension ) {
    return EdgeExtensionView<ImageT,ExtensionT>( v.impl(), extension );
  }

  /// This is an overloaded function provided for convenience; see vw::edge_extend.
  /// It does not alter the bounding box of the image, but still allows access
  /// beyond the bounding box using the default vw::ConstantEdgeExtension mode.
  template <class ImageT>
  EdgeExtensionView<ImageT,ConstantEdgeExtension> edge_extend( ImageViewBase<ImageT> const& v ) {
    return EdgeExtensionView<ImageT,ConstantEdgeExtension>( v.impl() );
  }

} // namespace vw

#include "EdgeExtension.tcc"

#endif // __VW_IMAGE_EDGEEXTENSION_H__
