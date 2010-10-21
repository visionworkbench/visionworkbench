// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
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
  struct NoEdgeExtension : EdgeExtensionBase {
    template <class ViewT>
    inline typename ViewT::pixel_type operator()( const ViewT &view, int32 i, int32 j, int32 p ) const {
      return view(i,j,p);
    }
    template <class ViewT>
    inline BBox2i source_bbox( ViewT const& /*view*/, BBox2i const& bbox ) const {
      return bbox;
    }
  };

  /// An edge extention type that extends the image with zeroes in
  /// all directions.
  struct ZeroEdgeExtension : EdgeExtensionBase {
    template <class ViewT>
    inline typename ViewT::pixel_type operator()( const ViewT &view, int32 i, int32 j, int32 p ) const {
      if( i>=0 && j>=0 && i<view.cols() && j<view.rows() )
        return view(i,j,p);
      else
        return typename ViewT::pixel_type();
    }
    template <class ViewT>
    inline BBox2i source_bbox( ViewT const& view, BBox2i const& bbox ) const {
      BBox2i result = bbox;
      result.crop( BBox2i( 0, 0, view.cols(), view.rows() ) );
      return result;
    }
  };

  /// An edge extention type that extends the image with a
  /// user-supplied value in all directions.
  template <class PixelT>
  struct ValueEdgeExtension : EdgeExtensionBase {
    PixelT m_pix;
    ValueEdgeExtension(PixelT pix) : m_pix(pix) {}

    template <class ViewT>
    inline PixelT operator()( const ViewT &view, int32 i, int32 j, int32 p ) const {
      if( i>=0 && j>=0 && i<view.cols() && j<view.rows() )
        return view(i,j,p);
      else
        return m_pix;
    }

    template <class ViewT>
    inline BBox2i source_bbox( ViewT const& view, BBox2i const& bbox ) const {
      BBox2i result = bbox;
      result.crop( BBox2i( 0, 0, view.cols(), view.rows() ) );
      return result;
    }
  };

  /// An edge extension type that extends the image using constant
  /// functions in all directions.  In other words, it returns the
  /// nearest valid pixel.
  struct ConstantEdgeExtension : EdgeExtensionBase {
    template <class ViewT>
    inline typename ViewT::pixel_type operator()( const ViewT &view, int32 i, int32 j, int32 p ) const {
      return view((i<0) ? 0 : (i>=view.cols()) ? (view.cols()-1) : i,
                  (j<0) ? 0 : (j>=view.rows()) ? (view.rows()-1) : j, p);
    }
    template <class ViewT>
    inline BBox2i source_bbox( ViewT const& view, BBox2i const& bbox ) const {
      BBox2i result = bbox;
      if( bbox.min().x() < 0 ) result.min().x() = 0;
      else if( bbox.min().x() >= view.cols() ) result.min().x() = view.cols()-1;
      if( bbox.min().y() < 0 ) result.min().y() = 0;
      else if( bbox.min().y() >= view.rows() ) result.min().y() = view.rows()-1;
      if( bbox.max().x() > view.cols() ) result.max().x() = view.cols();
      else if( bbox.max().x() <= 0 ) result.max().x() = 1;
      if( bbox.max().y() > view.rows() ) result.max().y() = view.rows();
      else if( bbox.max().y() <= 0 ) result.max().y() = 1;
      return result;
    }
  };

  /// A periodic edge extension type.
  struct PeriodicEdgeExtension : EdgeExtensionBase {
    template <class ViewT>
    inline typename ViewT::pixel_type operator()( const ViewT &view, int32 i, int32 j, int32 p ) const {
      int32 d_i=i, d_j=j;
      d_i %= int(view.cols());
      if( d_i < 0 ) d_i += view.cols();
      d_j %= int(view.rows());
      if( d_j < 0 ) d_j += view.rows();
      return view(d_i,d_j,p);
    }
    template <class ViewT>
    inline BBox2i source_bbox( ViewT const& view, BBox2i const& bbox ) const {
      BBox2i result;
      if( bbox.width() >= view.cols() ) {
        result.min().x() = 0;
        result.max().x() = view.cols();
      }
      else {
        result.min().x() = mod(bbox.min().x(), view.cols());
        result.max().x() = mod(bbox.max().x()-1, view.cols())+1;
        if( result.min().x() >= result.max().x() ) {
          result.min().x() = 0;
          result.max().x() = view.cols();
        }
      }
      if( bbox.height() >= view.rows() ) {
        result.min().y() = 0;
        result.max().y() = view.rows();
      }
      else {
        result.min().y() = mod(bbox.min().y(), view.rows());
        result.max().y() = mod(bbox.max().y()-1, view.rows())+1;
        if( result.min().y() >= result.max().y() ) {
          result.min().y() = 0;
          result.max().y() = view.rows();
        }
      }
      return result;
    }
  };

  /// A cylindrical edge extension type: periodic in the x axis, constant in the y axis.
  struct CylindricalEdgeExtension : EdgeExtensionBase {
    template <class ViewT>
    inline typename ViewT::pixel_type operator()( const ViewT &view, int32 i, int32 j, int32 p ) const {
      int32 d_i=i;
      d_i %= int(view.cols());
      if( d_i < 0 ) d_i += view.cols();
      int32 d_j = (j<0) ? 0 : (j>=view.rows()) ? (view.rows()-1) : j;
      return view(d_i,d_j,p);
    }
    template <class ViewT>
    inline BBox2i source_bbox( ViewT const& view, BBox2i const& bbox ) const {
      BBox2i result;
      if( bbox.width() >= view.cols() ) {
        result.min().x() = 0;
        result.max().x() = view.cols();
      }
      else {
        result.min().x() = mod(bbox.min().x(), view.cols());
        result.max().x() = mod(bbox.max().x()-1, view.cols())+1;
        if( result.min().x() >= result.max().x() ) {
          result.min().x() = 0;
          result.max().x() = view.cols();
        }
      }
      if( bbox.min().y() < 0 ) result.min().y() = 0;
      else if( bbox.min().y() >= view.rows() ) result.min().y() = view.rows()-1;
      else result.min().y() = bbox.min().y();
      if( bbox.max().y() > view.rows() ) result.max().y() = view.rows();
      else if( bbox.max().y() <= 0 ) result.max().y() = 1;
      else result.max().y() = bbox.max().y();
      return result;
    }
  };

  /// A reflection edge extension type.
  struct ReflectEdgeExtension : EdgeExtensionBase {
    template <class ViewT>
    inline typename ViewT::pixel_type operator()( const ViewT &view, int32 i, int32 j, int32 p ) const {
      int32 d_i=i, d_j=j;
      if( d_i < 0 ) d_i = -d_i;
      int32 vcm1 = view.cols() - 1;
      d_i %= 2*vcm1;
      if( d_i > vcm1 ) d_i = 2*vcm1 - d_i;
      if( d_j<0 ) d_j=-d_j;
      int32 vrm1 = view.rows() - 1;
      d_j %= 2*vrm1;
      if( d_j > vrm1 ) d_j = 2*vrm1 - d_j;
      return view(d_i,d_j,p);
    }
    template <class ViewT>
    inline BBox2i source_bbox( ViewT const& view, BBox2i const& bbox ) const {
      BBox2i result;
      if( bbox.width() >= 2*view.cols()-2 ) {
        result.min().x() = 0;
        result.max().x() = view.cols();
      }
      else {
        result.min().x() = mod(bbox.min().x(), 2*view.cols()-2);
        result.max().x() = mod(bbox.max().x()-1, 2*view.cols()-2)+1;
        if( result.min().x() >= result.max().x() ) {
          result.max().x() = (std::min)((std::max)(result.max().x(),2*view.cols()-1-result.min().x()),view.cols());
          result.min().x() = 0;
        }
        else if( result.min().x() >= view.cols() ) {
          int tmp = 2*view.cols()-1-result.min().x();
          result.min().x() = 2*view.cols()-1-result.max().x();
          result.max().x() = tmp;
        }
        else if( result.max().x() > view.cols() ) {
          result.min().x() = (std::min)(result.min().x(),2*view.cols()-1-result.max().x());
          result.max().x() = view.cols();
        }
      }
      if( bbox.height() >= 2*view.rows()-2 ) {
        result.min().y() = 0;
        result.max().y() = view.rows();
      }
      else {
        result.min().y() = mod(bbox.min().y(), 2*view.rows()-2);
        result.max().y() = mod(bbox.max().y()-1, 2*view.rows()-2)+1;
        if( result.min().y() >= result.max().y() ) {
          result.max().y() = (std::min)((std::max)(result.max().y(),2*view.rows()-1-result.min().y()),view.rows());
          result.min().y() = 0;
        }
        else if( result.min().y() >= view.rows() ) {
          int tmp = 2*view.rows()-1-result.min().y();
          result.min().y() = 2*view.rows()-1-result.max().y();
          result.max().y() = tmp;
        }
        else if( result.max().y() > view.rows() ) {
          result.min().y() = (std::min)(result.min().y(),2*view.rows()-1-result.max().y());
          result.max().y() = view.rows();
        }
      }
      return result;
    }
  };

  /// A linear extrapolation edge extension type.
  struct LinearEdgeExtension : EdgeExtensionBase {
    template <class ViewT>
    inline typename ViewT::pixel_type operator()( const ViewT &view, int32 i, int32 j, int32 p ) const {
      int32 vcm1 = view.cols() - 1;
      int32 vrm1 = view.rows() - 1;
      if( i < 0 ) {
        if( j < 0 ) return view(0,0,p) - i*(view(0,0,p)-view(1,0,p)) - j*(view(0,0,p)-view(0,1,p));
        else if( j > vrm1 ) return view(0,vrm1,p) - i*(view(0,vrm1,p)-view(1,vrm1,p)) + (j-vrm1)*(view(0,vrm1,p)-view(0,vrm1-1,p));
        else return view(0,j,p) - i*(view(0,j,p)-view(1,j,p));
      }
      else if( i > vcm1 ) {
        if( j < 0 ) return view(vcm1,0,p) + (i-vcm1)*(view(vcm1,0,p)-view(vcm1-1,0,p)) - j*(view(vcm1,0,p)-view(vcm1,1,p));
        else if( j > vrm1 ) return view(vcm1,vrm1,p) + (i-vcm1)*(view(vcm1,vrm1,p)-view(vcm1-1,vrm1,p)) + (j-vrm1)*(view(vcm1,vrm1,p)-view(vcm1,vrm1-1,p));
        else return view(vcm1,j,p) + (i-vcm1)*(view(vcm1,j,p)-view(vcm1-1,j,p));
      }
      else {
        if( j < 0 ) return view(i,0,p) - j*(view(i,0,p)-view(i,1,p));
        else if( j > vrm1 ) return view(i,vrm1,p) + (j-vrm1)*(view(i,vrm1,p)-view(i,vrm1-1,p));
        else return view(i,j,p);
      }
    }
    template <class ViewT>
    inline BBox2i source_bbox( ViewT const& view, BBox2i const& bbox ) const {
      BBox2i result = bbox;
      if( bbox.min().x() < 0 ) {
        result.min().x() = 0;
        result.max().x() = (std::max)( result.max().x(), 2 );
      }
      if( bbox.max().x() > view.cols() ) {
        result.max().x() = view.cols();
        result.min().x() = (std::min)( result.min().x(), view.cols()-2 );
      }
      if( bbox.min().y() < 0 ) {
        result.min().y() = 0;
        result.max().y() = (std::max)( result.max().y(), 2 );
      }
      if( bbox.max().y() > view.rows() ) {
        result.max().y() = view.rows();
        result.min().y() = (std::min)( result.min().y(), view.rows()-2 );
      }
      return result;
    }
  };

  // *******************************************************************
  // The edge extension view
  // *******************************************************************

  /// A wrapper view supporting procedural extension of an image beyond its
  /// boundaries using a user-specified edge extension mode.
  template <class ImageT, class ExtensionT>
  class EdgeExtensionView : public ImageViewBase<EdgeExtensionView<ImageT,ExtensionT> >
  {
  private:
    ImageT m_image;
    int32 m_xoffset, m_yoffset;
    int32 m_cols, m_rows;
    ExtensionT m_extension_func;
  public:

    typedef typename ImageT::pixel_type pixel_type;
    typedef pixel_type result_type;
    typedef ProceduralPixelAccessor<EdgeExtensionView<ImageT, ExtensionT> > pixel_accessor;

    EdgeExtensionView( ImageT const& image )
      : m_image(image), m_xoffset(0), m_yoffset(0), m_cols(image.cols()), m_rows(image.rows()), m_extension_func() {}

    EdgeExtensionView( ImageT const& image, ExtensionT const& extension )
      : m_image(image), m_xoffset(0), m_yoffset(0), m_cols(image.cols()), m_rows(image.rows()), m_extension_func(extension) {}

    EdgeExtensionView( ImageT const& image, int32 xoffset, int32 yoffset, int32 cols, int32 rows )
      : m_image(image), m_xoffset(xoffset), m_yoffset(yoffset), m_cols(cols), m_rows(rows), m_extension_func() {}

    EdgeExtensionView( ImageT const& image, int32 xoffset, int32 yoffset, int32 cols, int32 rows, ExtensionT const& extension )
      : m_image(image), m_xoffset(xoffset), m_yoffset(yoffset), m_cols(cols), m_rows(rows), m_extension_func(extension) {}

    inline int32 cols() const { return m_cols; }
    inline int32 rows() const { return m_rows; }
    inline int32 planes() const { return m_image.planes(); }

    inline pixel_accessor origin() const { return pixel_accessor(*this,0,0); }
    inline result_type operator()( int32 i, int32 j, int32 p = 0 ) const { return m_extension_func(m_image,i+m_xoffset,j+m_yoffset,p); }

    ImageT const& child() const { return m_image; }
    ExtensionT const& func() const { return m_extension_func; }
    BBox2i source_bbox( BBox2i const& bbox ) const {
      return m_extension_func.source_bbox( m_image, bbox + Vector2i( m_xoffset, m_yoffset ) );
    }

    typedef EdgeExtensionView<typename ImageT::prerasterize_type,ExtensionT> prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i const& bbox ) const {
      BBox2i src_bbox = source_bbox( bbox );
      // Make degenerate bboxes sane
      if( src_bbox.empty() ) src_bbox = BBox2i(0,0,0,0);
      vw_out(VerboseDebugMessage, "image") << "EdgeExtensionView: prerasterizing child view with bbox " << src_bbox << ".\n";
      return prerasterize_type(m_image.prerasterize(src_bbox), m_xoffset, m_yoffset, m_cols, m_rows, m_extension_func );
    }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i bbox ) const { vw::rasterize( prerasterize(bbox), dest, bbox ); }
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

#endif // __VW_IMAGE_EDGEEXTENSION_H__
