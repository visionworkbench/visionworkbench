// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file Manipulation.h
///
/// Simple image manipulation functions, such as flipping and
/// cropping.  All of the functions in this file except copy() return
/// <i>shallow</i> views of the ImageView.  That is, they do not copy
/// the data underneath but instead they refer to the same data,
/// indexing and accessing it in a different way.
///
/// The first collection of functions in this file perform basic
/// transformations to the domain of the image, such as transposition,
/// rotation by 90 degree increments, flipping, and cropping.
///
/// This file also provides views and functions that take simple
/// "slices" of images, returning a new view composed from individual
/// channels or planes of the source image. These include:
///
/// - select_col() : takes a single-column slice of an image
/// - select_row() : takes a single-row slice of an image
/// - select_plane() : takes a single-plane slice of an image
/// - select_channel() : takes a single-channel slice of an image
/// - channels_to_planes() : reinterprets a multi-channel image as a multi-plane image
/// - planes_to_channels() : reinterprets a multi-plane image as a multi-channel image
/// - pixel_cast() : casts the pixels of an image to a new pixel type
/// - pixel_cast_rescale() : same as above, but rescales pixel values appropriately
/// - channel_cast() : casts the channels of an image while retaining the pixel format
/// - channel_cast_rescale() : same as above, but rescales pixel values appropriately
///
#ifndef __VW_IMAGE_MANIPULATION_H__
#define __VW_IMAGE_MANIPULATION_H__

#include <boost/mpl/logical.hpp>

#include <vw/Image/ImageView.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Image/PerPixelViews.h>

namespace vw {

  // *******************************************************************
  // copy()
  // *******************************************************************

  // Class definition
  template <class ImageT>
  class CopyView : public ImageViewBase<CopyView<ImageT> >
  {
  private:
    ImageView<typename ImageT::pixel_type> m_child;
  public:
    typedef typename ImageView<typename ImageT::pixel_type>::pixel_type pixel_type;
    typedef pixel_type const& result_type;
    typedef typename ImageView<typename ImageT::pixel_type>::pixel_accessor pixel_accessor;

    CopyView( ImageT const& image ) : m_child(image.cols(),image.rows(),image.planes()) {
      image.rasterize( m_child, BBox2i(0,0,image.cols(),image.rows()) );
    }

    inline int32 cols() const { return m_child.cols(); }
    inline int32 rows() const { return m_child.rows(); }
    inline int32 planes() const { return m_child.planes(); }

    inline pixel_accessor origin() const { return m_child.origin(); }
    inline result_type operator()( int32 i, int32 j, int32 p=0 ) const { return m_child(i,j,p); }

    typedef CopyView prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i const& /*bbox*/ ) const { return *this; }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const { m_child.rasterize( dest, bbox ); }
  };

  template <class ImageT>
  struct IsMultiplyAccessible<CopyView<ImageT> > : public true_type {};

  /// Make a (deep) copy of an image.
  template <class ImageT>
  CopyView<ImageT> copy( ImageViewBase<ImageT> const& v ) {
    return CopyView<ImageT>( v.impl() );
  }


  // *******************************************************************
  // Transpose
  // *******************************************************************

  // Specialized pixel accessor
  template <class ChildT>
  class TransposePixelAccessor
  {
  private:
    ChildT m_child;
  public:
    typedef typename ChildT::pixel_type pixel_type;
    typedef typename ChildT::result_type result_type;
    TransposePixelAccessor( ChildT const& acc ) : m_child(acc) {}

    inline TransposePixelAccessor& next_col() { m_child.next_row(); return *this; }
    inline TransposePixelAccessor& prev_col() { m_child.prev_row(); return *this; }
    inline TransposePixelAccessor& next_row() { m_child.next_col(); return *this; }
    inline TransposePixelAccessor& prev_row() { m_child.prev_col(); return *this; }
    inline TransposePixelAccessor& next_plane() { m_child.next_plane(); return *this; }
    inline TransposePixelAccessor& prev_plane() { m_child.prev_plane(); return *this; }
    inline TransposePixelAccessor& advance( ptrdiff_t di, ptrdiff_t dj, ptrdiff_t dp=0 ) { m_child.advance(dj,di,dp); return *this; }

    inline result_type operator*() const { return *m_child; }
  };

  // Class definition
  template <class ImageT>
  class TransposeView : public ImageViewBase<TransposeView<ImageT> >
  {
  private:
    ImageT m_child;
  public:

    typedef typename ImageT::pixel_type pixel_type;
    typedef typename ImageT::result_type result_type;
    typedef TransposePixelAccessor<typename ImageT::pixel_accessor> pixel_accessor;

    TransposeView( ImageT const& image ) : m_child(image) {}

    inline int32 cols() const { return m_child.rows(); }
    inline int32 rows() const { return m_child.cols(); }
    inline int32 planes() const { return m_child.planes(); }

    inline pixel_accessor origin() const { return m_child.origin(); }

    inline result_type operator()( int32 i, int32 j, int32 p=0 ) const { return m_child(j,i,p); }

    template <class ViewT>
    TransposeView const& operator=( ImageViewBase<ViewT> const& view ) const {
      view.impl().rasterize( *this, BBox2i(0,0,view.impl().cols(),view.impl().rows()) );
      return *this;
    }

    template <class ViewT>
    TransposeView& operator=( ImageViewBase<ViewT> const& view ) {
      *const_cast<const TransposeView*>(this) = view.impl();
      return *this;
    }

    ImageT const& child() const {
      return m_child;
    }

    /// \cond INTERNAL
    typedef TransposeView<typename ImageT::prerasterize_type> prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i const& bbox ) const {
      BBox2i child_bbox( bbox.min().y(), bbox.min().x(), bbox.height(), bbox.width() );
      return prerasterize_type( m_child.prerasterize(child_bbox) );
    }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const { vw::rasterize( prerasterize(bbox), dest, bbox ); }
    /// \endcond
  };

  /// \cond INTERNAL
  template <class ImageT>
  struct IsMultiplyAccessible<TransposeView<ImageT> > : public IsMultiplyAccessible<ImageT> {};
  /// \endcond

  /// Transpose an image.
  template <class ImageT>
  TransposeView<ImageT> transpose( ImageViewBase<ImageT> const& v ) {
    return TransposeView<ImageT>( v.impl() );
  }


  // *******************************************************************
  // Rotate180
  // *******************************************************************

  // Specialized pixel accessor
  template <class ChildT>
  class Rotate180PixelAccessor
  {
  private:
    ChildT m_child;
  public:
    typedef typename ChildT::pixel_type pixel_type;
    typedef typename ChildT::result_type result_type;
    Rotate180PixelAccessor( ChildT const& image ) : m_child(image) {}

    inline Rotate180PixelAccessor& next_col() { m_child.prev_col(); return *this; }
    inline Rotate180PixelAccessor& prev_col() { m_child.next_col(); return *this; }
    inline Rotate180PixelAccessor& next_row() { m_child.prev_row(); return *this; }
    inline Rotate180PixelAccessor& prev_row() { m_child.next_row(); return *this; }
    inline Rotate180PixelAccessor& next_plane() { m_child.next_plane(); return *this; }
    inline Rotate180PixelAccessor& prev_plane() { m_child.prev_plane(); return *this; }
    inline Rotate180PixelAccessor& advance( ptrdiff_t di, ptrdiff_t dj, ptrdiff_t dp=0 ) { m_child.advance(-di,-dj,dp); return *this; }

    inline result_type operator*() const { return *m_child; }
  };

  // Image View Class
  template <class ImageT>
  class Rotate180View : public ImageViewBase<Rotate180View<ImageT> >
  {
  private:
    ImageT m_child;
  public:

    typedef typename ImageT::pixel_type pixel_type;
    typedef typename ImageT::result_type result_type;
    typedef Rotate180PixelAccessor<typename ImageT::pixel_accessor> pixel_accessor;

    Rotate180View( ImageT const& image ) : m_child(image) {}

    inline int32 cols() const { return m_child.cols(); }
    inline int32 rows() const { return m_child.rows(); }
    inline int32 planes() const { return m_child.planes(); }

    inline pixel_accessor origin() const { return m_child.origin().advance(cols()-1,rows()-1); }

    inline result_type operator()( int32 i, int32 j, int32 p=0 ) const { return m_child(cols()-1-i,rows()-1-j,p); }

    template <class ViewT>
    Rotate180View const& operator=( ImageViewBase<ViewT> const& view ) const {
      view.impl().rasterize( *this, BBox2i(0,0,view.impl().cols(),view.impl().rows()) );
      return *this;
    }

    template <class ViewT>
    Rotate180View& operator=( ImageViewBase<ViewT> const& view ) {
      *const_cast<const Rotate180View*>(this) = view.impl();
      return *this;
    }

    ImageT const& child() const {
      return m_child;
    }

    /// \cond INTERNAL
    typedef Rotate180View<typename ImageT::prerasterize_type> prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i const& bbox ) const {
      BBox2i child_bbox( cols()-bbox.max().x(), rows()-bbox.max().y(), bbox.width(), bbox.height() );
      return prerasterize_type( m_child.prerasterize(child_bbox) );
    }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const { vw::rasterize( prerasterize(bbox), dest, bbox ); }
    /// \endcond
  };

  /// \cond INTERNAL
  template <class ImageT>
  struct IsMultiplyAccessible<Rotate180View<ImageT> > : public IsMultiplyAccessible<ImageT> {};
  /// \endcond

  /// Rotate an image 180 degrees.
  template <class ImageT>
  Rotate180View<ImageT> rotate_180( ImageViewBase<ImageT> const& v ) {
    return Rotate180View<ImageT>( v.impl() );
  }


  // *******************************************************************
  // Rotate90CW
  // *******************************************************************

  // Specialized pixel accessor
  template <class ChildT>
  class Rotate90CWPixelAccessor
  {
  private:
    ChildT m_child;
  public:
    typedef typename ChildT::pixel_type pixel_type;
    typedef typename ChildT::result_type result_type;
    Rotate90CWPixelAccessor( ChildT const& acc ) : m_child(acc) {}
    inline Rotate90CWPixelAccessor& next_col() { m_child.prev_row(); return *this; }
    inline Rotate90CWPixelAccessor& prev_col() { m_child.next_row(); return *this; }
    inline Rotate90CWPixelAccessor& next_row() { m_child.next_col(); return *this; }
    inline Rotate90CWPixelAccessor& prev_row() { m_child.prev_col(); return *this; }
    inline Rotate90CWPixelAccessor& next_plane() { m_child.next_plane(); return *this; }
    inline Rotate90CWPixelAccessor& prev_plane() { m_child.prev_plane(); return *this; }
    inline Rotate90CWPixelAccessor& advance( ptrdiff_t di, ptrdiff_t dj, ptrdiff_t dp=0 ) { m_child.advance(dj,-di,dp); return *this; }

    inline result_type operator*() const { return *m_child; }
  };

  // Class definition
  template <class ImageT>
  class Rotate90CWView : public ImageViewBase<Rotate90CWView<ImageT> >
  {
  private:
    ImageT m_child;
  public:

    typedef typename ImageT::pixel_type pixel_type;
    typedef typename ImageT::result_type result_type;
    typedef Rotate90CWPixelAccessor<typename ImageT::pixel_accessor> pixel_accessor;

    Rotate90CWView( ImageT const& image ) : m_child(image) {}

    inline int32 cols() const { return m_child.rows(); }
    inline int32 rows() const { return m_child.cols(); }
    inline int32 planes() const { return m_child.planes(); }

    inline pixel_accessor origin() const { return m_child.origin().advance(0,cols()-1); }

    inline result_type operator()( int32 i, int32 j, int32 p=0 ) const { return m_child(j,cols()-1-i,p); }

    template <class ViewT>
    Rotate90CWView const& operator=( ImageViewBase<ViewT> const& view ) const {
      view.impl().rasterize( *this, BBox2i(0,0,view.impl().cols(),view.impl().rows()) );
      return *this;
    }

    template <class ViewT>
    Rotate90CWView& operator=( ImageViewBase<ViewT> const& view ) {
      *const_cast<const Rotate90CWView*>(this) = view.impl();
      return *this;
    }

    ImageT const& child() const {
      return m_child;
    }

    /// \cond INTERNAL
    typedef Rotate90CWView<typename ImageT::prerasterize_type> prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i const& bbox ) const {
      BBox2i child_bbox( bbox.min().y(), cols()-bbox.max().x(), bbox.height(), bbox.width() );
      return prerasterize_type( m_child.prerasterize(child_bbox) );
    }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const { vw::rasterize( prerasterize(bbox), dest, bbox ); }
    /// \endcond
  };

  /// \cond INTERNAL
  template <class ImageT>
  struct IsMultiplyAccessible<Rotate90CWView<ImageT> > : public IsMultiplyAccessible<ImageT> {};
  /// \endcond

  /// Rotate an image 90 degrees clockwise.
  template <class ImageT>
  Rotate90CWView<ImageT> rotate_90_cw( ImageViewBase<ImageT> const& v ) {
    return Rotate90CWView<ImageT>( v.impl() );
  }


  // *******************************************************************
  // Rotate90CCW
  // *******************************************************************

  // Specialized pixel accessor
  template <class ChildT>
  class Rotate90CCWPixelAccessor
  {
  private:
    ChildT m_child;
  public:
    typedef typename ChildT::pixel_type pixel_type;
    typedef typename ChildT::result_type result_type;
    Rotate90CCWPixelAccessor( ChildT const& acc ) : m_child(acc) {}

    inline Rotate90CCWPixelAccessor& next_col() { m_child.next_row(); return *this; }
    inline Rotate90CCWPixelAccessor& prev_col() { m_child.prev_row(); return *this; }
    inline Rotate90CCWPixelAccessor& next_row() { m_child.prev_col(); return *this; }
    inline Rotate90CCWPixelAccessor& prev_row() { m_child.next_col(); return *this; }
    inline Rotate90CCWPixelAccessor& next_plane() { m_child.next_plane(); return *this; }
    inline Rotate90CCWPixelAccessor& prev_plane() { m_child.prev_plane(); return *this; }
    inline Rotate90CCWPixelAccessor& advance( ptrdiff_t di, ptrdiff_t dj, ptrdiff_t dp=0 ) { m_child.advance(-dj,di,dp); return *this; }

    inline result_type operator*() const { return *m_child; }
  };

  // Class definition
  template <class ImageT>
  class Rotate90CCWView : public ImageViewBase<Rotate90CCWView<ImageT> >
  {
  private:
    ImageT m_child;
  public:

    typedef typename ImageT::pixel_type pixel_type;
    typedef typename ImageT::result_type result_type;
    typedef Rotate90CCWPixelAccessor<typename ImageT::pixel_accessor> pixel_accessor;

    Rotate90CCWView( ImageT const& image ) : m_child(image) {}

    inline int32 cols() const { return m_child.rows(); }
    inline int32 rows() const { return m_child.cols(); }
    inline int32 planes() const { return m_child.planes(); }

    inline pixel_accessor origin() const { return m_child.origin().advance(rows()-1,0); }

    inline result_type operator()( int32 i, int32 j, int32 p=0 ) const { return m_child(rows()-1-j,i,p); }

    template <class ViewT>
    Rotate90CCWView const& operator=( ImageViewBase<ViewT> const& view ) const {
      view.impl().rasterize( *this, BBox2i(0,0,view.impl().cols(),view.impl().rows()) );
      return *this;
    }

    template <class ViewT>
    Rotate90CCWView& operator=( ImageViewBase<ViewT> const& view ) {
      *const_cast<const Rotate90CCWView*>(this) = view.impl();
      return *this;
    }

    ImageT const& child() const {
      return m_child;
    }

    /// \cond INTERNAL
    typedef Rotate90CCWView<typename ImageT::prerasterize_type> prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i const& bbox ) const {
      BBox2i child_bbox( rows()-bbox.max().y(), bbox.min().x(), bbox.height(), bbox.width() );
      return prerasterize_type( m_child.prerasterize(child_bbox) );
    }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const { vw::rasterize( prerasterize(bbox), dest, bbox ); }
    /// \endcond
  };

  /// \cond INTERNAL
  template <class ImageT>
  struct IsMultiplyAccessible<Rotate90CCWView<ImageT> > : public IsMultiplyAccessible<ImageT> {};
  /// \endcond

  /// Rotate an image 90 degrees counter-clockwise.
  template <class ImageT>
  Rotate90CCWView<ImageT> rotate_90_ccw( ImageViewBase<ImageT> const& v ) {
    return Rotate90CCWView<ImageT>( v.impl() );
  }


  // *******************************************************************
  // FlipVertical
  // *******************************************************************

  // Specialized pixel accessor
  template <class ChildT>
  class FlipVerticalPixelAccessor
  {
  private:
    ChildT m_child;
  public:
    typedef typename ChildT::pixel_type pixel_type;
    typedef typename ChildT::result_type result_type;
    FlipVerticalPixelAccessor( ChildT const& acc ) : m_child(acc) {}

    inline FlipVerticalPixelAccessor& next_col() { m_child.next_col(); return *this; }
    inline FlipVerticalPixelAccessor& prev_col() { m_child.prev_col(); return *this; }
    inline FlipVerticalPixelAccessor& next_row() { m_child.prev_row(); return *this; }
    inline FlipVerticalPixelAccessor& prev_row() { m_child.next_row(); return *this; }
    inline FlipVerticalPixelAccessor& next_plane() { m_child.next_plane(); return *this; }
    inline FlipVerticalPixelAccessor& prev_plane() { m_child.prev_plane(); return *this; }
    inline FlipVerticalPixelAccessor& advance( ptrdiff_t di, ptrdiff_t dj, ptrdiff_t dp=0 ) { m_child.advance(di,-dj,dp); return *this; }

    inline result_type operator*() const { return *m_child; }
  };

  // Class definition
  template <class ImageT>
  class FlipVerticalView : public ImageViewBase<FlipVerticalView<ImageT> >
  {
  private:
    ImageT m_child;
  public:

    typedef typename ImageT::pixel_type pixel_type;
    typedef typename ImageT::result_type result_type;
    typedef FlipVerticalPixelAccessor<typename ImageT::pixel_accessor> pixel_accessor;

    FlipVerticalView( ImageT const& image ) : m_child(image) {}

    inline int32 cols() const { return m_child.cols(); }
    inline int32 rows() const { return m_child.rows(); }
    inline int32 planes() const { return m_child.planes(); }

    inline pixel_accessor origin() const { return m_child.origin().advance(0,rows()-1); }

    inline result_type operator()( int32 i, int32 j, int32 p=0 ) const { return m_child(i,rows()-1-j,p); }

    ImageT const& child() const {
      return m_child;
    }

    template <class ViewT>
    FlipVerticalView const& operator=( ImageViewBase<ViewT> const& view ) const {
      view.impl().rasterize( *this, BBox2i(0,0,view.impl().cols(),view.impl().rows()) );
      return *this;
    }

    template <class ViewT>
    FlipVerticalView& operator=( ImageViewBase<ViewT> const& view ) {
      *const_cast<const FlipVerticalView*>(this) = view.impl();
      return *this;
    }

    /// \cond INTERNAL
    typedef FlipVerticalView<typename ImageT::prerasterize_type> prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i const& bbox ) const {
      BBox2i child_bbox( bbox.min().x(), rows()-bbox.max().y(), bbox.width(), bbox.height() );
      return prerasterize_type( m_child.prerasterize(child_bbox) );
    }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const { vw::rasterize( prerasterize(bbox), dest, bbox ); }
    /// \endcond
  };

  /// \cond INTERNAL
  template <class ImageT>
  struct IsMultiplyAccessible<FlipVerticalView<ImageT> > : public IsMultiplyAccessible<ImageT> {};
  /// \endcond

  /// Flip an image vertically.
  template <class ImageT>
  FlipVerticalView<ImageT> flip_vertical( ImageViewBase<ImageT> const& v ) {
    return FlipVerticalView<ImageT>( v.impl() );
  }


  // *******************************************************************
  // FlipHorizontal
  // *******************************************************************

  // Specialized pixel accessor
  template <class ChildT>
  class FlipHorizontalPixelAccessor
  {
  private:
    ChildT m_child;
  public:
    typedef typename ChildT::pixel_type pixel_type;
    typedef typename ChildT::result_type result_type;
    FlipHorizontalPixelAccessor( ChildT const& acc ) : m_child(acc) {}

    inline FlipHorizontalPixelAccessor& next_col() { m_child.prev_col(); return *this; }
    inline FlipHorizontalPixelAccessor& prev_col() { m_child.next_col(); return *this; }
    inline FlipHorizontalPixelAccessor& next_row() { m_child.next_row(); return *this; }
    inline FlipHorizontalPixelAccessor& prev_row() { m_child.prev_row(); return *this; }
    inline FlipHorizontalPixelAccessor& next_plane() { m_child.next_plane(); return *this; }
    inline FlipHorizontalPixelAccessor& prev_plane() { m_child.prev_plane(); return *this; }
    inline FlipHorizontalPixelAccessor& advance( ptrdiff_t di, ptrdiff_t dj, ptrdiff_t dp=0 ) { m_child.advance(-di,dj,dp); return *this; }

    inline result_type operator*() const { return *m_child; }
  };

  // Class definition
  template <class ImageT>
  class FlipHorizontalView : public ImageViewBase<FlipHorizontalView<ImageT> >
  {
  private:
    ImageT m_child;
  public:

    typedef typename ImageT::pixel_type pixel_type;
    typedef typename ImageT::result_type result_type;
    typedef FlipHorizontalPixelAccessor<typename ImageT::pixel_accessor> pixel_accessor;

    FlipHorizontalView( ImageT const& image ) : m_child(image) {}

    inline int32 cols() const { return m_child.cols(); }
    inline int32 rows() const { return m_child.rows(); }
    inline int32 planes() const { return m_child.planes(); }

    inline pixel_accessor origin() const { return m_child.origin().advance(cols()-1,0); }

    inline result_type operator()( int32 i, int32 j, int32 p=0 ) const { return m_child(cols()-1-i,j,p); }

    ImageT const& child() const {
      return m_child;
    }

    template <class ViewT>
    FlipHorizontalView const& operator=( ImageViewBase<ViewT> const& view ) const {
      view.impl().rasterize( *this, BBox2i(0,0,view.impl().cols(),view.impl().rows()) );
      return *this;
    }

    template <class ViewT>
    FlipHorizontalView& operator=( ImageViewBase<ViewT> const& view ) {
      *const_cast<const FlipHorizontalView*>(this) = view.impl();
      return *this;
    }

    /// \cond INTERNAL
    typedef FlipHorizontalView<typename ImageT::prerasterize_type> prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i const& bbox ) const {
      BBox2i child_bbox( cols()-bbox.max().x(), bbox.min().y(), bbox.width(), bbox.height() );
      return prerasterize_type( m_child.prerasterize(child_bbox) );
    }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const { vw::rasterize( prerasterize(bbox), dest, bbox ); }
    /// \endcond
  };

  /// \cond INTERNAL
  template <class ImageT>
  struct IsMultiplyAccessible<FlipHorizontalView<ImageT> > : public IsMultiplyAccessible<ImageT> {};
  /// \endcond

  /// Flip an image horizontally.
  template <class ImageT>
  FlipHorizontalView<ImageT> flip_horizontal( ImageViewBase<ImageT> const& v ) {
    return FlipHorizontalView<ImageT>( v.impl() );
  }


  // *******************************************************************
  // crop()
  // *******************************************************************

  // Class definition
  template <class ImageT>
  class CropView : public ImageViewBase< CropView<ImageT> > {
  private:
    typedef typename boost::mpl::if_<IsFloatingPointIndexable<ImageT>, double, int32>::type offset_type;

    ImageT m_child;
    offset_type m_ci, m_cj;
    int32 m_di, m_dj;

  public:
    typedef typename ImageT::pixel_type pixel_type;
    typedef typename ImageT::result_type result_type;
    typedef typename ImageT::pixel_accessor pixel_accessor;

    CropView( ImageT const& image, offset_type const upper_left_i, offset_type const upper_left_j, int32 const width, int32 const height ) :
      m_child(image), m_ci(upper_left_i), m_cj(upper_left_j), m_di(width), m_dj(height) {}

    template<class RealT>
    CropView( ImageT const& image, BBox<RealT,2> const& bbox) :
      m_child(image),
      m_ci((offset_type)(bbox.min()[0])),
      m_cj((offset_type)(bbox.min()[1])),
      m_di(int32(.5+(bbox.width()))),
      m_dj(int32(.5+(bbox.height()))) {}

    inline int32 cols() const { return m_di; }
    inline int32 rows() const { return m_dj; }
    inline int32 planes() const { return m_child.planes(); }

    inline pixel_accessor origin() const { return m_child.origin().advance(m_ci, m_cj); }

    inline result_type operator()( offset_type i, offset_type j, int32 p=0 ) const { return m_child(m_ci + i, m_cj + j, p); }

    CropView const& operator=( CropView const& view ) const {
      view.rasterize( *this, BBox2i(0,0,view.impl().cols(),view.impl().rows()) );
      return *this;
    }

    template <class ViewT>
    CropView const& operator=( ImageViewBase<ViewT> const& view ) const {
      view.impl().rasterize( *this, BBox2i(0,0,view.impl().cols(),view.impl().rows()) );
      return *this;
    }

    ImageT const& child() const {
      return m_child;
    }

    /// \cond INTERNAL
    typedef CropView<typename ImageT::prerasterize_type> prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i const& bbox ) const {
      return prerasterize_type( m_child.prerasterize(bbox+Vector2i(m_ci,m_cj)), m_ci, m_cj, m_di, m_dj );
    }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const {
      // FIXME Warning: This does not respect floating-point offsets!
      m_child.rasterize( dest, bbox+Vector2i(int(m_ci),int(m_cj)) );
    }
    /// \endcond
  };

  /// \cond INTERNAL
  // Type traits
  template <class ImageT>
  struct IsFloatingPointIndexable<CropView<ImageT> >  : public IsFloatingPointIndexable<ImageT> {};

  template <class ImageT>
  struct IsMultiplyAccessible<CropView<ImageT> > : public IsMultiplyAccessible<ImageT> {};
  /// \endcond

  /// Crop an image.
  template <class ImageT>
  inline CropView<ImageT> crop( ImageViewBase<ImageT> const& v, int32 upper_left_x, int32 upper_left_y, int32 width, int32 height ) {
    return CropView<ImageT>( v.impl(), upper_left_x, upper_left_y, width, height );
  }

  /// Crop an image.
  template <class ImageT, class BBoxRealT>
  inline CropView<ImageT> crop( ImageViewBase<ImageT> const& v, BBox<BBoxRealT,2> const& bbox ) {
    return CropView<ImageT>( v.impl(), bbox.min().x(), bbox.min().y(), int(.5+(bbox.width())), int(.5+(bbox.height())) );
  }


  // *******************************************************************
  // subsample()
  // *******************************************************************

  // Specialized image accessor
  template <class ChildT>
  class SubsamplePixelAccessor {
    ChildT m_child;
    int32 m_xdelta, m_ydelta;
  public:
    typedef typename ChildT::pixel_type pixel_type;
    typedef typename ChildT::result_type result_type;
    SubsamplePixelAccessor( ChildT const& acc , int32 xdelta, int32 ydelta) : m_child(acc), m_xdelta(xdelta), m_ydelta(ydelta) {}

    inline SubsamplePixelAccessor& next_col() { m_child.advance(  m_xdelta, 0 ); return *this; }
    inline SubsamplePixelAccessor& prev_col() { m_child.advance( -m_xdelta, 0 ); return *this; }
    inline SubsamplePixelAccessor& next_row() { m_child.advance( 0,  m_ydelta ); return *this; }
    inline SubsamplePixelAccessor& prev_row() { m_child.advance( 0, -m_ydelta ); return *this; }
    inline SubsamplePixelAccessor& next_plane() { m_child.next_plane(); return *this; }
    inline SubsamplePixelAccessor& prev_plane() { m_child.prev_plane(); return *this; }
    inline SubsamplePixelAccessor& advance( ptrdiff_t di, ptrdiff_t dj, ptrdiff_t dp=0 ) { m_child.advance((ptrdiff_t)m_xdelta*di,(ptrdiff_t)m_ydelta*dj,dp); return *this; }

    inline result_type operator*() const { return *m_child; }
  };

  // Class definition
  template <class ImageT>
  class SubsampleView : public ImageViewBase<SubsampleView<ImageT> > {
    ImageT m_child;
    int32 m_xdelta, m_ydelta;
  public:
    typedef typename ImageT::pixel_type pixel_type;
    typedef typename ImageT::result_type result_type;
    typedef SubsamplePixelAccessor<typename ImageT::pixel_accessor> pixel_accessor;

    SubsampleView( ImageT const& image, int32 subsampling_factor ) : m_child(image), m_xdelta(subsampling_factor), m_ydelta(subsampling_factor) {}
    SubsampleView( ImageT const& image, int32 xfactor, int32 yfactor ) : m_child(image), m_xdelta(xfactor), m_ydelta(yfactor) {}

    inline int32 cols() const { return 1 + (m_child.cols()-1)/m_xdelta; }
    inline int32 rows() const { return 1 + (m_child.rows()-1)/m_ydelta; }
    inline int32 planes() const { return m_child.planes(); }

    inline pixel_accessor origin() const { return pixel_accessor(m_child.origin(), m_xdelta, m_ydelta); }
    inline result_type operator()( int32 i, int32 j, int32 p=0 ) const { return m_child(m_xdelta*i,m_ydelta*j,p); }

    ImageT const& child() const { return m_child; }

    /// \cond INTERNAL
    typedef SubsampleView<typename ImageT::prerasterize_type> prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i const& bbox ) const {
      return prerasterize_type( m_child.prerasterize(BBox2i(m_xdelta*bbox.min().x(),m_ydelta*bbox.min().y(),m_xdelta*bbox.width(),m_ydelta*bbox.height())), m_xdelta, m_ydelta );
    }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const { vw::rasterize( prerasterize(bbox), dest, bbox ); }
    /// \endcond
  };

  /// \cond INTERNAL
  template <class ImageT>
  struct IsMultiplyAccessible<SubsampleView<ImageT> > : public IsMultiplyAccessible<ImageT> {};
  /// \endcond

  /// Subsample an image by an integer factor.  Note that this
  /// function does not pre-smooth the image prior to subsampling: it
  /// simply selects every Nth pixel.  You will typically want to
  /// apply some sort of anti-aliasing filter prior to calling this
  /// function.
  template <class ImageT>
  inline SubsampleView<ImageT> subsample( ImageT const& v, int32 subsampling_factor ) {
    return SubsampleView<ImageT>( v, subsampling_factor );
  }

  /// Subsample an image by integer factors in x and y.  Note that
  /// this function does not pre-smooth the image prior to
  /// subsampling: it simply selects every Nth pixel.  You will
  /// typically want to apply some sort of anti-aliasing filter prior
  /// to calling this function.
  template <class ImageT>
  inline SubsampleView<ImageT> subsample( ImageT const& v, int32 xfactor, int32 yfactor ) {
    return SubsampleView<ImageT>( v, xfactor, yfactor );
  }


  // *******************************************************************
  // select_col()
  // *******************************************************************

  /// Return a single column from an image
  /// \see vw::select_col
  template <class ImageT>
  class SelectColView : public ImageViewBase<SelectColView<ImageT> > {
    ImageT m_child;
    int32 m_col;
  public:
    typedef typename ImageT::pixel_type pixel_type;
    typedef typename ImageT::result_type result_type;
    typedef typename ImageT::pixel_accessor pixel_accessor;

    SelectColView( ImageT const& image, int32 col ) : m_child(image), m_col(col) {}

    inline int32 cols() const { return 1; }
    inline int32 rows() const { return m_child.rows(); }
    inline int32 planes() const { return m_child.planes(); }

    inline pixel_accessor origin() const { return m_child.origin().advance(m_col,0,0); }
    inline result_type operator()( int32 /*i*/, int32 j, int32 p=0) const { return m_child(m_col,j,p); }

    SelectColView const& operator=( SelectColView const& view ) const {
      view.rasterize( *this, BBox2i(0,0,view.cols(),view.rows()) );
      return *this;
    }

    template <class ViewT>
    SelectColView const& operator=( ImageViewBase<ViewT> const& view ) const {
      view.impl().rasterize( *this, BBox2i(0,0,view.impl().cols(),view.impl().rows()) );
      return *this;
    }

    /// \cond INTERNAL
    typedef SelectColView<typename ImageT::prerasterize_type> prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i const& bbox ) const { return prerasterize_type( m_child.prerasterize(bbox), m_col ); }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const { vw::rasterize( prerasterize(bbox), dest, bbox ); }
    /// \endcond
  };

  /// \cond INTERNAL
  // View type Traits
  template <class ImageT>
  struct IsMultiplyAccessible<SelectColView<ImageT> > : public IsMultiplyAccessible<ImageT> {};
  /// \endcond

  /// Extracts a single column of an image.  This function returns a
  /// writeable view of a single column of a multi-column image.
  /// \see vw::SelectColView
  template <class ImageT>
  SelectColView<ImageT> select_col( ImageViewBase<ImageT> const& v, int32 col ) {
    return SelectColView<ImageT>( v.impl(), col );
  }


  // *******************************************************************
  // select_row()
  // *******************************************************************

  /// Return a single row from an image
  /// \see vw::select_row
  template <class ImageT>
  class SelectRowView : public ImageViewBase<SelectRowView<ImageT> > {
    ImageT m_child;
    int32 m_row;
  public:
    typedef typename ImageT::pixel_type pixel_type;
    typedef typename ImageT::result_type result_type;
    typedef typename ImageT::pixel_accessor pixel_accessor;

    SelectRowView( ImageT const& image, int32 row ) : m_child(image), m_row(row) {}

    inline int32 cols() const { return m_child.cols(); }
    inline int32 rows() const { return 1; }
    inline int32 planes() const { return m_child.planes(); }

    inline pixel_accessor origin() const { return m_child.origin().advance(0,m_row,0); }

    inline result_type operator()( int32 i, int32 /*j*/, int32 p=0) const { return m_child(i,m_row,p); }

    SelectRowView const& operator=( SelectRowView const& view ) const {
      view.rasterize( *this, BBox2i(0,0,view.cols(),view.rows()) );
      return *this;
    }

    template <class ViewT>
    SelectRowView const& operator=( ImageViewBase<ViewT> const& view ) const {
      view.impl().rasterize( *this, BBox2i(0,0,view.impl().cols(),view.impl().rows()) );
      return *this;
    }

    typedef SelectRowView<typename ImageT::prerasterize_type> prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i const& bbox ) const { return prerasterize_type( m_child.prerasterize(bbox), m_row ); }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const { vw::rasterize( prerasterize(bbox), dest, bbox ); }
  };

  // View type Traits
  template <class ImageT>
  struct IsMultiplyAccessible<SelectRowView<ImageT> > : public IsMultiplyAccessible<ImageT> {};

  /// Extracts a single row of an image.  This function returns a
  /// writeable view of a single row of a multi-row image.
  /// \see vw::SelectRowView
  template <class ImageT>
  SelectRowView<ImageT> select_row( ImageViewBase<ImageT> const& v, int32 row ) {
    return SelectRowView<ImageT>( v.impl(), row );
  }


  // *******************************************************************
  // select_plane()
  // *******************************************************************

  /// Return a single plane from a multi-plane image
  /// \see vw::select_plane
  template <class ImageT>
  class SelectPlaneView : public ImageViewBase<SelectPlaneView<ImageT> > {
    ImageT m_child;
    int32 m_plane;
  public:
    typedef typename ImageT::pixel_type pixel_type;
    typedef typename ImageT::result_type result_type;
    typedef typename ImageT::pixel_accessor pixel_accessor;

    SelectPlaneView( ImageT const& image, int32 plane ) : m_child(image), m_plane(plane) {}

    inline int32 cols() const { return m_child.cols(); }
    inline int32 rows() const { return m_child.rows(); }
    inline int32 planes() const { return 1; }

    inline pixel_accessor origin() const { return m_child.origin().advance(0,0,m_plane); }
    inline result_type operator()( int32 i, int32 j, int32 p=0 ) const { return m_child(i,j,m_plane+p); }

    template <class ViewT>
    SelectPlaneView const& operator=( ImageViewBase<ViewT> const& view ) const {
      view.impl().rasterize( *this, BBox2i(0,0,view.impl().cols(),view.impl().rows()) );
      return *this;
    }

    template <class ViewT>
    SelectPlaneView& operator=( ImageViewBase<ViewT> const& view ) {
      *const_cast<const SelectPlaneView*>(this) = view.impl();
      return *this;
    }

    /// \cond INTERNAL
    typedef SelectPlaneView<typename ImageT::prerasterize_type> prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i const& bbox ) const { return prerasterize_type( m_child.prerasterize(bbox), m_plane ); }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const { vw::rasterize( prerasterize(bbox), dest, bbox ); }
    /// \endcond
  };

  /// \cond INTERNAL
  // View type Traits
  template <class ImageT>
  struct IsMultiplyAccessible<SelectPlaneView<ImageT> > : public IsMultiplyAccessible<ImageT> {};
  /// \endcond

  /// Extracts a single plane of a multi-plane image.  This function
  /// returns a writeable view of a single plane of a multi-plane
  /// image.  \see vw::SelectPlaneView
  template <class ImageT>
  SelectPlaneView<ImageT> select_plane( ImageViewBase<ImageT> const& v, int32 plane ) {
    return SelectPlaneView<ImageT>( v.impl(), plane );
  }


  // *******************************************************************
  // select_channel()
  // *******************************************************************

  /// A channel selecting functor, used by \ref select_channel().
  template <class ImageT>
  struct SelectChannelFunctor {
    int32 m_channel;
  public:
    SelectChannelFunctor( int32 channel ) : m_channel(channel) {}

    // Computes an appropriate reference-to-channel type.
    typedef typename CompoundChannelType<typename ImageT::pixel_type>::type channel_type;
    typedef typename CopyCVR<typename ImageT::result_type,channel_type>::type result_type;

    result_type operator()( typename ImageT::result_type pixel ) const {
      return compound_select_channel<result_type>(pixel,m_channel);
    }
  };

  /// Extracts a single channel of a multi-channel image.  This function
  /// returns a writeable view of a single channel of a multi-channel
  /// image.
  template <class ImageT>
  UnaryPerPixelView<ImageT,SelectChannelFunctor<ImageT> >
  inline select_channel( ImageViewBase<ImageT>& image, int32 channel ) {
    return UnaryPerPixelView<ImageT,SelectChannelFunctor<ImageT> >( image.impl(), SelectChannelFunctor<ImageT>(channel) );
  }

  /// Extracts a single channel of a multi-channel image (const overload).
  /// This function returns a writeable view of a single channel of a
  /// multi-channel image.
  template <class ImageT>
  UnaryPerPixelView<ImageT,SelectChannelFunctor<const ImageT> >
  inline select_channel( ImageViewBase<ImageT> const& image, int32 channel ) {
    return UnaryPerPixelView<ImageT,SelectChannelFunctor<const ImageT> >( image.impl(), SelectChannelFunctor<const ImageT>(channel) );
  }

  /// A convenience function to select the alpha channel of an image.
  template <class ImageT>
  UnaryPerPixelView<ImageT,SelectChannelFunctor<ImageT> >
  inline select_alpha_channel( ImageViewBase<ImageT>& image ) {
    // FIXME: This should be a static assertion
    if( ! PixelHasAlpha<typename ImageT::pixel_type>::value )
      vw_throw( ArgumentErr() << "Image has no alpha channel in call to select_alpha_channel()" );
    return select_channel( image, PixelNumChannels<typename ImageT::pixel_type>::value - 1 );
  }

  /// A convenience function to select the alpha channel of an image
  /// (const overload).
  template <class ImageT>
  UnaryPerPixelView<ImageT,SelectChannelFunctor<const ImageT> >
  inline select_alpha_channel( ImageViewBase<ImageT> const& image ) {
    // FIXME: This should be a static assertion
    if( ! PixelHasAlpha<typename ImageT::pixel_type>::value )
      vw_throw( ArgumentErr() << "Image has no alpha channel in call to select_alpha_channel()" );
    return select_channel( image, PixelNumChannels<typename ImageT::pixel_type>::value - 1 );
  }


  // *******************************************************************
  // channels_to_planes()
  // *******************************************************************

  /// A channels-to-planes pixel accessor adaptor.
  ///
  /// This is a special wrapper pixel accessor type, used by
  /// \ref vw::ChannelsToPlanesView, that treats the channels in a
  /// multi-channel image as planes.
  template <class ChildT>
  class ChannelsToPlanesAccessor
  {
  private:
    ChildT m_child;
    int32 m_channel;
  public:
    typedef typename CompoundChannelType<typename ChildT::pixel_type>::type pixel_type;
    typedef typename CopyCVR<typename ChildT::result_type, pixel_type>::type result_type;

    ChannelsToPlanesAccessor( ChildT const& acc ) : m_child(acc), m_channel(0) {}
    inline ChannelsToPlanesAccessor& next_col() { m_child.next_col(); return *this; }
    inline ChannelsToPlanesAccessor& prev_col() { m_child.prev_col(); return *this; }
    inline ChannelsToPlanesAccessor& next_row() { m_child.next_row(); return *this; }
    inline ChannelsToPlanesAccessor& prev_row() { m_child.prev_row(); return *this; }
    inline ChannelsToPlanesAccessor& next_plane() { ++m_channel; return *this; }
    inline ChannelsToPlanesAccessor& prev_plane() { --m_channel; return *this; }
    inline ChannelsToPlanesAccessor& advance( ptrdiff_t di, ptrdiff_t dj, ptrdiff_t dp=0 ) { m_child.advance(di,dj); m_channel+=dp; return *this; }

    inline result_type operator*() const { return compound_select_channel<result_type>(*m_child,m_channel); }
  };

  /// A view that turns a one plane, multi-channel view into a mulit-plane, one channel view.
  /// \see vw::channels_to_planes
  template <class ImageT>
  class ChannelsToPlanesView : public ImageViewBase<ChannelsToPlanesView<ImageT> > {
    ImageT m_child;
  public:

    typedef typename CompoundChannelType<typename ImageT::pixel_type>::type pixel_type;
    typedef typename CopyCVR<typename ImageT::result_type, pixel_type>::type result_type;

    typedef typename boost::mpl::if_< IsCompound<typename ImageT::pixel_type>,
                                      ChannelsToPlanesAccessor<typename ImageT::pixel_accessor>,
                                      typename ImageT::pixel_accessor >::type pixel_accessor;

    ChannelsToPlanesView( ImageT const& image ) : m_child(image) {
      VW_ASSERT( m_child.planes()==1 , ArgumentErr() << "ChannelsToPlanesView: The image must be single plane.");
    }

    inline int32 cols() const { return m_child.cols(); }
    inline int32 rows() const { return m_child.rows(); }
    inline int32 planes() const { return m_child.channels(); }

    inline pixel_accessor origin() const { return pixel_accessor(m_child.origin()); }

    inline result_type operator()( int32 i, int32 j, int32 p=0 ) const {
      return compound_select_channel<result_type>(m_child(i,j),p);
    }

    template <class ViewT>
    ChannelsToPlanesView const& operator=( ImageViewBase<ViewT> const& view ) const {
      view.impl().rasterize( *this, BBox2i(0,0,view.impl().cols(),view.impl().rows()) );
      return *this;
    }

    template <class ViewT>
    ChannelsToPlanesView& operator=( ImageViewBase<ViewT> const& view ) {
      *const_cast<const ChannelsToPlanesView*>(this) = view.impl();
      return *this;
    }

    ImageT const& child() const {
      return m_child;
    }

    /// \cond INTERNAL
    typedef ChannelsToPlanesView<typename ImageT::prerasterize_type> prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i const& bbox ) const { return prerasterize_type( m_child.prerasterize(bbox) ); }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const { vw::rasterize( prerasterize(bbox), dest, bbox ); }
    /// \endcond
  };

  /// \cond INTERNAL
  // View type Traits
  template <class ImageT>
  struct IsMultiplyAccessible<ChannelsToPlanesView<ImageT> > : public IsMultiplyAccessible<ImageT> {};
  /// \endcond

  /// Adapts a multi-channel image view so the channels are treated as
  /// planes.  This function returns a writeable view of of a
  /// single-plane, multi-channel image that treats the channels as
  /// planes.  This is primarily intended to simplify interfacing with
  /// legacy non-Vision-Workbench code that can only support
  /// fundamental pixel types.  If you are thinking about using this
  /// function and you are not trying to interface to legacy code then
  /// you are almost certainly doing something wrong.
  /// \see vw::ChannelsToPlanesView
  template <class ImageT>
  ChannelsToPlanesView<ImageT> channels_to_planes( ImageViewBase<ImageT> const& v ) {
    return ChannelsToPlanesView<ImageT>( v.impl() );
  }


  // *******************************************************************
  // planes_to_channels()
  // *******************************************************************

  /// A view that turns a multi-plane, single-channel view into a
  /// one-plane, multi-channel view.
  /// \see vw::planes_to_channels
  template <class PixelT, class ImageT>
  class PlanesToChannelsView : public ImageViewBase<PlanesToChannelsView<PixelT,ImageT> > {
    ImageT m_child;
  public:

    typedef PixelT pixel_type;
    typedef PixelT result_type;
    typedef ProceduralPixelAccessor<PlanesToChannelsView> pixel_accessor;

    PlanesToChannelsView( ImageT const& image ) : m_child(image) {
      VW_ASSERT( m_child.channels()==1 &&
                 boost::numeric_cast<size_t>(m_child.planes())==CompoundNumChannels<PixelT>::value,
                 ArgumentErr() << "PlanesToChannelsView: The image must be multi-plane, single-channel.");
    }

    inline int32 cols() const { return m_child.cols(); }
    inline int32 rows() const { return m_child.rows(); }
    inline int32 planes() const { return 1; }

    inline pixel_accessor origin() const { return pixel_accessor(*this); }

    inline result_type operator()( int32 i, int32 j, int32 /*p*/=0 ) const {
      result_type result;
      typedef typename CompoundChannelType<result_type>::type channel_type;
      for ( size_t c=0; c<CompoundNumChannels<PixelT>::value; ++c )
        compound_select_channel<channel_type&>(result,c) = m_child(i,j,c);
      return result;
    }

    ImageT const& child() const {
      return m_child;
    }

    /// \cond INTERNAL
    typedef PlanesToChannelsView<PixelT, typename ImageT::prerasterize_type> prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i const& bbox ) const { return prerasterize_type( m_child.prerasterize(bbox) ); }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const { vw::rasterize( prerasterize(bbox), dest, bbox ); }
    /// \endcond
  };

  /// Adapts a multi-plane image view so the planes are treated as
  /// channels of the given pixel type.  This is primarily intended to
  /// simplify interfacing with legacy non-Vision-Workbench code that
  /// can only support fundamental pixel types.  If you are thinking
  /// about using this function and you are not trying to interface to
  /// legacy code then you are almost certainly doing something wrong.
  /// \see vw::PlanesToChannelsView
  /// \see vw::channels_to_planes
  template <class PixelT, class ImageT>
  PlanesToChannelsView<PixelT,ImageT> planes_to_channels( ImageViewBase<ImageT> const& v ) {
    return PlanesToChannelsView<PixelT,ImageT>( v.impl() );
  }


  // *******************************************************************
  // pixel_cast()
  // *******************************************************************

  /// A pixel casting functor, used by \ref pixel_cast().
  template <class PixelT>
  struct PixelCastFunctor : ReturnFixedType<PixelT> {
    template <class ArgT>
    inline PixelT operator()( ArgT pixel ) const {
      return pixel_cast<PixelT>(pixel);
    }
  };

  /// Create a new image view by statically casting the pixels to a
  /// new type.
  template <class PixelT, class ImageT>
  inline UnaryPerPixelView<ImageT,PixelCastFunctor<PixelT> > pixel_cast( ImageViewBase<ImageT> const& image ) {
    return UnaryPerPixelView<ImageT,PixelCastFunctor<PixelT> >( image.impl() );
  }


  // *******************************************************************
  // pixel_cast_rescale()
  // *******************************************************************

  /// A pixel casting functor, used by \ref pixel_cast_rescale().
  template <class PixelT>
  struct PixelCastRescaleFunctor : ReturnFixedType<PixelT> {
    template <class ArgT>
    inline PixelT operator()( ArgT pixel ) const {
      return pixel_cast_rescale<PixelT>(pixel);
    }
  };

  /// Create a new image view by casting and rescaling the pixels to a
  /// new type.
  template <class PixelT, class ImageT>
  inline UnaryPerPixelView<ImageT,PixelCastRescaleFunctor<PixelT> > pixel_cast_rescale( ImageViewBase<ImageT> const& image ) {
    return UnaryPerPixelView<ImageT,PixelCastRescaleFunctor<PixelT> >( image.impl() );
  }


  // *******************************************************************
  // channel_cast()
  // *******************************************************************

  /// A pixel channel casting functor, used by \ref channel_cast().
  template <class ChannelT>
  struct PixelChannelCastFunctor : UnaryReturnBinaryTemplateBind2nd<CompoundChannelCast,ChannelT> {
    template <class ArgT>
    inline typename CompoundChannelCast<ArgT,ChannelT>::type operator()( ArgT const& pixel ) const {
      return channel_cast<ChannelT>(pixel);
    }
  };

  /// Create a new image view by statically casting the channels of the pixels to a new type.
  template <class ChannelT, class ImageT>
  inline UnaryPerPixelView<ImageT,PixelChannelCastFunctor<ChannelT> > channel_cast( ImageViewBase<ImageT> const& image ) {
    return UnaryPerPixelView<ImageT,PixelChannelCastFunctor<ChannelT> >( image.impl() );
  }


  // *******************************************************************
  // channel_cast_rescale()
  // *******************************************************************

  /// A pixel channel casting and rescaling functor, used by
  /// \ref channel_cast_rescale().
  template <class ChannelT>
  struct PixelChannelCastRescaleFunctor : UnaryReturnBinaryTemplateBind2nd<CompoundChannelCast,ChannelT> {
    template <class ArgT>
    inline typename CompoundChannelCast<ArgT,ChannelT>::type operator()( ArgT const& pixel ) const {
      return channel_cast_rescale<ChannelT>(pixel);
    }
  };

  /// Create a new image view by casting and rescaling the channels of
  /// the pixels to a new type.
  template <class ChannelT, class ImageT>
  inline UnaryPerPixelView<ImageT,PixelChannelCastRescaleFunctor<ChannelT> > channel_cast_rescale( ImageViewBase<ImageT> const& image ) {
    return UnaryPerPixelView<ImageT,PixelChannelCastRescaleFunctor<ChannelT> >( image.impl() );
  }


  // *******************************************************************
  // weighted_rgb_to_gray()
  // *******************************************************************

  /// A weighted rgb-to-gray pixel conversion functor.
  class WeightedRGBToGrayFunctor {
    double m_rw, m_gw, m_bw;
  public:
    template <class ArgsT> struct result {};
    template <class FuncT, class ChannelT> struct result<FuncT(PixelRGB<ChannelT>)> { typedef PixelGray<ChannelT> type; };
    template <class FuncT, class ChannelT> struct result<FuncT(PixelRGBA<ChannelT>)> { typedef PixelGrayA<ChannelT> type; };
    WeightedRGBToGrayFunctor( double rw, double gw, double bw ) : m_rw(rw), m_gw(gw), m_bw(bw) {}
    template <class ChannelT> inline PixelGrayA<ChannelT> operator()( PixelRGBA<ChannelT> const& rgb ) const {
      return weighted_rgb_to_gray( rgb, m_rw, m_gw, m_bw );
    }
    template <class ChannelT> inline PixelGray<ChannelT> operator()( PixelRGB<ChannelT> const& rgb ) const {
      return weighted_rgb_to_gray( rgb, m_rw, m_gw, m_bw );
    }
  };


  /// Weighted conversion from PixelRGBA to PixelGrayA using user-specified weights.
  template <class ImageT>
  inline UnaryPerPixelView<ImageT,WeightedRGBToGrayFunctor> weighted_rgb_to_gray( ImageViewBase<ImageT> const& image, double rw, double gw, double bw ) {
    return UnaryPerPixelView<ImageT,WeightedRGBToGrayFunctor>( image.impl(), WeightedRGBToGrayFunctor(rw,gw,bw) );
  }

  /// Weighted conversion from PixelRGBA to PixelGrayA using the default weights.
  template <class ImageT>
  inline UnaryPerPixelView<ImageT,WeightedRGBToGrayFunctor> weighted_rgb_to_gray( ImageViewBase<ImageT> const& image ) {
    WeightedRGBToGrayFunctor func( VW_RGB_TO_GRAY_R_WEIGHT, VW_RGB_TO_GRAY_G_WEIGHT, VW_RGB_TO_GRAY_B_WEIGHT );
    return UnaryPerPixelView<ImageT,WeightedRGBToGrayFunctor>( image.impl(), func );
  }

} // namespace vw

#endif // __VW_IMAGE_MANIPULATION_H__
