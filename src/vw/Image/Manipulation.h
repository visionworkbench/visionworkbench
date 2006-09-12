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

/// \file Manipulation.h
///
/// Simple image manipulation functions, such as flipping and
/// cropping.  All of the functions in this file return <i>shallow</i>
/// views of the ImageView.  That is, they do not copy the data underneath 
/// but instead they refer to the same data, indexing and accessing it in 
/// a different way.
/// 
/// The first collection of functions in this file perform basic
/// transformations to the domain of the image, such as transposition,
/// rotation by 90 degree increments, flipping, and cropping.
/// 
/// This file also provides views and functions that take simple 
/// "slices" of images, returning a new view composed from individual
/// channels or planes of the source image. These include:
/// 
/// - select_plane() : takes a single-plane slice of an image
/// - select_channel() : takes a single-channel slice of an image
/// - channels_to_planes() : reinterprets a multi-channel image as a multi-plane image
///
#ifndef __VW_IMAGE__MANIPULATION_H__
#define __VW_IMAGE__MANIPULATION_H__

#include <boost/mpl/logical.hpp>

#include <vw/Core/CompoundTypes.h>
#include <vw/Image/ImageViewBase.h>
//#include <vw/BBox.h>

namespace vw {

  // *******************************************************************
  // Copy
  // *******************************************************************

  // Class definition
  template <class ImageT>
  class CopyView : public ImageViewBase<CopyView<ImageT> >
  {
  private:
    ImageView<typename ImageT::pixel_type> m_image;
  public:

    typedef typename ImageView<typename ImageT::pixel_type>::pixel_type pixel_type;
    typedef typename ImageView<typename ImageT::pixel_type>::pixel_accessor pixel_accessor;

    CopyView( ImageT const& image ) : m_image(image.cols(),image.rows(),image.planes()) {
      vw::rasterize( image, m_image );
    }

    inline unsigned cols() const { return m_image.cols(); }
    inline unsigned rows() const { return m_image.rows(); }
    inline unsigned planes() const { return m_image.planes(); }

    inline pixel_accessor origin() const { return m_image.origin(); }

    typename boost::mpl::if_< IsReferenceable<ImageT>, pixel_type&, pixel_type >::type
    inline operator()( int i, int j ) const { return m_image(i,j); }

    typename boost::mpl::if_< IsReferenceable<ImageT>, pixel_type&, pixel_type >::type
    inline operator()( int i, int j, int p ) const { return m_image(i,j,p); }

    typedef CopyView prerasterize_type;
    inline prerasterize_type prerasterize() const { return *this; }
    template <class DestT> inline void rasterize( DestT const& dest ) const { vw::rasterize( m_image, dest ); }
  };

  // Type Traits
  template <class ImageT>
  struct IsReferenceable<CopyView<ImageT> > : public boost::true_type {};

  template <class ImageT>
  struct IsMultiplyAccessible<CopyView<ImageT> > : public boost::true_type {};

  /// \endcond

  /// Make a (deep) copy of an image.
  template <class ImageT>
  CopyView<ImageT> copy( ImageViewBase<ImageT> const& v ) {
    return CopyView<ImageT>( v.impl() );
  }


  // *******************************************************************
  // Transpose
  // *******************************************************************

  // Specialized pixel accessor
  template <class ImageAccT>
  class TransposePixelAccessor
  {
  private:
    ImageAccT m_acc;
  public:
    typedef typename ImageAccT::pixel_type pixel_type;
    TransposePixelAccessor( ImageAccT const& iter ) : m_acc(iter) {}
    inline TransposePixelAccessor& next_col() { m_acc.next_row(); return *this; }
    inline TransposePixelAccessor& prev_col() { m_acc.prev_row(); return *this; }
    inline TransposePixelAccessor& next_row() { m_acc.next_col(); return *this; }
    inline TransposePixelAccessor& prev_row() { m_acc.prev_col(); return *this; }
    inline TransposePixelAccessor& next_plane() { m_acc.next_plane(); return *this; }
    inline TransposePixelAccessor& prev_plane() { m_acc.prev_plane(); return *this; }
    inline TransposePixelAccessor& advance( ptrdiff_t di, ptrdiff_t dj, ptrdiff_t dp=0 ) { m_acc.advance(dj,di,dp); return *this; }

    typename boost::mpl::if_< IsReferenceable<ImageAccT>, pixel_type&, pixel_type >::type
    inline operator*() const { return *m_acc; }
  };

  /// \cond INTERNAL
  // Accessor type traits
  template <class ImageAccT>
  struct IsReferenceable<TransposePixelAccessor<ImageAccT> > : public IsReferenceable<ImageAccT> {};
  /// \endcond

  // Class definition
  template <class ImageT>
  class TransposeView : public ImageViewBase<TransposeView<ImageT> >
  {
  private:
    ImageT m_image;
  public:

    typedef typename ImageT::pixel_type pixel_type;
    typedef TransposePixelAccessor<typename ImageT::pixel_accessor> pixel_accessor;

    TransposeView( ImageT const& image ) : m_image(image) {}

    inline unsigned cols() const { return m_image.rows(); }
    inline unsigned rows() const { return m_image.cols(); }
    inline unsigned planes() const { return m_image.planes(); }

    inline pixel_accessor origin() const { return m_image.origin(); }

    typename boost::mpl::if_< IsReferenceable<ImageT>, pixel_type&, pixel_type >::type
    inline operator()( int i, int j ) const { return m_image(j,i); }

    typename boost::mpl::if_< IsReferenceable<ImageT>, pixel_type&, pixel_type >::type
    inline operator()( int i, int j, int p ) const { return m_image(j,i,p); }

    template <class ViewT>
    TransposeView& operator=( ImageViewBase<ViewT> const& view ) {
      view.impl().rasterize( *this );
      return *this;
    }

    ImageT const& child() const {
      return m_image;
    }

    /// \cond INTERNAL
    typedef TransposeView<typename ImageT::prerasterize_type> prerasterize_type;
    inline prerasterize_type prerasterize() const { return prerasterize_type( m_image.prerasterize() ); }
    template <class DestT> inline void rasterize( DestT const& dest ) const { vw::rasterize( prerasterize(), dest ); }
    /// \endcond
  };

  /// \cond INTERNAL
  // Type Traits
  template <class ImageT>
  struct IsReferenceable<TransposeView<ImageT> > : public IsReferenceable<ImageT> {};

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
  template <class ImageAccT>
  class Rotate180PixelAccessor
  {
  private:
    ImageAccT m_image;
  public:
    typedef typename ImageAccT::pixel_type pixel_type;
    Rotate180PixelAccessor( ImageAccT const& image ) : m_image(image) {}

    inline Rotate180PixelAccessor& next_col() { m_image.prev_col(); return *this; }
    inline Rotate180PixelAccessor& prev_col() { m_image.next_col(); return *this; }
    inline Rotate180PixelAccessor& next_row() { m_image.prev_row(); return *this; }
    inline Rotate180PixelAccessor& prev_row() { m_image.next_row(); return *this; }
    inline Rotate180PixelAccessor& next_plane() { m_image.prev_plane(); return *this; }
    inline Rotate180PixelAccessor& prev_plane() { m_image.next_plane(); return *this; }
    inline Rotate180PixelAccessor& advance( ptrdiff_t di, ptrdiff_t dj, ptrdiff_t dp=0 ) { m_image.advance(-di,-dj,dp); return *this; }

    typename boost::mpl::if_< IsReferenceable<ImageAccT>, pixel_type&, pixel_type >::type
    inline operator*() const { return *m_image; }
  };

  /// \cond INTERNAL
  template <class ImageAccT>
  struct IsReferenceable<Rotate180PixelAccessor<ImageAccT> > : public IsReferenceable<ImageAccT> {};
  /// \endcond

  // Image View Class
  template <class ImageT>
  class Rotate180View : public ImageViewBase<Rotate180View<ImageT> >
  {
  private:
    ImageT m_image;
  public:

    typedef typename ImageT::pixel_type pixel_type;
    typedef Rotate180PixelAccessor<typename ImageT::pixel_accessor> pixel_accessor;

    Rotate180View( ImageT const& image ) : m_image(image) {}

    inline unsigned cols() const { return m_image.cols(); }
    inline unsigned rows() const { return m_image.rows(); }
    inline unsigned planes() const { return m_image.planes(); }

    inline pixel_accessor origin() const { return m_image.origin().advance(cols()-1,rows()-1); }

    typename boost::mpl::if_< IsReferenceable<ImageT>, pixel_type&, pixel_type >::type
    inline operator()( int i, int j ) const { return m_image(cols()-1-i,rows()-1-j); }

    typename boost::mpl::if_< IsReferenceable<ImageT>, pixel_type&, pixel_type >::type
    inline operator()( int i, int j, int p ) const { return m_image(cols()-1-i,rows()-1-j,p); }

    template <class ViewT>
    Rotate180View& operator=( ImageViewBase<ViewT> const& view ) {
      view.impl().rasterize( *this );
      return *this;
    }

    ImageT const& child() const {
      return m_image;
    }

    /// \cond INTERNAL
    typedef Rotate180View<typename ImageT::prerasterize_type> prerasterize_type;
    inline prerasterize_type prerasterize() const { return prerasterize_type( m_image.prerasterize() ); }
    template <class DestT> inline void rasterize( DestT const& dest ) const { vw::rasterize( prerasterize(), dest ); }
    /// \endcond
  };

  /// \cond INTERNAL
  // Type traits
  template <class ImageT>
  struct IsReferenceable<Rotate180View<ImageT> > : public IsReferenceable<ImageT> {};

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
  template <class ImageAccT>
  class Rotate90CWPixelAccessor
  {
  private:
    ImageAccT m_acc;
  public:
    typedef typename ImageAccT::pixel_type pixel_type;
    Rotate90CWPixelAccessor( ImageAccT const& iter ) : m_acc(iter) {}
    inline Rotate90CWPixelAccessor& next_col() { m_acc.prev_row(); return *this; }
    inline Rotate90CWPixelAccessor& prev_col() { m_acc.next_row(); return *this; }
    inline Rotate90CWPixelAccessor& next_row() { m_acc.next_col(); return *this; }
    inline Rotate90CWPixelAccessor& prev_row() { m_acc.prev_col(); return *this; }
    inline Rotate90CWPixelAccessor& next_plane() { m_acc.next_plane(); return *this; }
    inline Rotate90CWPixelAccessor& prev_plane() { m_acc.prev_plane(); return *this; }
    inline Rotate90CWPixelAccessor& advance( ptrdiff_t di, ptrdiff_t dj, ptrdiff_t dp=0 ) { m_acc.advance(dj,-di,dp); return *this; }

    typename boost::mpl::if_< IsReferenceable<ImageAccT>, pixel_type&, pixel_type >::type
    inline operator*() const { return *m_acc; }
  };

  /// \cond INTERNAL
  // Accessor type traits
  template <class ImageAccT>
  struct IsReferenceable<Rotate90CWPixelAccessor<ImageAccT> > : public IsReferenceable<ImageAccT> {};
  /// \endcond

  // Class definition
  template <class ImageT>
  class Rotate90CWView : public ImageViewBase<Rotate90CWView<ImageT> >
  {
  private:
    ImageT m_image;
  public:

    typedef typename ImageT::pixel_type pixel_type;
    typedef Rotate90CWPixelAccessor<typename ImageT::pixel_accessor> pixel_accessor;

    Rotate90CWView( ImageT const& image ) : m_image(image) {}

    inline unsigned cols() const { return m_image.rows(); }
    inline unsigned rows() const { return m_image.cols(); }
    inline unsigned planes() const { return m_image.planes(); }

    inline pixel_accessor origin() const { return m_image.origin().advance(0,cols()-1); }

    typename boost::mpl::if_< IsReferenceable<ImageT>, pixel_type&, pixel_type >::type
    inline operator()( int i, int j ) const { return m_image(j,cols()-1-i); }

    typename boost::mpl::if_< IsReferenceable<ImageT>, pixel_type&, pixel_type >::type
    inline operator()( int i, int j, int p ) const { return m_image(j,cols()-1-i,p); }

    template <class ViewT>
    Rotate90CWView& operator=( ImageViewBase<ViewT> const& view ) {
      view.impl().rasterize( *this );
      return *this;
    }

    ImageT const& child() const {
      return m_image;
    }

    /// \cond INTERNAL
    typedef Rotate90CWView<typename ImageT::prerasterize_type> prerasterize_type;
    inline prerasterize_type prerasterize() const { return prerasterize_type( m_image.prerasterize() ); }
    template <class DestT> inline void rasterize( DestT const& dest ) const { vw::rasterize( prerasterize(), dest ); }
    /// \endcond
  };

  /// \cond INTERNAL
  // Type tTraits
  template <class ImageT>
  struct IsReferenceable<Rotate90CWView<ImageT> > : public IsReferenceable<ImageT> {};

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
  template <class ImageAccT>
  class Rotate90CCWPixelAccessor
  {
  private:
    ImageAccT m_acc;
  public:
    typedef typename ImageAccT::pixel_type pixel_type;
    Rotate90CCWPixelAccessor( ImageAccT const& iter ) : m_acc(iter) {}
    inline Rotate90CCWPixelAccessor& next_col() { m_acc.next_row(); return *this; }
    inline Rotate90CCWPixelAccessor& prev_col() { m_acc.prev_row(); return *this; }
    inline Rotate90CCWPixelAccessor& next_row() { m_acc.prev_col(); return *this; }
    inline Rotate90CCWPixelAccessor& prev_row() { m_acc.next_col(); return *this; }
    inline Rotate90CCWPixelAccessor& next_plane() { m_acc.next_plane(); return *this; }
    inline Rotate90CCWPixelAccessor& prev_plane() { m_acc.prev_plane(); return *this; }
    inline Rotate90CCWPixelAccessor& advance( ptrdiff_t di, ptrdiff_t dj, ptrdiff_t dp=0 ) { m_acc.advance(-dj,di,dp); return *this; }

    typename boost::mpl::if_< IsReferenceable<ImageAccT>, pixel_type&, pixel_type >::type
    inline operator*() const { return *m_acc; }
  };

  /// \cond INTERNAL
  // Accessor type traits
  template <class ImageAccT>
  struct IsReferenceable<Rotate90CCWPixelAccessor<ImageAccT> > : public IsReferenceable<ImageAccT> {};
  /// \endcond

  // Class definition
  template <class ImageT>
  class Rotate90CCWView : public ImageViewBase<Rotate90CCWView<ImageT> >
  {
  private:
    ImageT m_image;
  public:

    typedef typename ImageT::pixel_type pixel_type;
    typedef Rotate90CCWPixelAccessor<typename ImageT::pixel_accessor> pixel_accessor;

    Rotate90CCWView( ImageT const& image ) : m_image(image) {}

    inline unsigned cols() const { return m_image.rows(); }
    inline unsigned rows() const { return m_image.cols(); }
    inline unsigned planes() const { return m_image.planes(); }

    inline pixel_accessor origin() const { return m_image.origin().advance(rows()-1,0); }

    typename boost::mpl::if_< IsReferenceable<ImageT>, pixel_type&, pixel_type >::type
    inline operator()( int i, int j ) const { return m_image(rows()-1-j,i); }

    typename boost::mpl::if_< IsReferenceable<ImageT>, pixel_type&, pixel_type >::type
    inline operator()( int i, int j, int p ) const { return m_image(rows()-1-j,i,p); }

    template <class ViewT>
    Rotate90CCWView& operator=( ImageViewBase<ViewT> const& view ) {
      view.impl().rasterize( *this );
      return *this;
    }

    ImageT const& child() const {
      return m_image;
    }

    /// \cond INTERNAL
    typedef Rotate90CCWView<typename ImageT::prerasterize_type> prerasterize_type;
    inline prerasterize_type prerasterize() const { return prerasterize_type( m_image.prerasterize() ); }
    template <class DestT> inline void rasterize( DestT const& dest ) const { vw::rasterize( prerasterize(), dest ); }
    /// \endcond
  };

  /// \cond INTERNAL
  // Type traits
  template <class ImageT>
  struct IsReferenceable<Rotate90CCWView<ImageT> > : public IsReferenceable<ImageT> {};

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
  template <class ImageAccT>
  class FlipVerticalPixelAccessor
  {
  private:
    ImageAccT m_acc;
  public:
    typedef typename ImageAccT::pixel_type pixel_type;
    FlipVerticalPixelAccessor( ImageAccT const& iter ) : m_acc(iter) {}
    inline FlipVerticalPixelAccessor& next_col() { m_acc.next_col(); return *this; }
    inline FlipVerticalPixelAccessor& prev_col() { m_acc.prev_col(); return *this; }
    inline FlipVerticalPixelAccessor& next_row() { m_acc.prev_row(); return *this; }
    inline FlipVerticalPixelAccessor& prev_row() { m_acc.next_row(); return *this; }
    inline FlipVerticalPixelAccessor& next_plane() { m_acc.next_plane(); return *this; }
    inline FlipVerticalPixelAccessor& prev_plane() { m_acc.prev_plane(); return *this; }
    inline FlipVerticalPixelAccessor& advance( ptrdiff_t di, ptrdiff_t dj, ptrdiff_t dp=0 ) { m_acc.advance(di,-dj,dp); return *this; }

    typename boost::mpl::if_< IsReferenceable<ImageAccT>, pixel_type&, pixel_type >::type
    inline operator*() const { return *m_acc; }
  };

  /// \cond INTERNAL
  // Accessor type traits
  template <class ImageAccT>
  struct IsReferenceable<FlipVerticalPixelAccessor<ImageAccT> > : public IsReferenceable<ImageAccT> {};
  /// \endcond

  // Class definition
  template <class ImageT>
  class FlipVerticalView : public ImageViewBase<FlipVerticalView<ImageT> >
  {
  private:
    ImageT m_image;
  public:

    typedef typename ImageT::pixel_type pixel_type;
    typedef FlipVerticalPixelAccessor<typename ImageT::pixel_accessor> pixel_accessor;

    FlipVerticalView( ImageT const& image ) : m_image(image) {}

    inline unsigned cols() const { return m_image.cols(); }
    inline unsigned rows() const { return m_image.rows(); }
    inline unsigned planes() const { return m_image.planes(); }

    inline pixel_accessor origin() const { return m_image.origin().advance(0,rows()-1); }

    typename boost::mpl::if_< IsReferenceable<ImageT>, pixel_type&, pixel_type >::type
    inline operator()( int i, int j ) const { return m_image(i,rows()-1-j); }

    typename boost::mpl::if_< IsReferenceable<ImageT>, pixel_type&, pixel_type >::type
    inline operator()( int i, int j, int p ) const { return m_image(i,rows()-1-j,p); }

    ImageT const& child() const {
      return m_image;
    }

    template <class ViewT>
    FlipVerticalView& operator=( ImageViewBase<ViewT> const& view ) {
      view.impl().rasterize( *this );
      return *this;
    }

    /// \cond INTERNAL
    typedef FlipVerticalView<typename ImageT::prerasterize_type> prerasterize_type;
    inline prerasterize_type prerasterize() const { return prerasterize_type( m_image.prerasterize() ); }
    template <class DestT> inline void rasterize( DestT const& dest ) const { vw::rasterize( prerasterize(), dest ); }
    /// \endcond
  };

  /// \cond INTERNAL
  // Type traits
  template <class ImageT>
  struct IsReferenceable<FlipVerticalView<ImageT> > : public IsReferenceable<ImageT> {};

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
  template <class ImageAccT>
  class FlipHorizontalPixelAccessor
  {
  private:
    ImageAccT m_acc;
  public:
    typedef typename ImageAccT::pixel_type pixel_type;
    FlipHorizontalPixelAccessor( ImageAccT const& iter ) : m_acc(iter) {}
    inline FlipHorizontalPixelAccessor& next_col() { m_acc.prev_col(); return *this; }
    inline FlipHorizontalPixelAccessor& prev_col() { m_acc.next_col(); return *this; }
    inline FlipHorizontalPixelAccessor& next_row() { m_acc.next_row(); return *this; }
    inline FlipHorizontalPixelAccessor& prev_row() { m_acc.prev_row(); return *this; }
    inline FlipHorizontalPixelAccessor& next_plane() { m_acc.next_plane(); return *this; }
    inline FlipHorizontalPixelAccessor& prev_plane() { m_acc.prev_plane(); return *this; }
    inline FlipHorizontalPixelAccessor& advance( ptrdiff_t di, ptrdiff_t dj, ptrdiff_t dp=0 ) { m_acc.advance(-di,dj,dp); return *this; }

    typename boost::mpl::if_< IsReferenceable<ImageAccT>, pixel_type&, pixel_type >::type
    inline operator*() const { return *m_acc; }
  };

  /// \cond INTERNAL
  // Accessor type traits
  template <class ImageAccT>
  struct IsReferenceable<FlipHorizontalPixelAccessor<ImageAccT> > : public IsReferenceable<ImageAccT> {};
  /// \endcond

  // Class definition
  template <class ImageT>
  class FlipHorizontalView : public ImageViewBase<FlipHorizontalView<ImageT> >
  {
  private:
    ImageT m_image;
  public:

    typedef typename ImageT::pixel_type pixel_type;
    typedef FlipHorizontalPixelAccessor<typename ImageT::pixel_accessor> pixel_accessor;

    FlipHorizontalView( ImageT const& image ) : m_image(image) {}

    inline unsigned cols() const { return m_image.cols(); }
    inline unsigned rows() const { return m_image.rows(); }
    inline unsigned planes() const { return m_image.planes(); }

    inline pixel_accessor origin() const { return m_image.origin().advance(cols()-1,0); }

    typename boost::mpl::if_< IsReferenceable<ImageT>, pixel_type&, pixel_type >::type
    inline operator()( int i, int j ) const { return m_image(cols()-1-i,j); }

    typename boost::mpl::if_< IsReferenceable<ImageT>, pixel_type&, pixel_type >::type
    inline operator()( int i, int j, int p ) const { return m_image(cols()-1-i,j,p); }

    template <class ViewT>
    FlipHorizontalView& operator=( ImageViewBase<ViewT> const& view ) {
      view.impl().rasterize( *this );
      return *this;
    }

    ImageT const& child() const {
      return m_image;
    }

    /// \cond INTERNAL
    typedef FlipHorizontalView<typename ImageT::prerasterize_type> prerasterize_type;
    inline prerasterize_type prerasterize() const { return prerasterize_type( m_image.prerasterize() ); }
    template <class DestT> inline void rasterize( DestT const& dest ) const { vw::rasterize( prerasterize(), dest ); }
    /// \endcond
  };

  /// \cond INTERNAL
  // Type traits
  template <class ImageT>
  struct IsReferenceable<FlipHorizontalView<ImageT> > : public IsReferenceable<ImageT> {};

  template <class ImageT>
  struct IsMultiplyAccessible<FlipHorizontalView<ImageT> > : public IsMultiplyAccessible<ImageT> {};
  /// \endcond

  /// Flip an image horizontally.
  template <class ImageT>
  FlipHorizontalView<ImageT> flip_horizontal( ImageViewBase<ImageT> const& v ) {
    return FlipHorizontalView<ImageT>( v.impl() );
  }


  // *******************************************************************
  // Crop
  // *******************************************************************

  // Class definition
  template <class ImageT>
  class CropView : public ImageViewBase< CropView<ImageT> >
  {
  private:
    ImageT m_image;
    unsigned m_ci, m_cj, m_di, m_dj;
  public:

    typedef typename ImageT::pixel_type pixel_type;
    typedef typename ImageT::pixel_accessor pixel_accessor;

    CropView( ImageT const& image, unsigned const &upper_left_i, unsigned const &upper_left_j, unsigned const &width, unsigned const &height ) : 
      m_image(image), m_ci(upper_left_i), m_cj(upper_left_j), m_di(width), m_dj(height) {      
      VW_ASSERT( m_ci + m_di <= m_image.cols() && m_cj + m_dj <= m_image.rows(),
                 ArgumentErr() << "CropView: Crop dimensions exceed image boundary." );
    }

    inline unsigned cols() const { return m_di; }
    inline unsigned rows() const { return m_dj; }
    inline unsigned planes() const { return m_image.planes(); }

    inline pixel_accessor origin() const { return m_image.origin().advance(m_ci, m_cj); }

    typename boost::mpl::if_< IsReferenceable<ImageT>, pixel_type&, pixel_type >::type
    inline operator()( int i, int j ) const { return m_image(m_ci + i, m_cj + j); }

    typename boost::mpl::if_< IsReferenceable<ImageT>, pixel_type&, pixel_type >::type
    inline operator()( int i, int j, int p ) const { return m_image(m_ci + i, m_cj + j, p); }

    CropView& operator=( CropView const& view ) {
      view.rasterize( *this );
      return *this;
    }

    template <class ViewT>
    CropView& operator=( ImageViewBase<ViewT> const& view ) {
      view.impl().rasterize( *this );
      return *this;
    }

    ImageT const& child() const {
      return m_image;
    }

    /// \cond INTERNAL
    typedef CropView<typename ImageT::prerasterize_type> prerasterize_type;
    inline prerasterize_type prerasterize() const { return prerasterize_type( m_image.prerasterize(), m_ci, m_cj, m_di, m_dj ); }
    template <class DestT> inline void rasterize( DestT const& dest ) const { vw::rasterize( prerasterize(), dest ); }
    /// \endcond
  };

  /// \cond INTERNAL
  // Type traits
  template <class ImageT>
  struct IsReferenceable<CropView<ImageT> >  : public IsReferenceable<ImageT> {}; 

  template <class ImageT>
  struct IsMultiplyAccessible<CropView<ImageT> > : public IsMultiplyAccessible<ImageT> {};
  /// \endcond

  /// Crop an image.
  template <class ImageT>
  inline CropView<ImageT> crop( ImageT const& v, int upper_left_x, int upper_left_y, int width, int height ) {
    return CropView<ImageT>( v, upper_left_x, upper_left_y, width, height );
  }

  // /// Crop an image.
  // template <class ImageT>
  // inline CropView<ImageT> crop( ImageT const& v, BBox<int,2> const& bbox ) {
  //   return CropView<ImageT>( v, bbox.min().x(), bbox.min().y(), bbox.width(), bbox.height() );
  // }


  // *******************************************************************
  // Subsample
  // *******************************************************************

  // Specialized image accessor
  template <class ImageAccT>
  class SubsamplePixelAccessor
  {
  private:
    ImageAccT m_acc;
    unsigned m_xdelta, m_ydelta;
  public:
    typedef typename ImageAccT::pixel_type pixel_type;
    SubsamplePixelAccessor( ImageAccT const& iter , unsigned xdelta, unsigned ydelta) : m_acc(iter), m_xdelta(xdelta), m_ydelta(ydelta) {}

    inline SubsamplePixelAccessor& next_col() { m_acc.advance(  m_xdelta, 0 ); return *this; }
    inline SubsamplePixelAccessor& prev_col() { m_acc.advance( -m_xdelta, 0 ); return *this; }
    inline SubsamplePixelAccessor& next_row() { m_acc.advance( 0,  m_ydelta ); return *this; }
    inline SubsamplePixelAccessor& prev_row() { m_acc.advance( 0, -m_ydelta ); return *this; }
    inline SubsamplePixelAccessor& next_plane() { m_acc.next_plane(); return *this; }
    inline SubsamplePixelAccessor& prev_plane() { m_acc.prev_plane(); return *this; }
    inline SubsamplePixelAccessor& advance( ptrdiff_t di, ptrdiff_t dj, ptrdiff_t dp=0 ) { m_acc.advance((ptrdiff_t)m_xdelta*di,(ptrdiff_t)m_ydelta*dj,dp); return *this; }

    typename boost::mpl::if_< IsReferenceable<ImageAccT>, pixel_type&, pixel_type >::type
    inline operator*() const { return *m_acc; }
  };

  /// \cond INTERNAL
  template <class ImageAccT>
  struct IsReferenceable<SubsamplePixelAccessor<ImageAccT> > : public IsReferenceable<ImageAccT> {};
  /// \endcond

  // Class definition
  template <class ImageT>
  class SubsampleView : public ImageViewBase<SubsampleView<ImageT> >
  {
  private:
    ImageT m_image;
    unsigned m_xdelta, m_ydelta;
  public:

    typedef typename ImageT::pixel_type pixel_type;
    typedef SubsamplePixelAccessor<typename ImageT::pixel_accessor> pixel_accessor;

    SubsampleView( ImageT const& image, unsigned subsampling_factor ) : m_image(image), m_xdelta(subsampling_factor), m_ydelta(subsampling_factor) {}
    SubsampleView( ImageT const& image, unsigned xfactor, unsigned yfactor ) : m_image(image), m_xdelta(xfactor), m_ydelta(yfactor) {}

    inline unsigned cols() const { return 1 + (m_image.cols()-1)/m_xdelta; }
    inline unsigned rows() const { return 1 + (m_image.rows()-1)/m_ydelta; }
    inline unsigned planes() const { return m_image.planes(); }

    inline pixel_accessor origin() const { return pixel_accessor(m_image.origin(), m_xdelta, m_ydelta); }

    typename boost::mpl::if_< IsReferenceable<ImageT>, pixel_type&, pixel_type >::type
    inline operator()( int i, int j ) const { return m_image(m_xdelta*i,m_ydelta*j); }

    typename boost::mpl::if_< IsReferenceable<ImageT>, pixel_type&, pixel_type >::type
    inline operator()( int i, int j, int p ) const { return m_image(m_xdelta*i,m_ydelta*j,p); }

    ImageT const& child() const {
      return m_image;
    }

    /// \cond INTERNAL
    typedef SubsampleView<typename ImageT::prerasterize_type> prerasterize_type;
    inline prerasterize_type prerasterize() const { return prerasterize_type( m_image.prerasterize(), m_xdelta, m_ydelta ); }
    template <class DestT> inline void rasterize( DestT const& dest ) const { vw::rasterize( prerasterize(), dest ); }
    /// \endcond
  };

  /// \cond INTERNAL
  // Type Traits
  template <class ImageT>
  struct IsReferenceable<SubsampleView<ImageT> > : public IsReferenceable<ImageT> {};

  template <class ImageT>
  struct IsMultiplyAccessible<SubsampleView<ImageT> > : public IsMultiplyAccessible<ImageT> {};
  /// \endcond

  /// Subsample an image by an integer factor.
  template <class ImageT>
  inline SubsampleView<ImageT> subsample( ImageT const& v, unsigned subsampling_factor ) {
    return SubsampleView<ImageT>( v, subsampling_factor );
  }

  /// Subsample an image by integer factors in x and y.
  template <class ImageT>
  inline SubsampleView<ImageT> subsample( ImageT const& v, unsigned xfactor, unsigned yfactor ) {
    return SubsampleView<ImageT>( v, xfactor, yfactor );
  }


  // *******************************************************************
  // SelectPlane 
  // *******************************************************************

  /// Return a single plane from a multi-plane image
  /// \see vw::select_plane
  template <class ImageT>
  class SelectPlaneView : public ImageViewBase<SelectPlaneView<ImageT> >
  {
  private:
    ImageT m_image;
    unsigned m_plane;
  public:

    typedef typename ImageT::pixel_type pixel_type;
    typedef typename ImageT::pixel_accessor pixel_accessor;

    SelectPlaneView( ImageT const& image, unsigned plane ) : m_image(image), m_plane(plane) {}

    inline unsigned cols() const { return m_image.cols(); }
    inline unsigned rows() const { return m_image.rows(); }
    inline unsigned planes() const { return 1; }

    inline pixel_accessor origin() const { return m_image.origin().advance(0,0,m_plane); }

    typename boost::mpl::if_< IsReferenceable<ImageT>, pixel_type&, pixel_type >::type
    inline operator()( int i, int j ) const { return m_image(i,j,m_plane); }

    typename boost::mpl::if_< IsReferenceable<ImageT>, pixel_type&, pixel_type >::type
    inline operator()( int i, int j, int p ) const { return m_image(i,j,m_plane+p); }

    template <class ViewT>
    SelectPlaneView& operator=( ImageViewBase<ViewT> const& view ) {
      view.impl().rasterize( *this );
      return *this;
    }

    /// \cond INTERNAL
    typedef SelectPlaneView<typename ImageT::prerasterize_type> prerasterize_type;
    inline prerasterize_type prerasterize() const { return prerasterize_type( m_image.prerasterize(), m_plane ); }
    template <class DestT> inline void rasterize( DestT const& dest ) const { vw::rasterize( prerasterize(), dest ); }
    /// \endcond
  };

  /// \cond INTERNAL
  // View type Traits
  template <class ImageT>
  struct IsReferenceable<SelectPlaneView<ImageT> > : public IsReferenceable<ImageT> {};

  template <class ImageT>
  struct IsMultiplyAccessible<SelectPlaneView<ImageT> > : public IsMultiplyAccessible<ImageT> {};
  /// \endcond

  /// Extracts a single plane of a multi-plane image.  This function
  /// returns a writeable view of a single plane of a multi-plane
  /// image.  \see vw::SelectPlaneView
  template <class ImageT>
  SelectPlaneView<ImageT> select_plane( ImageViewBase<ImageT> const& v, unsigned plane ) {
    return SelectPlaneView<ImageT>( v.impl(), plane );
  }


  // *******************************************************************
  // SelectChannel
  // *******************************************************************

  /// A single-channel pixel accessor adaptor.
  ///
  /// This is a special wrapper pixel accessor type, used by 
  /// \ref vw::SelectChannelView, that refers to a single channel 
  /// in a multi-channel image.
  template <class ImageAccT>
  class SelectChannelAccessor
  {
  private:
    ImageAccT m_acc;
    unsigned m_channel;
  public:
    typedef typename CompoundChannelType<typename ImageAccT::pixel_type>::type base_pixel_type;
    typedef typename boost::mpl::if_<boost::is_const<typename ImageAccT::pixel_type>,typename boost::add_const<base_pixel_type>::type,base_pixel_type>::type pixel_type;
    typedef typename boost::mpl::if_<IsReferenceable<ImageAccT>,pixel_type&,pixel_type>::type result_type;

    SelectChannelAccessor( ImageAccT const& iter, unsigned channel ) : m_acc(iter), m_channel(channel) {}
    inline SelectChannelAccessor& next_col() { m_acc.next_col(); return *this; }
    inline SelectChannelAccessor& prev_col() { m_acc.prev_col(); return *this; }
    inline SelectChannelAccessor& next_row() { m_acc.next_row(); return *this; }
    inline SelectChannelAccessor& prev_row() { m_acc.prev_row(); return *this; }
    inline SelectChannelAccessor& next_plane() { m_acc.next_plane(); return *this; }
    inline SelectChannelAccessor& prev_plane() { m_acc.prev_plane(); return *this; }
    inline SelectChannelAccessor& advance( ptrdiff_t di, ptrdiff_t dj, ptrdiff_t dp=0 ) { m_acc.advance(di,dj,dp); return *this; }

    inline result_type operator*() const { return compound_select_channel<result_type>(*m_acc,m_channel); }
  };

  /// \cond INTERNAL
  // Accessor type traits
  template <class ImageAccT>
  struct IsReferenceable<SelectChannelAccessor<ImageAccT> > : public IsReferenceable<ImageAccT> {};
  /// \endcond

  /// Returns a single channel of a multi-channel image.
  /// \see vw::select_channel
  template <class ImageT>
  class SelectChannelView : public ImageViewBase<SelectChannelView<ImageT> >
  {
  private:
    ImageT m_image;
    unsigned m_channel;
  public:

    typedef typename CompoundChannelType<typename ImageT::pixel_type>::type base_pixel_type;
    typedef typename boost::mpl::if_<boost::is_const<typename ImageT::pixel_type>,typename boost::add_const<base_pixel_type>::type,base_pixel_type>::type pixel_type;
    typedef typename boost::mpl::if_<IsReferenceable<ImageT>,pixel_type&,pixel_type>::type result_type;

    typedef SelectChannelAccessor<typename ImageT::pixel_accessor> pixel_accessor;

    SelectChannelView( ImageT const& image, unsigned channel ) : m_image(image), m_channel(channel) {}

    inline unsigned cols() const { return m_image.cols(); }
    inline unsigned rows() const { return m_image.rows(); }
    inline unsigned planes() const { return m_image.planes(); }

    inline pixel_accessor origin() const { return pixel_accessor(m_image.origin(), m_channel); }

    inline result_type operator()( int i, int j, int p=0 ) const {
      return compound_select_channel<result_type>(m_image(i,j,p),m_channel);
    }

    template <class ViewT>
    SelectChannelView& operator=( ImageViewBase<ViewT> const& view ) {
      view.impl().rasterize( *this );
      return *this;
    }

    ImageT const& child() const {
      return m_image;
    }

    /// \cond INTERNAL
    typedef SelectChannelView<typename ImageT::prerasterize_type> prerasterize_type;
    inline prerasterize_type prerasterize() const { return prerasterize_type( m_image.prerasterize(), m_channel ); }
    template <class DestT> inline void rasterize( DestT const& dest ) const { vw::rasterize( prerasterize(), dest ); }
    /// \endcond
  };

  /// \cond INTERNAL
  // View type traits
  template <class ImageT>
  struct IsReferenceable<SelectChannelView<ImageT> > : public IsReferenceable<ImageT> {};

  template <class ImageT>
  struct IsMultiplyAccessible<SelectChannelView<ImageT> > : public IsMultiplyAccessible<ImageT> {};
  /// \endcond

  /// Extracts a single channel of a multi-channel image.
  /// This function returns a writeable view of a single channel 
  /// of a multi-channel image.
  /// \see vw::SelectChannelView
  template <class ImageT>
  SelectChannelView<ImageT> select_channel( ImageViewBase<ImageT> const& v, unsigned channel ) {
    return SelectChannelView<ImageT>( v.impl(), channel );
  }


  // *******************************************************************
  // ChannelsToPlanes
  // *******************************************************************

  /// A channels-to-planes pixel accessor adaptor.
  ///
  /// This is a special wrapper pixel accessor type, used by 
  /// \ref vw::ChannelsToPlanesView, that treats the channels in a
  /// multi-channel image as planes.
  template <class ImageAccT>
  class ChannelsToPlanesAccessor
  {
  private:
    ImageAccT m_acc;
    unsigned m_channel;
  public:
    typedef typename CompoundChannelType<typename ImageAccT::pixel_type>::type base_pixel_type;
    typedef typename boost::mpl::if_<boost::is_const<typename ImageAccT::pixel_type>,typename boost::add_const<base_pixel_type>::type,base_pixel_type>::type pixel_type;
    typedef typename boost::mpl::if_<IsReferenceable<ImageAccT>,pixel_type&,pixel_type>::type result_type;

    ChannelsToPlanesAccessor( ImageAccT const& iter ) : m_acc(iter), m_channel(0) {}
    inline ChannelsToPlanesAccessor& next_col() { m_acc.next_col(); return *this; }
    inline ChannelsToPlanesAccessor& prev_col() { m_acc.prev_col(); return *this; }
    inline ChannelsToPlanesAccessor& next_row() { m_acc.next_row(); return *this; }
    inline ChannelsToPlanesAccessor& prev_row() { m_acc.prev_row(); return *this; }
    inline ChannelsToPlanesAccessor& next_plane() { ++m_channel; return *this; }
    inline ChannelsToPlanesAccessor& prev_plane() { --m_channel; return *this; }
    inline ChannelsToPlanesAccessor& advance( ptrdiff_t di, ptrdiff_t dj, ptrdiff_t dp=0 ) { m_acc.advance(di,dj); m_channel+=dp; return *this; }

    inline result_type operator*() const { return compound_select_channel<result_type>(*m_acc,m_channel); }
  };

  /// \cond INTERNAL
  // Accessor type traits
  template <class ImageAccT>
  struct IsReferenceable<ChannelsToPlanesAccessor<ImageAccT> > : public IsReferenceable<ImageAccT> {};
  /// \endcond

  /// A view that turns a one plane, multi-channel view into a mulit-plane, one channel view.
  /// \see vw::channels_to_planes
  template <class ImageT>
  class ChannelsToPlanesView : public ImageViewBase<ChannelsToPlanesView<ImageT> >
  {
  private:
    ImageT m_image;
  public:

    typedef typename CompoundChannelType<typename ImageT::pixel_type>::type base_pixel_type;
    typedef typename boost::mpl::if_<boost::is_const<typename ImageT::pixel_type>,typename boost::add_const<base_pixel_type>::type,base_pixel_type>::type pixel_type;
    typedef typename boost::mpl::if_< IsReferenceable<ImageT>, pixel_type&, pixel_type >::type result_type;

    typedef typename boost::mpl::if_< IsCompound<typename ImageT::pixel_type>, 
                                      ChannelsToPlanesAccessor<typename ImageT::pixel_accessor>,
                                      typename ImageT::pixel_accessor >::type pixel_accessor;
    
    ChannelsToPlanesView( ImageT const& image ) : m_image(image) {
      VW_ASSERT( m_image.planes()==1 || m_image.channels()==1, ArgumentErr() <<
                 "ChannelsToPlanesView: The image must be multi-channel, single plane or single-channel, multi-plane.");
    }
    
    inline unsigned cols() const { return m_image.cols(); }
    inline unsigned rows() const { return m_image.rows(); }
    inline unsigned planes() const { 
      if( IsCompound<typename ImageT::pixel_type>::value ) {
        return m_image.channels();
      } else {
        return m_image.planes();
      }
    }

    inline pixel_accessor origin() const { return pixel_accessor(m_image.origin()); }

    inline result_type operator()( int i, int j, int p=0 ) const {
      return compound_select_channel<result_type>(m_image(i,j),p);
    }

    template <class ViewT>
    ChannelsToPlanesView& operator=( ImageViewBase<ViewT> const& view ) {
      view.impl().rasterize( *this );
      return *this;
    }

    ImageT const& child() const {
      return m_image;
    }

    /// \cond INTERNAL
    typedef ChannelsToPlanesView<typename ImageT::prerasterize_type> prerasterize_type;
    inline prerasterize_type prerasterize() const { return prerasterize_type( m_image.prerasterize() ); }
    template <class DestT> inline void rasterize( DestT const& dest ) const { vw::rasterize( prerasterize(), dest ); }
    /// \endcond
  };

  /// \cond INTERNAL
  // View type Traits
  template <class ImageT>
  struct IsReferenceable<ChannelsToPlanesView<ImageT> > : public IsReferenceable<ImageT> {};

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
  // PlanesToChannels
  // *******************************************************************

  /// A view that turns a multi-plane, single-channel view into a
  /// one-plane, multi-channel view.
  /// \see vw::planes_to_channels
  template <class PixelT,class ImageT>
  class PlanesToChannelsView : public ImageViewBase<PlanesToChannelsView<PixelT,ImageT> >
  {
  private:
    ImageT m_image;
  public:

    typedef PixelT pixel_type;
    typedef PixelT result_type;

    typedef ProceduralPixelAccessor<PlanesToChannelsView> pixel_accessor;
    
    PlanesToChannelsView( ImageT const& image ) : m_image(image) {
      VW_ASSERT( m_image.channels()==1 && m_image.planes()==CompoundNumChannels<PixelT>::value, 
                 ArgumentErr() << "PlanesToChannelsView: The image must be multi-plane, single-channel.");
    }
    
    inline unsigned cols() const { return m_image.cols(); }
    inline unsigned rows() const { return m_image.rows(); }
    inline unsigned planes() const { return 1; }

    inline pixel_accessor origin() const { return pixel_accessor(*this); }

    inline result_type operator()( int i, int j, int p=0 ) const {
      result_type result;
      typedef typename CompoundChannelType<result_type>::type channel_type;
      for( int c=0; c<CompoundNumChannels<PixelT>::value; ++c )
        compound_select_channel<channel_type&>(result,c) = m_image(i,j,c);
      return result;
    }

    ImageT const& child() const {
      return m_image;
    }

    /// \cond INTERNAL
    typedef PlanesToChannelsView<PixelT, typename ImageT::prerasterize_type> prerasterize_type;
    inline prerasterize_type prerasterize() const { return prerasterize_type( m_image.prerasterize() ); }
    template <class DestT> inline void rasterize( DestT const& dest ) const { vw::rasterize( prerasterize(), dest ); }
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

} // namespace vw

#endif // __VW_IMAGE__MANIPULATION_H__
