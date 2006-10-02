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
/// - channel_cast() : casts the channels of an image while retaining the pixel format
///
#ifndef __VW_IMAGE_MANIPULATION_H__
#define __VW_IMAGE_MANIPULATION_H__

#include <boost/mpl/logical.hpp>

#include <vw/Core/CompoundTypes.h>
#include <vw/Image/ImageViewBase.h>
#include <vw/Image/PerPixelViews.h>
#include <vw/Math/BBox.h>

namespace vw {

  // *******************************************************************
  // copy()
  // *******************************************************************

  // Class definition
  template <class ImageT>
  class CopyView : public ImageViewBase<CopyView<ImageT> >
  {
  private:
    ImageView<typename ImageT::pixel_type> m_image;
  public:
    typedef typename ImageView<typename ImageT::pixel_type>::pixel_type pixel_type;
    typedef pixel_type const& result_type;
    typedef typename ImageView<typename ImageT::pixel_type>::pixel_accessor pixel_accessor;

    CopyView( ImageT const& image ) : m_image(image.cols(),image.rows(),image.planes()) {
      image.rasterize( m_image, BBox2i(0,0,image.cols(),image.rows()) );
    }

    inline unsigned cols() const { return m_image.cols(); }
    inline unsigned rows() const { return m_image.rows(); }
    inline unsigned planes() const { return m_image.planes(); }

    inline pixel_accessor origin() const { return m_image.origin(); }
    inline result_type operator()( int i, int j ) const { return m_image(i,j); }
    inline result_type operator()( int i, int j, int p ) const { return m_image(i,j,p); }

    typedef CopyView prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i bbox ) const { return *this; }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i bbox ) const { m_image.rasterize( dest, bbox ); }
  };

  template <class ImageT>
  struct IsMultiplyAccessible<CopyView<ImageT> > : public boost::true_type {};

  /// Make a (deep) copy of an image.
  template <class ImageT>
  CopyView<ImageT> copy( ImageViewBase<ImageT> const& v ) {
    return CopyView<ImageT>( v.impl() );
  }


  // *******************************************************************
  // Axis remapping views: transpose, flipping and rotation.
  // *******************************************************************

  template <int Mode>
  struct AxisTraits {};

  template <> struct AxisTraits<1> {
    template <class AccessT> static inline void iterate( AccessT& acc ) { acc.next_col(); }
    static inline int advance( int dc, int dr ) { return dc; }
    template <class ViewT> static inline int index( int c, int r, ViewT const& view ) { return c; }
    template <class ViewT> static inline int size( ViewT const& view ) { return view.cols(); }
  };

  template <> struct AxisTraits<-1> {
    template <class AccessT> static inline void iterate( AccessT& acc ) { acc.prev_col(); }
    static inline int advance( int dc, int dr ) { return -dc; }
    template <class ViewT> static inline int index( int c, int r, ViewT const& view ) { return view.cols()-1-c; }
    template <class ViewT> static inline int size( ViewT const& view ) { return view.cols(); }
  };

  template <> struct AxisTraits<2> {
    template <class AccessT> static inline void iterate( AccessT& acc ) { acc.next_row(); }
    static inline int advance( int dc, int dr ) { return dr; }
    template <class ViewT> static inline int index( int c, int r, ViewT const& view ) { return r; }
    template <class ViewT> static inline int size( ViewT const& view ) { return view.rows(); }
  };

  template <> struct AxisTraits<-2> {
    template <class AccessT> static inline void iterate( AccessT& acc ) { acc.prev_row(); }
    static inline int advance( int dc, int dr ) { return -dr; }
    template <class ViewT> static inline int index( int c, int r, ViewT const& view ) { return view.rows()-1-r; }
    template <class ViewT> static inline int size( ViewT const& view ) { return view.rows(); }
  };

  template <class ChildT, int FwdColMode, int FwdRowMode, int RevColMode, int RevRowMode>
  class RemapPixelAccessor {
    ChildT m_child;
  public:
    typedef typename ChildT::pixel_type pixel_type;
    typedef typename ChildT::result_type result_type;
    RemapPixelAccessor( ChildT const& child ) : m_child(child) {}
    
    inline RemapPixelAccessor& next_col() { AxisTraits<FwdColMode>::iterate(m_child); return *this; }
    inline RemapPixelAccessor& prev_col() { AxisTraits<-FwdColMode>::iterate(m_child); return *this; }
    inline RemapPixelAccessor& next_row() { AxisTraits<FwdRowMode>::iterate(m_child); return *this; }
    inline RemapPixelAccessor& prev_row() { AxisTraits<-FwdRowMode>::iterate(m_child); return *this; }
    inline RemapPixelAccessor& next_plane() { m_child.next_plane(); return *this; }
    inline RemapPixelAccessor& prev_plane() { m_child.prev_plane(); return *this; }
    inline RemapPixelAccessor& advance( ptrdiff_t di, ptrdiff_t dj, ptrdiff_t dp=0 ) { 
      m_child.advance( AxisTraits<RevColMode>::advance(di,dj), AxisTraits<RevRowMode>::advance(di,dj), dp );
      return *this;
    }

    inline result_type operator*() const { return *m_child; }
  };

  template <class ChildT, int FwdColMode, int FwdRowMode, int RevColMode, int RevRowMode>
  class RemapView : public ImageViewBase<RemapView<ChildT,FwdColMode,FwdRowMode,RevColMode,RevRowMode> > {
    ChildT m_child;
  public:

    typedef typename ChildT::pixel_type pixel_type;
    typedef typename ChildT::result_type result_type;
    typedef RemapPixelAccessor<typename ChildT::pixel_accessor,FwdColMode,FwdRowMode,RevColMode,RevRowMode> pixel_accessor;

    RemapView( ChildT const& child ) : m_child(child) {}

    inline unsigned cols() const { return AxisTraits<FwdColMode>::size( m_child ); }
    inline unsigned rows() const { return AxisTraits<FwdRowMode>::size( m_child ); }
    inline unsigned planes() const { return m_child.planes(); }

    inline pixel_accessor origin() const {
      return m_child.origin().advance( AxisTraits<RevColMode>::index(0,0,*this),
                                       AxisTraits<RevRowMode>::index(0,0,*this) );
    }

    inline result_type operator()( int c, int r, int p=0 ) const {
      return m_child( AxisTraits<RevColMode>::index(c,r,*this),
                      AxisTraits<RevRowMode>::index(c,r,*this), p );
    }

    template <class ViewT>
    RemapView& operator=( ImageViewBase<ViewT> const& view ) {
      view.impl().rasterize( *this, BBox2i(0,0,view.impl().cols(),view.impl().rows()) );
      return *this;
    }

    ChildT const& child() const {
      return m_child;
    }

    /// \cond INTERNAL
    typedef RemapView<typename ChildT::prerasterize_type,FwdColMode,FwdRowMode,RevColMode,RevRowMode> prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i bbox ) const { return prerasterize_type( m_child.prerasterize(bbox) ); }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i bbox ) const { vw::rasterize( prerasterize(bbox), dest, bbox ); }
    /// \endcond
  };

  // Type Traits
  template <class ChildT, int FwdColMode, int FwdRowMode, int RevColMode, int RevRowMode>
  struct IsMultiplyAccessible<RemapView<ChildT,FwdColMode,FwdRowMode,RevColMode,RevRowMode> > : public IsMultiplyAccessible<ChildT> {};

  /// Transpose an image.
  template <class ImageT>
  RemapView<ImageT,2,1,2,1> transpose( ImageViewBase<ImageT> const& v ) {
    return RemapView<ImageT,2,1,2,1>( v.impl() );
  }

  /// Rotate an image 180 degrees.
  template <class ImageT>
  RemapView<ImageT,-1,-2,-1,-2> rotate_180( ImageViewBase<ImageT> const& v ) {
    return RemapView<ImageT,-1,-2,-1,-2>( v.impl() );
  }

  /// Rotate an image 90 degrees clockwise.
  template <class ImageT>
  RemapView<ImageT,-2,1,2,-1> rotate_90_cw( ImageViewBase<ImageT> const& v ) {
    return RemapView<ImageT,-2,1,2,-1>( v.impl() );
  }

  /// Rotate an image 90 degrees counter-clockwise.
  template <class ImageT>
  RemapView<ImageT,2,-1,-2,1> rotate_90_ccw( ImageViewBase<ImageT> const& v ) {
    return RemapView<ImageT,2,-1,-2,1>( v.impl() );
  }

  /// Flip an image vertically.
  template <class ImageT>
  RemapView<ImageT,1,-2,1,-2> flip_vertical( ImageViewBase<ImageT> const& v ) {
    return RemapView<ImageT,1,-2,1,-2>( v.impl() );
  }

  /// Flip an image horizontally.
  template <class ImageT>
  RemapView<ImageT,-1,2,-1,2> flip_horizontal( ImageViewBase<ImageT> const& v ) {
    return RemapView<ImageT,-1,2,-1,2>( v.impl() );
  }


  // *******************************************************************
  // crop()
  // *******************************************************************

  // Class definition
  template <class ImageT>
  class CropView : public ImageViewBase< CropView<ImageT> > {
  private:
    typedef typename boost::mpl::if_<IsFloatingPointIndexable<ImageT>, float, int>::type offset_type;

    ImageT m_image;
    offset_type m_ci, m_cj;
    unsigned m_di, m_dj;

  public:
    typedef typename ImageT::pixel_type pixel_type;
    typedef typename ImageT::result_type result_type;
    typedef typename ImageT::pixel_accessor pixel_accessor;

    CropView( ImageT const& image, offset_type const upper_left_i, offset_type const upper_left_j, unsigned const width, unsigned const height ) : 
      m_image(image), m_ci(upper_left_i), m_cj(upper_left_j), m_di(width), m_dj(height) {}

    template<class RealT>
    CropView( ImageT const& image, BBox<RealT,2> const& bbox) :
      m_image(image), 
      m_ci((offset_type)(bbox.min()[0])), 
      m_cj((offset_type)(bbox.min()[1])), 
      m_di((unsigned)round(bbox.width())), 
      m_dj((unsigned)round(bbox.height())) {}

    inline unsigned cols() const { return m_di; }
    inline unsigned rows() const { return m_dj; }
    inline unsigned planes() const { return m_image.planes(); }

    inline pixel_accessor origin() const { return m_image.origin().advance(m_ci, m_cj); }

    inline result_type operator()( offset_type i, offset_type j, int p=0 ) const { return m_image(m_ci + i, m_cj + j, p); }

    CropView& operator=( CropView const& view ) {
      view.rasterize( *this, BBox2i(0,0,view.impl().cols(),view.impl().rows()) );
      return *this;
    }

    template <class ViewT>
    CropView& operator=( ImageViewBase<ViewT> const& view ) {
      view.impl().rasterize( *this, BBox2i(0,0,view.impl().cols(),view.impl().rows()) );
      return *this;
    }

    ImageT const& child() const {
      return m_image;
    }

    /// \cond INTERNAL
    typedef CropView<typename ImageT::prerasterize_type> prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i bbox ) const {
      return prerasterize_type( m_image.prerasterize(bbox+Vector2i(m_ci,m_cj)), m_ci, m_cj, m_di, m_dj );
    }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i bbox ) const {
      m_image.rasterize( dest, bbox+Vector2i(m_ci,m_cj) );
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
  inline CropView<ImageT> crop( ImageT const& v, int upper_left_x, int upper_left_y, int width, int height ) {
    return CropView<ImageT>( v, upper_left_x, upper_left_y, width, height );
  }

  /// Crop an image.
  template <class ImageT>
  inline CropView<ImageT> crop( ImageT const& v, BBox<int,2> const& bbox ) {
    return CropView<ImageT>( v, bbox.min().x(), bbox.min().y(), bbox.width(), bbox.height() );
  }


  // *******************************************************************
  // subsample()
  // *******************************************************************

  // Specialized image accessor
  template <class ImageAccT>
  class SubsamplePixelAccessor {
    ImageAccT m_acc;
    unsigned m_xdelta, m_ydelta;
  public:
    typedef typename ImageAccT::pixel_type pixel_type;
    typedef typename ImageAccT::result_type result_type;
    SubsamplePixelAccessor( ImageAccT const& iter , unsigned xdelta, unsigned ydelta) : m_acc(iter), m_xdelta(xdelta), m_ydelta(ydelta) {}

    inline SubsamplePixelAccessor& next_col() { m_acc.advance(  m_xdelta, 0 ); return *this; }
    inline SubsamplePixelAccessor& prev_col() { m_acc.advance( -m_xdelta, 0 ); return *this; }
    inline SubsamplePixelAccessor& next_row() { m_acc.advance( 0,  m_ydelta ); return *this; }
    inline SubsamplePixelAccessor& prev_row() { m_acc.advance( 0, -m_ydelta ); return *this; }
    inline SubsamplePixelAccessor& next_plane() { m_acc.next_plane(); return *this; }
    inline SubsamplePixelAccessor& prev_plane() { m_acc.prev_plane(); return *this; }
    inline SubsamplePixelAccessor& advance( ptrdiff_t di, ptrdiff_t dj, ptrdiff_t dp=0 ) { m_acc.advance((ptrdiff_t)m_xdelta*di,(ptrdiff_t)m_ydelta*dj,dp); return *this; }

    inline result_type operator*() const { return *m_acc; }
  };

  // Class definition
  template <class ImageT>
  class SubsampleView : public ImageViewBase<SubsampleView<ImageT> > {
    ImageT m_image;
    unsigned m_xdelta, m_ydelta;
  public:
    typedef typename ImageT::pixel_type pixel_type;
    typedef typename ImageT::result_type result_type;
    typedef SubsamplePixelAccessor<typename ImageT::pixel_accessor> pixel_accessor;

    SubsampleView( ImageT const& image, unsigned subsampling_factor ) : m_image(image), m_xdelta(subsampling_factor), m_ydelta(subsampling_factor) {}
    SubsampleView( ImageT const& image, unsigned xfactor, unsigned yfactor ) : m_image(image), m_xdelta(xfactor), m_ydelta(yfactor) {}

    inline unsigned cols() const { return 1 + (m_image.cols()-1)/m_xdelta; }
    inline unsigned rows() const { return 1 + (m_image.rows()-1)/m_ydelta; }
    inline unsigned planes() const { return m_image.planes(); }

    inline pixel_accessor origin() const { return pixel_accessor(m_image.origin(), m_xdelta, m_ydelta); }
    inline result_type operator()( int i, int j ) const { return m_image(m_xdelta*i,m_ydelta*j); }
    inline result_type operator()( int i, int j, int p ) const { return m_image(m_xdelta*i,m_ydelta*j,p); }

    ImageT const& child() const { return m_image; }

    /// \cond INTERNAL
    typedef SubsampleView<typename ImageT::prerasterize_type> prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i bbox ) const { return prerasterize_type( m_image.prerasterize(bbox), m_xdelta, m_ydelta ); }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i bbox ) const { vw::rasterize( prerasterize(bbox), dest, bbox ); }
    /// \endcond
  };

  /// \cond INTERNAL
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
  // select_col()
  // *******************************************************************

  /// Return a single column from an image
  /// \see vw::select_col
  template <class ImageT>
  class SelectColView : public ImageViewBase<SelectColView<ImageT> > {
    ImageT m_image;
    unsigned m_col;
  public:
    typedef typename ImageT::pixel_type pixel_type;
    typedef typename ImageT::result_type result_type;
    typedef typename ImageT::pixel_accessor pixel_accessor;

    SelectColView( ImageT const& image, unsigned col ) : m_image(image), m_col(col) {}

    inline unsigned cols() const { return 1; }
    inline unsigned rows() const { return m_image.rows(); }
    inline unsigned planes() const { return m_image.planes(); }

    inline pixel_accessor origin() const { return m_image.origin().advance(m_col,0,0); }
    inline result_type operator()( int i, int j, int p=0) const { return m_image(m_col,j,p); }

    template <class ViewT>
    SelectColView& operator=( ImageViewBase<ViewT> const& view ) {
      view.impl().rasterize( *this, BBox2i(0,0,view.impl().cols(),view.impl().rows()) );
      return *this;
    }

    /// \cond INTERNAL
    typedef SelectColView<typename ImageT::prerasterize_type> prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i bbox ) const { return prerasterize_type( m_image.prerasterize(bbox), m_col ); }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i bbox ) const { vw::rasterize( prerasterize(bbox), dest, bbox ); }
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
  SelectColView<ImageT> select_col( ImageViewBase<ImageT> const& v, unsigned col ) {
    return SelectColView<ImageT>( v.impl(), col );
  }


  // *******************************************************************
  // select_row()
  // *******************************************************************

  /// Return a single row from an image
  /// \see vw::select_row
  template <class ImageT>
  class SelectRowView : public ImageViewBase<SelectRowView<ImageT> > {
    ImageT m_image;
    unsigned m_row;
  public:
    typedef typename ImageT::pixel_type pixel_type;
    typedef typename ImageT::result_type result_type;
    typedef typename ImageT::pixel_accessor pixel_accessor;

    SelectRowView( ImageT const& image, unsigned row ) : m_image(image), m_row(row) {}

    inline unsigned cols() const { return m_image.cols(); }
    inline unsigned rows() const { return 1; }
    inline unsigned planes() const { return m_image.planes(); }

    inline pixel_accessor origin() const { return m_image.origin().advance(0,m_row,0); }

    inline result_type operator()( int i, int j, int p=0) const { return m_image(i,m_row,p); }

    template <class ViewT>
    SelectRowView& operator=( ImageViewBase<ViewT> const& view ) {
      view.impl().rasterize( *this, BBox2i(0,0,view.impl().cols(),view.impl().rows()) );
      return *this;
    }

    typedef SelectRowView<typename ImageT::prerasterize_type> prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i bbox ) const { return prerasterize_type( m_image.prerasterize(bbox), m_row ); }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i bbox ) const { vw::rasterize( prerasterize(bbox), dest, bbox ); }
  };

  // View type Traits
  template <class ImageT>
  struct IsMultiplyAccessible<SelectRowView<ImageT> > : public IsMultiplyAccessible<ImageT> {};

  /// Extracts a single row of an image.  This function returns a
  /// writeable view of a single row of a multi-row image.  
  /// \see vw::SelectRowView
  template <class ImageT>
  SelectRowView<ImageT> select_row( ImageViewBase<ImageT> const& v, unsigned row ) {
    return SelectRowView<ImageT>( v.impl(), row );
  }


  // *******************************************************************
  // select_plane() 
  // *******************************************************************

  /// Return a single plane from a multi-plane image
  /// \see vw::select_plane
  template <class ImageT>
  class SelectPlaneView : public ImageViewBase<SelectPlaneView<ImageT> > {
    ImageT m_image;
    unsigned m_plane;
  public:
    typedef typename ImageT::pixel_type pixel_type;
    typedef typename ImageT::result_type result_type;
    typedef typename ImageT::pixel_accessor pixel_accessor;

    SelectPlaneView( ImageT const& image, unsigned plane ) : m_image(image), m_plane(plane) {}

    inline unsigned cols() const { return m_image.cols(); }
    inline unsigned rows() const { return m_image.rows(); }
    inline unsigned planes() const { return 1; }

    inline pixel_accessor origin() const { return m_image.origin().advance(0,0,m_plane); }
    inline result_type operator()( int i, int j ) const { return m_image(i,j,m_plane); }
    inline result_type operator()( int i, int j, int p ) const { return m_image(i,j,m_plane+p); }

    template <class ViewT>
    SelectPlaneView& operator=( ImageViewBase<ViewT> const& view ) {
      view.impl().rasterize( *this, BBox2i(0,0,view.impl().cols(),view.impl().rows()) );
      return *this;
    }

    /// \cond INTERNAL
    typedef SelectPlaneView<typename ImageT::prerasterize_type> prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i bbox ) const { return prerasterize_type( m_image.prerasterize(bbox), m_plane ); }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i bbox ) const { vw::rasterize( prerasterize(bbox), dest, bbox ); }
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
  SelectPlaneView<ImageT> select_plane( ImageViewBase<ImageT> const& v, unsigned plane ) {
    return SelectPlaneView<ImageT>( v.impl(), plane );
  }


  // *******************************************************************
  // select_channel()
  // *******************************************************************

  /// A channel selecting functor, used by \ref select_channel().
  template <class ImageT>
  struct SelectChannelFunctor {
    int m_channel;
  public:
    SelectChannelFunctor( int channel ) : m_channel(channel) {}

    // Computes an appropriate reference-to-channel type.
    typedef typename CompoundChannelType<typename ImageT::pixel_type>::type base_channel_type;
    typedef typename boost::mpl::if_<boost::is_const<ImageT>,typename boost::add_const<base_channel_type>::type,base_channel_type>::type channel_type;
    typedef typename boost::mpl::if_<boost::is_reference<typename ImageT::result_type>,channel_type&,channel_type>::type result_type;

    template <class ArgT>
    result_type operator()( ArgT& pixel ) const {
      return compound_select_channel<result_type>(pixel,m_channel);
    }
  };

  /// Extracts a single channel of a multi-channel image.  This function
  /// returns a writeable view of a single channel of a multi-channel
  /// image.
  template <class ImageT>
  UnaryPerPixelView<ImageT,SelectChannelFunctor<ImageT> >
  inline select_channel( ImageViewBase<ImageT>& image, int channel ) {
    return UnaryPerPixelView<ImageT,SelectChannelFunctor<ImageT> >( image.impl(), SelectChannelFunctor<ImageT>(channel) );
  }

  /// Extracts a single channel of a multi-channel image (const overload).
  /// This function returns a writeable view of a single channel of a
  /// multi-channel image.
  template <class ImageT>
  UnaryPerPixelView<ImageT,SelectChannelFunctor<const ImageT> >
  inline select_channel( ImageViewBase<ImageT> const& image, int channel ) {
    return UnaryPerPixelView<ImageT,SelectChannelFunctor<const ImageT> >( image.impl(), SelectChannelFunctor<const ImageT>(channel) );
  }


  // *******************************************************************
  // channels_to_planes()
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
    typedef typename CompoundChannelType<typename ImageAccT::pixel_type>::type pixel_type;
    typedef typename CopyCVR<typename ImageAccT::result_type, pixel_type>::type result_type;

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

  /// A view that turns a one plane, multi-channel view into a mulit-plane, one channel view.
  /// \see vw::channels_to_planes
  template <class ImageT>
  class ChannelsToPlanesView : public ImageViewBase<ChannelsToPlanesView<ImageT> > {
    ImageT m_image;
  public:

    typedef typename CompoundChannelType<typename ImageT::pixel_type>::type pixel_type;
    typedef typename CopyCVR<typename ImageT::result_type, pixel_type>::type result_type;

    typedef typename boost::mpl::if_< IsCompound<typename ImageT::pixel_type>, 
                                      ChannelsToPlanesAccessor<typename ImageT::pixel_accessor>,
                                      typename ImageT::pixel_accessor >::type pixel_accessor;
    
    ChannelsToPlanesView( ImageT const& image ) : m_image(image) {
      VW_ASSERT( m_image.planes()==1 , ArgumentErr() << "ChannelsToPlanesView: The image must be single plane.");
    }
    
    inline unsigned cols() const { return m_image.cols(); }
    inline unsigned rows() const { return m_image.rows(); }
    inline unsigned planes() const { return m_image.channels(); }

    inline pixel_accessor origin() const { return pixel_accessor(m_image.origin()); }

    inline result_type operator()( int i, int j, int p=0 ) const {
      return compound_select_channel<result_type>(m_image(i,j),p);
    }

    template <class ViewT>
    ChannelsToPlanesView& operator=( ImageViewBase<ViewT> const& view ) {
      view.impl().rasterize( *this, BBox2i(0,0,view.impl().cols(),view.impl().rows()) );
      return *this;
    }

    ImageT const& child() const {
      return m_image;
    }

    /// \cond INTERNAL
    typedef ChannelsToPlanesView<typename ImageT::prerasterize_type> prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i bbox ) const { return prerasterize_type( m_image.prerasterize(bbox) ); }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i bbox ) const { vw::rasterize( prerasterize(bbox), dest, bbox ); }
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
    inline prerasterize_type prerasterize( BBox2i bbox ) const { return prerasterize_type( m_image.prerasterize(bbox) ); }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i bbox ) const { vw::rasterize( prerasterize(bbox), dest, bbox ); }
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
    PixelT operator()( ArgT const& pixel ) const {
      return (PixelT)(pixel);
    }
  };

  /// Create a new image view by statically casting the pixels to a new type.
  template <class PixelT, class ImageT>
  inline UnaryPerPixelView<ImageT,PixelCastFunctor<PixelT> > pixel_cast( ImageViewBase<ImageT> const& image ) {
    return UnaryPerPixelView<ImageT,PixelCastFunctor<PixelT> >( image.impl() );
  }


  // *******************************************************************
  // channel_cast()
  // *******************************************************************

  /// A pixel channel casting functor, used by \ref channel_cast().
  template <class ChannelT>
  struct PixelChannelCastFunctor : UnaryReturnBinaryTemplateBind2nd<CompoundChannelCast,ChannelT> {
    template <class ArgT>
    typename CompoundChannelCast<ArgT,ChannelT>::type operator()( ArgT const& pixel ) const {
      return compound_channel_cast<ChannelT>(pixel);
    }
  };

  /// Create a new image view by statically casting the channels of the pixels to a new type.
  template <class ChannelT, class ImageT>
  inline UnaryPerPixelView<ImageT,PixelChannelCastFunctor<ChannelT> > channel_cast( ImageViewBase<ImageT> const& image ) {
    return UnaryPerPixelView<ImageT,PixelChannelCastFunctor<ChannelT> >( image.impl() );
  }

} // namespace vw

#endif // __VW_IMAGE_MANIPULATION_H__
