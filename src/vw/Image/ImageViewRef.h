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

/// \file ImageViewRef.h
///
/// A generic image view reference class.
///
/// Most Vision Workbench image processing functions simply return 
/// image view objects that lazily represent the processed data.
/// Under some circumstances it is helpful to be able to hold onto 
/// such a processed view without rasterizing it.  Ordinarily this 
/// requires knowing the full type of the view.  When this is not 
/// acceptable, the \ref vw::ImageViewRef class allows you to hold 
/// a virtualized reference to an arbitrary image view with a given 
/// pixel type.
///
#ifndef __VW_IMAGE_IMAGE_VIEW_REF_H__
#define __VW_IMAGE_IMAGE_VIEW_REF_H__

#include <boost/type_traits.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>

#include <vw/Image/ImageView.h>
#include <vw/Image/Manipulation.h>

namespace vw {

  /// \cond INTERNAL
  template <class PixelT>
  class ImageViewRefAccessorBase {
  public:
    virtual ~ImageViewRefAccessorBase() {}
    virtual ImageViewRefAccessorBase* copy() const = 0;
    virtual void next_col() = 0;
    virtual void prev_col() = 0;
    virtual void next_row() = 0;
    virtual void prev_row() = 0;
    virtual void next_plane() = 0;
    virtual void prev_plane() = 0;
    virtual void advance( ptrdiff_t di, ptrdiff_t dj, ptrdiff_t dp=0 ) = 0;
    virtual PixelT operator*() const = 0;
  };

  template <class IterT>
  class ImageViewRefAccessorImpl : public ImageViewRefAccessorBase<typename IterT::pixel_type> {
  private:
    IterT m_iter;
  public:
    typedef typename IterT::pixel_type pixel_type;

    ImageViewRefAccessorImpl( IterT const& iter ) : m_iter(iter) {}
    virtual ~ImageViewRefAccessorImpl() {}

    virtual ImageViewRefAccessorBase<pixel_type>* copy() const { return new ImageViewRefAccessorImpl(m_iter); }

    virtual void next_col() { m_iter.next_col(); }
    virtual void prev_col() { m_iter.prev_col(); }
    virtual void next_row() { m_iter.next_row(); }
    virtual void prev_row() { m_iter.prev_row(); }
    virtual void next_plane() { m_iter.next_plane(); }
    virtual void prev_plane() { m_iter.prev_plane(); }
    virtual void advance( ptrdiff_t di, ptrdiff_t dj, ptrdiff_t dp=0 ) { m_iter.advance(di,dj,dp); }
    virtual pixel_type operator*() const { return *m_iter; }
  };
  /// \endcond

  /// A special virtualized accessor adaptor.
  ///
  /// This accessor adaptor is used by the \ref vw::ImageViewRef class.
  template <class PixelT>
  class ImageViewRefAccessor {
  private:
    boost::scoped_ptr< ImageViewRefAccessorBase<PixelT> > m_iter;
  public:
    typedef PixelT pixel_type;
    typedef PixelT result_type;

    template <class IterT> ImageViewRefAccessor( IterT const& iter ) : m_iter( new ImageViewRefAccessorImpl<IterT>(iter) ) {}
    ~ImageViewRefAccessor() {}

    ImageViewRefAccessor( ImageViewRefAccessor const& other ) : m_iter( other.m_iter->copy() ) {}
    ImageViewRefAccessor& operator=( ImageViewRefAccessor const& other ) { m_iter = other.m_iter->copy(); }

    inline ImageViewRefAccessor& next_col() { m_iter->next_col(); return *this; }
    inline ImageViewRefAccessor& prev_col() { m_iter->prev_col(); return *this; }
    inline ImageViewRefAccessor& next_row() { m_iter->next_row(); return *this; }
    inline ImageViewRefAccessor& prev_row() { m_iter->prev_row(); return *this; }
    inline ImageViewRefAccessor& next_plane() { m_iter->next_plane(); return *this; }
    inline ImageViewRefAccessor& prev_plane() { m_iter->prev_plane(); return *this; }
    inline ImageViewRefAccessor& advance( ptrdiff_t di, ptrdiff_t dj, ptrdiff_t dp=0 ) { m_iter->advance(di,dj,dp=0); return *this; }
    inline pixel_type operator*() const { return *(*m_iter); }
  };


  /// \cond INTERNAL
  // Base class definition
  template <class PixelT>
  class ImageViewRefBase {
  public:
    typedef PixelT pixel_type;
    typedef ImageViewRefAccessor<PixelT> pixel_accessor;

    virtual ~ImageViewRefBase() {}
    
    virtual unsigned cols() const = 0;
    virtual unsigned rows() const = 0;
    virtual unsigned planes() const = 0;
    virtual pixel_type operator()( int i, int j ) const = 0;
    virtual pixel_type operator()( int i, int j, int p ) const = 0;
    virtual pixel_accessor origin() const = 0;

    virtual void rasterize( ImageView<pixel_type> const& dest, BBox2i bbox ) const = 0;
  };

  // ImageViewRef class implementation
  template <class ViewT>
  class ImageViewRefImpl : public ImageViewRefBase<typename ViewT::pixel_type> {
  private:
    ViewT m_view;
  public:
    typedef typename ViewT::pixel_type pixel_type;
    typedef ImageViewRefAccessor<typename ViewT::pixel_type> pixel_accessor;

    ImageViewRefImpl( ImageViewBase<ViewT> const& view ) : m_view(view.impl()) {}
    virtual ~ImageViewRefImpl() {}

    virtual unsigned cols() const { return m_view.cols(); }
    virtual unsigned rows() const { return m_view.rows(); }
    virtual unsigned planes() const { return m_view.planes(); }
    virtual pixel_type operator()( int i, int j ) const { return m_view(i,j); }
    virtual pixel_type operator()( int i, int j, int p ) const { return m_view(i,j,p); }
    virtual pixel_accessor origin() const { return m_view.origin(); }

    virtual void rasterize( ImageView<pixel_type> const& dest, BBox2i bbox ) const { m_view.rasterize( dest, bbox ); }

    ViewT const& child() const { return m_view; }
  };
  /// \endcond


  /// A virtualized image view reference object.
  ///
  /// This class behaves as a reference to an arbitrary image view 
  /// with the given pixel type.  The purpose of this class is to 
  /// hide the full type of a view behind a veil of abstraction, 
  /// making things like run-time polymorphic behavior possible. 
  /// The inevitable cost of this flexibility is one virtual 
  /// function call per method invocation.  In many cases there 
  /// are additional costs associated with not being able to 
  /// perform template-based optimizations at compile time.
  ///
  /// Like any C++ reference, you bind an ImageViewRef to a view 
  /// using a constructor and future operations act on the bound 
  /// object instead of the reference object itself.
  ///
  /// The current implementation of ImageViewRef is read-only.
  template <class PixelT>
  class ImageViewRef : public ImageViewBase<ImageViewRef<PixelT> > {
  private:
    boost::shared_ptr< ImageViewRefBase<PixelT> > m_view;
  public:
    typedef PixelT pixel_type;
    typedef PixelT result_type;
    typedef ImageViewRefAccessor<PixelT> pixel_accessor;

    template <class ViewT> ImageViewRef( ImageViewBase<ViewT> const& view ) : m_view( new ImageViewRefImpl<ViewT>(view) ) {}
    ~ImageViewRef() {}

    inline unsigned cols() const { return m_view->cols(); }
    inline unsigned rows() const { return m_view->rows(); }
    inline unsigned planes() const { return m_view->planes(); }
    inline pixel_type operator()( int i, int j ) const { return m_view->operator()(i,j); }
    inline pixel_type operator()( int i, int j, int p ) const { return m_view->operator()(i,j,p); }
    inline pixel_accessor origin() const { return m_view->origin(); }

    /// \cond INTERNAL
    typedef CropView<ImageView<PixelT> > prerasterize_type;

    inline prerasterize_type prerasterize( BBox2i bbox ) const {
      // If we're wrapping a plain ImageView, we can avoid copying the data.
      ImageViewRefImpl<ImageView<PixelT> > *image_ptr = dynamic_cast<ImageViewRefImpl<ImageView<PixelT> >*>( m_view.get() );
      if( image_ptr ) return CropView<ImageView<PixelT> >( image_ptr->child(), 0, 0, cols(), rows() );
      // Otherwise, we must rasterize ourselves....
      ImageView<PixelT> buf( bbox.width(), bbox.height(), planes() );
      m_view->rasterize( buf, bbox );
      return CropView<ImageView<PixelT> >( buf, BBox2i(-bbox.min().x(),-bbox.min().y(),bbox.width(),bbox.height()) );
    }

    template <class DestT> inline void rasterize( DestT const& dest, BBox2i bbox ) const {
      vw::rasterize( prerasterize(bbox), dest, bbox );
    }

    // A special performance-enhancing overload for rasterizing directly into 
    // an ImageView with the proper pixel type.  This cannot be templatized 
    // or otherwise generalized because it calls m_view's virtual rasterize 
    // method.
    inline void rasterize( ImageView<PixelT> const& dest, BBox2i bbox ) const {
      m_view->rasterize( dest, bbox );
    }
    /// \endcond
  };

} // namespace vw

#endif // __VW_IMAGE_IMAGE_VIEW_REF_H__
