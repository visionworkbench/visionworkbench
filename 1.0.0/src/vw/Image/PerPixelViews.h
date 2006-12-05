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

/// \file PerPixelViews.h
/// 
/// Per-pixel and per-pixel-channel image views, including pixel cast views.
/// 
#ifndef __VW_IMAGE_PERPIXELVIEWS_H__
#define __VW_IMAGE_PERPIXELVIEWS_H__

#include <boost/type_traits.hpp>
#include <boost/mpl/logical.hpp>
#include <boost/mpl/if.hpp>
#include <boost/utility/result_of.hpp>

#include <vw/Core/Functors.h>
#include <vw/Image/ImageViewBase.h>

namespace vw {

  // *******************************************************************
  // UnaryPerPixelView
  // *******************************************************************
  
  // Specialized Accessor
  template <class ImageIterT, class FuncT>
  class UnaryPerPixelAccessor {
    ImageIterT m_iter;
    FuncT const& m_func;
  public:
    typedef typename boost::result_of<FuncT(typename ImageIterT::pixel_type)>::type result_type;
    typedef typename boost::remove_reference<result_type>::type pixel_type;
    UnaryPerPixelAccessor( ImageIterT const& iter, FuncT const& func ) : m_iter(iter), m_func(func) {}
    inline UnaryPerPixelAccessor& next_col() { m_iter.next_col(); return *this; }
    inline UnaryPerPixelAccessor& prev_col() { m_iter.prev_col(); return *this; }
    inline UnaryPerPixelAccessor& next_row() { m_iter.next_row(); return *this; }
    inline UnaryPerPixelAccessor& prev_row() { m_iter.prev_row(); return *this; }
    inline UnaryPerPixelAccessor& next_plane() { m_iter.next_plane(); return *this; }
    inline UnaryPerPixelAccessor& prev_plane() { m_iter.prev_plane(); return *this; }
    inline UnaryPerPixelAccessor& advance( ptrdiff_t di, ptrdiff_t dj, ptrdiff_t dp=0 ) { m_iter.advance(di,dj,dp); return *this; }
    inline result_type operator*() const { return m_func(*m_iter); }
  };

  // Image View Class Declaration
  template <class ImageT, class FuncT>
  class UnaryPerPixelView : public ImageViewBase<UnaryPerPixelView<ImageT,FuncT> > {
    ImageT m_image;
    FuncT m_func;
  public:
    typedef typename boost::result_of<FuncT(typename ImageT::pixel_type)>::type result_type;
    typedef typename boost::remove_reference<result_type>::type pixel_type;
    typedef UnaryPerPixelAccessor<typename ImageT::pixel_accessor, FuncT> pixel_accessor;

    UnaryPerPixelView( ImageT const& image ) : m_image(image), m_func() {}
    UnaryPerPixelView( ImageT const& image, FuncT const& func ) : m_image(image), m_func(func) {}

    inline unsigned cols() const { return m_image.cols(); }
    inline unsigned rows() const { return m_image.rows(); }
    inline unsigned planes() const { return m_image.planes(); }

    inline pixel_accessor origin() const { return pixel_accessor(m_image.origin(),m_func); }
    inline result_type operator()( int i, int j, int p=0 ) const { return m_func(m_image(i,j,p)); }

    template <class ViewT>
    UnaryPerPixelView& operator=( ImageViewBase<ViewT> const& view ) {
      view.impl().rasterize( *this, BBox2i(0,0,view.impl().cols(),view.impl().rows()) );
      return *this;
    }

    /// \cond INTERNAL
    typedef UnaryPerPixelView<typename ImageT::prerasterize_type, FuncT> prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i bbox ) const { return prerasterize_type( m_image.prerasterize(bbox), m_func ); }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i bbox ) const { vw::rasterize( prerasterize(bbox), dest, bbox ); }
    /// \endcond
  };
  
  /// \cond INTERNAL
  // View type Traits.  This exists mainly for select_channel(), and may not 
  // be correct in all cases.  Perhaps it should be specialized there instead?
  template <class ImageT, class FuncT>
  struct IsMultiplyAccessible<UnaryPerPixelView<ImageT,FuncT> > : boost::is_reference<typename UnaryPerPixelView<ImageT,FuncT>::result_type>::type {};
  /// \endcond


  // *******************************************************************
  // BinaryPerPixelView
  // *******************************************************************

  // Specialized Accessor
  template <class Image1IterT, class Image2IterT, class FuncT>
  class BinaryPerPixelAccessor {
    Image1IterT m_iter1;
    Image2IterT m_iter2;
    FuncT const& m_func;
  public:
    typedef typename boost::result_of<FuncT(typename Image1IterT::pixel_type,typename Image2IterT::pixel_type)>::type result_type;
    typedef typename boost::remove_reference<result_type>::type pixel_type;

    BinaryPerPixelAccessor( Image1IterT const& iter1, Image2IterT const& iter2, FuncT const& func ) : m_iter1(iter1), m_iter2(iter2), m_func(func) {}
    inline BinaryPerPixelAccessor& next_col() { m_iter1.next_col(); m_iter2.next_col(); return *this; }
    inline BinaryPerPixelAccessor& prev_col() { m_iter1.prev_col(); m_iter2.prev_col(); return *this; }
    inline BinaryPerPixelAccessor& next_row() { m_iter1.next_row(); m_iter2.next_row(); return *this; }
    inline BinaryPerPixelAccessor& prev_row() { m_iter1.prev_row(); m_iter2.prev_row(); return *this; }
    inline BinaryPerPixelAccessor& next_plane() { m_iter1.next_plane(); m_iter2.next_plane(); return *this; }
    inline BinaryPerPixelAccessor& prev_plane() { m_iter1.prev_plane(); m_iter2.prev_plane(); return *this; }
    inline BinaryPerPixelAccessor& advance( ptrdiff_t di, ptrdiff_t dj, ptrdiff_t dp=0 ) 
      { m_iter1.advance(di,dj,dp); m_iter2.advance(di,dj,dp); return *this; }
    inline result_type operator*() const { return m_func(*m_iter1,*m_iter2); }
  };

  // Image View Class Definition
  template <class Image1T, class Image2T, class FuncT>
  class BinaryPerPixelView : public ImageViewBase<BinaryPerPixelView<Image1T,Image2T,FuncT> >
  {
  private:
    Image1T m_image1;
    Image2T m_image2;
    FuncT m_func;
  public:
    typedef typename boost::result_of<FuncT(typename Image1T::pixel_type, typename Image2T::pixel_type)>::type result_type;
    typedef typename boost::remove_reference<result_type>::type pixel_type;

    typedef BinaryPerPixelAccessor<typename Image1T::pixel_accessor,
                                   typename Image2T::pixel_accessor,
                                   FuncT> pixel_accessor;

    BinaryPerPixelView( Image1T const& image1, Image2T const& image2 )
      : m_image1(image1), m_image2(image2), m_func()
    {
      VW_ASSERT( m_image1.cols()==m_image2.cols() &&
		 m_image1.rows()==m_image1.rows() &&
		 m_image1.planes()==m_image2.planes(),
                 ArgumentErr() << "BinaryPerPixelView: Images must have same dimensions in binary image operation." );
    }

    BinaryPerPixelView( Image1T const& image1, Image2T const& image2, FuncT const& func )
      : m_image1(image1), m_image2(image2), m_func(func)
    {
      VW_ASSERT( m_image1.cols()==m_image2.cols() &&
		 m_image1.rows()==m_image1.rows() &&
		 m_image1.planes()==m_image2.planes(),
                 ArgumentErr() << "BinaryPerPixelView: Images must have same dimensions in binary image operation." );
    }

    inline unsigned cols() const { return m_image1.cols(); }
    inline unsigned rows() const { return m_image1.rows(); }
    inline unsigned planes() const { return m_image1.planes(); }

    inline pixel_accessor origin() const { return pixel_accessor(m_image1.origin(),m_image2.origin(),m_func); }
    inline result_type operator()( int i, int j, int p=0 ) const { return m_func(m_image1(i,j,p),m_image2(i,j,p)); }

    /// \cond INTERNAL
    typedef BinaryPerPixelView<typename Image1T::prerasterize_type, typename Image2T::prerasterize_type, FuncT> prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i bbox ) const { return prerasterize_type( m_image1.prerasterize(bbox), m_image2.prerasterize(bbox), m_func ); }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i bbox ) const { vw::rasterize( prerasterize(bbox), dest, bbox ); }
    /// \endcond
  };

};

#endif // __VW_IMAGE_PERPIXELVIEWS_H__
