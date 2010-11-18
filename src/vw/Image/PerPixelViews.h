// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
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

#include <vw/Image/ImageViewBase.h>
#include <vw/Image/PixelAccessors.h>

namespace vw {

  // *******************************************************************
  // PerPixelIndexView
  // *******************************************************************

  template<class FuncT>
  class PerPixelIndexView : public ImageViewBase<PerPixelIndexView<FuncT> > {
    FuncT m_func;
    int32 m_cols, m_rows, m_planes;
  public:
    typedef typename FuncT::result_type pixel_type;
    typedef pixel_type const result_type;
    typedef ProceduralPixelAccessor<PerPixelIndexView> pixel_accessor;

    PerPixelIndexView( FuncT func, int32 cols, int32 rows, int32 planes = 1 )
      : m_func(func), m_cols(cols), m_rows(rows), m_planes(planes) {}

    inline int32 cols() const { return m_cols; }
    inline int32 rows() const { return m_rows; }
    inline int32 planes() const { return m_planes; }

    inline pixel_accessor origin() const { return pixel_accessor(*this); }

    inline result_type operator()( double i, double j, int32 p = 0 ) const {
      return m_func( i, j, p );
    }

    typedef PerPixelIndexView prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i /*bbox*/ ) const { return *this; }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i bbox ) const {
      vw::rasterize( prerasterize(bbox), dest, bbox );
    }
  };

  template <class FuncT>
  struct IsMultiplyAccessible<PerPixelIndexView<FuncT> > : public true_type {};

  template <class FuncT>
  struct IsFloatingPointIndexable<PerPixelIndexView<FuncT> > : public true_type {};

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
    typedef typename boost::remove_cv<typename boost::remove_reference<result_type>::type>::type pixel_type;
    typedef typename ImageIterT::offset_type offset_type;

    UnaryPerPixelAccessor( ImageIterT const& iter, FuncT const& func ) : m_iter(iter), m_func(func) {}
    inline UnaryPerPixelAccessor& next_col() { m_iter.next_col(); return *this; }
    inline UnaryPerPixelAccessor& prev_col() { m_iter.prev_col(); return *this; }
    inline UnaryPerPixelAccessor& next_row() { m_iter.next_row(); return *this; }
    inline UnaryPerPixelAccessor& prev_row() { m_iter.prev_row(); return *this; }
    inline UnaryPerPixelAccessor& next_plane() { m_iter.next_plane(); return *this; }
    inline UnaryPerPixelAccessor& prev_plane() { m_iter.prev_plane(); return *this; }
    inline UnaryPerPixelAccessor& advance( offset_type di, offset_type dj, ssize_t dp=0 ) { m_iter.advance(di,dj,dp); return *this; }
    inline result_type operator*() const { return m_func(*m_iter); }
  };

  // Image View Class Declaration
  template <class ImageT, class FuncT>
  class UnaryPerPixelView : public ImageViewBase<UnaryPerPixelView<ImageT,FuncT> > {
    ImageT m_image;
    FuncT m_func;
  public:
    typedef typename boost::result_of<FuncT(typename ImageT::pixel_type)>::type result_type;
    typedef typename boost::remove_cv<typename boost::remove_reference<result_type>::type>::type pixel_type;
    typedef UnaryPerPixelAccessor<typename ImageT::pixel_accessor, FuncT> pixel_accessor;

    UnaryPerPixelView( ImageT const& image ) : m_image(image), m_func() {}
    UnaryPerPixelView( ImageT const& image, FuncT const& func ) : m_image(image), m_func(func) {}

    inline int32 cols() const { return m_image.cols(); }
    inline int32 rows() const { return m_image.rows(); }
    inline int32 planes() const { return m_image.planes(); }

    inline pixel_accessor origin() const { return pixel_accessor(m_image.origin(),m_func); }
    inline result_type operator()( int32 i, int32 j, int32 p=0 ) const { return m_func(m_image(i,j,p)); }

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
    typedef typename boost::remove_cv<typename boost::remove_reference<result_type>::type>::type pixel_type;
    typedef typename boost::mpl::if_<boost::is_same<typename Image1IterT::offset_type, typename Image2IterT::offset_type>,
                                     typename Image1IterT::offset_type, int32>::type offset_type;

    BinaryPerPixelAccessor( Image1IterT const& iter1, Image2IterT const& iter2, FuncT const& func ) : m_iter1(iter1), m_iter2(iter2), m_func(func) {}
    inline BinaryPerPixelAccessor& next_col() { m_iter1.next_col(); m_iter2.next_col(); return *this; }
    inline BinaryPerPixelAccessor& prev_col() { m_iter1.prev_col(); m_iter2.prev_col(); return *this; }
    inline BinaryPerPixelAccessor& next_row() { m_iter1.next_row(); m_iter2.next_row(); return *this; }
    inline BinaryPerPixelAccessor& prev_row() { m_iter1.prev_row(); m_iter2.prev_row(); return *this; }
    inline BinaryPerPixelAccessor& next_plane() { m_iter1.next_plane(); m_iter2.next_plane(); return *this; }
    inline BinaryPerPixelAccessor& prev_plane() { m_iter1.prev_plane(); m_iter2.prev_plane(); return *this; }
    inline BinaryPerPixelAccessor& advance( offset_type di, offset_type dj, ssize_t dp=0 )
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
    typedef typename boost::remove_cv<typename boost::remove_reference<result_type>::type>::type pixel_type;

    typedef BinaryPerPixelAccessor<typename Image1T::pixel_accessor,
                                   typename Image2T::pixel_accessor,
                                   FuncT> pixel_accessor;

    BinaryPerPixelView( Image1T const& image1, Image2T const& image2 )
      : m_image1(image1), m_image2(image2), m_func()
    {
      VW_ASSERT( m_image1.cols()==m_image2.cols() &&
                 m_image1.rows()==m_image2.rows() &&
                 m_image1.planes()==m_image2.planes(),
                 ArgumentErr() << "BinaryPerPixelView: Images must have same dimensions in binary image operation." );
    }

    BinaryPerPixelView( Image1T const& image1, Image2T const& image2, FuncT const& func )
      : m_image1(image1), m_image2(image2), m_func(func)
    {
      VW_ASSERT( m_image1.cols()==m_image2.cols() &&
                 m_image1.rows()==m_image2.rows() &&
                 m_image1.planes()==m_image2.planes(),
                 ArgumentErr() << "BinaryPerPixelView: Images must have same dimensions in binary image operation." );
    }

    inline int32 cols() const { return m_image1.cols(); }
    inline int32 rows() const { return m_image1.rows(); }
    inline int32 planes() const { return m_image1.planes(); }

    inline pixel_accessor origin() const { return pixel_accessor(m_image1.origin(),m_image2.origin(),m_func); }
    inline result_type operator()( int32 i, int32 j, int32 p=0 ) const { return m_func(m_image1(i,j,p),m_image2(i,j,p)); }

    /// \cond INTERNAL
    typedef BinaryPerPixelView<typename Image1T::prerasterize_type, typename Image2T::prerasterize_type, FuncT> prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i bbox ) const { return prerasterize_type( m_image1.prerasterize(bbox), m_image2.prerasterize(bbox), m_func ); }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i bbox ) const { vw::rasterize( prerasterize(bbox), dest, bbox ); }
    /// \endcond
  };

  // *******************************************************************
  // TrinaryPerPixelView
  // *******************************************************************

  // Specialized Accessor
  template <class Image1IterT, class Image2IterT, class Image3IterT, class FuncT>
  class TrinaryPerPixelAccessor {
    Image1IterT m_iter1;
    Image2IterT m_iter2;
    Image3IterT m_iter3;
    FuncT const& m_func;
  public:
    typedef typename boost::result_of<FuncT(typename Image1IterT::pixel_type,typename Image2IterT::pixel_type,typename Image3IterT::pixel_type)>::type result_type;
    typedef typename boost::remove_cv<typename boost::remove_reference<result_type>::type>::type pixel_type;
    typedef typename boost::mpl::if_<boost::mpl::and_<boost::is_same<typename Image1IterT::offset_type, typename Image2IterT::offset_type>,
                                                      boost::is_same<typename Image1IterT::offset_type, typename Image3IterT::offset_type> >,
                                     typename Image1IterT::offset_type, int32>::type offset_type;

    TrinaryPerPixelAccessor( Image1IterT const& iter1, Image2IterT const& iter2, Image2IterT const& iter3, FuncT const& func ) : m_iter1(iter1), m_iter2(iter2), m_iter3(iter3), m_func(func) {}
    inline TrinaryPerPixelAccessor& next_col() { m_iter1.next_col(); m_iter2.next_col(); m_iter3.next_col(); return *this; }
    inline TrinaryPerPixelAccessor& prev_col() { m_iter1.prev_col(); m_iter2.prev_col(); m_iter3.prev_col(); return *this; }
    inline TrinaryPerPixelAccessor& next_row() { m_iter1.next_row(); m_iter2.next_row(); m_iter3.next_row(); return *this; }
    inline TrinaryPerPixelAccessor& prev_row() { m_iter1.prev_row(); m_iter2.prev_row(); m_iter3.prev_row(); return *this; }
    inline TrinaryPerPixelAccessor& next_plane() { m_iter1.next_plane(); m_iter2.next_plane(); m_iter3.next_plane(); return *this; }
    inline TrinaryPerPixelAccessor& prev_plane() { m_iter1.prev_plane(); m_iter2.prev_plane(); m_iter3.prev_plane(); return *this; }
    inline TrinaryPerPixelAccessor& advance( offset_type di, offset_type dj, ssize_t dp=0 )
      { m_iter1.advance(di,dj,dp); m_iter2.advance(di,dj,dp); m_iter3.advance(di,dj,dp); return *this; }
    inline result_type operator*() const { return m_func(*m_iter1,*m_iter2,*m_iter3); }
  };

  // Image View Class Definition
  template <class Image1T, class Image2T, class Image3T, class FuncT>
  class TrinaryPerPixelView : public ImageViewBase<TrinaryPerPixelView<Image1T,Image2T,Image3T,FuncT> >
  {
  private:
    Image1T m_image1;
    Image2T m_image2;
    Image3T m_image3;
    FuncT m_func;
  public:
    typedef typename boost::result_of<FuncT(typename Image1T::pixel_type, typename Image2T::pixel_type, typename Image3T::pixel_type)>::type result_type;
    typedef typename boost::remove_cv<typename boost::remove_reference<result_type>::type>::type pixel_type;

    typedef TrinaryPerPixelAccessor<typename Image1T::pixel_accessor,
                                    typename Image2T::pixel_accessor,
                                    typename Image3T::pixel_accessor,
                                    FuncT> pixel_accessor;

    TrinaryPerPixelView( Image1T const& image1, Image2T const& image2, Image3T const& image3 )
      : m_image1(image1), m_image2(image2), m_image3(image3), m_func()
    {
      VW_ASSERT( m_image1.cols()==m_image2.cols() &&
                 m_image1.rows()==m_image2.rows() &&
                 m_image1.planes()==m_image2.planes() &&
                 m_image1.cols()==m_image3.cols() &&
                 m_image1.rows()==m_image3.rows() &&
                 m_image1.planes()==m_image3.planes(),
                 ArgumentErr() << "TrinaryPerPixelView: Images must have same dimensions in trinary image operation." );
    }

    TrinaryPerPixelView( Image1T const& image1, Image2T const& image2, Image3T const& image3, FuncT const& func )
      : m_image1(image1), m_image2(image2), m_image3(image3), m_func(func)
    {
      VW_ASSERT( m_image1.cols()==m_image2.cols() &&
                 m_image1.rows()==m_image2.rows() &&
                 m_image1.planes()==m_image2.planes() &&
                 m_image1.cols()==m_image2.cols() &&
                 m_image1.rows()==m_image2.rows() &&
                 m_image1.planes()==m_image2.planes(),
                 ArgumentErr() << "TrinaryPerPixelView: Images must have same dimensions in trinary image operation." );
    }

    inline int32 cols() const { return m_image1.cols(); }
    inline int32 rows() const { return m_image1.rows(); }
    inline int32 planes() const { return m_image1.planes(); }

    inline pixel_accessor origin() const { return pixel_accessor(m_image1.origin(),m_image2.origin(),m_image3.origin(),m_func); }
    inline result_type operator()( int32 i, int32 j, int32 p=0 ) const { return m_func(m_image1(i,j,p),m_image2(i,j,p),m_image3(i,j,p)); }

    /// \cond INTERNAL
    typedef TrinaryPerPixelView<typename Image1T::prerasterize_type, typename Image2T::prerasterize_type, typename Image3T::prerasterize_type, FuncT> prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i bbox ) const { return prerasterize_type( m_image1.prerasterize(bbox), m_image2.prerasterize(bbox), m_image3.prerasterize(bbox), m_func ); }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i bbox ) const { vw::rasterize( prerasterize(bbox), dest, bbox ); }
    /// \endcond
  };

};

#endif // __VW_IMAGE_PERPIXELVIEWS_H__
