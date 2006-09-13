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

/// \file PerPixelViews.h
/// 
/// Per-pixel and per-pixel-channel image views, including pixel cast views.
/// 
#ifndef __VW_IMAGE__PER_PIXEL_VIEWS_H__
#define __VW_IMAGE__PER_PIXEL_VIEWS_H__

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
  class UnaryPerPixelAccessor
  {
  private:
    ImageIterT m_iter;
    FuncT const& m_func;
  public:
    typedef typename boost::result_of<FuncT(typename ImageIterT::pixel_type)>::type reference_type;
    typename boost::remove_reference<reference_type>::type pixel_type;
    UnaryPerPixelAccessor( ImageIterT const& iter, FuncT const& func ) : m_iter(iter), m_func(func) {}
    inline UnaryPerPixelAccessor& next_col() { m_iter.next_col(); return *this; }
    inline UnaryPerPixelAccessor& prev_col() { m_iter.prev_col(); return *this; }
    inline UnaryPerPixelAccessor& next_row() { m_iter.next_row(); return *this; }
    inline UnaryPerPixelAccessor& prev_row() { m_iter.prev_row(); return *this; }
    inline UnaryPerPixelAccessor& next_plane() { m_iter.next_plane(); return *this; }
    inline UnaryPerPixelAccessor& prev_plane() { m_iter.prev_plane(); return *this; }
    inline UnaryPerPixelAccessor& advance( ptrdiff_t di, ptrdiff_t dj, ptrdiff_t dp=0 ) { m_iter.advance(di,dj,dp); return *this; }
    inline reference_type operator*() const { return m_func(*m_iter); }
  };

  /// \cond INTERNAL
  // Accessor type Traits
  template <class ImageIterT, class FuncT>
  struct IsReferenceable<UnaryPerPixelAccessor<ImageIterT,FuncT> > : boost::is_reference<typename UnaryPerPixelAccessor<ImageIterT,FuncT>::reference_type>::type {};
  /// \endcond

  // Image View Class Declaration
  template <class ImageT, class FuncT>
  class UnaryPerPixelView : public ImageViewBase<UnaryPerPixelView<ImageT,FuncT> > {
  private:
    ImageT m_image;
    FuncT m_func;
  public:
    typedef typename boost::result_of<FuncT(typename ImageT::pixel_type)>::type reference_type;
    typedef typename boost::remove_reference<reference_type>::type pixel_type;
    typedef UnaryPerPixelAccessor<typename ImageT::pixel_accessor, FuncT> pixel_accessor;

    UnaryPerPixelView( ImageT const& image ) : m_image(image), m_func() {}
    UnaryPerPixelView( ImageT const& image, FuncT const& func ) : m_image(image), m_func(func) {}

    inline unsigned cols() const { return m_image.cols(); }
    inline unsigned rows() const { return m_image.rows(); }
    inline unsigned planes() const { return m_image.planes(); }

    inline pixel_accessor origin() const { return pixel_accessor(m_image.origin(),m_func); }

    inline reference_type operator()( int i, int j ) const { return m_func(m_image(i,j)); }
    inline reference_type operator()( int i, int j, int p ) const { return m_func(m_image(i,j,p)); }

    /// \cond INTERNAL
    typedef UnaryPerPixelView<typename ImageT::prerasterize_type, FuncT> prerasterize_type;
    inline prerasterize_type prerasterize() const
    { return prerasterize_type( m_image.prerasterize(), m_func ); }
    template <class DestT> inline void rasterize( DestT const& dest ) const { vw::rasterize( prerasterize(), dest ); }
    /// \endcond
  };
  
  /// \cond INTERNAL
  // View type Traits
  template <class ImageT, class FuncT>
  struct IsReferenceable<UnaryPerPixelView<ImageT,FuncT> > : boost::is_reference<typename UnaryPerPixelView<ImageT,FuncT>::reference_type>::type {};

  template <class ImageT, class FuncT>
  struct IsMultiplyAccessible<UnaryPerPixelView<ImageT,FuncT> > : boost::is_reference<typename UnaryPerPixelView<ImageT,FuncT>::reference_type>::type {};
  /// \endcond


  // *******************************************************************
  // UnaryPerPixelChannelView
  // *******************************************************************

  // Specialized Accessor
  template <class ImageIterT, class FuncT> 
  class UnaryPerPixelChannelAccessor
  {
  private:
    ImageIterT m_iter;
    FuncT const& m_func;
    typedef typename boost::result_of<FuncT(typename CompoundChannelType<typename ImageIterT::pixel_type>::type)>::type result_type;
  public:
    typedef typename CompoundChannelCast<typename ImageIterT::pixel_type, result_type>::type pixel_type;
    UnaryPerPixelChannelAccessor( ImageIterT const& iter, FuncT const& func ) : m_iter(iter), m_func(func) {}
    inline UnaryPerPixelChannelAccessor& next_col() { m_iter.next_col(); return *this; }
    inline UnaryPerPixelChannelAccessor& prev_col() { m_iter.prev_col(); return *this; }
    inline UnaryPerPixelChannelAccessor& next_row() { m_iter.next_row(); return *this; }
    inline UnaryPerPixelChannelAccessor& prev_row() { m_iter.prev_row(); return *this; }
    inline UnaryPerPixelChannelAccessor& next_plane() { m_iter.next_plane(); return *this; }
    inline UnaryPerPixelChannelAccessor& prev_plane() { m_iter.prev_plane(); return *this; }
    inline UnaryPerPixelChannelAccessor& advance( ptrdiff_t di, ptrdiff_t dj, ptrdiff_t dp=0 ) { m_iter.advance(di,dj,dp); return *this; }
    inline pixel_type operator*() const { return apply_per_pixel_channel(m_func,*m_iter); }
  };

  // Image View Class Declaration
  template <class ImageT, class FuncT>
  class UnaryPerPixelChannelView : public ImageViewBase<UnaryPerPixelChannelView<ImageT,FuncT> >
  {
  private:
    ImageT m_image;
    FuncT m_func;
  public:
    typedef typename CompoundResult<FuncT, typename ImageT::pixel_type>::type pixel_type;
    typedef UnaryPerPixelChannelAccessor<typename ImageT::pixel_accessor, FuncT> pixel_accessor;

    UnaryPerPixelChannelView( ImageT const& image, FuncT const& func ) : m_image(image), m_func(func) {}

    inline unsigned cols() const { return m_image.cols(); }
    inline unsigned rows() const { return m_image.rows(); }
    inline unsigned planes() const { return m_image.planes(); }

    inline pixel_accessor origin() const { return pixel_accessor(m_image.origin(),m_func); }

    inline pixel_type operator()( int i, int j ) const { return apply_per_pixel_channel(m_func,m_image(i,j)); }
    inline pixel_type operator()( int i, int j, int p ) const { return apply_per_pixel_channel(m_func,m_image(i,j,p)); }

    /// \cond INTERNAL
    typedef UnaryPerPixelChannelView<typename ImageT::prerasterize_type, FuncT> prerasterize_type;
    inline prerasterize_type prerasterize() const { return prerasterize_type( m_image.prerasterize(), m_func ); }
    template <class DestT> inline void rasterize( DestT const& dest ) const { vw::rasterize( prerasterize(), dest ); }
    /// \endcond
  };


  // *******************************************************************
  // BinaryPerPixelView
  // *******************************************************************

  // Specialized Accessor
  template <class Image1IterT, class Image2IterT, class FuncT>
  class BinaryPerPixelAccessor
  {
  private:
    Image1IterT m_iter1;
    Image2IterT m_iter2;
    FuncT const& m_func;
  public:
    typedef typename boost::result_of<FuncT(typename Image1IterT::pixel_type,typename Image2IterT::pixel_type)>::type reference_type;
    typedef typename boost::remove_reference<reference_type>::type pixel_type;

    BinaryPerPixelAccessor( Image1IterT const& iter1, Image2IterT const& iter2, FuncT const& func ) : m_iter1(iter1), m_iter2(iter2), m_func(func) {}
    inline BinaryPerPixelAccessor& next_col() { m_iter1.next_col(); m_iter2.next_col(); return *this; }
    inline BinaryPerPixelAccessor& prev_col() { m_iter1.prev_col(); m_iter2.prev_col(); return *this; }
    inline BinaryPerPixelAccessor& next_row() { m_iter1.next_row(); m_iter2.next_row(); return *this; }
    inline BinaryPerPixelAccessor& prev_row() { m_iter1.prev_row(); m_iter2.prev_row(); return *this; }
    inline BinaryPerPixelAccessor& next_plane() { m_iter1.next_plane(); m_iter2.next_plane(); return *this; }
    inline BinaryPerPixelAccessor& prev_plane() { m_iter1.prev_plane(); m_iter2.prev_plane(); return *this; }
    inline BinaryPerPixelAccessor& advance( ptrdiff_t di, ptrdiff_t dj, ptrdiff_t dp=0 ) 
      { m_iter1.advance(di,dj,dp); m_iter2.advance(di,dj,dp); return *this; }
    inline reference_type operator*() const { return m_func(*m_iter1,*m_iter2); }
  };

  /// \cond INTERNAL
  // Accessor type Traits
  template <class Image1IterT, class Image2IterT, class FuncT>
  struct IsReferenceable<BinaryPerPixelAccessor<Image1IterT,Image2IterT,FuncT> > : boost::is_reference<typename BinaryPerPixelAccessor<Image1IterT,Image2IterT,FuncT>::reference_type>::type {};
  /// \endcond

  // Image View Class Definition
  template <class Image1T, class Image2T, class FuncT>
  class BinaryPerPixelView : public ImageViewBase<BinaryPerPixelView<Image1T,Image2T,FuncT> >
  {
  private:
    Image1T m_image1;
    Image2T m_image2;
    FuncT m_func;
  public:
    typedef typename boost::result_of<FuncT(typename Image1T::pixel_type, typename Image2T::pixel_type)>::type reference_type;
    typedef typename boost::remove_reference<reference_type>::type pixel_type;

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

    inline reference_type operator()( int i, int j ) const
    { return m_func(m_image1(i,j),m_image2(i,j)); }

    inline reference_type operator()( int i, int j, int p ) const
    { return m_func(m_image1(i,j,p),m_image2(i,j,p)); }

    /// \cond INTERNAL
    typedef BinaryPerPixelView<typename Image1T::prerasterize_type, typename Image2T::prerasterize_type, FuncT> prerasterize_type;

    inline prerasterize_type prerasterize() const
    { return prerasterize_type( m_image1.prerasterize(), m_image2.prerasterize(), m_func ); }

    template <class DestT> inline void rasterize( DestT const& dest )
      const { vw::rasterize( prerasterize(), dest ); }
    /// \endcond
  };

  /// \cond INTERNAL
  // View type Traits
  template <class Image1T, class Image2T, class FuncT>
  struct IsReferenceable<BinaryPerPixelView<Image1T,Image2T,FuncT> > : boost::is_reference<typename BinaryPerPixelView<Image1T,Image2T,FuncT>::reference_type>::type {};

  template <class Image1T, class Image2T, class FuncT>
  struct IsMultiplyAccessible<BinaryPerPixelView<Image1T,Image2T,FuncT> > : boost::is_reference<typename BinaryPerPixelView<Image1T,Image2T,FuncT>::reference_type>::type {};
  /// \endcond


  // *******************************************************************
  // BinaryPerPixelChannelView
  // *******************************************************************

  // Specialized Accessor
  template <class Image1IterT, class Image2IterT, class FuncT>
  class BinaryPerPixelChannelAccessor
  {
    Image1IterT m_iter1;
    Image1IterT m_iter2;
    FuncT const& m_func;
  public:
    typedef typename CompoundResult<FuncT, typename Image1IterT::pixel_type, typename Image2IterT::pixel_type>::type pixel_type;
    BinaryPerPixelChannelAccessor( Image1IterT const& iter1, Image2IterT const& iter2, FuncT const& func ) : m_iter1(iter1), m_iter2(iter2), m_func(func) {}
    inline BinaryPerPixelChannelAccessor& next_col() { m_iter1.next_col(); m_iter2.next_col(); return *this; }
    inline BinaryPerPixelChannelAccessor& prev_col() { m_iter1.prev_col(); m_iter2.prev_col(); return *this; }
    inline BinaryPerPixelChannelAccessor& next_row() { m_iter1.next_row(); m_iter2.next_row(); return *this; }
    inline BinaryPerPixelChannelAccessor& prev_row() { m_iter1.prev_row(); m_iter2.prev_row(); return *this; }
    inline BinaryPerPixelChannelAccessor& next_plane() { m_iter1.next_plane(); m_iter2.next_plane(); return *this; }
    inline BinaryPerPixelChannelAccessor& prev_plane() { m_iter1.prev_plane(); m_iter2.prev_plane(); return *this; }
    inline BinaryPerPixelChannelAccessor& advance( ptrdiff_t di, ptrdiff_t dj, ptrdiff_t dp=0 ) { m_iter1.advance(di,dj,dp); m_iter2.advance(di,dj,dp); return *this; }
    inline pixel_type operator*() const { return apply_per_pixel_channel(m_func,*m_iter1,*m_iter2); }
  };

  // Image View Class Definition
  template <class Image1T, class Image2T, class FuncT>
  class BinaryPerPixelChannelView : public ImageViewBase<BinaryPerPixelChannelView<Image1T,Image2T,FuncT> > {
    Image1T m_image1;
    Image2T m_image2;
    FuncT m_func;
  public:
    typedef typename CompoundResult<FuncT, typename Image1T::pixel_type, typename Image2T::pixel_type>::type pixel_type;
    typedef BinaryPerPixelChannelAccessor<typename Image1T::pixel_accessor, typename Image2T::pixel_accessor, FuncT> pixel_accessor;

    BinaryPerPixelChannelView( Image1T const& image1, Image2T const& image2, FuncT const& func )
      : m_image1(image1), m_image2(image2), m_func(func)
    {
      VW_ASSERT( m_image1.cols()==m_image2.cols() &&
		 m_image1.rows()==m_image1.rows() &&
		 m_image1.planes()==m_image2.planes(),
                 ArgumentErr() << "BinaryPerPixelChannelView: Images must have same dimensions in binary image operation." );
    }
    
    inline unsigned cols() const { return m_image1.cols(); }
    inline unsigned rows() const { return m_image1.rows(); }
    inline unsigned planes() const { return m_image1.planes(); }
    
    inline pixel_accessor origin() const
    { return pixel_accessor(m_image1.origin(),m_image2.origin(),m_func); }
    
    inline pixel_type operator()( int i, int j ) const
    { return apply_per_pixel_channel(m_func,m_image1(i,j),m_image2(i,j)); }
    
    inline pixel_type operator()( int i, int j, int p ) const
    { return apply_per_pixel_channel(m_func,m_image1(i,j,p),m_image2(i,j,p)); }
    
    /// \cond INTERNAL
    typedef BinaryPerPixelChannelView<typename Image1T::prerasterize_type, typename Image2T::prerasterize_type, FuncT> prerasterize_type;
    
    inline prerasterize_type prerasterize() const
    { return prerasterize_type( m_image1.prerasterize(),
				m_image2.prerasterize(), m_func ); }
    
    template <class DestT> inline void rasterize( DestT const& dest ) const
    { vw::rasterize( prerasterize(), dest ); }
    /// \endcond
  };


  // *******************************************************************
  // Some simple per-pixel functions that don't belong anywhere else.
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

  /// A channel selecting functor, used by \ref select_channel().
  template <class ImageT>
  struct SelectChannelFunctor {
  private:
    int m_channel;
  public:
    SelectChannelFunctor( int channel ) : m_channel(channel) {}

    // Computes an appropriate reference-to-channel type.
    typedef typename CompoundChannelType<typename ImageT::pixel_type>::type base_channel_type;
    typedef typename boost::mpl::if_<boost::is_const<ImageT>,typename boost::add_const<base_channel_type>::type,base_channel_type>::type channel_type;
    typedef typename boost::mpl::if_<IsReferenceable<ImageT>,channel_type&,channel_type>::type result_type;

    template <class ArgT>
    result_type operator()( ArgT& pixel ) const {
      return compound_select_channel<result_type>(pixel,m_channel);
    }
  };

  /// Create a new image view by statically casting the channels of
  /// the pixels to a new type.
  template <class ImageT>
  UnaryPerPixelView<ImageT,SelectChannelFunctor<ImageT> >
  inline select_channel( ImageViewBase<ImageT>& image, int channel ) {
    return UnaryPerPixelView<ImageT,SelectChannelFunctor<ImageT> >( image.impl(), SelectChannelFunctor<ImageT>(channel) );
  }

  /// Create a new image view by statically casting the channels of
  /// the pixels to a new type (const overload).
  template <class ImageT>
  UnaryPerPixelView<ImageT,SelectChannelFunctor<const ImageT> >
  inline select_channel( ImageViewBase<ImageT> const& image, int channel ) {
    return UnaryPerPixelView<ImageT,SelectChannelFunctor<const ImageT> >( image.impl(), SelectChannelFunctor<const ImageT>(channel) );
  }

};

#endif // __VW_IMAGE__PER_PIXEL_VIEWS_H__
