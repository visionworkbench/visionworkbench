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

  /// Calls the passed in functor on every pixel LOCATION
  template<class FuncT>
  class PerPixelIndexView : public ImageViewBase<PerPixelIndexView<FuncT> > {
    FuncT m_func;
    int32 m_cols, m_rows, m_planes;
  public:
    typedef typename FuncT::result_type pixel_type;
    typedef pixel_type const            result_type;
    typedef ProceduralPixelAccessor<PerPixelIndexView> pixel_accessor;

    PerPixelIndexView( FuncT func, int32 cols, int32 rows, int32 planes = 1 )
      : m_func(func), m_cols(cols), m_rows(rows), m_planes(planes) {}

    inline int32 cols  () const { return m_cols;   }
    inline int32 rows  () const { return m_rows;   }
    inline int32 planes() const { return m_planes; }

    inline pixel_accessor origin() const { return pixel_accessor(*this); }

    inline result_type operator()( double i, double j, int32 p = 0 ) const {
      return m_func( i, j, p );
    }

    typedef PerPixelIndexView prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i const& /*bbox*/ ) const { return *this; }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const {
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
    ImageIterT   m_iter;
    FuncT const& m_func;
  public:
    typedef typename boost::result_of<FuncT(typename ImageIterT::pixel_type)>::type              result_type;
    typedef typename boost::remove_cv<typename boost::remove_reference<result_type>::type>::type pixel_type;
    typedef typename ImageIterT::offset_type offset_type;

    UnaryPerPixelAccessor( ImageIterT const& iter, FuncT const& func ) : m_iter(iter), m_func(func) {}
    inline UnaryPerPixelAccessor& next_col  () { m_iter.next_col();   return *this; }
    inline UnaryPerPixelAccessor& prev_col  () { m_iter.prev_col();   return *this; }
    inline UnaryPerPixelAccessor& next_row  () { m_iter.next_row();   return *this; }
    inline UnaryPerPixelAccessor& prev_row  () { m_iter.prev_row();   return *this; }
    inline UnaryPerPixelAccessor& next_plane() { m_iter.next_plane(); return *this; }
    inline UnaryPerPixelAccessor& prev_plane() { m_iter.prev_plane(); return *this; }
    inline UnaryPerPixelAccessor& advance( offset_type di, offset_type dj, ssize_t dp=0 ) { m_iter.advance(di,dj,dp); return *this; }
    inline result_type operator*() const { return m_func(*m_iter); }
  };

  /// Calls the passed in functor on every pixel VALUE
  template <class ImageT, class FuncT>
  class UnaryPerPixelView : public ImageViewBase<UnaryPerPixelView<ImageT,FuncT> > {
    ImageT m_image;
    FuncT  m_func;
  public:
    typedef typename boost::result_of<FuncT(typename ImageT::pixel_type)>::type                  result_type;
    typedef typename boost::remove_cv<typename boost::remove_reference<result_type>::type>::type pixel_type;
    typedef UnaryPerPixelAccessor<typename ImageT::pixel_accessor, FuncT>                        pixel_accessor;

    UnaryPerPixelView( ImageT const& image                    ) : m_image(image), m_func()     {}
    UnaryPerPixelView( ImageT const& image, FuncT const& func ) : m_image(image), m_func(func) {}

    inline int32 cols  () const { return m_image.cols();   }
    inline int32 rows  () const { return m_image.rows();   }
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
    inline prerasterize_type prerasterize( BBox2i const& bbox ) const { return prerasterize_type( m_image.prerasterize(bbox), m_func ); }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const { vw::rasterize( prerasterize(bbox), dest, bbox ); }
    /// \endcond
  };

  /// \cond INTERNAL
  // View type Traits.  This exists mainly for select_channel(), and may not
  // be correct in all cases.  Perhaps it should be specialized there instead?
  template <class ImageT, class FuncT>
  struct IsMultiplyAccessible<UnaryPerPixelView<ImageT,FuncT> > : 
      boost::is_reference<typename UnaryPerPixelView<ImageT,FuncT>::result_type>::type {};
  /// \endcond


  // *******************************************************************
  // BinaryPerPixelView
  // *******************************************************************

  // Specialized Accessor
  template <class Image1IterT, class Image2IterT, class FuncT>
  class BinaryPerPixelAccessor;

  // Image View Class Definition
  template <class Image1T, class Image2T, class FuncT>
  class BinaryPerPixelView;

  // *******************************************************************
  // TrinaryPerPixelView
  // *******************************************************************

  // Specialized Accessor
  template <class Image1IterT, class Image2IterT, class Image3IterT, class FuncT>
  class TrinaryPerPixelAccessor;

  // Image View Class Definition
  template <class Image1T, class Image2T, class Image3T, class FuncT>
  class TrinaryPerPixelView;


  // *******************************************************************
  // QuaternaryPerPixelView
  // *******************************************************************

  // Specialized Accessor
  template <class Image1IterT, class Image2IterT, class Image3IterT, class Image4IterT, class FuncT>
  class QuaternaryPerPixelAccessor;

  // Image View Class Definition
  template <class Image1T, class Image2T, class Image3T, class Image4T, class FuncT>
  class QuaternaryPerPixelView;



  // *******************************************************************
  // UnaryPerPixelIndexView
  // *******************************************************************

  /// Class for applying a functor to every pixel in an image.
  /// - This is nearly identical to UnaryPerPixelView but it passes the
  ///   pixel coordinates to the function in addition to the pixel value.
  ///   Basically a combination of UnaryPerPixelView and PerPixelIndexView.
  template <class ImageT, class FuncT>
  class UnaryPerPixelIndexView : public ImageViewBase<UnaryPerPixelIndexView<ImageT,FuncT> > {
    ImageT m_image;
    FuncT  m_func;
  public:
    typedef typename boost::result_of<FuncT(typename ImageT::pixel_type, int32, int32, int32)>::type result_type;
    typedef typename boost::remove_cv<typename boost::remove_reference<result_type>::type>::type     pixel_type;
    typedef ProceduralPixelAccessor<UnaryPerPixelIndexView>                                          pixel_accessor;

    UnaryPerPixelIndexView( ImageT const& image                    ) : m_image(image), m_func() {}
    UnaryPerPixelIndexView( ImageT const& image, FuncT const& func ) : m_image(image), m_func(func) {}

    inline int32 cols  () const { return m_image.cols  (); }
    inline int32 rows  () const { return m_image.rows  (); }
    inline int32 planes() const { return m_image.planes(); }

    //inline pixel_accessor origin() const { return pixel_accessor(m_image.origin(), m_func); }
    inline pixel_accessor origin() const { return pixel_accessor(*this); }
    inline result_type operator()( int32 i, int32 j, int32 p=0 ) const { return m_func(m_image(i,j,p), i, j, p); }

    template <class ViewT>
    UnaryPerPixelIndexView& operator=( ImageViewBase<ViewT> const& view ) {
      view.impl().rasterize( *this, BBox2i(0,0,view.impl().cols(),view.impl().rows()) );
      return *this;
    }

    /// \cond INTERNAL
    typedef UnaryPerPixelIndexView<typename ImageT::prerasterize_type, FuncT> prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i const& bbox ) const { return prerasterize_type( m_image.prerasterize(bbox), m_func ); }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const { vw::rasterize( prerasterize(bbox), dest, bbox ); }
    /// \endcond
  }; // End class UnaryPerPixelIndexView

  /// \cond INTERNAL
  // View type Traits.  This exists mainly for select_channel(), and may not
  // be correct in all cases.  Perhaps it should be specialized there instead?
  template <class ImageT, class FuncT>
  struct IsMultiplyAccessible<UnaryPerPixelIndexView<ImageT,FuncT> > : 
      boost::is_reference<typename UnaryPerPixelIndexView<ImageT,FuncT>::result_type>::type {};
  /// \endcond

//=================================================================================================
// Here are some convenience functions to use



  template <class Image1T, class Functor>
  UnaryPerPixelView<Image1T, Functor>
  per_pixel_view(Image1T const& image1, Functor const& func) {
    typedef UnaryPerPixelView<Image1T, Functor> result_type;
    return result_type(image1, func);
  }

  template <class Image1T, class Image2T, class Functor>
  BinaryPerPixelView<Image1T, Image2T, Functor>
  per_pixel_view(Image1T const& image1, Image2T const& image2, Functor const& func) {
    typedef BinaryPerPixelView<Image1T, Image2T, Functor> result_type;
    return result_type(image1, image2, func);
  }

  template <class Image1T, class Image2T, class Image3T, class Functor>
  TrinaryPerPixelView<Image1T, Image2T, Image3T, Functor>
  per_pixel_view(Image1T const& image1, Image2T const& image2, 
                 Image3T const& image3, Functor const& func) {
    typedef TrinaryPerPixelView<Image1T, Image2T, Image3T, Functor> result_type;
    return result_type(image1, image2, image3, func);
  }

  template <class Image1T, class Image2T, class Image3T, class Image4T, class Functor>
  QuaternaryPerPixelView<Image1T, Image2T, Image3T, Image4T, Functor>
  per_pixel_view(Image1T const& image1, Image2T const& image2, 
                 Image3T const& image3, Image4T const& image4, Functor const& func) {
    typedef QuaternaryPerPixelView<Image1T, Image2T, Image3T, Image4T, Functor> result_type;
    return result_type(image1, image2, image3, image4, func);
  }



} // End namespace vw

#include "PerPixelViews.tcc"

#endif // __VW_IMAGE_PERPIXELVIEWS_H__
