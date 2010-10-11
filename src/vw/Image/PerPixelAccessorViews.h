// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file PerPixelAccessorViews.h
///
/// Per-pixel accessor image views.
///
/// To use the PerPixelAccessorView, you write a unary functor that
/// takes a PixelAccessor as the argument to its call operator, and
/// returns the PixelAccessor's pixel_type as the result.  The
/// returned value can be a function computed by moving the accessor
/// around this region of the image.
///
/// In order to inform the PerPixelAccessorView of the correct level
/// of padding for the edges, your functor must also provide a method
/// called work_space() that provides a bounding box indicating the
/// extent of the region around this pixel that may be accessed by the
/// pixel accessor.
///
/// An example functor would look something like this:
///
///  class ExampleFunc : public UnaryReturnTemplateType<PixelTypeFromPixelAccessor>
///  {
///    BBox2i work_area() const { ... }
///
///    template <class PixelAccessorT>
///    typename PixelAccessorT::pixel_type
///    operator() (PixelAccessorT const& acc) const { ... }
///  };
///
#ifndef __VW_IMAGE_PERPIXELACCESSORVIEWS_H__
#define __VW_IMAGE_PERPIXELACCESSORVIEWS_H__

#include <boost/type_traits.hpp>
#include <boost/mpl/logical.hpp>
#include <boost/mpl/if.hpp>
#include <boost/utility/result_of.hpp>

#include <vw/Image/ImageViewBase.h>
#include <vw/Image/Manipulation.h>
#include <vw/Image/EdgeExtension.h>

namespace vw {

  // *******************************************************************
  // UnaryPerPixelAccessorView
  // *******************************************************************

  // Image View Class Declaration
  template <class ImageT, class FuncT>
  class UnaryPerPixelAccessorView : public ImageViewBase<UnaryPerPixelAccessorView<ImageT,FuncT> > {
    ImageT m_image;
    FuncT m_func;
  public:

    typedef ProceduralPixelAccessor<UnaryPerPixelAccessorView> pixel_accessor;
    typedef typename boost::result_of<FuncT(typename ImageT::pixel_accessor)>::type result_type;
    typedef typename boost::remove_reference<result_type>::type pixel_type;

    UnaryPerPixelAccessorView( ImageT const& image ) : m_image(image), m_func() {}
    UnaryPerPixelAccessorView( ImageT const& image, FuncT const& func ) : m_image(image), m_func(func) {}

    inline int32 cols() const { return m_image.cols(); }
    inline int32 rows() const { return m_image.rows(); }
    inline int32 planes() const { return m_image.planes(); }

    inline pixel_accessor origin() const { return pixel_accessor(*this); }
    inline result_type operator()( int32 i, int32 j, int32 p=0 ) const { return m_func(m_image.origin().advance(i,j,p)); }

    /// \cond INTERNAL
    // We can make an optimization here.  If the pixels in the child
    // view cannot be repeatedly accessed without incurring any
    // additional overhead  then we should rasterize the child
    // before we proceed to rasterize ourself.
    BBox2i pad_bbox(BBox2i const& bbox) const {
      int32 padded_width = bbox.width() + m_func.work_area().width();
      int32 padded_height = bbox.height() + m_func.work_area().height();
      return BBox2i (bbox.min().x() + m_func.work_area().min().x(),
                     bbox.min().y() + m_func.work_area().min().y(),
                     padded_width, padded_height);
    }

    typedef typename boost::mpl::if_< IsMultiplyAccessible<ImageT>,
                                      UnaryPerPixelAccessorView<typename ImageT::prerasterize_type, FuncT>,
                                      UnaryPerPixelAccessorView<CropView<ImageView<typename ImageT::pixel_type> >, FuncT > >::type prerasterize_type;

    template <class PreRastImageT>
    prerasterize_type prerasterize_helper( BBox2i bbox, PreRastImageT const& image, true_type ) const {
      return prerasterize_type( image.prerasterize(this->pad_bbox(bbox)), m_func );
    }

    template <class PreRastImageT>
    prerasterize_type prerasterize_helper( BBox2i bbox, PreRastImageT const& image, false_type ) const {
      BBox2i adjusted_bbox = this->pad_bbox(bbox);
      ImageView<typename ImageT::pixel_type> buf( adjusted_bbox.width(), adjusted_bbox.height(), m_image.planes() );
      m_image.rasterize( buf, adjusted_bbox );
      return prerasterize_type( CropView<ImageView<typename ImageT::pixel_type> >( buf, BBox2i(-adjusted_bbox.min().x(), -adjusted_bbox.min().y(),
                                                                                               image.cols(), image.rows())), m_func);
    }

    inline prerasterize_type prerasterize( BBox2i bbox ) const {
      return prerasterize_helper(bbox, m_image, typename IsMultiplyAccessible<ImageT>::type() );
    }

    template <class DestT> inline void rasterize( DestT const& dest, BBox2i bbox ) const { vw::rasterize( prerasterize(bbox), dest, bbox ); }
    /// \endcond
  };

  template <class ViewT, class FuncT, class EdgeT>
  UnaryPerPixelAccessorView<EdgeExtensionView<ViewT,EdgeT>, FuncT>
  per_pixel_accessor_filter(ImageViewBase<ViewT> const& image, FuncT const& func, EdgeT edge) {
    return UnaryPerPixelAccessorView<EdgeExtensionView<ViewT,EdgeT>, FuncT>(edge_extend(image.impl(), edge), func);
  }

  template <class ViewT, class FuncT>
  UnaryPerPixelAccessorView<EdgeExtensionView<ViewT,ZeroEdgeExtension>, FuncT>
  per_pixel_accessor_filter(ImageViewBase<ViewT> const& image, FuncT const& func) {
    return UnaryPerPixelAccessorView<EdgeExtensionView<ViewT,ZeroEdgeExtension>, FuncT>(edge_extend(image.impl(), ZeroEdgeExtension()), func);
  }

} // namespace vw

#endif // __VW_IMAGE_PERPIXELACCESSORVIEWS_H__
