// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file Interpolation.h
///
/// Image views that can be accessed with real (floating point) pixel
/// indices.  These image views will perform one of:
///
/// - bilinear interpolation       ( BilinearInterpolation()      )
/// - bicubic interpolation        ( BicubicInterpolation()       )
/// - nearest pixel interpolation  ( NearestPixelInterpolaiton()  )
///   (unless the underlying image view can also be accessed using real pixel values)
///
#ifndef __VW_IMAGE_INTERPOLATION_H__
#define __VW_IMAGE_INTERPOLATION_H__

#include <boost/type_traits.hpp>
#include <boost/mpl/logical.hpp>
#include <boost/utility/enable_if.hpp>

#include <vw/Image/ImageView.h>
#include <vw/Image/EdgeExtension.h>
#include <vw/Image/Manipulation.h>
#include <vw/Image/SparseImageCheck.h>

namespace vw {

  /// \cond INTERNAL
  // Stub classes defining common interpolation modes.  You may define
  // your own class similar to those that appear below.  The stub class
  // must functor must have an Interpolator type function that returns
  // the type of the interpolation functor, and a static interpolator()
  // function that returns the interpolator itself.  For simple
  // functions, the stub class and the interpolation function may be
  // the same class (e.g. see NearestPixelInterpolation).  This extra
  // level of indirection makes it possible to do more flexible
  // template-based special-casing for optimization.

  /// A base class for interpolation functors that provides the
  /// common return type deduction logic in case users want to use
  /// these types in a more general manner.
  ///
  /// pixel_buffer is the number of additional pixels to prerasterize
  /// along the edge of the child image (in case we do need to
  /// rasterize the child when pre-rasterize is called). The subclass
  /// of InterpolationBase _must_ override this value and set it to
  /// the number of pixels that the interpolation algorithm will need
  /// to search outside the boundaries of the image on each side.
  struct InterpolationBase {
    static const int32 pixel_buffer = 0;
    template <class ArgsT> struct result {};
    template <class FuncT, class ViewT, class IT, class JT, class PT>
    struct result<FuncT(ViewT,IT,JT,PT)> {
      typedef typename boost::remove_reference<ViewT>::type::pixel_type type;
    };
  };

  // This is broken out so that the implementation can be overridden
  // by pixel type.  Optimized versions go at the bottom of the file
  // for clarity.
  template <class ViewT, class PixelT = typename ViewT::pixel_type>
  struct BilinearInterpolationImpl : InterpolationBase {
    PixelT operator()( const ViewT &view, double i, double j, int32 p ) const {
      typedef typename ViewT::pixel_type pixel_type;
      typedef typename CompoundChannelType<pixel_type>::type channel_type;
      typedef typename FloatType<channel_type>::type real_type;
      typedef typename CompoundChannelCast<pixel_type,real_type>::type result_type;

      int32 x = math::impl::_floor(i), y = math::impl::_floor(j);
      real_type normx = real_type(i)-real_type(x), normy = real_type(j)-real_type(y), norm1mx = 1-normx, norm1my = 1-normy;

      typename ViewT::pixel_accessor acc = view.origin().advance(x,y,p);
      result_type result = (*acc) * norm1mx;
      acc.next_col();
      result += (*acc) * normx;
      result *= norm1my;
      acc.advance(-1,1);
      result_type row = (*acc) * norm1mx;
      acc.next_col();
      row += (*acc) * normx;
      result += row * normy;

      // Linear interpolation is, well, linear, so there's no need to clamp.
      return channel_cast_round_if_int<channel_type>(result);
    }
  };

  // Bilinear interpolation operator
  struct BilinearInterpolation {
    static const int32 pixel_buffer = 1;
    template <class ViewT>
    struct Interpolator {
      typedef BilinearInterpolationImpl<ViewT> type;
    };
    template <class ViewT>
    static typename Interpolator<ViewT>::type interpolator( const ViewT& /*view*/ ) {
      return typename Interpolator<ViewT>::type();
    }
    // This function is here for backwards-compatibility and is deprecated.
    template <class ViewT>
    inline typename ViewT::pixel_type operator()( const ViewT &view, double i, double j, int32 p ) const VW_DEPRECATED;
  };

  template <class ViewT>
  inline typename ViewT::pixel_type BilinearInterpolation::operator()( const ViewT &view, double i, double j, int32 p ) const {
    return interpolator(view)( view, i, j, p );
  }


  // This is broken out so that the implementation can be overridden
  // by pixel type.  Optimized versions go at the bottom of the file
  // for clarity.
  template <class ViewT, class PixelT = typename ViewT::pixel_type>
  struct BicubicInterpolationImpl {
    PixelT operator()( const ViewT &view, double i, double j, int32 p ) const {
      typedef typename CompoundChannelType<PixelT>::type channel_type;
      typedef typename CompoundChannelCast<PixelT,double>::type result_type;

      int32 x = math::impl::_floor(i), y = math::impl::_floor(j);
      double normx = i-x, normy = j-y;

      double s0 = ((2-normx)*normx-1)*normx;      double t0 = ((2-normy)*normy-1)*normy;
      double s1 = (3*normx-5)*normx*normx+2;      double t1 = (3*normy-5)*normy*normy+2;
      double s2 = ((4-3*normx)*normx+1)*normx;    double t2 = ((4-3*normy)*normy+1)*normy;
      double s3 = (normx-1)*normx*normx;          double t3 = (normy-1)*normy*normy;

      typename ViewT::pixel_accessor acc = view.origin().advance(x-1,y-1,p);
      result_type row =         s0*(*acc);
      acc.next_col();    row += s1*(*acc);
      acc.next_col();    row += s2*(*acc);
      acc.next_col();    row += s3*(*acc);
      result_type result =      t0*row;
      acc.advance(-3,1); row =  s0*(*acc);
      acc.next_col();    row += s1*(*acc);
      acc.next_col();    row += s2*(*acc);
      acc.next_col();    row += s3*(*acc);
      result +=                 t1*row;
      acc.advance(-3,1); row =  s0*(*acc);
      acc.next_col();    row += s1*(*acc);
      acc.next_col();    row += s2*(*acc);
      acc.next_col();    row += s3*(*acc);
      result +=                 t2*row;
      acc.advance(-3,1); row =  s0*(*acc);
      acc.next_col();    row += s1*(*acc);
      acc.next_col();    row += s2*(*acc);
      acc.next_col();    row += s3*(*acc);
      result +=                 t3*row;
      result *= 0.25;

      // Bicubic interpolation is non-convex, so we must clamp if integer
      return channel_cast_round_and_clamp_if_int<channel_type>( result );
    }
  };

  // Bicubic interpolation operator
  struct BicubicInterpolation {
    static const int32 pixel_buffer = 2;
    template <class ViewT>
    struct Interpolator {
      typedef BicubicInterpolationImpl<ViewT> type;
    };
    template <class ViewT>
    static typename Interpolator<ViewT>::type interpolator( const ViewT & /*view*/ ) {
      return typename Interpolator<ViewT>::type();
    }
    // This function is here for backwards-compatibility and is deprecated.
    template <class ViewT>
    inline typename ViewT::pixel_type operator()( const ViewT &view, double i, double j, int32 p ) const VW_DEPRECATED;
  };

  template <class ViewT>
  inline typename ViewT::pixel_type BicubicInterpolation::operator()( const ViewT &view, double i, double j, int32 p ) const {
    return interpolator(view)( view, i, j, p );
  }

  // NearestPixel interpolation operator.
  struct NearestPixelInterpolation {
    static const int32 pixel_buffer = 1;
    template <class ViewT>
    struct Interpolator {
      typedef NearestPixelInterpolation type;
    };
    template <class ViewT>
    static NearestPixelInterpolation interpolator( const ViewT & /*view*/ ) {
      return NearestPixelInterpolation();
    }
    template <class ViewT>
    typename ViewT::pixel_type operator()( const ViewT &view, double i, double j, int32 p ) const {
      int32 x = math::impl::_round(i), y = math::impl::_round(j);
      return view(x,y,p);
    }
  };


  /// Interpolation View Class
  ///
  /// An image view that excepts real numbers as pixel coordinates and
  /// interpolates the value of the image at these coordinates using
  /// some interpolation method.  For pixels that fall outside the range
  /// of the wrapped image view, an newly constructed (empty) pixeltype
  /// is returned.
  template <class ImageT, class InterpT>
  class InterpolationView : public ImageViewBase<InterpolationView<ImageT, InterpT> >
  {
  private:
    typedef typename InterpT::template Interpolator<ImageT>::type interp_type;
    ImageT m_image;
    interp_type m_interp_func;
  public:

    typedef typename ImageT::pixel_type pixel_type;
    typedef pixel_type result_type;
    typedef ProceduralPixelAccessor<InterpolationView<ImageT, InterpT> > pixel_accessor;

    InterpolationView( ImageT const& image,
                       InterpT const& /*interp_stub*/ = InterpT()) :
      m_image(image), m_interp_func(InterpT::interpolator(image)) {}

    inline int32 cols() const { return m_image.cols(); }
    inline int32 rows() const { return m_image.rows(); }
    inline int32 planes() const { return m_image.planes(); }

    inline pixel_accessor origin() const { return pixel_accessor(*this, 0, 0); }

    inline result_type operator() (double i, double j, int32 p = 0) const { return m_interp_func(m_image,i,j,p); }

    ImageT const& child() const { return m_image; }
    InterpT const& func() const { return m_interp_func; }

    /// \cond INTERNAL
    // We can make an optimization here.  If the pixels in the child
    // view cannot be repeatedly accessed without incurring any
    // additional overhead  then we should rasterize the child
    // before we proceed to rasterize ourself.
    typedef typename boost::mpl::if_< IsMultiplyAccessible<ImageT>,
                                      InterpolationView<typename ImageT::prerasterize_type, InterpT>,
                                      InterpolationView<CropView<ImageView<pixel_type> >, InterpT> >::type prerasterize_type;

    template <class PreRastImageT>
    prerasterize_type prerasterize_helper( BBox2i bbox, PreRastImageT const& image, true_type ) const {
      return prerasterize_type( image.prerasterize(bbox) );
    }

    template <class PreRastImageT>
    prerasterize_type prerasterize_helper( BBox2i bbox, PreRastImageT const& image, false_type ) const {
      ImageView<pixel_type> buf( bbox.width(), bbox.height(), m_image.planes() );
      image.rasterize( buf, bbox );
      return prerasterize_type( CropView<ImageView<pixel_type> >( buf, BBox2i(-bbox.min().x(),-bbox.min().y(),
                                                                              image.cols(), image.rows())));
    }

    inline prerasterize_type prerasterize( BBox2i bbox ) const {
      int32 padded_width = bbox.width() + 2 * InterpT::pixel_buffer;
      int32 padded_height = bbox.height() + 2 * InterpT::pixel_buffer;
      BBox2i adjusted_bbox(bbox.min().x() - InterpT::pixel_buffer,
                           bbox.min().y() - InterpT::pixel_buffer,
                           padded_width, padded_height);
      return prerasterize_helper(adjusted_bbox, m_image, typename IsMultiplyAccessible<ImageT>::type() );
    }

    template <class DestT> inline void rasterize( DestT const& dest, BBox2i bbox ) const { vw::rasterize( prerasterize(bbox), dest, bbox ); }
    /// \endcond
  };

  /// \cond INTERNAL
  // Type traits
  template <class ImageT, class InterpT>
  struct IsFloatingPointIndexable<InterpolationView<ImageT, InterpT> > : public true_type {};
  /// \endcond

  template <class ImageT, class InterpT>
  class SparseImageCheck<InterpolationView<ImageT, InterpT> > {
    InterpolationView<ImageT, InterpT> const& m_view;
  public:
    SparseImageCheck(InterpolationView<ImageT, InterpT> const& view)
      : m_view(view) {}
    bool operator()( BBox2i const& bbox ) const {
      if( bbox.empty() ) return false;
      BBox2i src_bbox = bbox;
      src_bbox.expand( InterpT::pixel_buffer );
      return SparseImageCheck<ImageT>(m_view.child())( src_bbox );
    }
  };

  // -------------------------------------------------------------------------------
  // Functional API
  // -------------------------------------------------------------------------------

  /// Use this free function to pass in an arbitrary interpolation
  /// functor.  You can use of the predefined functors at the top of
  /// this file or use one of your own devising.
  ///
  /// This version of interpolate takes an extra argument, the edge
  /// extension functor, and it automatically edge extends the image
  /// before interpolating.  See EdgeExtension.h for a list of built-in
  /// functors.
  template <class ImageT, class InterpT, class EdgeExtensionT>
  InterpolationView<EdgeExtensionView<ImageT,EdgeExtensionT>, InterpT> interpolate( ImageViewBase<ImageT> const& v,
                                                                                    InterpT const& interp_func,
                                                                                    EdgeExtensionT const& edge_extend_func) {
    return InterpolationView<EdgeExtensionView<ImageT, EdgeExtensionT>, InterpT>( edge_extend(v, edge_extend_func) , interp_func );
  }

  /// Use this free function to pass in an arbitrary interpolation
  /// functor.  You can use of the predefined functors at the top of
  /// this file or even one of your own devising.
  ///
  /// This version of the interpolation function uses Constant edge
  /// extension by default.
  template <class ImageT, class InterpT> InterpolationView<EdgeExtensionView<ImageT, ConstantEdgeExtension>, InterpT>
  interpolate( ImageViewBase<ImageT> const& v, InterpT const& interp_func) {
    return InterpolationView<EdgeExtensionView<ImageT, ConstantEdgeExtension>, InterpT>( edge_extend(v, ConstantEdgeExtension()), interp_func );
  }

  /// Use this free function to pass in an arbitrary interpolation
  /// functor.  You can use of the predefined functors at the top of
  /// this file or even one of your own devising.
  ///
  /// This version of the interpolation function uses Constant edge
  /// extension by default.
  template <class ImageT> InterpolationView<EdgeExtensionView<ImageT, ConstantEdgeExtension>, BilinearInterpolation>
  interpolate( ImageViewBase<ImageT> const& v ) {
    return InterpolationView<EdgeExtensionView<ImageT, ConstantEdgeExtension>, BilinearInterpolation>( edge_extend(v, ConstantEdgeExtension()), BilinearInterpolation() );
  }

} // namespace vw


// -------------------------------------------------------------------------------
// SSE Optimizations
// -------------------------------------------------------------------------------

#if defined(VW_ENABLE_SSE) && (VW_ENABLE_SSE==1)
#include <xmmintrin.h>

namespace vw {

  extern const float bicubic_coeffs[16] __attribute__ ((aligned (16)));

  template <class ViewT>
  struct BicubicInterpolationImpl<ViewT, PixelRGBA<uint8> > {
    PixelRGBA<uint8> operator()( const ViewT &view, double i, double j, int32 p ) const {
      // Front matter: 1.3% of samples.
      typedef typename ViewT::pixel_accessor acc_type;
      typedef float v4f[4] __attribute__ ((aligned (16)));

      // Compute integers and normalized offsets.  7.5% of samples.
      int32 x = math::impl::_floor(i), y = math::impl::_floor(j);
      float normx = i-x, normy = j-y;

      // Compute the bicubic coefficients.  The next three bits get
      // blurred together a bit by GCC, but total 16.9% of samples.
      __m128 a, b, c, d, s0, s1, s2, s3;
      a = _mm_set_ps1( normx );
      b = _mm_set_ps1( normy );
      s0 = _mm_load_ps( bicubic_coeffs );
      s1 = _mm_load_ps( bicubic_coeffs+4 );
      s2 = _mm_load_ps( bicubic_coeffs+8 );
      s3 = _mm_load_ps( bicubic_coeffs+12 );
      c = _mm_mul_ps( a, s0 );
      d = _mm_mul_ps( b, s0 );
      c = _mm_add_ps( c, s1 );
      d = _mm_add_ps( d, s1 );
      c = _mm_mul_ps( c, a );
      d = _mm_mul_ps( d, b );
      c = _mm_add_ps( c, s2 );
      d = _mm_add_ps( d, s2 );
      c = _mm_mul_ps( c, a );
      d = _mm_mul_ps( d, b );
      c = _mm_add_ps( c, s3 );
      d = _mm_add_ps( d, s3 );

      // Move the coefficients into place.
      v4f tmp;
      _mm_store_ps(tmp,c);
      s0 = _mm_set_ps1(tmp[0]);
      s1 = _mm_set_ps1(tmp[1]);
      s2 = _mm_set_ps1(tmp[2]);
      s3 = _mm_set_ps1(tmp[3]);
      _mm_store_ps(tmp,d);

      // Get ready to loop over the source pixels.
      acc_type acc = view.origin().advance(x-1,y-1);
      v4f pixel1, pixel2, pixel3, pixel4;

      // Loop over the source rows, accumulating the result;
      d = _mm_set_ps1( 0.0f );
      for( int i=0; i<4; ++i ) {
        // GCC does a decent job with this.  17.8% of samples.
        *(PixelRGBA<float>*)(pixel1) = *acc;  acc.next_col();
        *(PixelRGBA<float>*)(pixel2) = *acc;  acc.next_col();
        *(PixelRGBA<float>*)(pixel3) = *acc;  acc.next_col();
        *(PixelRGBA<float>*)(pixel4) = *acc;  acc.advance(-3,1);

        // Multiply-and-add one row of pixels. 45.7% of samples.
        a = _mm_load_ps(pixel1);
        b = _mm_load_ps(pixel2);
        c = _mm_mul_ps(a, s0);
        b = _mm_mul_ps(b, s1);
        c = _mm_add_ps(c, b);
        a = _mm_load_ps(pixel3);
        b = _mm_load_ps(pixel4);
        a = _mm_mul_ps(a, s2);
        b = _mm_mul_ps(b, s3);
        c = _mm_add_ps(c, a);
        c = _mm_add_ps(c, b);
        a = _mm_set_ps1(tmp[i]);
        c = _mm_mul_ps(c, a);
        d = _mm_add_ps(d, c);
      }

      // Clamp the values to the range of uint8, and add an offset of 0.5
      // so that truncation results in rounding.  2.0% of samples.
      a = _mm_set_ps1( 0.0f );
      b = _mm_set_ps1( 255.0f );
      c = _mm_set_ps1( 0.5f );
      d = _mm_max_ps( d, a );
      d = _mm_min_ps( d, b );
      d = _mm_add_ps( d, c );

      // Move the result out of SSE-land.  1.0% of samples.
      _mm_store_ps(tmp,d);

      // Truncate to packed uint8.  6.2% of samples.
      return PixelRGBA<uint8>( *(PixelRGBA<float>*)(tmp) );

      // Wrap-up.  1.6% of samples.
    }
  };

} // namespace vw

#endif // VW_ENABLE_SSE

#endif // __VW_IMAGE_INTERPOLATION_H__
