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

/// \file Interpolation.h
/// 
/// Image views that can be accessed with real (floating point) pixel
/// indices.  These image views will perform one of:
/// 
/// - bilinear interpolation       ( vw::interpolation::Bilinear()      )
/// - bicubic interpolation        ( vw::interpolation::Bicubic()       )
/// - nearest pixel interpolation  ( vw::interpolation::NearestPixel()  )
///   (unless the underlying image view can also be accessed using real pixel values)
///
#ifndef __VW_INTERPOLATION_H__
#define __VW_INTERPOLATION_H__

#include <boost/type_traits.hpp>
#include <boost/mpl/logical.hpp>
#include <boost/utility/enable_if.hpp>

#include <vw/Image/ImageView.h>
#include <vw/Image/PixelAccessors.h>
#include <vw/Image/EdgeExtend.h>

namespace vw {

  /// The namespace for interpolation modes.  We define the
  /// interpolation modes as placeholder classes that may be passed to
  /// the interpolation routines as dummy variables.

  /// Bilinear interpolation between the four nearest pixels
  class BilinearInterpolation {};
  /// Bicubic interpolation between the nine nearest pixels
  class BicubicInterpolation {};
  /// Interpolate by using the value of the nearest pixel
  class NearestPixelInterpolation {};

	/// \cond INTERNAL
	// Abstract "Base" template for interpolation methods
	//
  // The logic for determining the type of the interpolation to use in
  // the InterpolationView class activates the appropriate method
  // using the dummy classes above.  You can easily extend the list of
  // supported interpolation modes by creating a new dummy type in the
  // namespace interpolation, and then define a specialization of
  // InterpolationImplementation similar to those below with a method
  // called interpolate(i,j,p,view).  Be certain that the view
  // argument is passed by reference in this function.
  template <class InterpT, class ViewT>
  struct InterpolationImplementation {};

  // Bilinear interpolation operator
  template <class ViewT>
  struct InterpolationImplementation<BilinearInterpolation,ViewT> {
    static inline typename ViewT::pixel_type interpolate(const ViewT &view, float i, float j, unsigned p ) { 
      typedef typename ViewT::pixel_type pixel_type;

      int x = int(floor(i));       int y = int(floor(j));
      double normx = i-x;          double normy = j - y;

      pixel_type i1( view(x,y,p)   + (view(x,y+1,p)   - view(x,y,p)  ) * normy );
      pixel_type i2( view(x+1,y,p) + (view(x+1,y+1,p) - view(x+1,y,p)) * normy );
      return pixel_type(i1 + (i2 - i1) * normx); 
    }
  };

  // Bicubic interpolation operator
  template <class ViewT>
  struct InterpolationImplementation<BicubicInterpolation,ViewT> {
    static inline typename ViewT::pixel_type interpolate( const ViewT &view, float i, float j, unsigned p ) { 
      typedef typename ViewT::pixel_type pixel_type;
      
      int x = int(floor(i));       int y = int(floor(j));
      double normx = i-x;          double normy = j - y;
  
      // This code was taken almost verbatim from vil_bicub_interp.txx
      // in the VXL source tree.
      double s0 = ((2-normx)*normx-1)*normx;      double t0 = ((2-normy)*normy-1)*normy;
      double s1 = (3*normx-5)*normx*normx+2;      double t1 = (3*normy-5)*normy*normy+2;
      double s2 = ((4-3*normx)*normx+1)*normx;    double t2 = ((4-3*normy)*normy+1)*normy;
      double s3 = (normx-1)*normx*normx;          double t3 = (normy-1)*normy*normy;

      pixel_type xi0( s0*view(x-1,y-1,p) + s1*view(x+0,y-1,p) + s2*view(x+1,y-1,p) + s3*view(x+2,y-1,p) );
      pixel_type xi1( s0*view(x-1,y+0,p) + s1*view(x+0,y+0,p) + s2*view(x+1,y+0,p) + s3*view(x+2,y+0,p) );
      pixel_type xi2( s0*view(x-1,y+1,p) + s1*view(x+0,y+1,p) + s2*view(x+1,y+1,p) + s3*view(x+2,y+1,p) );
      pixel_type xi3( s0*view(x-1,y+2,p) + s1*view(x+0,y+2,p) + s2*view(x+1,y+2,p) + s3*view(x+2,y+2,p) );
      return pixel_type(0.25 * ( xi0*t0 + xi1*t1 + xi2*t2 + xi3*t3 ));
    }
  };

  // NearestPixel interpolation operator.  
  template <class ViewT>
  struct InterpolationImplementation<NearestPixelInterpolation, ViewT> {
    static inline typename ViewT::pixel_type interpolate( const ViewT &view, float i, float j, unsigned p ) {
      int x = int(lroundf(i));       int y = int(lroundf(j));
      return view(x,y,p);
    }
  };
  /// \endcond


  /// Interpolation View Class
  ///
  /// An image view that excepts real numbers as pixel coordinates and
  /// interpolates the value of the image at these coordinates using
  /// some interpolation method.  For pixels that fall outside the range
  /// of the wrapped image view, an newly constructed (empty) pixeltype
  /// is returned.
  template <class ImageT, class InterpT = BilinearInterpolation, class EdgeT = ZeroEdgeExtend>
  class InterpolationView : public ImageViewBase<InterpolationView<ImageT, InterpT, EdgeT> >
  {
  private:
    ImageT m_image;
    EdgeExtendView<ImageT,EdgeT> m_extend;
  public:

    typedef typename boost::remove_cv<typename ImageT::pixel_type>::type base_pixel_type;
    typedef const base_pixel_type pixel_type;
    typedef ProceduralPixelAccessor<InterpolationView<ImageT, InterpT, EdgeT> > pixel_accessor;
    
    InterpolationView( ImageT const& image ) : m_image(image), m_extend(image) {}

    inline unsigned cols() const { return m_image.cols(); }
    inline unsigned rows() const { return m_image.rows(); }
    inline unsigned planes() const { return m_image.planes(); }

    inline pixel_accessor origin() const { return pixel_accessor(*this, 0, 0); }

    inline pixel_type operator()(float i, float j, int p = 0) const { return InterpolationImplementation<InterpT, EdgeExtendView<ImageT,EdgeT> >::interpolate(m_extend,i,j,p); }

    /// \cond INTERNAL
    // We can make an optimization here.  If the pixels in the child
    // view can be repeatedly accessed without incurring any
    // additional overhead (e.g. a TransposeView of an ImageView), we
    // do not need to rasterize the child before we proceed to
    // rasterize ourself.
    typedef typename boost::mpl::if_< IsMultiplyAccessible<ImageT>, 
 				      InterpolationView<typename ImageT::prerasterize_type, InterpT, EdgeT>,
 				      InterpolationView<ImageView<base_pixel_type>, InterpT, EdgeT> >::type prerasterize_type;

    inline prerasterize_type prerasterize() const {
      if (IsMultiplyAccessible<ImageT>::value) {
				return prerasterize_type( m_image.prerasterize() );
      } else {
				ImageView<base_pixel_type> buf( m_image.cols(), m_image.rows() );
				m_image.rasterize( buf );
				return prerasterize_type( buf );
      }
    }

    template <class DestT> inline void rasterize( DestT const& dest ) const { vw::rasterize( prerasterize(), dest ); }
    /// \endcond
  };
  
  /// \cond INTERNAL
  // Type traits 
  template <class ImageT, class InterpT, class EdgeT>
  struct IsFloatingPointIndexable<InterpolationView<ImageT, InterpT, EdgeT> > : public boost::true_type {};
  /// \endcond
	

  // -------------------------------------------------------------------------------
  // Functional API
  // -------------------------------------------------------------------------------

	/// Helper function for quickly creating a bilinear interpolation view of an image.
  template <class ImageT>
  InterpolationView<ImageT, BilinearInterpolation> bilinear_interpolation( ImageViewBase<ImageT> const& v ) {
    return InterpolationView<ImageT, BilinearInterpolation>( v.impl() );
  }

	/// Helper function for quickly creating a bicubic interpolation view of an image.
  template <class ImageT>
  InterpolationView<ImageT, BicubicInterpolation> bicubic_interpolation( ImageViewBase<ImageT> const& v ) {
    return InterpolationView<ImageT, BicubicInterpolation>( v.impl() );
  }
	
	/// Helper function for quickly creating a nearest neighbor interpolation view of an image.
  template <class ImageT>
  InterpolationView<ImageT, NearestPixelInterpolation> nearest_pixel_interpolation( ImageViewBase<ImageT> const& v ) {
    return InterpolationView<ImageT, NearestPixelInterpolation>( v.impl() );
  }
	
} // namespace vw

#endif // __VW_ITERPOLATION_H__
