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

/// \file Convolution.h
/// 
/// One- and two-dimensional convolution numerical functions, and 
/// standard and separable two-dimensional image convolution view 
/// classes used by the filtering functions in \ref Filter.h.
///
#ifndef __VW_IMAGE_CONVOLUTION_H__
#define __VW_IMAGE_CONVOLUTION_H__

#include <vector>
#include <iterator>

#include <vw/Core/Exception.h>
#include <vw/Core/TypeDeduction.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/Manipulation.h>
#include <vw/Image/EdgeExtend.h>

namespace vw {

  // *******************************************************************
  // The core correlation functions
  // *******************************************************************

  /// \cond INTERNAL

  // These primitive operations are currently only used by the
  // convolution views, but I kept them free functions because I
  // thought they might be useful down the road in other
  // correlation-type functions....

  template <class SrcAccessT, class KernelIterT>
  typename ProductType<typename SrcAccessT::pixel_type, typename std::iterator_traits<KernelIterT>::value_type>::type
  inline correlate_1d_at_point( SrcAccessT const& src, KernelIterT const& kernel, unsigned n ) {
    typedef typename ProductType<typename SrcAccessT::pixel_type, typename std::iterator_traits<KernelIterT>::value_type>::type result_type;
    result_type result = result_type();
    SrcAccessT s = src;
    KernelIterT k = kernel;
    for( unsigned i=n; i; --i ) {
      result += (*k)*(*s);
      s.next_col();
      ++k;
    }
    return result;
  }

  template <class SrcAccessT, class KernelAccessT>
  typename ProductType<typename SrcAccessT::pixel_type, typename KernelAccessT::pixel_type>::type
  inline correlate_2d_at_point( SrcAccessT const& src, KernelAccessT const& kernel, unsigned cols, unsigned rows ) {
    typedef typename ProductType<typename SrcAccessT::pixel_type, typename KernelAccessT::pixel_type>::type result_type;
    result_type result = result_type();
    SrcAccessT srow = src;
    KernelAccessT krow = kernel;
    for( unsigned j=rows; j; --j ) {
      SrcAccessT scol = srow;
      KernelAccessT kcol = krow;
      for( unsigned i=cols; i; --i ) {
        result += (*kcol)*(*scol);
        scol.next_col();
        kcol.next_col();
      }
      srow.next_row();
      krow.next_row();
    }
    return result;
  }

  /// \endcond


  // *******************************************************************
  // The standard 2D convolution view type
  // *******************************************************************

  /// A standard 2D convolution image view.
  ///
  /// Represents the convolution of an image with a 2D kernel.
  ///
  /// \see ConvolutionFilter
  template <class ImageT, class KernelT, class EdgeT>
  class ConvolutionView : public ImageViewBase<ConvolutionView<ImageT,KernelT,EdgeT> >
  {
  private:
    EdgeExtendView<ImageT,EdgeT> m_image;
    Rotate180View<KernelT> m_kernel;
    int m_ci, m_cj;

  public:
    /// The pixel type of the image view.
    typedef typename ProductType<typename ImageT::pixel_type, typename KernelT::pixel_type>::type pixel_type;

    /// The view's %pixel_accessor type.
    typedef ProceduralPixelAccessor<ConvolutionView<ImageT, KernelT, EdgeT> > pixel_accessor;

    /// Constructs a ConvolutionView with the given image and kernel and with the origin of the kernel located at the point (ci,cj).
    ConvolutionView( ImageT const& image, KernelT const& kernel, unsigned ci, unsigned cj, EdgeT const& edge = EdgeT() )
      : m_image(image,edge), m_kernel(kernel), m_ci(ci), m_cj(cj) {}

    /// Constructs a ConvolutionView with the given image and kernel and with the origin of the kernel located at the center.
    ConvolutionView( ImageT const& image, KernelT const& kernel, EdgeT const& edge = EdgeT() )
      : m_image(image,edge), m_kernel(kernel), m_ci(kernel.cols()/2), m_cj(kernel.rows()/2) {}

    /// Returns the number of columns in the image.
    inline unsigned cols() const { return m_image.cols(); }

    /// Returns the number of rows in the image.
    inline unsigned rows() const { return m_image.rows(); }

    /// Returns the number of planes in the image.
    inline unsigned planes() const { return m_image.planes(); }

    /// Returns a pixel_accessor pointing to the origin.
    inline pixel_accessor origin() const { return pixel_accessor( *this ); }

    /// Returns the pixel at the given position in the given plane.
    inline pixel_type operator()( int x, int y, int p=0 ) const {
      if( (x >= m_ci) && (y >= m_cj) &&
	  (x <= int(m_image.cols())-int(m_kernel.cols())+m_ci) &&
	  (y <= int(m_image.rows())-int(m_kernel.rows())+m_cj) ) {
        return correlate_2d_at_point( m_image.child().origin().advance(x-m_ci,y-m_cj,p),
                                      m_kernel.origin(), m_kernel.cols(), m_kernel.rows() );
      }
      else {
	return correlate_2d_at_point( m_image.origin().advance(x-m_ci,y-m_cj,p),
                                      m_kernel.origin(), m_kernel.cols(), m_kernel.rows() );
      }
    }

    /// \cond INTERNAL

    // If the pixels in the child view can be repeatedly accessed
    // without incurring any additional overhead (e.g. an ImageView)
    // we do not need to rasterize the child before we proceed to
    // rasterize ourself.
    typedef typename boost::mpl::if_< IsMultiplyAccessible<ImageT>, 
 				      ConvolutionView<typename ImageT::prerasterize_type, KernelT, EdgeT>,
 				      ConvolutionView<ImageView<typename ImageT::pixel_type>, KernelT, EdgeT> >::type prerasterize_type;
    inline prerasterize_type prerasterize() const {
      if (IsMultiplyAccessible<ImageT>::value) {
	return prerasterize_type( m_image.child().prerasterize(), m_kernel.child(), m_ci, m_cj, m_image.func() );
      } else {
	return prerasterize_type( copy(m_image.child()), m_kernel.child(), m_ci, m_cj, m_image.func() );
      }
    }

    template <class DestT>
    void rasterize( DestT const& dest ) const {
      prerasterize().rasterize_helper( dest );
    }

    template <class DestT>
    void rasterize_helper( DestT const& dest ) const {
      typedef typename ImageT::pixel_accessor SrcAccessT;
      typedef typename EdgeExtendView<ImageT,EdgeT>::pixel_accessor EdgeAccessT;
      typedef typename DestT::pixel_accessor DestAccessT;
      typedef typename DestT::pixel_type DestPixelT;
      SrcAccessT splane = m_image.child().origin();
      EdgeAccessT eplane = m_image.origin().advance(-int(m_ci),-int(m_cj));
      DestAccessT dplane = dest.origin();

      for( int p=0; p<int(m_image.planes()); ++p ) {
        SrcAccessT srow = splane;
        EdgeAccessT erow = eplane;
        DestAccessT drow = dplane;
        for( int r=0; r<int(m_image.rows()); ++r ) {
          EdgeAccessT ecol = erow;
          DestAccessT dcol = drow;
          if( r<m_cj || r>int(m_image.rows())-int(m_kernel.rows())+m_cj ) {
            for( int i=0; i<int(m_image.cols()); ++i ) {
              *dcol = correlate_2d_at_point( ecol, m_kernel.origin(), m_kernel.cols(), m_kernel.rows() );
              dcol.next_col();
              ecol.next_col();
            }
          }
          else {
            SrcAccessT scol = srow;
            int i=0;
            for( ; i<m_ci; ++i ) {
              *dcol = DestPixelT( correlate_2d_at_point( ecol, m_kernel.origin(), m_kernel.cols(), m_kernel.rows() ) );
              ecol.next_col();
              dcol.next_col();
            }
            for( ; i<=int(m_image.cols())-int(m_kernel.cols())+m_ci ; ++i ) {
              *dcol = DestPixelT( correlate_2d_at_point( scol, m_kernel.origin(), m_kernel.cols(), m_kernel.rows() ) );
              scol.next_col();
              dcol.next_col();
            }
            ecol.advance( m_image.cols()-m_kernel.cols()+1, 0 );
            for( ; i<int(m_image.cols()); ++i ) {
              *dcol = DestPixelT( correlate_2d_at_point( ecol, m_kernel.origin(), m_kernel.cols(), m_kernel.rows() ) );
              ecol.next_col();
              dcol.next_col();
            }
            srow.next_row();
          }
          erow.next_row();
          drow.next_row();
        }
        splane.next_plane();
        eplane.next_plane();
        dplane.next_plane();
      }
    }
    /// \endcond
  };


  // *******************************************************************
  // The separable 2D convolution view type
  // *******************************************************************

  /// A separable 2D convolution view.
  ///
  /// Represents the convolution of an image with a separable 2D
  /// kernel, specified in terms of its x and y factors.
  template <class ImageT, class KernelT, class EdgeT >
  class SeparableConvolutionView : public ImageViewBase<SeparableConvolutionView<ImageT,KernelT,EdgeT> >
  {
  private:
    ImageT m_image;
    std::vector<KernelT> m_i_kernel, m_j_kernel;
    unsigned m_ci, m_cj;
    EdgeT m_edge;
    mutable ImageView<KernelT> m_kernel2d;

    void generate2DKernel() const {
      unsigned ni = m_i_kernel.size() ? m_i_kernel.size() : 1;
      unsigned nj = m_j_kernel.size() ? m_j_kernel.size() : 1;
      m_kernel2d.set_size( ni, nj );
      for( unsigned i=0; i<ni; ++i )
        for( unsigned j=0; j<nj; ++j )
          m_kernel2d(ni-1-i,nj-1-j) = (m_i_kernel.size()?m_i_kernel[i]:1)*(m_j_kernel.size()?m_j_kernel[j]:1);
    }

  public:
    /// The pixel type of the view.
    typedef typename ProductType<typename ImageT::pixel_type, KernelT>::type pixel_type;

    /// The view's %pixel_accessor type.
    typedef ProceduralPixelAccessor<SeparableConvolutionView<ImageT, KernelT, EdgeT> > pixel_accessor;

    /// Constructs a SeparableConvolutionView with the given image and kernels and with the origin of the kernel located at the center.
    template <class KRangeT>
    SeparableConvolutionView( ImageT const& image, KRangeT const& ik, KRangeT const& jk, unsigned ci, unsigned cj, EdgeT const& edge = EdgeT() ) : 
      m_image(image), m_i_kernel(ik.begin(),ik.end()), m_j_kernel(jk.begin(),jk.end()), m_ci(ci), m_cj(cj), m_edge(edge) {}

    /// Constructs a SeparableConvolutionView with the given image and kernels and with the origin of the kernel located at the point (ci,cj).
    template <class KRangeT>
    SeparableConvolutionView( ImageT const& image, KRangeT const& ik, KRangeT const& jk, EdgeT const& edge = EdgeT() ) :
      m_image(image), m_i_kernel(ik.begin(),ik.end()), m_j_kernel(jk.begin(),jk.end()), m_ci(m_i_kernel.size()/2), m_cj(m_j_kernel.size()/2), m_edge(edge) {}

    /// Returns the number of columns in the image.
    inline unsigned cols() const { return m_image.cols(); }

    /// Returns the number of rows in the image.
    inline unsigned rows() const { return m_image.rows(); }

    /// Returns the number of planes in the image.
    inline unsigned planes() const { return m_image.planes(); }

    /// Returns a pixel_accessor pointing to the top-left corner of the first plane.
    inline pixel_accessor origin() const { return pixel_accessor( *this ); }

    /// Returns the pixel at the given position in the given plane.
    inline pixel_type operator()( int x, int y, int p = 0 ) const
    {
      if( m_kernel2d.cols()==0 ) generate2DKernel();

      if( (x >= int(m_ci)) && (y >= int(m_cj)) &&
	  (x <= int(m_image.cols())-int(m_kernel2d.cols())+int(m_ci)) &&
	  (y <= int(m_image.rows())-int(m_kernel2d.rows())+int(m_cj)) ) {
        return correlate_2d_at_point( m_image.origin().advance(x-m_ci,y-m_cj,p),
                                      m_kernel2d.origin(), m_kernel2d.cols(), m_kernel2d.rows() );
      }
      else {
	return correlate_2d_at_point( edge_extend(m_image,m_edge).origin().advance(x-m_ci,y-m_cj,p),
                                      m_kernel2d.origin(), m_kernel2d.cols(), m_kernel2d.rows() );
      }
    }

    /// \cond INTERNAL

    // The seperable convolution view knows that it is fastest to
    // fully rasterize itself if it is about to be part of a
    // rasterization operation with nested views.
    typedef ImageView<pixel_type> prerasterize_type;
    inline prerasterize_type prerasterize() const {
      prerasterize_type dest( m_image.cols(), m_image.rows(), m_image.planes() );
      rasterize( dest );
      return dest;
    }
    
    // If the pixels in the child view can be repeatedly accessed without 
    // incurring any additional overhead (e.g. a TransposeView of an 
    // ImageView), we do not need to rasterize the child before we proceed 
    // to rasterize ourself.
    template <class DestT>
    void rasterize( DestT const& dest ) const {
      if( IsMultiplyAccessible<ImageT>::value ) {
	rasterize_helper( m_image.prerasterize(), dest );
      } else {
	rasterize_helper( copy(m_image), dest );
      }
    }

    template <class SrcT, class DestT>
    void rasterize_helper( SrcT const& src, DestT const& dest ) const {	
      unsigned ni=m_i_kernel.size(), nj=m_j_kernel.size();
      if( ni>0 && nj>0 ) {
        ImageView<pixel_type> work( src.cols(), src.rows(), src.planes() );
        convolve_1d( src, work, m_i_kernel, m_ci );
        convolve_1d( transpose(work), transpose(dest), m_j_kernel, m_cj );
      }
      else if( ni>0 ) {
        convolve_1d( src, dest, m_i_kernel, m_ci );
      }
      else if( nj>0 ) {
        convolve_1d( transpose(src), transpose(dest), m_j_kernel, m_cj );
      }
      else {
        src.rasterize( dest );
      }
    }

    template <class SrcT, class DestT>
    void convolve_1d( SrcT const& src, DestT const& dest, std::vector<KernelT> const& kernel, unsigned c ) const {
      typedef typename SrcT::pixel_accessor SrcAccessT;
      EdgeExtendView<SrcT,EdgeT> edge_view( src, m_edge );
      typedef typename EdgeExtendView<SrcT,EdgeT>::pixel_accessor EdgeAccessT;
      typedef typename DestT::pixel_accessor DestAccessT;
      typedef typename DestT::pixel_type DestPixelT;
      typedef typename ProductType<typename SrcT::pixel_type, KernelT>::type AccumT;

      SrcAccessT splane = src.origin();
      EdgeAccessT eplane = edge_view.origin().advance(-c,0);
      DestAccessT dplane = dest.origin();
      for( unsigned p=0; p<src.planes(); ++p ) {
        SrcAccessT srow = splane;
        EdgeAccessT erow = eplane;
        DestAccessT drow = dplane;
        for( unsigned j=0; j<src.rows(); ++j ) {
          SrcAccessT scol = srow;
          EdgeAccessT ecol = erow;
          DestAccessT dcol = drow;
          unsigned i=0;
          if( kernel.size() <= src.cols() ) {
            for( ; i<c; ++i ) {
              *dcol = DestPixelT( correlate_1d_at_point( ecol, kernel.rbegin(), kernel.size() ) );
              ecol.next_col();
              dcol.next_col();
            }
            for( ; i<=src.cols()-kernel.size()+c ; ++i ) {
              *dcol = DestPixelT( correlate_1d_at_point( scol, kernel.rbegin(), kernel.size() ) );
              scol.next_col();
              dcol.next_col();
            }
            ecol.advance( src.cols()-kernel.size()+1, 0 );
          }
          for( ; i<src.cols(); ++i ) {
            *dcol = DestPixelT( correlate_1d_at_point( ecol, kernel.rbegin(), kernel.size() ) );
            ecol.next_col();
            dcol.next_col();
          }
          srow.next_row();
          erow.next_row();
          drow.next_row();
        }
        splane.next_plane();
        eplane.next_plane();
        dplane.next_plane();
      }
    }
    
    /// \endcond
  };

} // namespace vw

#endif // __VW_IMAGE_CONVOLUTION_H__
