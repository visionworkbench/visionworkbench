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

/// \file ImageAlgorithms.h
/// 
/// Basic algorithms operating on images.
/// 
#ifndef __VW_IMAGE_ALGORITHMS_H__
#define __VW_IMAGE_ALGORITHMS_H__

#include <vw/Core/CompoundTypes.h>
#include <vw/Image/ImageViewBase.h>
#include <vw/Image/PerPixelViews.h>
#include <vw/Image/Statistics.h>

namespace vw {

  // *******************************************************************
  // bounding_box()
  // *******************************************************************

  template <class ViewT>
  BBox2i bounding_box( ImageViewBase<ViewT> const& image_ ) {
    const typename ViewT::pixel_type zero = typename ViewT::pixel_type();
    ViewT const& image = static_cast<ViewT const&>(image_);
    int x=0, y=0, cols=0, rows=0;
    int i, j, icols = image.cols(), irows = image.rows();
    for( j=0; j<irows; ++j ) {
      for( i=0; i<icols; ++i ) {
        if( image(i,j) != zero ) break;
      }
      if( i != icols ) break;
    }
    if( j != irows ) {
      y = j;
      for( j=irows-1; j; --j ) {
        for( i=0; i<icols; ++i ) {
          if( image(i,j) != zero ) break;
        }
        if( i != icols ) break;
      }
      rows = j - y + 1;
      for( i=0; i<icols; ++i ) {
        for( j=y; j<y+rows; ++j ) {
          if( image(i,j) != zero ) break;
        }
        if( j != y+rows ) break;
      }
      x = i;
      for( i=icols-1; i; --i ) {
        for( j=y; j<y+rows; ++j ) {
          if( image(i,j) != zero ) break;
        }
        if( j != y+rows ) break;
      }
      cols = i - x + 1;
    }
    return BBox2i( x, y, cols, rows );
  }


  // *******************************************************************
  // fill()
  // *******************************************************************

  /// Fill an image with a constant pixel value
  template <class ImageT, class ValueT>
  void fill( ImageViewBase<ImageT> const& image, ValueT value_ ) {
    unsigned planes=image.impl().planes(), rows=image.impl().rows(), cols=image.impl().cols();
    typename ImageT::pixel_type value( value_ );
    typename ImageT::pixel_accessor plane = image.impl().origin();
    for( unsigned p=planes; p; --p ) {
      typename ImageT::pixel_accessor row = plane;
      for( unsigned r=rows; r; --r ) {
        typename ImageT::pixel_accessor col = row;
        for( unsigned c=cols; c; --c ) {
          *col = value;
          col.next_col();
        }
        row.next_row();
      }
      plane.next_plane();
    }
  }


  // *******************************************************************
  // clamp()
  // *******************************************************************

  template <class PixelT>
  class ChannelClampFunctor: public UnaryReturnSameType {
    typedef typename CompoundChannelType<PixelT>::type channel_type;
    channel_type m_low, m_high;
  public:
    ChannelClampFunctor( channel_type low, channel_type high ) :
      m_low(low), m_high(high) {}

    channel_type operator()( channel_type value ) const {
      if      (value > m_high) { return m_high; }
      else if (value < m_low ) { return m_low;  }
      else                     { return value;  }
    }
  };

  /// Clamp the values in an image to fall within the range [low,high).
  template <class ImageT, class ArgT>
  UnaryPerPixelView<ImageT,UnaryCompoundFunctor<ChannelClampFunctor<typename ImageT::pixel_type> > >
  inline clamp( ImageViewBase<ImageT> const& image, ArgT low, ArgT high ) {
    typedef UnaryCompoundFunctor<ChannelClampFunctor<typename ImageT::pixel_type> > func_type;
    func_type func( ChannelClampFunctor<typename ImageT::pixel_type>(low,high) );
    return UnaryPerPixelView<ImageT,func_type>( image.impl(), func );
  }


  // *******************************************************************
  // normalize()
  // *******************************************************************

  /// \cond INTERNAL
  class ChannelNormalizeFunctor: public UnaryReturnSameType {
  private:
    double m_old_min, m_new_min, m_old_to_new_ratio;
  public:
    ChannelNormalizeFunctor( double old_min, double old_max, double new_min, double new_max )
      : m_old_min(old_min), m_new_min(new_min)
    {
      if( old_max == old_min ) { m_old_to_new_ratio = 0.0; }
      else { m_old_to_new_ratio = (new_max - new_min)/(old_max - old_min); }
    }

    template <class ChannelT>
    ChannelT operator() (ChannelT value) const {
      return (ChannelT)((value - m_old_min) * m_old_to_new_ratio + m_new_min);
    }
  };
  /// \endcond

  /// Renormalize the values in an image to fall within the range [low,high).
  template <class ImageT>
  UnaryPerPixelView<ImageT, UnaryCompoundFunctor<ChannelNormalizeFunctor> >
  inline normalize( ImageViewBase<ImageT> const& image, double low = 0.0, double high = 1.0 ) {
    typedef UnaryCompoundFunctor<ChannelNormalizeFunctor> func_type;
    typename CompoundChannelType<typename ImageT::pixel_type>::type old_min, old_max;
    min_max_channel_values( image, old_min, old_max );
    func_type func( ChannelNormalizeFunctor( old_min, old_max, low, high ) );
    return UnaryPerPixelView<ImageT, func_type >( image.impl(), func );
  }


  // *******************************************************************
  // threshold()
  // *******************************************************************

  // A per-pixel thresholding filter with adjustable threshold and 
  // high and low values.
  template <class PixelT>
  class ChannelThresholdFunctor {
    typedef typename CompoundChannelType<PixelT>::type channel_type;
    channel_type m_thresh, m_low, m_high;
  public:

    ChannelThresholdFunctor( channel_type thresh, channel_type low, channel_type high )
      : m_thresh(thresh), m_low(low), m_high(high) {}
    
    template <class Args> struct result {
      typedef channel_type type;
    };

    inline channel_type operator()( channel_type const& val ) const {
      return (val > m_thresh) ? m_high : m_low;
    }
  };

  /// Threshold the values in the image and set the range to [low,high).
  template <class ImageT>
  UnaryPerPixelView<ImageT,UnaryCompoundFunctor<ChannelThresholdFunctor<typename ImageT::pixel_type> > >
  inline threshold( ImageViewBase<ImageT> const& image ) {
    typedef UnaryCompoundFunctor<ChannelThresholdFunctor<typename ImageT::pixel_type> > func_type;
    func_type func( ChannelThresholdFunctor<typename ImageT::pixel_type>(0,0,1) );
    return UnaryPerPixelView<ImageT,func_type>( image.impl(), func );
  }

  template <class ImageT, class ThreshT>
  UnaryPerPixelView<ImageT,UnaryCompoundFunctor<ChannelThresholdFunctor<typename ImageT::pixel_type> > >
  inline threshold( ImageViewBase<ImageT> const& image, ThreshT thresh ) {
    typedef UnaryCompoundFunctor<ChannelThresholdFunctor<typename ImageT::pixel_type> > func_type;
    func_type func( ChannelThresholdFunctor<typename ImageT::pixel_type>(thresh,0,1) );
    return UnaryPerPixelView<ImageT,func_type>( image.impl(), func );
  }

  /* This overload disagrees with the optional-argument semantics of clamp() and normalize()!
  template <class ImageT, class ThreshT, class HighT>
  UnaryPerPixelView<ImageT,UnaryCompoundFunctor<ChannelThresholdFunctor<typename ImageT::pixel_type> > >
  inline threshold( ImageViewBase<ImageT> const& image, ThreshT thresh, HighT high ) {
    typedef UnaryCompoundFunctor<ChannelThresholdFunctor<typename ImageT::pixel_type> > func_type;
    func_type func( ChannelThresholdFunctor<typename ImageT::pixel_type>(thresh,0,high) );
    return UnaryPerPixelView<ImageT,func_type>( image.impl(), func );
  }
  */

  template <class ImageT, class ThreshT, class LowT, class HighT>
  UnaryPerPixelView<ImageT,UnaryCompoundFunctor<ChannelThresholdFunctor<typename ImageT::pixel_type> > >
  inline threshold( ImageViewBase<ImageT> const& image, ThreshT thresh, LowT low, HighT high ) {
    typedef UnaryCompoundFunctor<ChannelThresholdFunctor<typename ImageT::pixel_type> > func_type;
    func_type func( ChannelThresholdFunctor<typename ImageT::pixel_type>(thresh,low,high) );
    return UnaryPerPixelView<ImageT,func_type>( image.impl(), func );
  }

} // namespace vw

#endif // __VW_IMAGE_ALGORITHMS_H__
