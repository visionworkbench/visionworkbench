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

/// \file Algorithms.h
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
  // fill()
  // *******************************************************************

  /// Fill an image with a constant pixel value
  template <class ImageT, class ValueT>
  void fill( ImageViewBase<ImageT> const& image, ValueT value_ ) {
    int32 planes=image.impl().planes(), rows=image.impl().rows(), cols=image.impl().cols();
    typename ImageT::pixel_type value( value_ );
    typename ImageT::pixel_accessor plane = image.impl().origin();
    for( int32 p=planes; p; --p ) {
      typename ImageT::pixel_accessor row = plane;
      for( int32 r=rows; r; --r ) {
        typename ImageT::pixel_accessor col = row;
        for( int32 c=cols; c; --c ) {
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
      m_low(low), m_high(high) {
    }

    channel_type operator()( channel_type value ) const {
      if      (value > m_high) { return m_high; }
      else if (value < m_low ) { return m_low;  }
      else                     { return value;  }
    }
  };

  /// Clamp the values in an image to fall within the range [low,high].
  template <class ImageT, class LowT, class HighT>
  UnaryPerPixelView<ImageT,UnaryCompoundFunctor<ChannelClampFunctor<typename ImageT::pixel_type> > >
  inline clamp( ImageViewBase<ImageT> const& image, LowT low, HighT high ) {
    typedef UnaryCompoundFunctor<ChannelClampFunctor<typename ImageT::pixel_type> > func_type;
    func_type func( ChannelClampFunctor<typename ImageT::pixel_type>(low,high) );
    return UnaryPerPixelView<ImageT,func_type>( image.impl(), func );
  }

  /// Clamp the values in an image to fall within the range [0,high].
  /// The low end of the range is actually determined by the
  /// ChannelRange type trait but is generally zero.
  template <class ImageT, class HighT>
  UnaryPerPixelView<ImageT,UnaryCompoundFunctor<ChannelClampFunctor<typename ImageT::pixel_type> > >
  inline clamp( ImageViewBase<ImageT> const& image, HighT high ) {
    typedef UnaryCompoundFunctor<ChannelClampFunctor<typename ImageT::pixel_type> > func_type;
    typedef ChannelRange<typename CompoundChannelType<typename ImageT::pixel_type>::type> range_type;
    typename CompoundChannelType<typename ImageT::pixel_type>::type min_val = range_type::min();
    func_type func( ChannelClampFunctor<typename ImageT::pixel_type>(min_val,high) );
    return UnaryPerPixelView<ImageT,func_type>( image.impl(), func );
  }

  /// Clamp the values in an image to fall within the range [min,max],
  /// where min and max are determined by the ChannelRange type trait
  /// and are generally equal to 0.0 and 1.0 for floating point types
  /// and 0 and the largest positve value for integral types.
  template <class ImageT>
  UnaryPerPixelView<ImageT,UnaryCompoundFunctor<ChannelClampFunctor<typename ImageT::pixel_type> > >
  inline clamp( ImageViewBase<ImageT> const& image ) {
    typedef UnaryCompoundFunctor<ChannelClampFunctor<typename ImageT::pixel_type> > func_type;
    typedef ChannelRange<typename CompoundChannelType<typename ImageT::pixel_type>::type> range_type;
    typename CompoundChannelType<typename ImageT::pixel_type>::type min_val = range_type::min();
    typename CompoundChannelType<typename ImageT::pixel_type>::type max_val = range_type::max();
    func_type func( ChannelClampFunctor<typename ImageT::pixel_type>(min_val,max_val) );
    return UnaryPerPixelView<ImageT,func_type>( image.impl(), func );
  }

  // *******************************************************************
  // normalize()
  // *******************************************************************

  /// \cond INTERNAL
  template <class PixelT>
  class ChannelNormalizeFunctor: public UnaryReturnSameType {
    typedef typename CompoundChannelType<PixelT>::type channel_type;
    channel_type m_old_min, m_new_min;
    double m_old_to_new_ratio;
  public:
    ChannelNormalizeFunctor( channel_type old_min, channel_type old_max, 
                             channel_type new_min, channel_type new_max )
      : m_old_min(old_min), m_new_min(new_min)
    {
      if( old_max == old_min ) { m_old_to_new_ratio = 0.0; }
      else { m_old_to_new_ratio = (new_max - new_min)/(double)(old_max - old_min); }
    }

    template <class ChannelT>
    ChannelT operator()( ChannelT value ) const {
      return (ChannelT)((value - m_old_min) * m_old_to_new_ratio + m_new_min);
    }
  };

  template <class PixelT>
  class ChannelNormalizeRetainAlphaFunctor: public UnaryReturnSameType {
    typedef typename CompoundChannelType<PixelT>::type channel_type;
    channel_type m_old_min, m_new_min;
    double m_old_to_new_ratio;
  public:
    ChannelNormalizeRetainAlphaFunctor( channel_type old_min, channel_type old_max, 
                                        channel_type new_min, channel_type new_max )
      : m_old_min(old_min), m_new_min(new_min)
    {
      if( old_max == old_min ) { m_old_to_new_ratio = 0.0; }
      else { m_old_to_new_ratio = (new_max - new_min)/(double)(old_max - old_min); }
    }

    PixelT operator()( PixelT value ) const {
      if (is_transparent(value)) return value;
      else {
        PixelT result;
        non_alpha_channels(result) = (non_alpha_channels(value) - m_old_min) * m_old_to_new_ratio + m_new_min;
        alpha_channel(result) = alpha_channel(value);
        return result;
      }
    }
  };
  /// \endcond

  /// Renormalize the values in an image to fall within the range [low,high).
  template <class ImageT, class OldLowT, class OldHighT, class NewLowT, class NewHighT>
  UnaryPerPixelView<ImageT, ChannelNormalizeRetainAlphaFunctor<typename ImageT::pixel_type> >
  inline normalize_retain_alpha( ImageViewBase<ImageT> const& image, OldLowT old_low, OldHighT old_high, NewLowT new_low, NewHighT new_high  ) {
    typedef ChannelNormalizeRetainAlphaFunctor<typename ImageT::pixel_type> func_type;
    func_type func ( old_low, old_high, new_low, new_high );
    return UnaryPerPixelView<ImageT, func_type >( image.impl(), func );
  }

  /// Renormalize the values in an image to fall within the range [low,high).
  template <class ImageT, class OldLowT, class OldHighT, class NewLowT, class NewHighT>
  UnaryPerPixelView<ImageT, UnaryCompoundFunctor<ChannelNormalizeFunctor<typename ImageT::pixel_type> > >
  inline normalize( ImageViewBase<ImageT> const& image, OldLowT old_low, OldHighT old_high, NewLowT new_low, NewHighT new_high  ) {
    typedef UnaryCompoundFunctor<ChannelNormalizeFunctor<typename ImageT::pixel_type> > func_type;
    func_type func( ChannelNormalizeFunctor<typename ImageT::pixel_type>( old_low, old_high, new_low, new_high ) );
    return UnaryPerPixelView<ImageT, func_type >( image.impl(), func );
  }

  /// Renormalize the values in an image to fall within the range [low,high).
  template <class ImageT, class LowT, class HighT>
  UnaryPerPixelView<ImageT, UnaryCompoundFunctor<ChannelNormalizeFunctor<typename ImageT::pixel_type> > >
  inline normalize( ImageViewBase<ImageT> const& image, LowT low, HighT high ) {
    typedef UnaryCompoundFunctor<ChannelNormalizeFunctor<typename ImageT::pixel_type> > func_type;
    typename CompoundChannelType<typename ImageT::pixel_type>::type old_min, old_max;
    min_max_channel_values( image, old_min, old_max );
    func_type func( ChannelNormalizeFunctor<typename ImageT::pixel_type>( old_min, old_max, low, high ) );
    return UnaryPerPixelView<ImageT, func_type >( image.impl(), func );
  }

  /// Renormalize the values in an image to fall within the range
  /// [0,high).  The low end of the range is actually determined by
  /// the ChannelRange type trait but is generally zero.
  template <class ImageT, class HighT>
  UnaryPerPixelView<ImageT, UnaryCompoundFunctor<ChannelNormalizeFunctor<typename ImageT::pixel_type> > >
  inline normalize( ImageViewBase<ImageT> const& image, HighT high ) {
    typedef UnaryCompoundFunctor<ChannelNormalizeFunctor<typename ImageT::pixel_type> > func_type;
    typedef ChannelRange<typename CompoundChannelType<typename ImageT::pixel_type>::type> range_type;
    typename CompoundChannelType<typename ImageT::pixel_type>::type old_min, old_max;
    min_max_channel_values( image, old_min, old_max );
    func_type func( ChannelNormalizeFunctor<typename ImageT::pixel_type>( old_min, old_max, range_type::min(), high ) );
    return UnaryPerPixelView<ImageT, func_type >( image.impl(), func );
  }

  /// Renormalize the values in an image to fall within the range
  /// [min,max), where min and max are determined by the ChannelRange
  /// type trait and are generally equal to 0.0 and 1.0 for floating
  /// point types and 0 and the largest positve value for integral
  /// types.
  template <class ImageT>
  UnaryPerPixelView<ImageT, UnaryCompoundFunctor<ChannelNormalizeFunctor<typename ImageT::pixel_type> > >
  inline normalize( ImageViewBase<ImageT> const& image ) {
    typedef UnaryCompoundFunctor<ChannelNormalizeFunctor<typename ImageT::pixel_type> > func_type;
    typedef ChannelRange<typename CompoundChannelType<typename ImageT::pixel_type>::type> range_type;
    typename CompoundChannelType<typename ImageT::pixel_type>::type old_min, old_max;
    min_max_channel_values( image, old_min, old_max );
    func_type func( ChannelNormalizeFunctor<typename ImageT::pixel_type>( old_min, old_max, range_type::min(), range_type::max() ) );
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

  /// Threshold the values in an image, generating a two-valued output
  /// image with values low and high.
  template <class ImageT, class ThreshT, class LowT, class HighT>
  UnaryPerPixelView<ImageT, UnaryCompoundFunctor<ChannelThresholdFunctor<typename ImageT::pixel_type> > >
  inline threshold( ImageViewBase<ImageT> const& image, ThreshT thresh, LowT low, HighT high ) {
    typedef UnaryCompoundFunctor<ChannelThresholdFunctor<typename ImageT::pixel_type> > func_type;
    func_type func( ChannelThresholdFunctor<typename ImageT::pixel_type>(thresh,low,high) );
    return UnaryPerPixelView<ImageT,func_type>( image.impl(), func );
  }

  /// Threshold the values in an image, generating a two-valued output
  /// image with values 0 and high.  The low value is actually
  /// determined by the ChannelRange type trait but is generally zero.
  template <class ImageT, class ThreshT, class HighT>
  UnaryPerPixelView<ImageT,UnaryCompoundFunctor<ChannelThresholdFunctor<typename ImageT::pixel_type> > >
  inline threshold( ImageViewBase<ImageT> const& image, ThreshT thresh, HighT high ) {
    typedef UnaryCompoundFunctor<ChannelThresholdFunctor<typename ImageT::pixel_type> > func_type;
    typedef ChannelRange<typename CompoundChannelType<typename ImageT::pixel_type>::type> range_type;
    func_type func( ChannelThresholdFunctor<typename ImageT::pixel_type>(thresh,range_type::min(),high) );
    return UnaryPerPixelView<ImageT,func_type>( image.impl(), func );
  }

  /// Threshold the values in an image, generating a two-valued output
  /// where the values are determined by the ChannelRange type trait
  /// and are generally equal to 0.0 and 1.0 for floating point types
  /// and 0 and the largest positve value for integral types.
  template <class ImageT, class ThreshT>
  UnaryPerPixelView<ImageT,UnaryCompoundFunctor<ChannelThresholdFunctor<typename ImageT::pixel_type> > >
  inline threshold( ImageViewBase<ImageT> const& image, ThreshT thresh ) {
    typedef UnaryCompoundFunctor<ChannelThresholdFunctor<typename ImageT::pixel_type> > func_type;
    typedef ChannelRange<typename CompoundChannelType<typename ImageT::pixel_type>::type> range_type;
    func_type func( ChannelThresholdFunctor<typename ImageT::pixel_type>(thresh,range_type::min(),range_type::max()) );
    return UnaryPerPixelView<ImageT,func_type>( image.impl(), func );
  }

  /// Threshold the values in an image against zero, generating a
  /// two-valued output where the values are determined by the
  /// ChannelRange type trait and are generally equal to 0.0 and 1.0
  /// for floating point types and 0 and the largest positve value for
  /// integral types.
  template <class ImageT>
  UnaryPerPixelView<ImageT,UnaryCompoundFunctor<ChannelThresholdFunctor<typename ImageT::pixel_type> > >
  inline threshold( ImageViewBase<ImageT> const& image ) {
    typedef UnaryCompoundFunctor<ChannelThresholdFunctor<typename ImageT::pixel_type> > func_type;
    typedef ChannelRange<typename CompoundChannelType<typename ImageT::pixel_type>::type> range_type;
    func_type func( ChannelThresholdFunctor<typename ImageT::pixel_type>(0,range_type::min(),range_type::max()) );
    return UnaryPerPixelView<ImageT,func_type>( image.impl(), func );
  }


  // *******************************************************************
  // grassfire()
  // *******************************************************************

  // Computes the 4-connected grassfire image of an image.
  // (Specifically, computes the Manhattan distance from each pixel to
  // the nearest pixel with zero value, assuming the borders of the
  // image are zero.)
  template <class SourceT, class OutputT>
  void grassfire( ImageViewBase<SourceT> const& src, ImageView<OutputT>& dst ) {
    int32 cols = src.impl().cols(), rows = src.impl().rows();
    dst.set_size( cols, rows );
    typename SourceT::pixel_accessor srow = src.impl().origin();
    typename ImageView<OutputT>::pixel_accessor drow = dst.origin();
    const typename SourceT::pixel_type zero = typename SourceT::pixel_type();
    { // First row
      typename SourceT::pixel_accessor scol = srow;
      typename ImageView<OutputT>::pixel_accessor dcol = drow;
      for( int32 col=cols; col; --col ) {
        *dcol = ((*scol)==zero)?0:1;
        scol.next_col();
        dcol.next_col();
      }
      srow.next_row();
      drow.next_row();
    }
    for( int32 row=rows-2; row; --row ) {
      typename SourceT::pixel_accessor scol = srow;
      typename ImageView<OutputT>::pixel_accessor dcol = drow;
      *dcol = ((*scol)==zero)?0:1;
      scol.next_col();
      dcol.next_col();
      for( int32 col=cols-2; col; --col ) {
        if( (*scol)==zero ) (*dcol)=0;
        else {
          typename ImageView<OutputT>::pixel_accessor s1 = dcol, s2 = dcol;
          (*dcol) = 1 + std::min( *(s1.prev_col()), *(s2.prev_row()) );
        }
        scol.next_col();
        dcol.next_col();
      }
      *dcol = ((*scol)==zero)?0:1;
      srow.next_row();
      drow.next_row();
    }
    { // Last row
      typename SourceT::pixel_accessor scol = srow;
      typename ImageView<OutputT>::pixel_accessor dcol = drow;
      for( int32 col=cols; col; --col ) {
        *dcol = ((*scol)==zero)?0:1;
        scol.next_col();
        dcol.next_col();
      }
    }
    drow.advance(cols-2,-1);
    for( int32 row=rows-2; row; --row ) {
      typename ImageView<OutputT>::pixel_accessor dcol = drow;
      for( int32 col=cols-2; col; --col ) {
        if( (*dcol)!=0 ) {
          typename ImageView<OutputT>::pixel_accessor s1 = dcol, s2 = dcol;
          int32 m = std::min( *(s1.next_col()), *(s2.next_row()) );
          if( m < *dcol ) *dcol = m + 1;
        }
        dcol.prev_col();
      }
      drow.prev_row();
    }
  }

  // Without destination given, return in a newly-created ImageView<int32>
  template <class SourceT>
  ImageView<int32> grassfire( ImageViewBase<SourceT> const& src ) {
    int32 cols = src.impl().cols(), rows = src.impl().rows();
    ImageView<int32> result( cols, rows );
    grassfire( src, result );
    return result;
  }

  // *******************************************************************
  // bounding_box()
  // *******************************************************************

  template <class ViewT>
  BBox2i bounding_box( ImageViewBase<ViewT> const& image_ ) {
    return BBox2i( 0, 0, image_.impl().cols(), image_.impl().rows() );
  }


  // *******************************************************************
  // nonzero_data_bounding_box()
  // *******************************************************************

  template <class ViewT>
  BBox2i nonzero_data_bounding_box( ImageViewBase<ViewT> const& image_ ) {
    const typename ViewT::pixel_type zero = typename ViewT::pixel_type();
    ViewT const& image = static_cast<ViewT const&>(image_);
    int32 x=0, y=0, cols=0, rows=0;
    int32 i, j, icols = image.cols(), irows = image.rows();
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
  // is_opaque()
  // *******************************************************************

  template <class ImageT>
  bool is_opaque_helper( ImageT const& image, true_type ) {
    typename PixelChannelType<typename ImageT::pixel_type>::type maxval = ChannelRange<typename ImageT::pixel_type>::max();
    for( int32 y=0; y<image.rows(); ++y )
      for( int32 x=0; x<image.cols(); ++x )
        if( image(x,y)[PixelNumChannels<typename ImageT::pixel_type>::value-1] < maxval ) return false;
    return true;
  }
  
  template <class ImageT>
  bool is_opaque_helper( ImageT const& image, false_type ) {
    return true;
  }

  /// Returns true if the given image is entirely opaque, or false if 
  /// it is at least partially transparent.
  template <class ImageT>
  bool is_opaque( ImageViewBase<ImageT> const& image ) {
    return is_opaque_helper( image.impl(), typename PixelHasAlpha<typename ImageT::pixel_type>::type() );
  }


  // *******************************************************************
  // image_blocks()
  // *******************************************************************

  /// A utility routine that, given an image, returns a vector of
  /// bounding boxes for sub-regions of the image of the specified
  /// size.  Note that bounding boxes along the right and bottom edges
  /// of the image will not have the specified dimension unless the
  /// image width and height are perfectly divisible by the bounding
  /// box width and height, respectively. This routine is useful if you
  /// want to apply an operation to a large image one region at a time.
  template <class ViewT>
  std::vector<BBox2i> image_blocks(ImageViewBase<ViewT> const& image, 
                                   int32 block_width, int32 block_height) {  
    std::vector<BBox2i> bboxes;

    int32 j_offset = 0;
    while ( j_offset < image.impl().rows() ) {
      int32 j_dim = (image.impl().rows() - j_offset) < block_height ? (image.impl().rows() - j_offset) : block_height;
      int32 i_offset = 0;
      while ( i_offset < image.impl().cols() ) {
        int32 i_dim = (image.impl().cols() - i_offset) < block_width ? (image.impl().cols() - i_offset) : block_width;      
        bboxes.push_back(BBox2i(i_offset,j_offset,i_dim,j_dim));
        i_offset += i_dim;
      }
      j_offset += j_dim;
    }
    return bboxes;
  }


} // namespace vw

#endif // __VW_IMAGE_ALGORITHMS_H__
