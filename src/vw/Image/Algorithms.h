// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file Algorithms.h
///
/// Basic algorithms operating on images.
///
#ifndef __VW_IMAGE_ALGORITHMS_H__
#define __VW_IMAGE_ALGORITHMS_H__

#include <vw/Image/ImageView.h>
#include <vw/Image/PerPixelViews.h>
#include <vw/Image/Statistics.h>

/// Used in blobindex
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

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
  UnaryPerPixelView<ImageT,UnaryCompoundFunctor<ChannelClampFunctor<typename ImageT::pixel_type>, typename ImageT::pixel_type > >
  inline clamp( ImageViewBase<ImageT> const& image, LowT low, HighT high ) {
    typedef UnaryCompoundFunctor<ChannelClampFunctor<typename ImageT::pixel_type>, typename ImageT::pixel_type> func_type;
    func_type func( ChannelClampFunctor<typename ImageT::pixel_type>(low,high) );
    return UnaryPerPixelView<ImageT,func_type>( image.impl(), func );
  }

  /// Clamp the values in an image to fall within the range [0,high].
  /// The low end of the range is actually determined by the
  /// ChannelRange type trait but is generally zero.
  template <class ImageT, class HighT>
  UnaryPerPixelView<ImageT,UnaryCompoundFunctor<ChannelClampFunctor<typename ImageT::pixel_type>, typename ImageT::pixel_type> >
  inline clamp( ImageViewBase<ImageT> const& image, HighT high ) {
    typedef UnaryCompoundFunctor<ChannelClampFunctor<typename ImageT::pixel_type>, typename ImageT::pixel_type> func_type;
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
  UnaryPerPixelView<ImageT,UnaryCompoundFunctor<ChannelClampFunctor<typename ImageT::pixel_type>, typename ImageT::pixel_type> >
  inline clamp( ImageViewBase<ImageT> const& image ) {
    typedef UnaryCompoundFunctor<ChannelClampFunctor<typename ImageT::pixel_type>, typename ImageT::pixel_type> func_type;
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

  /// Renormalize the values in an image to fall within the range
  /// [low,high), but leave the values in the alpha channel untouched.
  template <class ImageT>
  UnaryPerPixelView<ImageT, ChannelNormalizeRetainAlphaFunctor<typename ImageT::pixel_type> >
  inline normalize_retain_alpha( ImageViewBase<ImageT> const& image,
                                 typename ImageChannelType<ImageT>::type old_low,
                                 typename ImageChannelType<ImageT>::type old_high,
                                 typename ImageChannelType<ImageT>::type new_low,
                                 typename ImageChannelType<ImageT>::type new_high  ) {
    typedef ChannelNormalizeRetainAlphaFunctor<typename ImageT::pixel_type> func_type;
    func_type func ( old_low, old_high, new_low, new_high );
    return UnaryPerPixelView<ImageT, func_type >( image.impl(), func );
  }

  /// Renormalize the values in an image to fall within the range [low,high).
  template <class ImageT>
  UnaryPerPixelView<ImageT, UnaryCompoundFunctor<ChannelNormalizeFunctor<typename ImageT::pixel_type>, typename ImageT::pixel_type> >
  inline normalize( ImageViewBase<ImageT> const& image,
                    typename ImageChannelType<ImageT>::type old_low,
                    typename ImageChannelType<ImageT>::type old_high,
                    typename ImageChannelType<ImageT>::type new_low,
                    typename ImageChannelType<ImageT>::type new_high  ) {
    typedef UnaryCompoundFunctor<ChannelNormalizeFunctor<typename ImageT::pixel_type>, typename ImageT::pixel_type> func_type;
    func_type func( ChannelNormalizeFunctor<typename ImageT::pixel_type>( old_low, old_high, new_low, new_high ) );
    return UnaryPerPixelView<ImageT, func_type >( image.impl(), func );
  }

  /// Renormalize the values in an image to fall within the range [low,high).
  template <class ImageT>
  UnaryPerPixelView<ImageT, UnaryCompoundFunctor<ChannelNormalizeFunctor<typename ImageT::pixel_type>, typename ImageT::pixel_type> >
  inline normalize( ImageViewBase<ImageT> const& image,
                    typename ImageChannelType<ImageT>::type low, typename ImageChannelType<ImageT>::type high ) {
    typedef UnaryCompoundFunctor<ChannelNormalizeFunctor<typename ImageT::pixel_type>, typename ImageT::pixel_type> func_type;
    typename ImageChannelType<ImageT>::type old_min, old_max;
    min_max_channel_values( image, old_min, old_max );
    func_type func( ChannelNormalizeFunctor<typename ImageT::pixel_type>( old_min, old_max, low, high ) );
    return UnaryPerPixelView<ImageT, func_type >( image.impl(), func );
  }

  /// Renormalize the values in an image to fall within the range
  /// [0,high).  The low end of the range is actually determined by
  /// the ChannelRange type trait but is generally zero.
  template <class ImageT>
  UnaryPerPixelView<ImageT, UnaryCompoundFunctor<ChannelNormalizeFunctor<typename ImageT::pixel_type>, typename ImageT::pixel_type> >
  inline normalize( ImageViewBase<ImageT> const& image, typename ImageChannelType<ImageT>::type high ) {
    typedef UnaryCompoundFunctor<ChannelNormalizeFunctor<typename ImageT::pixel_type>, typename ImageT::pixel_type> func_type;
    typedef ChannelRange<typename ImageChannelType<ImageT>::type> range_type;
    typename ImageChannelType<ImageT>::type old_min, old_max;
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
  UnaryPerPixelView<ImageT, UnaryCompoundFunctor<ChannelNormalizeFunctor<typename ImageT::pixel_type>, typename ImageT::pixel_type> >
  inline normalize( ImageViewBase<ImageT> const& image ) {
    typedef UnaryCompoundFunctor<ChannelNormalizeFunctor<typename ImageT::pixel_type>, typename ImageT::pixel_type> func_type;
    typedef ChannelRange<typename ImageChannelType<ImageT>::type> range_type;
    typename ImageChannelType<ImageT>::type old_min, old_max;
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
  UnaryPerPixelView<ImageT, UnaryCompoundFunctor<ChannelThresholdFunctor<typename ImageT::pixel_type>, typename ImageT::pixel_type> >
  inline threshold( ImageViewBase<ImageT> const& image, ThreshT thresh, LowT low, HighT high ) {
    typedef UnaryCompoundFunctor<ChannelThresholdFunctor<typename ImageT::pixel_type>, typename ImageT::pixel_type> func_type;
    func_type func( ChannelThresholdFunctor<typename ImageT::pixel_type>(thresh,low,high) );
    return UnaryPerPixelView<ImageT,func_type>( image.impl(), func );
  }

  /// Threshold the values in an image, generating a two-valued output
  /// image with values 0 and high.  The low value is actually
  /// determined by the ChannelRange type trait but is generally zero.
  template <class ImageT, class ThreshT, class HighT>
  UnaryPerPixelView<ImageT,UnaryCompoundFunctor<ChannelThresholdFunctor<typename ImageT::pixel_type>, typename ImageT::pixel_type> >
  inline threshold( ImageViewBase<ImageT> const& image, ThreshT thresh, HighT high ) {
    typedef UnaryCompoundFunctor<ChannelThresholdFunctor<typename ImageT::pixel_type>, typename ImageT::pixel_type> func_type;
    typedef ChannelRange<typename ImageChannelType<ImageT>::type> range_type;
    func_type func( ChannelThresholdFunctor<typename ImageT::pixel_type>(thresh,range_type::min(),high) );
    return UnaryPerPixelView<ImageT,func_type>( image.impl(), func );
  }

  /// Threshold the values in an image, generating a two-valued output
  /// where the values are determined by the ChannelRange type trait
  /// and are generally equal to 0.0 and 1.0 for floating point types
  /// and 0 and the largest positve value for integral types.
  template <class ImageT, class ThreshT>
  UnaryPerPixelView<ImageT,UnaryCompoundFunctor<ChannelThresholdFunctor<typename ImageT::pixel_type>, typename ImageT::pixel_type> >
  inline threshold( ImageViewBase<ImageT> const& image, ThreshT thresh ) {
    typedef UnaryCompoundFunctor<ChannelThresholdFunctor<typename ImageT::pixel_type>, typename ImageT::pixel_type> func_type;
    typedef ChannelRange<typename ImageChannelType<ImageT>::type> range_type;
    func_type func( ChannelThresholdFunctor<typename ImageT::pixel_type>(thresh,range_type::min(),range_type::max()) );
    return UnaryPerPixelView<ImageT,func_type>( image.impl(), func );
  }

  /// Threshold the values in an image against zero, generating a
  /// two-valued output where the values are determined by the
  /// ChannelRange type trait and are generally equal to 0.0 and 1.0
  /// for floating point types and 0 and the largest positve value for
  /// integral types.
  template <class ImageT>
  UnaryPerPixelView<ImageT,UnaryCompoundFunctor<ChannelThresholdFunctor<typename ImageT::pixel_type>, typename ImageT::pixel_type> >
  inline threshold( ImageViewBase<ImageT> const& image ) {
    typedef UnaryCompoundFunctor<ChannelThresholdFunctor<typename ImageT::pixel_type>, typename ImageT::pixel_type> func_type;
    typedef ChannelRange<typename ImageChannelType<ImageT>::type> range_type;
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

    typedef typename SourceT::pixel_accessor src_accessor;
    typedef typename ImageView<OutputT>::pixel_accessor dst_accessor;

    src_accessor srow = src.impl().origin();
    dst_accessor drow = dst.origin();
    const typename SourceT::pixel_type zero = typename SourceT::pixel_type();

    { // First row
      src_accessor scol = srow;
      dst_accessor dcol = drow;
      for( int32 col=cols; col; --col ) {
        *dcol = ((*scol)==zero)?0:1;
        scol.next_col();
        dcol.next_col();
      }
      srow.next_row();
      drow.next_row();
    }
    for( int32 row=rows-2; row; --row ) {
      src_accessor scol = srow;
      dst_accessor dcol = drow;
      *dcol = ((*scol)==zero)?0:1;
      scol.next_col();
      dcol.next_col();
      for( int32 col=cols-2; col; --col ) {
        if( (*scol)==zero ) (*dcol)=0;
        else {
          dst_accessor s1 = dcol, s2 = dcol;
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
      src_accessor scol = srow;
      dst_accessor dcol = drow;
      for( int32 col=cols; col; --col ) {
        *dcol = ((*scol)==zero)?0:1;
        scol.next_col();
        dcol.next_col();
      }
    }
    drow.advance(cols-2,-1);
    for( int32 row=rows-2; row; --row ) {
      dst_accessor dcol = drow;
      for( int32 col=cols-2; col; --col ) {
        if( (*dcol)!=0 ) {
          dst_accessor s1 = dcol, s2 = dcol;
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
    for( int32 y=0; y<image.rows(); ++y )
      for( int32 x=0; x<image.cols(); ++x )
        if( ! (is_opaque( image(x,y) ) ) )
          return false;
    return true;
  }

  template <class ImageT>
  bool is_opaque_helper( ImageT const& /*image*/, false_type ) {
    return true;
  }

  /// Returns true if the given image is entirely opaque, or false if
  /// it is at least partially transparent.
  template <class ImageT>
  bool is_opaque( ImageViewBase<ImageT> const& image ) {
    return is_opaque_helper( image.impl(), typename PixelHasAlpha<typename ImageT::pixel_type>::type() );
  }


  // *******************************************************************
  // is_transparent()
  // *******************************************************************

  template <class ImageT>
  bool is_transparent_helper( ImageT const& image, true_type ) {
    for( int32 y=0; y<image.rows(); ++y )
      for( int32 x=0; x<image.cols(); ++x )
        if( ! is_transparent(image(x,y)) ) return false;
    return true;
  }

  template <class ImageT>
  bool is_transparent_helper( ImageT const& /*image*/, false_type ) {
    return false;
  }

  /// Returns true if the given image is entirely transparent, or false if
  /// it is opaque or only partially transparent.
  template <class ImageT>
  bool is_transparent( ImageViewBase<ImageT> const& image ) {
    return is_transparent_helper( image.impl(), typename PixelHasAlpha<typename ImageT::pixel_type>::type() );
  }


  // *******************************************************************
  // clear_nonopaque_pixels()
  //
  // This filter is useful for eliminating fringe effects along the
  // edges of images with some transparent or nodata values that have
  // be transformed with bilinear or bicubic interpolation.
  // *******************************************************************
  template <class PixelT>
  class ClearNonOpaqueFunctor: public UnaryReturnSameType {
  public:
    ClearNonOpaqueFunctor() {}

    PixelT operator()( PixelT const& value ) const {
      if (is_opaque(value)) return value;
      else return PixelT();
    }
  };

  /// Zero out any pixels that aren't completely opaque.
  template <class ImageT>
  UnaryPerPixelView<ImageT,ClearNonOpaqueFunctor<typename ImageT::pixel_type> >
  inline clear_nonopaque_pixels( ImageViewBase<ImageT> const& image ) {
    typedef ClearNonOpaqueFunctor<typename ImageT::pixel_type> func_type;
    return UnaryPerPixelView<ImageT,func_type>( image.impl(), func_type() );
  }

  // *******************************************************************
  // remap_pixel_value()
  //
  // This filter can be used to map one pixel value to another.  This
  // can be useful in many situations, for example when you need to
  // remap the nodata value used in a DEM.
  // *******************************************************************
  template <class PixelT>
  class RemapPixelFunctor: public UnaryReturnSameType {
    typename PixelChannelType<PixelT>::type m_src_val, m_dst_val;
  public:
    RemapPixelFunctor(typename PixelChannelType<PixelT>::type src_val,
                      typename PixelChannelType<PixelT>::type dst_val) :
      m_src_val(src_val), m_dst_val(dst_val) {}

    PixelT operator()( PixelT const& value ) const {
      if (value == m_src_val) return m_dst_val;
      else return value;
    }
  };

  /// Zero out any pixels that aren't completely opaque.
  template <class ImageT>
  UnaryPerPixelView<ImageT,RemapPixelFunctor<typename ImageT::pixel_type> >
  inline remap_pixel_value( ImageViewBase<ImageT> const& image,
                            typename PixelChannelType<typename ImageT::pixel_type>::type src_val,
                            typename PixelChannelType<typename ImageT::pixel_type>::type dst_val) {
    typedef RemapPixelFunctor<typename ImageT::pixel_type> func_type;
    return UnaryPerPixelView<ImageT,func_type>( image.impl(), func_type(src_val, dst_val) );
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
  /// It will operate on any object that has cols() and rows() methods.

  template <class T>
  inline std::vector<BBox2i> image_blocks(T const& object, int32 block_width, int32 block_height) {
    std::vector<BBox2i> bboxes;

    int32 j_offset = 0;
    while ( j_offset < object.rows() ) {
      int32 j_dim = (object.rows() - j_offset) < block_height ? (object.rows() - j_offset) : block_height;
      int32 i_offset = 0;
      while ( i_offset < object.cols() ) {
        int32 i_dim = (object.cols() - i_offset) < block_width ? (object.cols() - i_offset) : block_width;
        bboxes.push_back(BBox2i(i_offset,j_offset,i_dim,j_dim));
        i_offset += i_dim;
      }
      j_offset += j_dim;
    }
    return bboxes;
  }

  // ********************************************************************
  // blob_index()
  // ********************************************************************

  /// A utility that numbers off blobs of valid pixels. So it's
  /// important that you use pixel mask, or the entire image will be
  /// numbered one. This is for the most part a clone of bwlabel.
  ///

  template <class SourceT>
  void blob_index( ImageViewBase<SourceT> const& src,
                   ImageView<uint32>& dst ) {

    if ( src.impl().planes() > 1 )
      vw_throw( NoImplErr() << "Blob index currently only works with 2D images." );
    dst.set_size( src.impl().cols(),
                  src.impl().rows() );
    fill(dst,0);

    // Initialize Graph used to pair blobs
    typedef boost::adjacency_list<boost::vecS,boost::vecS,boost::undirectedS> Graph;
    Graph connections;
    uint32 m_blob_count=1;

    { // Initial Pass
      typename SourceT::pixel_accessor s_acc = src.impl().origin();
      typename SourceT::pixel_accessor p_s_acc = src.impl().origin(); // previous
      typename ImageView<uint32>::pixel_accessor d_acc = dst.origin();
      typename ImageView<uint32>::pixel_accessor p_d_acc = dst.origin(); // previous

      // Top Corner
      if ( is_valid(*s_acc) ) {
        *d_acc = m_blob_count;
        m_blob_count++;
      }

      // Top Row
      s_acc.next_col();
      d_acc.next_col();
      for ( int32 i = 1; i < dst.cols(); i++ ) {
        if ( is_valid(*s_acc) ) {
          if ( is_valid(*p_s_acc) ) {
            *d_acc = *p_d_acc;
          } else {
            *d_acc = m_blob_count;
            m_blob_count++;
          }
          }
        s_acc.next_col();
        d_acc.next_col();
        p_s_acc.next_col();
        p_d_acc.next_col();
      }
    }

    { // Everything else (9 connected)
      typename SourceT::pixel_accessor s_acc_row = src.impl().origin();
      typename ImageView<uint32>::pixel_accessor d_acc_row = dst.origin();
      s_acc_row.advance(0,1);
      d_acc_row.advance(0,1);

      for (int j = dst.rows()-1; j; --j ) { // Not for indexing
        typename SourceT::pixel_accessor s_acc = s_acc_row;
        typename SourceT::pixel_accessor p_s_acc = s_acc_row;
        typename ImageView<uint32>::pixel_accessor d_acc = d_acc_row;
        typename ImageView<uint32>::pixel_accessor p_d_acc = d_acc_row;

        // Process
        for ( int i = dst.cols(); i; --i ) {
          if ( is_valid(*s_acc) ) {
            if ( i != dst.cols() ) {
              // Left
              p_s_acc.advance(-1,0);
              p_d_acc.advance(-1,0);
              if ( is_valid(*p_s_acc) ) {
                if ( (*d_acc != 0) && (*d_acc != *p_d_acc) ) {
                  boost::add_edge(*p_d_acc,*d_acc,connections);
                } else
                  *d_acc = *p_d_acc;
              }
              // Upper Left
              p_s_acc.advance(0,-1);
              p_d_acc.advance(0,-1);
              if ( is_valid(*p_s_acc) ) {
                if ( (*d_acc != 0) && (*d_acc != *p_d_acc) ) {
                  boost::add_edge(*p_d_acc,*d_acc,connections);
                } else
                  *d_acc = *p_d_acc;
              }
            } else {
              p_s_acc.advance(-1,-1);
              p_d_acc.advance(-1,-1);
            }
            // Upper
            p_s_acc.advance(1,0);
            p_d_acc.advance(1,0);
            if ( is_valid(*p_s_acc) ) {
              if ( (*d_acc != 0) && (*d_acc != *p_d_acc) ) {
                boost::add_edge(*p_d_acc,*d_acc,connections);
              } else
                *d_acc = *p_d_acc;
            }
            // Upper Right
            p_s_acc.advance(1,0);
            p_d_acc.advance(1,0);
            if ( i != 1 )
              if ( is_valid(*p_s_acc) ) {
                if ( (*d_acc != 0) && (*d_acc != *p_d_acc) ) {
                  boost::add_edge(*p_d_acc,*d_acc,connections);
                } else
                  *d_acc = *p_d_acc;
              }
            // Setting if not
            p_s_acc.advance(-1,1);
            p_d_acc.advance(-1,1);
            if ( *d_acc == 0 ) {
              *d_acc = m_blob_count;
              m_blob_count++;
            }
          }
          s_acc.next_col();
          p_s_acc.next_col();
          d_acc.next_col();
          p_d_acc.next_col();
        } // end row process

        s_acc_row.next_row();
        d_acc_row.next_row();
      }
    }

    // Making sure connections has vertices for all indexes made
    add_edge(m_blob_count-1,m_blob_count-1,connections);
    std::vector<uint32> component(boost::num_vertices(connections));
    m_blob_count = boost::connected_components(connections, &component[0])-1;

    { // Update index map to optimal numbering
      ImageView<uint32>::pixel_accessor p_d_acc = dst.origin(); // previous

      for ( int32 r = 0; r < dst.rows(); r++ ) {
        ImageView<uint32>::pixel_accessor d_acc = p_d_acc;
        for ( int32 c = 0; c < dst.cols(); c++ ) {
          if ( (*d_acc) != 0 ) {

            *d_acc = component[*d_acc];
          }
          d_acc.next_col();
        }
        p_d_acc.next_row();
      }
    }
  }

  // Simple interface
  template <class SourceT>
  ImageView<uint32> blob_index( ImageViewBase<SourceT> const& src ) {
    ImageView<uint32> result( src.impl().cols(), src.impl().rows() );
    blob_index( src, result );
    return result;
  }

} // namespace vw

#endif // __VW_IMAGE_ALGORITHMS_H__
