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


/// \file Algorithms.h
///
/// Basic algorithms operating on images. This includes only lazy view
/// implementations.
///
#ifndef __VW_IMAGE_ALGORITHMS_H__
#define __VW_IMAGE_ALGORITHMS_H__

#include <vw/Image/AlgorithmFunctions.h>
#include <vw/Image/PerPixelViews.h>
#include <vw/Image/Statistics.h>
#include <vw/Image/Manipulation.h>
#include <vw/Image/EdgeExtension.h>

namespace vw {

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
    typedef typename PixelWithoutAlpha<PixelT>::type non_alpha_type;
    typedef ChannelNormalizeFunctor<non_alpha_type> norm_func_type;
    UnaryCompoundFunctor<norm_func_type, non_alpha_type> m_compound_func;
  public:
    ChannelNormalizeRetainAlphaFunctor( channel_type old_min, channel_type old_max,
                                        channel_type new_min, channel_type new_max )
      : m_compound_func( norm_func_type( old_min, old_max, new_min, new_max ) ) {}

    PixelT operator()( PixelT value ) const {
      if (is_transparent(value)) return value;
      else {
        PixelT result;
        non_alpha_channels(result) = m_compound_func( non_alpha_channels( value ) );
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

  // ******************************************************************
  // MeanFillTransparent
  // ******************************************************************

  // This is a preprocess step that set the value of transparent
  // pixels to the mean of the nearby opaque pixels. This will not
  // produce a visible difference to the image as it only modifies
  // completely transparent pixels. The reason for this is to remove a
  // "bath tub ring" that happens when interpolating/resampling an
  // image with transparent sections.

  template <class ImageT>
  class MeanFillTransparent : public ImageViewBase<MeanFillTransparent<ImageT> > {
    ImageT m_image;

    template <class SrcAccessT>
    typename SrcAccessT::pixel_type
    inline accumulate_mean( SrcAccessT const& src ) const {
      typedef typename SrcAccessT::pixel_type result_type;
      typedef typename CompoundChannelType<result_type>::type channel_type;
      typedef typename PixelWithoutAlpha<result_type>::type non_a_type;
      typedef typename AccumulatorType<channel_type>::type acc_type;
      typedef typename PixelChannelCast<non_a_type,acc_type>::type non_a_acc_type;
      non_a_acc_type sum_value;
      acc_type weight = 0;

      SrcAccessT px = src;
      px.next_col();
      sum_value += non_a_acc_type(non_alpha_channels(*px))*acc_type(alpha_channel(*px));
      weight += acc_type(alpha_channel(*px));
      px.next_row();
      sum_value += non_a_acc_type(non_alpha_channels(*px))*acc_type(alpha_channel(*px));
      weight += acc_type(alpha_channel(*px));
      px.prev_col();
      sum_value += non_a_acc_type(non_alpha_channels(*px))*acc_type(alpha_channel(*px));
      weight += acc_type(alpha_channel(*px));
      px.prev_col();
      sum_value += non_a_acc_type(non_alpha_channels(*px))*acc_type(alpha_channel(*px));
      weight += acc_type(alpha_channel(*px));
      px.prev_row();
      sum_value += non_a_acc_type(non_alpha_channels(*px))*acc_type(alpha_channel(*px));
      weight += acc_type(alpha_channel(*px));
      px.prev_row();
      sum_value += non_a_acc_type(non_alpha_channels(*px))*acc_type(alpha_channel(*px));
      weight += acc_type(alpha_channel(*px));
      px.next_col();
      sum_value += non_a_acc_type(non_alpha_channels(*px))*acc_type(alpha_channel(*px));
      weight += acc_type(alpha_channel(*px));
      px.next_col();
      sum_value += non_a_acc_type(non_alpha_channels(*px))*acc_type(alpha_channel(*px));
      weight += acc_type(alpha_channel(*px));

      if ( weight <= 0 )
        return result_type();

      result_type result(sum_value / weight);
      alpha_channel( result ) = ChannelRange<channel_type>::min();
      return result;
    }

  public:
    typedef typename ImageT::pixel_type pixel_type;
    typedef pixel_type result_type;
    typedef ProceduralPixelAccessor<MeanFillTransparent > pixel_accessor;

    MeanFillTransparent( ImageT const& image ) : m_image( image ) {}

    inline int32 cols() const { return m_image.cols(); }
    inline int32 rows() const { return m_image.rows(); }
    inline int32 planes() const { return m_image.planes(); }
    inline pixel_accessor origin() const { return pixel_accessor(*this); }

    inline result_type helper( int32 x, int32 y, int32 p, true_type ) const {
      if ( is_transparent(m_image(x,y,p)) ) {
        if ( x > 1 && y > 1 && x + 1 < cols() && y + 1 < rows() )
          return accumulate_mean( m_image.origin().advance(x, y, p ) );
        else
          return accumulate_mean( edge_extend(m_image, ConstantEdgeExtension()).origin().advance(x,y,p) );
      }
      return m_image(x,y,p);
    }

    inline result_type helper( int32 x, int32 y, int32 p, false_type ) const {
      return m_image(x,y,p);
    }

    inline result_type operator()( int32 x, int32 y, int32 p=0 ) const {
      return helper( x, y, p, typename PixelHasAlpha<pixel_type>::type() );
    }

    typedef MeanFillTransparent<CropView<ImageView<result_type> > > prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i const& bbox ) const {
      BBox2i actual = bbox;
      actual.expand(1);
      ImageView<result_type> src =
        edge_extend( m_image, actual, ConstantEdgeExtension() );
      return prerasterize_type( crop( src, -actual.min()[0], -actual.min()[1],
                                      cols(), rows() ) );
    }

    template <class DestT>
    inline void rasterize( DestT const& dest, BBox2i const& bbox ) const {
      vw::rasterize( prerasterize(bbox), dest, bbox );
    }
  };

  template <class SourceT>
  MeanFillTransparent<SourceT>
  inline mean_fill_transparent( ImageViewBase<SourceT> const& src ) {
    return MeanFillTransparent<SourceT>( src.impl() );
  }

} // namespace vw

#endif // __VW_IMAGE_ALGORITHMS_H__
