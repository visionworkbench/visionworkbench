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
clamp( ImageViewBase<ImageT> const& image, LowT low, HighT high ) {
  typedef UnaryCompoundFunctor<ChannelClampFunctor<typename ImageT::pixel_type>, typename ImageT::pixel_type> func_type;
  func_type func( ChannelClampFunctor<typename ImageT::pixel_type>(low,high) );
  return UnaryPerPixelView<ImageT,func_type>( image.impl(), func );
}

/// Clamp the values in an image to fall within the range [0,high].
/// The low end of the range is actually determined by the
/// ChannelRange type trait but is generally zero.
template <class ImageT, class HighT>
UnaryPerPixelView<ImageT,UnaryCompoundFunctor<ChannelClampFunctor<typename ImageT::pixel_type>, typename ImageT::pixel_type> >
clamp( ImageViewBase<ImageT> const& image, HighT high ) {
  typedef UnaryCompoundFunctor<ChannelClampFunctor<typename ImageT::pixel_type>, typename ImageT::pixel_type> func_type;
  typedef ChannelRange<typename CompoundChannelType<typename ImageT::pixel_type>::type> range_type;
  typename CompoundChannelType<typename ImageT::pixel_type>::type min_val = range_type::min();
  func_type func( ChannelClampFunctor<typename ImageT::pixel_type>(min_val,high) );
  return UnaryPerPixelView<ImageT,func_type>( image.impl(), func );
}

/// Clamp the values in an image to fall within the range [min,max],
/// where min and max are determined by the ChannelRange type trait
/// and are generally equal to 0.0 and 1.0 for floating point types
/// and 0 and the largest positive value for integral types.
template <class ImageT>
UnaryPerPixelView<ImageT,UnaryCompoundFunctor<ChannelClampFunctor<typename ImageT::pixel_type>, typename ImageT::pixel_type> >
clamp( ImageViewBase<ImageT> const& image ) {
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
normalize_retain_alpha( ImageViewBase<ImageT> const& image,
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
normalize( ImageViewBase<ImageT> const& image,
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
normalize( ImageViewBase<ImageT> const& image,
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
normalize( ImageViewBase<ImageT> const& image, typename ImageChannelType<ImageT>::type high ) {
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
/// point types and 0 and the largest positive value for integral types.
template <class ImageT>
UnaryPerPixelView<ImageT, UnaryCompoundFunctor<ChannelNormalizeFunctor<typename ImageT::pixel_type>, typename ImageT::pixel_type> >
normalize( ImageViewBase<ImageT> const& image ) {
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
threshold( ImageViewBase<ImageT> const& image, ThreshT thresh, LowT low, HighT high ) {
  typedef UnaryCompoundFunctor<ChannelThresholdFunctor<typename ImageT::pixel_type>, typename ImageT::pixel_type> func_type;
  func_type func( ChannelThresholdFunctor<typename ImageT::pixel_type>(thresh,low,high) );
  return UnaryPerPixelView<ImageT,func_type>( image.impl(), func );
}

/// Threshold the values in an image, generating a two-valued output
/// image with values 0 and high.  The low value is actually
/// determined by the ChannelRange type trait but is generally zero.
template <class ImageT, class ThreshT, class HighT>
UnaryPerPixelView<ImageT,UnaryCompoundFunctor<ChannelThresholdFunctor<typename ImageT::pixel_type>, typename ImageT::pixel_type> >
threshold( ImageViewBase<ImageT> const& image, ThreshT thresh, HighT high ) {
  typedef UnaryCompoundFunctor<ChannelThresholdFunctor<typename ImageT::pixel_type>, typename ImageT::pixel_type> func_type;
  typedef ChannelRange<typename ImageChannelType<ImageT>::type> range_type;
  func_type func( ChannelThresholdFunctor<typename ImageT::pixel_type>(thresh,range_type::min(),high) );
  return UnaryPerPixelView<ImageT,func_type>( image.impl(), func );
}

/// Threshold the values in an image, generating a two-valued output
/// where the values are determined by the ChannelRange type trait
/// and are generally equal to 0.0 and 1.0 for floating point types
/// and 0 and the largest positive value for integral types.
template <class ImageT, class ThreshT>
UnaryPerPixelView<ImageT,UnaryCompoundFunctor<ChannelThresholdFunctor<typename ImageT::pixel_type>, typename ImageT::pixel_type> >
threshold( ImageViewBase<ImageT> const& image, ThreshT thresh ) {
  typedef UnaryCompoundFunctor<ChannelThresholdFunctor<typename ImageT::pixel_type>, typename ImageT::pixel_type> func_type;
  typedef ChannelRange<typename ImageChannelType<ImageT>::type> range_type;
  func_type func( ChannelThresholdFunctor<typename ImageT::pixel_type>(thresh,range_type::min(),range_type::max()) );
  return UnaryPerPixelView<ImageT,func_type>( image.impl(), func );
}

/// Threshold the values in an image against zero, generating a
/// two-valued output where the values are determined by the
/// ChannelRange type trait and are generally equal to 0.0 and 1.0
/// for floating point types and 0 and the largest positive value for
/// integral types.
template <class ImageT>
UnaryPerPixelView<ImageT,UnaryCompoundFunctor<ChannelThresholdFunctor<typename ImageT::pixel_type>, typename ImageT::pixel_type> >
threshold( ImageViewBase<ImageT> const& image ) {
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
clear_nonopaque_pixels( ImageViewBase<ImageT> const& image ) {
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
remap_pixel_value( ImageViewBase<ImageT> const& image,
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

// ******************************************************************
// FillHoles
// ******************************************************************

// A class to fill holes in images. Given the number
// hole_fill_len, which gives the size of holes, look left, right,
// up, and down, as far as hole_fill_len/2. Must have valid pixels
// both left and right, otherwise both up and down, else do
// nothing. The motivation here is to fill in only pixels
// "surrounded" by valid pixels.

// Use one of the two approaches:
//
// 1. Interpolate using the left and right values. Interpolate using
// the up and down values. Average the results. Fast.
//
// 2. Find the weighted average of the points in the window
// of size hole_fill_len centered at the current point. Slow.

// After holes are filled, do several passes to average the results.

// The image being passed in should be a PixelMask.
template <class ImageT>
class FillHoles: public ImageViewBase< FillHoles<ImageT> > {
  ImageT m_img;
  int m_hole_fill_mode, m_hole_fill_num_smooth_iter;
  int m_hole_fill_half, m_smooth_half;
  ImageView<double> m_dist_kernel, m_gauss_kernel;
  typedef typename ImageT::pixel_type PixelT;

public:

  typedef PixelT pixel_type;
  typedef PixelT result_type;
  typedef ProceduralPixelAccessor<FillHoles> pixel_accessor;

  FillHoles( ImageViewBase<ImageT> const& img, int hole_fill_mode,
             int hole_fill_num_smooth_iter,
             int hole_fill_len) :
    m_img(img.impl()), m_hole_fill_mode(hole_fill_mode),
    m_hole_fill_num_smooth_iter(hole_fill_num_smooth_iter),
    m_hole_fill_half((hole_fill_len+1)/2), m_smooth_half(2) {

    if (m_hole_fill_mode == 2){
      // Use the inverse square distance kernel
      int h = m_hole_fill_half;
      m_dist_kernel.set_size(2*h+1, 2*h+1);
      for (int c = 0; c < m_dist_kernel.cols(); c++){
        for (int r = 0; r < m_dist_kernel.rows(); r++){
          double r2 = double(c-h)*(c-h) + double(r-h)*(r-h);
          m_dist_kernel(c, r) = 1.0/std::max(r2, 1.0); 
        }
      }
    }

    if (m_hole_fill_num_smooth_iter > 0){
      // Use a gaussian kernel
      int h = m_smooth_half;
      m_gauss_kernel.set_size(2*h+1, 2*h+1);
      double val = 0.25; // value to reach at kernel edge
      double sigma = -log(val)/double(h*h);
      for (int c = 0; c < m_gauss_kernel.cols(); c++){
        for (int r = 0; r < m_gauss_kernel.rows(); r++){
          double r2 = double(c-h)*(c-h) + double(r-h)*(r-h);
          m_gauss_kernel(c, r) = exp(-sigma*r2);
        }
      }
    }
  }
  
  inline int32 cols  () const { return m_img.cols(); }
  inline int32 rows  () const { return m_img.rows(); }
  inline int32 planes() const { return 1; }

  inline pixel_accessor origin() const { return pixel_accessor(*this); }

  inline result_type operator()( size_t i, size_t j, size_t p=0 ) const {
    vw_throw( NoImplErr() << "FillHoles: operator() not implemented.\n" );
  }
  
  typedef CropView< ImageView<PixelT> > prerasterize_type;
  inline prerasterize_type prerasterize( BBox2i const& bbox ) const {

    // Crop into an expanded box as to have enough pixels to do
    // averaging with given window at every pixel in the current box.
    int h = m_hole_fill_half; // shorten
    BBox2i biased_box = bbox;
    biased_box.expand(h+1);
    biased_box.crop(bounding_box(m_img));
    ImageView<PixelT> img( crop( m_img, biased_box ) );
    ImageView<PixelT> filled_img = copy(img);
    int nc = img.cols(), nr = img.rows(); // shorten

    for (int row = 0; row < nr; row++){
      for (int col = 0; col < nc; col++){

        if (is_valid(img(col, row))) continue; // skip valid
      
        // Look left, right, up, down, and find the closest valid pixels
        double r0 = -1, r1 = -1, c0 = -1, c1 = -1; // no good indices yet
        for (int k = row-1; k >= std::max(row-h, 0); k--)
          if (is_valid(img(col, k))){ r0 = k; break; }
        if (r0 >=0){
          // Found a point to the left, try to also find one to the right
          for (int k = row+1; k <= std::min(row + h, nr-1); k++)
            if (is_valid(img(col, k))){ r1 = k; break; }
        }
        for (int k = col-1; k >= std::max(col-h, 0); k--)
          if (is_valid(img(k, row))){ c0 = k; break; }
        if (c0 >=0){
          // Found a point up, try to also find one down
          for (int k = col+1; k <= std::min(col + h, nc-1); k++)
            if (is_valid(img(k, row))){ c1 = k; break; }
        }

        // skip if no good neighbors
        if ( (r0 < 0 || r1 < 0) && (c0 < 0 || c1 < 0) ) continue;

        double sum = 0.0;
        PixelT V; V.validate();
        if (m_hole_fill_mode == 1){
          // Interpolate between left and right, then between top and
          // bottom.  Average the results.
          if (r0 >= 0 && r1 >= 0){
            V += ((r1-row)*img(col, r0) + (row-r0)*img(col, r1))/double(r1-r0);
            sum++;
          }
          if (c0 >= 0 && c1 >= 0){
            V += ((c1-col)*img(c0, row) + (col-c0)*img(c1, row))/double(c1-c0);
            sum++;
          }
          
        }else{
          // Weighted average with given kernel
          for (int c = std::max(col-h, 0); c <= std::min(col+h, nc-1); c++){
            for (int r = std::max(row-h, 0); r <= std::min(row+h, nr-1); r++){
              if (!is_valid(img(c, r))) continue;
              double wt = m_dist_kernel(c-col+h, r-row+h);
              V   += wt*img(c, r);
              sum += wt;
            }
          }
        }
        
        if (sum > 0) filled_img(col, row) = V/sum;
        
      }
    }

    // Smooth the resulting image by repeated convolutions
    // with a small gaussian kernel
    for (int i = 0; i < m_hole_fill_num_smooth_iter; i++){
      
      ImageView<PixelT> curr_img = copy(filled_img);

      for (int row = 0; row < nr; row++){
        for (int col = 0; col < nc; col++){
          if (is_valid(img(col, row))) continue; // skip valid
          if (!is_valid(filled_img(col, row))) continue; // don't add more valid

          int nh = m_smooth_half;
          double sum = 0.0;
          PixelT V; V.validate();
          for (int c = std::max(col-nh, 0); c <= std::min(col+nh, nc-1); c++){
            for (int r = std::max(row-nh, 0); r <= std::min(row+nh, nr-1); r++){
              if (!is_valid(curr_img(c, r))) continue;
              double wt = m_gauss_kernel(c-col+nh, r-row+nh);
              V   += wt*curr_img(c, r);
              sum += wt;
            }
          }
          if (sum > 0) filled_img(col, row) = V/sum;
          
        }
      }
      
    } // end smoothing iterations
    
  return prerasterize_type(filled_img,
                           -biased_box.min().x(), -biased_box.min().y(),
                           cols(), rows());
  
  }
  template <class ImgT>
  inline void rasterize( ImgT const& img, BBox2i const& bbox ) const {
    vw::rasterize( prerasterize(bbox), img, bbox );
  }

};


// ----------------------------------------------------------------------------
/*
//  compute_normals()
//
// Compute a vector normal to the surface of a DEM for each given
// pixel.  The normal is computed by forming a plane with three points
// in the vicinity of the requested pixel, and then finding the vector
// normal to that plane.  The user must specify the scale in the [u,v]
// directions so that the direction of the vector in physical space
// can be properly ascertained.  This is often contained in the (0,0)
// and (1,1) entry of the georeference transform.
class ComputeNormalsFunc : public ReturnFixedType<PixelMask<Vector3f> >
{
  float m_u_scale, m_v_scale;

public:
  ComputeNormalsFunc(float u_scale, float v_scale) :
    m_u_scale(u_scale), m_v_scale(v_scale) {}

  BBox2i work_area() const { return BBox2i(Vector2i(0, 0), Vector2i(1, 1)); }

  template <class PixelAccessorT>
  PixelMask<Vector3f> operator() (PixelAccessorT const& accessor_loc) const {
    PixelAccessorT acc = accessor_loc;

    // Pick out the three altitude values.
    if (is_transparent(*acc))
      return PixelMask<Vector3f>();
    float alt1 = *acc;

    acc.advance(1,0);
    if (is_transparent(*acc))
      return PixelMask<Vector3f>();
    float alt2 = *acc;

    acc.advance(-1,1);
    if (is_transparent(*acc))
      return PixelMask<Vector3f>();
    float alt3 = *acc;

    // Form two orthogonal vectors in the plane containing the three
    // altitude points
    Vector3f n1(m_u_scale, 0, alt2-alt1);
    Vector3f n2(0, m_v_scale, alt3-alt1);

    // Return the vector normal to the local plane.
    return normalize(cross_prod(n1,n2));
  }
}; // End class ComputeNormalsFunc


/// Perform the dot product between each pixel and a constant vector.
class DotProdFunc : public ReturnFixedType<PixelMask<PixelGray<float> > > {
  Vector3f m_vec;
public:
  DotProdFunc(Vector3f const& vec) : m_vec(normalize(vec)) {}
  PixelMask<PixelGray<float> > operator() (PixelMask<Vector3f> const& pix) const {
    if (is_transparent(pix))
      return PixelMask<PixelGray<float> >();
    else
      return dot_prod(pix.child(),m_vec);
  }
};
*/
/// Apply a double threshold to an image.
/// - Pixels are set if they are above the high threshold.  In addition, a flood-fill is performed
///   from the high threshold pixels using pixels above the low threshold.
/// - No built-in masked pixel handling.
template <class ImageT>
class TwoThresholdFill: public ImageViewBase<TwoThresholdFill<ImageT> >{

  ImageT const& m_image;
  int    m_expand_size; ///< Tile expansion used to more accurately size flood.
  double m_low_threshold;   ///
  double m_high_threshold;  ///
  uint8  m_output_false;    /// Above thresholds
  uint8  m_output_true;     /// Below thresholds
public:
  /// Constructor
  /// - Pixels below thresholds will be set to output_true, others to output_false.
  TwoThresholdFill(ImageViewBase<ImageT> const& image, int expand_size, 
                   double low_threshold, double high_threshold,
                   uint8 output_false = 0, uint8 output_true = 1):
    m_image(image.impl()), m_expand_size(expand_size), 
    m_low_threshold(low_threshold), m_high_threshold(high_threshold),
    m_output_false(output_false), m_output_true(output_true) {}

  // Image View interface
  typedef uint8      pixel_type;
  typedef pixel_type result_type;
  typedef ProceduralPixelAccessor<TwoThresholdFill> pixel_accessor;

  inline int32 cols  () const { return m_image.cols(); }
  inline int32 rows  () const { return m_image.rows(); }
  inline int32 planes() const { return 1; }

  inline pixel_accessor origin() const { return pixel_accessor( *this, 0, 0 ); }

  inline pixel_type operator()( double i, double j, int32 p = 0 ) const {
    vw_throw(NoImplErr() << "operator()(...) is not implemented");
    return pixel_type();
  }

  typedef CropView<ImageView<pixel_type> > prerasterize_type;
  inline prerasterize_type prerasterize(BBox2i const& bbox) const {

    // Rasterize the region with an expanded tile size so that we can count
    //  flood fill some distance outside the tile borders.
    BBox2i big_bbox = bbox;
    big_bbox.expand(m_expand_size);
    big_bbox.crop(bounding_box(m_image));
    
    ImageView<pixel_type> output_tile(big_bbox.width(), big_bbox.height());
    ValueEdgeExtension<pixel_type> edge_wrapper(pixel_type(0)); // Edge handling for the output image
    
    ImageView<typename ImageT::pixel_type> input_tile = crop(m_image, big_bbox);   
    

    // Doing the flood fill is basically a blob labeling problem and can be done in two passes through the tile.
    // First pass is top-left to bottom-right.
    for (int r=0; r<output_tile.rows(); ++r) {
      for (int c=0; c<output_tile.cols(); ++c) {
        if ((input_tile(c,r) > m_high_threshold) || // Check the lower threshold
            ((input_tile(c,r) > m_low_threshold) && // Check neighboring inputs if below the higher threshold
             ((edge_wrapper(output_tile, c-1,r-1) > 0) || // Check top left
              (edge_wrapper(output_tile, c,  r-1) > 0) || // Check top
              (edge_wrapper(output_tile, c+1,r-1) > 0) || // Check top right
              (edge_wrapper(output_tile, c-1,r  ) > 0)    // Check left
             )
            )
           )
          output_tile(c,r) = m_output_true;
        else
          output_tile(c,r) = m_output_false;
      }
    } // Done with first pass through the tile
    
    // Second pass is bottom_right to top_left
    for (int r=output_tile.rows()-1; r>=0; --r) {
      for (int c=output_tile.cols()-1; c>=0; --c) {
        if (output_tile(c,r) == m_output_true) // Skip set pixels
          continue;
        // No need to check the lower threshold on the second pass.
        if ((input_tile(c,r) > m_low_threshold) && // Check neighboring inputs only if this pixel is low enough
            ((edge_wrapper(output_tile, c+1,r+1) > 0) || // Check bottom right
             (edge_wrapper(output_tile, c,  r+1) > 0) || // Check bottom
             (edge_wrapper(output_tile, c-1,r+1) > 0) || // Check bottom left
             (edge_wrapper(output_tile, c+1,r  ) > 0)    // Check right
            )
           )
          output_tile(c,r) = m_output_true;
      }
    } // Done with second pass through the tile    

    // Perform tile size faking trick to make small tile look the size of the entire image
    return prerasterize_type(output_tile,
                             -big_bbox.min().x(), -big_bbox.min().y(),
                             cols(), rows() );
  }

  template <class DestT>
  inline void rasterize(DestT const& dest, BBox2i bbox) const {
    vw::rasterize(prerasterize(bbox), dest, bbox);
  }
}; // End class FloodFill



} // namespace vw
