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

// ******************************************************************
// MeanFillTransparent
// ******************************************************************

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
