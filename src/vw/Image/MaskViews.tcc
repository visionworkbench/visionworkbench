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
  /// create_mask( view, value )
  ///
  /// Given a view with pixels of type PixelT and a pixel value to
  /// consider as the "no data" or masked value, returns a view with
  /// pixels that are of the PixelMask<PixelT>, with the appropriate
  /// pixels masked.
  /// - Should safely accept inputs which already contain a mask, in 
  ///   which case the input mask is ignored.
  
  // *******************************************************************
  /// apply_mask( view, value )
  ///
  /// Given a view with pixels of the type PixelMask<T>, this view
  /// returns an image with pixels of type T where any pixel that was
  /// marked as "invalid" in the mask is replaced with the constant
  /// pixel value passed in as value.  The value is T() by default.
  ///
  template <class PixelT>
  class ApplyPixelMask : public ReturnFixedType<PixelT> {
    PixelT m_nodata_value;
  public:
    ApplyPixelMask( PixelT const& nodata_value ) : m_nodata_value(nodata_value) {}
    inline PixelT operator()( PixelMask<PixelT> const& value ) const {
      return value.valid() ? value.child() : m_nodata_value;
    }
  };

  // *******************************************************************
  /// copy_mask(view, mask)
  ///
  /// Copies a mask from one image to another.
  ///
  template <class PixelT>
  class CopyPixelMask : public ReturnFixedType<typename MaskedPixelType<PixelT>::type> {
  public:
    template <class MaskPixelT>
    inline typename MaskedPixelType<PixelT>::type
    operator()( PixelT const& value, MaskPixelT const& mask ) const {
      typename MaskedPixelType<PixelT>::type result = value;
      if (is_transparent(mask)) {
        result.invalidate();
      }
      return result;
    }
  };

  // *******************************************************************
  /// mask_to_alpha(view)
  ///
  /// Converts a mask channel to an alpha channel, generating an image that
  /// is transparent wherever the data is masked.
  ///
  template <class PixelT>
  class MaskToAlpha : public ReturnFixedType<typename PixelWithAlpha<typename UnmaskedPixelType<PixelT>::type>::type> {
  public:
    typedef typename PixelWithAlpha<typename UnmaskedPixelType<PixelT>::type>::type result_type;
    inline result_type operator()( PixelT const& pixel ) const {
      if (is_transparent(pixel)) {
        return result_type();
      }
      else return result_type(pixel.child());
    }
  };

  // *******************************************************************
  /// alpha_to_mask(view)
  ///
  /// Converts an channel to a mask channel, generating an image that
  /// is transparent wherever the data has an alpha value of 0.
  ///
  template <class PixelT>
  class AlphaToMask : public ReturnFixedType<typename MaskedPixelType<typename PixelWithoutAlpha<PixelT>::type>::type> {
  public:
    typedef typename MaskedPixelType<typename PixelWithoutAlpha<PixelT>::type>::type result_type;
    inline result_type operator()( PixelT const& pixel ) const {
      if (is_transparent(pixel)) {
        return result_type();
      }
      return result_type(non_alpha_channels(pixel));
    }
  };

  // *******************************************************************
  /// EdgeMaskView
  ///
  /// Create an image with zero-valued (i.e. default contructor)
  /// pixels around the edges masked out.
  template <class ViewT>
  class EdgeMaskView : public ImageViewBase<EdgeMaskView<ViewT> >
  {
    ViewT m_view;
    //BlockCacheView<typename ViewT::pixel_type> m_view;

    // These vectors contain the indices of the first good pixel from
    // the edge of the image on each side.
    Vector<int> m_left, m_right;
    Vector<int> m_top,  m_bottom;

    // Use the edge vectors to determine if a pixel is valid.  Note:
    // this check fails for non convex edge masks!
    inline bool valid(int32 i, int32 j) const {
      if (i > m_left[j] && i < m_right[j] && j > m_top[i] && j < m_bottom[i])
        return true;
      else
        return false;
    }

  public:
    typedef typename ViewT::pixel_type            orig_pixel_type;
    typedef typename boost::remove_cv<typename boost::remove_reference<orig_pixel_type>::type>::type unmasked_pixel_type;
    typedef PixelMask<unmasked_pixel_type>        pixel_type;
    typedef PixelMask<unmasked_pixel_type>        result_type;
    typedef ProceduralPixelAccessor<EdgeMaskView> pixel_accessor;

    // EdgeMaskView( ViewT const& view,
    //               const ProgressCallback &progress_callback = ProgressCallback::dummy_instance() ) :
    //   m_view(view, Vector2i(512,512) ) {

    EdgeMaskView( ViewT const& view,
                  unmasked_pixel_type const& mask_value,
                  int32 mask_buffer,
                  const ProgressCallback &progress_callback = ProgressCallback::dummy_instance() ) :
      m_view(view) {

      m_left.set_size(view.rows());
      m_right.set_size(view.rows());

      for (int i = 0; i < view.rows(); ++i) {
        m_left[i] = 0;
        m_right[i] = view.cols();
      }

      m_top.set_size(view.cols());
      m_bottom.set_size(view.cols());

      for (int j = 0; j < view.cols(); ++j) {
        m_top[j] = 0;
        m_bottom[j] = view.rows();
      }

      // Scan over the image
      for (int j = 0; j < m_view.impl().rows(); ++j) {
        progress_callback.report_progress(float(j)/m_view.impl().rows()*0.5);

        // Search from the left side
        int i = 0;
        while ( i < m_view.impl().cols() && m_view.impl()(i,j) == mask_value )
          i++;
        m_left[j] = i + mask_buffer;

        // Search from the right side
        i = m_view.impl().cols() - 1;
        while ( i >= 0 && m_view.impl()(i,j) == mask_value )
          --i;
        m_right[j] = i - mask_buffer;
      }

      for (int i = 0; i < m_view.impl().cols(); ++i) {
        progress_callback.report_progress(0.5 + float(i)/m_view.impl().cols()*0.5);

        // Search from the top side of the image for black pixels
        int j = 0;
        while ( j < m_view.impl().rows() && m_view.impl()(i,j) == mask_value )
          ++j;
        m_top[i] = j + mask_buffer;

        // Search from the right side of the image for black pixels
        j = m_view.impl().rows() - 1;
        while ( j >= 0 && m_view.impl()(i,j) == mask_value )
          --j;
        m_bottom[i] = j - mask_buffer;
      }

      progress_callback.report_finished();
    }

    inline int32 cols  () const { return m_view.cols();   }
    inline int32 rows  () const { return m_view.rows();   }
    inline int32 planes() const { return m_view.planes(); }

    inline pixel_accessor origin() const { return pixel_accessor(*this); }

    inline result_type operator()( int32 i, int32 j, int32 p=0 ) const {
      if ( this->valid(i,j) )
        return pixel_type(m_view(i,j,p));
      else
        return pixel_type();
    }

    /// \cond INTERNAL
    typedef EdgeMaskView<ViewT> prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i const& /*bbox*/ ) const { return *this; }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const {
      vw::rasterize( prerasterize(bbox), dest, bbox );
    }
    /// \endcond
  };

  //******************************************************************
  /// invert_mask(view)
  ///
  /// Given a view with pixels of type PixelMask<T>, this will toggle
  /// all valids to invalid, and invalid to valids.
  template <class PixelT>
  class InvertPixelMask : public ReturnFixedType<PixelT> {
  public:
    inline PixelT operator()( PixelT value ) const {
      toggle(value);
      return value;
    }
  };

  //*****************************************************************
  /// validate_mask(view)
  ///
  /// Given an image of PixelMasks, this will make all pixels valid
  template <class PixelT>
  class ValidatePixelMask : public ReturnFixedType<PixelT> {
  public:
    inline PixelT operator()( PixelT value ) const {
      validate(value);
      return value;
    }
  };

  //*****************************************************************
  /// invalidate_mask(view)
  ///
  /// Given an image of PixelMasks, this will make all pixels invalid
  template <class PixelT>
  class InvalidatePixelMask : public ReturnFixedType<PixelT> {
  public:
    inline PixelT operator()( PixelT value ) const {
      invalidate(value);
      return value;
    }
  };

  //******************************************************************
  /// union_mask(view, mask)
  ///
  /// Unions 'mask' w/ view. View's data is returned
  ///
  template <class PixelT>
  class UnionPixelMask : public ReturnFixedType<typename MaskedPixelType<PixelT>::type> {
    typedef typename MaskedPixelType<PixelT>::type return_type;
  public:
    template <class MaskedPixelT>
    inline return_type operator()( PixelT const& value, MaskedPixelT const& mask ) const {
      return_type result = value;
      if ( is_valid(value) || is_valid(mask) )
        validate(result);
      else
        invalidate(result);
      return result;
    }
  };

  //******************************************************************
  /// intersect_mask(view, mask)
  ///
  /// Intersects 'mask' w/ view. View's data is returned
  ///
  template <class PixelT>
  class IntersectPixelMask : public ReturnFixedType<typename MaskedPixelType<PixelT>::type> {
    typedef typename MaskedPixelType<PixelT>::type return_type;
  public:
    template <class MaskedPixelT>
    inline return_type operator()( PixelT const& value, MaskedPixelT const& mask ) const {
      return_type result = value;
      if ( is_valid(value) && is_valid(mask) )
        validate(result);
      else
        invalidate(result);
      return result;
    }
  };

} // namespace vw

