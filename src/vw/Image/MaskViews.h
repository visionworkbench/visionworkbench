// __BEGIN_LICENSE__
//  Copyright (c) 2006-2025, United States Government as represented by the
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


/// \file PixelMask.h
///
/// Defines the useful pixel utility type that can wrap any existing
/// pixel type and add mask semantics.  Any operations with an
/// "invalid" pixel returns an invalid pixel as a result.
///
#ifndef __VW_IMAGE_MASK_VIEWS_H__
#define __VW_IMAGE_MASK_VIEWS_H__

#include <vw/Image/PixelMask.h>
#include <vw/Image/ImageViewBase.h>
#include <vw/Image/PerPixelViews.h>
#include <vw/Image/PixelAccessors.h>

#include <boost/math/special_functions/next.hpp>

namespace vw {

  // *******************************************************************
  /// create_mask(view, value)
  ///
  /// Given a view with pixels of type PixelT and a pixel value to
  /// consider as the "no data" or masked value, returns a view with
  /// pixels that are of the PixelMask<PixelT>, with the appropriate
  /// pixels masked.
  /// - Should safely accept inputs which already contain a mask, in
  ///   which case the input mask is ignored.

  /// Values are valid if they are different than nodata_val
  template <class PixelT>
  class CreatePixelMask: public ReturnFixedType<typename MaskedPixelType<PixelT>::type > {
    PixelT m_nodata_value;
  public:
    CreatePixelMask(PixelT const& nodata_value): m_nodata_value(nodata_value) {}

    inline typename MaskedPixelType<PixelT>::type operator()(PixelT const& value) const {
      typedef typename MaskedPixelType<PixelT>::type MPixelT;
      if (value != m_nodata_value && value == value)
        return MPixelT(value);

      if (value != value)
        return  MPixelT();  // Mask NaN values

      if (m_nodata_value != m_nodata_value)
        return MPixelT(value); // If value is non-NaN, but m_nodata_value is NaN, return good

      // We arrive here only if both value and m_nodata_value are not NaN,
      // and value == m_nodata_value.
      return MPixelT();
    }

  };

  /// Values are valid if nodata_val < val
  /// Mask values less than or equal to the nodata value.
  template <class PixelT>
  class CreatePixelMaskLE: public ReturnFixedType<typename MaskedPixelType<PixelT>::type > {
    PixelT m_nodata_value;
  public:
    CreatePixelMaskLE(PixelT const& nodata_value): m_nodata_value(nodata_value) {}

    inline typename MaskedPixelType<PixelT>::type operator()(PixelT const& value) const {
      typedef typename MaskedPixelType<PixelT>::type MPixelT;
      if (value > m_nodata_value)
        return MPixelT(value);

      if (value != value)
        return  MPixelT();  // Mask NaN values

      if (m_nodata_value != m_nodata_value)
        return MPixelT(value); // If value is non-NaN, but m_nodata_value is NaN, return good

      // We arrive here only if both value and m_nodata_value are not NaN,
      // and value <= m_nodata_value.
      return MPixelT();
    }

  };

  /// Values are valid if they min_val <= val <= max_val
  /// Mask values fall within a range.
  template <class PixelT>
  class CreatePixelRangeMask: public ReturnFixedType<typename MaskedPixelType<PixelT>::type > {
    PixelT m_valid_min;
    PixelT m_valid_max;
  public:
    CreatePixelRangeMask(PixelT const& valid_min, PixelT const& valid_max): m_valid_min(valid_min), m_valid_max(valid_max) {}

    // Helper to access only specific types of pixels
    template <bool CompoundB, class Arg1T, class Arg2T>
    struct Helper {
      static inline bool greater_than(Arg1T const& /*arg1*/, Arg2T const& /*arg2*/) {
        return true;
      }
      static inline bool less_than(Arg1T const& /*arg1*/, Arg2T const& /*arg2*/) {
        return true;
      }
    };

    // Specialization only for scalars
    template <class Arg1T, class Arg2T>
    struct Helper<false,Arg1T,Arg2T> {
      static inline bool greater_than(Arg1T const& arg1, Arg2T const& arg2) {
        return arg1 > arg2;
      }
      static inline bool less_than(Arg1T const& arg1, Arg2T const& arg2) {
        return arg1 < arg2;
      }
    };

    // Specialization for compounds
    template <class Arg1T, class Arg2T>
    struct Helper<true,Arg1T,Arg2T> {
      static inline bool greater_than(Arg1T const& arg1, Arg2T const& arg2) {
        return arg1[0] > arg2[0];
      }
      static inline bool less_than(Arg1T const& arg1, Arg2T const& arg2) {
        return arg1[0] < arg2[0];
      }
    };

    inline typename MaskedPixelType<PixelT>::type operator()(PixelT const& value) const {
      // Create Pixel Mask doesn't support thresholding of pixels with multiple channels
      BOOST_STATIC_ASSERT(CompoundNumChannels<PixelT>::value == 1);
      typedef typename MaskedPixelType<PixelT>::type MPixelT;

      typedef Helper<IsCompound<PixelT>::value,PixelT,PixelT> help_func;
      if (help_func::greater_than(value,m_valid_max) ||
          help_func::less_than(value,m_valid_min)    ||
          value != value) {
        return MPixelT();
      }

      return MPixelT(value);
    }
  }; // End class CreatePixelRangeMask

  /// Values are valid if nodata_val < val <= max_val
  /// Mask values less than or equal to the nodata value and greater than the given value.
  /// Use only with scalar types.
  template <class PixelT>
  class CreatePixelRangeMask2: public ReturnFixedType<typename MaskedPixelType<PixelT>::type > {
    PixelT m_nodata_value;
    PixelT m_max_valid_value;
  public:
    CreatePixelRangeMask2(PixelT const& nodata_value, PixelT const& max_valid_value):
      m_nodata_value(nodata_value), m_max_valid_value(max_valid_value) {}

    inline typename MaskedPixelType<PixelT>::type operator()(PixelT const& value) const {
      typedef typename MaskedPixelType<PixelT>::type MPixelT;

      if (value > m_nodata_value && value <= m_max_valid_value)
        return MPixelT(value);

      if (value != value)
        return  MPixelT();  // Mask NaN values

      // This code was not tested if nodata_value or max_valid_value is NaN.
      return MPixelT();
    }

  };

  /// Masks out pixels which are equal to NaN
  /// - Only use this with floats and doubles!
  /// Masks out pixels which are equal to NaN
  /// - Only use this with floats and doubles!
  template <class PixelT>
  class CreatePixelMaskNan: public ReturnFixedType<typename MaskedPixelType<PixelT>::type > {
  public:
    CreatePixelMaskNan() {}
    inline typename MaskedPixelType<PixelT>::type operator()(PixelT const& value) const {
      typedef typename MaskedPixelType<PixelT>::type MPixelT;
      if (boost::math::isnan(value)) {
        return MPixelT();
      }
      return MPixelT(value);
    }
  };

  /// Simple single value nodata
  template <class ViewT>
  UnaryPerPixelView<ViewT,CreatePixelMask<typename ViewT::pixel_type> >
  create_mask(ImageViewBase<ViewT> const& view, typename ViewT::pixel_type const& value) {
    typedef UnaryPerPixelView<ViewT,CreatePixelMask<typename ViewT::pixel_type> > view_type;
    return view_type(view.impl(), CreatePixelMask<typename ViewT::pixel_type>(value));
  }

  /// Valid if data falls within a range
  template <class ViewT>
  UnaryPerPixelView<ViewT,CreatePixelRangeMask<typename ViewT::pixel_type> >
  create_mask(ImageViewBase<ViewT> const& view,
               typename ViewT::pixel_type const& valid_min,
               typename ViewT::pixel_type const& valid_max) {
    typedef UnaryPerPixelView<ViewT,CreatePixelRangeMask<typename ViewT::pixel_type> > view_type;
    return view_type(view.impl(), CreatePixelRangeMask<typename ViewT::pixel_type>(valid_min, valid_max));
  }

  /// Default mask zero
  template <class ViewT>
  UnaryPerPixelView<ViewT,CreatePixelMask<typename ViewT::pixel_type> >
  create_mask(ImageViewBase<ViewT> const& view) {
    return create_mask(view.impl(), typename ViewT::pixel_type());
  }

  /// Mask values unless nodata_val < val <= max_val
  template <class ViewT>
  UnaryPerPixelView<ViewT,CreatePixelRangeMask2<typename ViewT::pixel_type> >
  create_pixel_range_mask2(ImageViewBase<ViewT> const& view,
                           typename ViewT::pixel_type const& nodata_val,
                           typename ViewT::pixel_type const& max_val) {
    typedef UnaryPerPixelView<ViewT,CreatePixelRangeMask2<typename ViewT::pixel_type> > view_type;
    return view_type(view.impl(),
                     CreatePixelRangeMask2<typename ViewT::pixel_type>(nodata_val, max_val));
  }

  /// Mask values less than or equal to the nodata value.
  template <class ViewT>
  UnaryPerPixelView<ViewT,CreatePixelMaskLE<typename ViewT::pixel_type> >
  create_mask_less_or_equal(ImageViewBase<ViewT> const& view, typename ViewT::pixel_type const& value) {
    typedef UnaryPerPixelView<ViewT,CreatePixelMaskLE<typename ViewT::pixel_type> > view_type;
    return view_type(view.impl(), CreatePixelMaskLE<typename ViewT::pixel_type>(value));
  }

  /// Mask out values which are NaN
  template <class ViewT>
  UnaryPerPixelView<ViewT,CreatePixelMaskNan<typename ViewT::pixel_type> >
  create_mask_nan(ImageViewBase<ViewT> const& view) {
    typedef UnaryPerPixelView<ViewT,CreatePixelMaskNan<typename ViewT::pixel_type> > view_type;
    return view_type(view.impl(), CreatePixelMaskNan<typename ViewT::pixel_type>());
  }

  // Indicate that create_mask is "reasonably fast" and should never
  // induce an extra rasterization step during prerasterization.
  template <class ViewT>
  struct IsMultiplyAccessible<UnaryPerPixelView<ViewT,CreatePixelMask<typename ViewT::pixel_type> > >: public IsMultiplyAccessible<ViewT> {};
  template <class ViewT>
  struct IsMultiplyAccessible<UnaryPerPixelView<ViewT,CreatePixelMaskLE<typename ViewT::pixel_type> > >: public IsMultiplyAccessible<ViewT> {};
  template <class ViewT>
  struct IsMultiplyAccessible<UnaryPerPixelView<ViewT,CreatePixelRangeMask<typename ViewT::pixel_type> > >: public IsMultiplyAccessible<ViewT> {};

  // *******************************************************************
  /// apply_mask(view, value)
  ///
  /// Given a view with pixels of the type PixelMask<T>, this view
  /// returns an image with pixels of type T where any pixel that was
  /// marked as "invalid" in the mask is replaced with the constant
  /// pixel value passed in as value.  The value is T() by default.
  ///
  template <class PixelT>
  class ApplyPixelMask: public ReturnFixedType<PixelT> {
    PixelT m_nodata_value;
  public:
    ApplyPixelMask(PixelT const& nodata_value): m_nodata_value(nodata_value) {}
    inline PixelT operator()(PixelMask<PixelT> const& value) const {
      return value.valid() ? value.child(): m_nodata_value;
    }
  };

  template <class ViewT>
  UnaryPerPixelView<ViewT,ApplyPixelMask<typename UnmaskedPixelType<typename ViewT::pixel_type>::type> >
  apply_mask(ImageViewBase<ViewT> const& view,
              typename UnmaskedPixelType<typename ViewT::pixel_type>::type const& value) {
    typedef UnaryPerPixelView<ViewT,ApplyPixelMask<typename UnmaskedPixelType<typename ViewT::pixel_type>::type> > view_type;
    return view_type(view.impl(), ApplyPixelMask<typename UnmaskedPixelType<typename ViewT::pixel_type>::type>(value));
  }

  // We overload the function rather than defaulting the value
  // argument to work around a compiler issue in MSVC 2005.
  template <class ViewT>
  UnaryPerPixelView<ViewT,ApplyPixelMask<typename UnmaskedPixelType<typename ViewT::pixel_type>::type> >
  apply_mask(ImageViewBase<ViewT> const& view) {
    return apply_mask(view.impl(), typename UnmaskedPixelType<typename ViewT::pixel_type>::type());
  }

  // Indicate that apply_mask is "reasonably fast" and should never
  // induce an extra rasterization step during prerasterization.
  template <class ViewT>
  struct IsMultiplyAccessible<UnaryPerPixelView<ViewT,ApplyPixelMask<typename UnmaskedPixelType<typename ViewT::pixel_type>::type> > >: public IsMultiplyAccessible<ViewT> {};

  // *******************************************************************
  /// copy_mask(view, mask)
  ///
  /// Copies a mask from one image to another.
  ///
  template <class PixelT>
  class CopyPixelMask: public ReturnFixedType<typename MaskedPixelType<PixelT>::type> {
  public:
    template <class MaskPixelT>
    inline typename MaskedPixelType<PixelT>::type
    operator()(PixelT const& value, MaskPixelT const& mask) const {
      typename MaskedPixelType<PixelT>::type result = value;
      if (is_transparent(mask)) {
        result.invalidate();
      }
      return result;
    }
  };

  /// Return a copy of the first argument with a mask copied from the second argument.
  template <class ViewT, class MaskViewT>
  BinaryPerPixelView<ViewT,MaskViewT,CopyPixelMask<typename ViewT::pixel_type> >
  copy_mask(ImageViewBase<    ViewT> const&      view,
             ImageViewBase<MaskViewT> const& mask_view) {
    typedef BinaryPerPixelView<ViewT,MaskViewT,CopyPixelMask<typename ViewT::pixel_type> > view_type;
    return view_type(view.impl(), mask_view.impl(), CopyPixelMask<typename ViewT::pixel_type>());
  }

  // Indicate that copy_mask is "reasonably fast" and should never
  // induce an extra rasterization step during prerasterization.
  template <class ViewT, class MaskViewT>
  struct IsMultiplyAccessible<BinaryPerPixelView<ViewT,MaskViewT,CopyPixelMask<typename ViewT::pixel_type> > >: public boost::mpl::and_<IsMultiplyAccessible<ViewT>,IsMultiplyAccessible<MaskViewT> >::type {};

  // *******************************************************************
  /// mask_to_alpha(view)
  ///
  /// Converts a mask channel to an alpha channel, generating an image that
  /// is transparent wherever the data is masked.
  ///
  template <class PixelT>
  class MaskToAlpha: public ReturnFixedType<typename PixelWithAlpha<typename UnmaskedPixelType<PixelT>::type>::type> {
  public:
    typedef typename PixelWithAlpha<typename UnmaskedPixelType<PixelT>::type>::type result_type;
    inline result_type operator()(PixelT const& pixel) const {
      if (is_transparent(pixel)) {
        return result_type();
      } else return result_type(pixel.child());
    }
  };

  template <class ViewT>
  UnaryPerPixelView<ViewT,MaskToAlpha<typename ViewT::pixel_type> >
  mask_to_alpha(ImageViewBase<ViewT> const& view) {
    typedef UnaryPerPixelView<ViewT,MaskToAlpha<typename ViewT::pixel_type> > view_type;
    return view_type(view.impl(), MaskToAlpha<typename ViewT::pixel_type>());
  }

  // Indicate that mask_to_alpha is "reasonably fast" and should never
  // induce an extra rasterization step during prerasterization.
  template <class ViewT>
  struct IsMultiplyAccessible<UnaryPerPixelView<ViewT,MaskToAlpha<typename ViewT::pixel_type> > >: public IsMultiplyAccessible<ViewT> {};

  // *******************************************************************
  /// alpha_to_mask(view)
  ///
  /// Converts an channel to a mask channel, generating an image that
  /// is transparent wherever the data has an alpha value of 0.
  ///
  template <class PixelT>
  class AlphaToMask: public ReturnFixedType<typename MaskedPixelType<typename PixelWithoutAlpha<PixelT>::type>::type> {
  public:
    typedef typename MaskedPixelType<typename PixelWithoutAlpha<PixelT>::type>::type result_type;
    inline result_type operator()(PixelT const& pixel) const {
      if (is_transparent(pixel)) {
        return result_type();
      }
      return result_type(non_alpha_channels(pixel));
    }
  };

  template <class ViewT>
  UnaryPerPixelView<ViewT,AlphaToMask<typename ViewT::pixel_type> >
  alpha_to_mask(ImageViewBase<ViewT> const& view) {
    typedef UnaryPerPixelView<ViewT,AlphaToMask<typename ViewT::pixel_type> > view_type;
    return view_type(view.impl(), AlphaToMask<typename ViewT::pixel_type>());
  }

  // Indicate that alpha_to_mask is "reasonably fast" and should never
  // induce an extra rasterization step during prerasterization.
  template <class ViewT>
  struct IsMultiplyAccessible<UnaryPerPixelView<ViewT,AlphaToMask<typename ViewT::pixel_type> > >: public IsMultiplyAccessible<ViewT> {};

  // *******************************************************************
  /// EdgeMaskView
  ///
  /// Create an image with zero-valued (i.e. default constructor)
  /// pixels around the edges masked out.
  template <class ViewT>
  class EdgeMaskView: public ImageViewBase<EdgeMaskView<ViewT> >
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

    // EdgeMaskView(ViewT const& view,
    //               const ProgressCallback &progress_callback = ProgressCallback::dummy_instance()):
    //   m_view(view, Vector2i(512,512)) {

    EdgeMaskView(ViewT const& view,
                  unmasked_pixel_type const& mask_value,
                  int32 mask_buffer,
                  const ProgressCallback &progress_callback = ProgressCallback::dummy_instance()):
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
        while (i < m_view.impl().cols() && m_view.impl()(i,j) == mask_value)
          i++;
        m_left[j] = i + mask_buffer;

        // Search from the right side
        i = m_view.impl().cols() - 1;
        while (i >= 0 && m_view.impl()(i,j) == mask_value)
          --i;
        m_right[j] = i - mask_buffer;
      }

      for (int i = 0; i < m_view.impl().cols(); ++i) {
        progress_callback.report_progress(0.5 + float(i)/m_view.impl().cols()*0.5);

        // Search from the top side of the image for black pixels
        int j = 0;
        while (j < m_view.impl().rows() && m_view.impl()(i,j) == mask_value)
          ++j;
        m_top[i] = j + mask_buffer;

        // Search from the right side of the image for black pixels
        j = m_view.impl().rows() - 1;
        while (j >= 0 && m_view.impl()(i,j) == mask_value)
          --j;
        m_bottom[i] = j - mask_buffer;
      }

      progress_callback.report_finished();
    }

    inline int32 cols  () const { return m_view.cols();   }
    inline int32 rows  () const { return m_view.rows();   }
    inline int32 planes() const { return m_view.planes(); }

    inline pixel_accessor origin() const { return pixel_accessor(*this); }

    inline result_type operator()(int32 i, int32 j, int32 p=0) const {
      if (this->valid(i,j))
        return pixel_type(m_view(i,j,p));
      else
        return pixel_type();
    }

    /// \cond INTERNAL
    typedef EdgeMaskView<ViewT> prerasterize_type;
    inline prerasterize_type prerasterize(BBox2i const& /*bbox*/) const { return *this; }
    template <class DestT> inline void rasterize(DestT const& dest, BBox2i const& bbox) const {
      vw::rasterize(prerasterize(bbox), dest, bbox);
    }
    /// \endcond
  };

  /// \cond INTERNAL
  template <class ViewT>
  struct IsMultiplyAccessible<EdgeMaskView<ViewT> >: public true_type {};
  /// \endcond

  /// edge_mask(view)
  ///
  /// Search from the edges of an image for the first "valid" pixels,
  /// masking invalid pixels along the way.  Unlike some other image
  /// views, the EdgeMaskView does a good portion of its work in its
  /// constructor, where it searches for valid/invalid pixels in the
  /// source view.  The results are efficiently stored in four "edge
  /// location" vectors (one for each side).
  ///
  /// XXX The following note currently appears to be false:
  /// Performance note: this algorithm stores the input view in an
  /// additional BlockCacheView<> since it scans over every pixel in
  /// the image both horizontally and vertically.  Be sure that your
  /// cache is large enough to store a full row or column of blocks!!
  template <class ViewT>
  EdgeMaskView<ViewT> edge_mask(ImageViewBase<ViewT> const& v,
                                 const ProgressCallback &progress_callback = ProgressCallback::dummy_instance()) {
    return EdgeMaskView<ViewT>(v.impl(), typename ViewT::pixel_type(), 0, progress_callback);
  }

  template <class ViewT>
  EdgeMaskView<ViewT> edge_mask(ImageViewBase<ViewT> const& v,
                                 typename ViewT::pixel_type value,
                                 int32 buffer = 0,
                                 const ProgressCallback &progress_callback = ProgressCallback::dummy_instance()) {
    return EdgeMaskView<ViewT>(v.impl(), value, buffer, progress_callback);
  }

  //******************************************************************
  /// invert_mask(view)
  ///
  /// Given a view with pixels of type PixelMask<T>, this will toggle
  /// all valid to invalid, and invalid to valid.
  template <class PixelT>
  class InvertPixelMask: public ReturnFixedType<PixelT> {
  public:
    inline PixelT operator()(PixelT value) const {
      toggle(value);
      return value;
    }
  };

  template <class ViewT>
  UnaryPerPixelView<ViewT,InvertPixelMask<typename ViewT::pixel_type> >
  invert_mask(ImageViewBase<ViewT> const& view) {
    typedef UnaryPerPixelView<ViewT,InvertPixelMask<typename ViewT::pixel_type> > view_type;
    return view_type(view.impl(), InvertPixelMask<typename ViewT::pixel_type>());
  }

  // Invert Pixel Mask is "reasonably fast"
  template <class ViewT>
  struct IsMultiplyAccessible<UnaryPerPixelView<ViewT,InvertPixelMask<typename ViewT::pixel_type> > >: public IsMultiplyAccessible<ViewT> {};

  //*****************************************************************
  /// validate_mask(view)
  ///
  /// Given an image of PixelMasks, this will make all pixels valid
  template <class PixelT>
  class ValidatePixelMask;

  template <class ViewT>
  UnaryPerPixelView<ViewT,ValidatePixelMask<typename ViewT::pixel_type> >
  validate_mask(ImageViewBase<ViewT> const& view) {
    typedef UnaryPerPixelView<ViewT,ValidatePixelMask<typename ViewT::pixel_type> > view_type;
    return view_type(view.impl(), ValidatePixelMask<typename ViewT::pixel_type>());
  }

  // Validate Pixel Mask is "reasonably fast"
  template <class ViewT>
  struct IsMultiplyAccessible<UnaryPerPixelView<ViewT,ValidatePixelMask<typename ViewT::pixel_type> > >: public IsMultiplyAccessible<ViewT> {};

  //*****************************************************************
  /// invalidate_mask(view)
  ///
  /// Given an image of PixelMasks, this will make all pixels invalid
  template <class PixelT>
  class InvalidatePixelMask: public ReturnFixedType<PixelT> {
  public:
    inline PixelT operator()(PixelT value) const {
      invalidate(value);
      return value;
    }
  };

  template <class ViewT>
  UnaryPerPixelView<ViewT,InvalidatePixelMask<typename ViewT::pixel_type> >
  invalidate_mask(ImageViewBase<ViewT> const& view) {
    typedef UnaryPerPixelView<ViewT,InvalidatePixelMask<typename ViewT::pixel_type> > view_type;
    return view_type(view.impl(), InvalidatePixelMask<typename ViewT::pixel_type>());
  }

  // Invalidate Pixel Mask is "reasonably fast"
  template <class ViewT>
  struct IsMultiplyAccessible<UnaryPerPixelView<ViewT,InvalidatePixelMask<typename ViewT::pixel_type> > >: public IsMultiplyAccessible<ViewT> {};

  //******************************************************************
  /// union_mask(view, mask)
  ///
  /// Unions 'mask' w/ view. View's data is returned
  ///
  template <class PixelT>
  class UnionPixelMask: public ReturnFixedType<typename MaskedPixelType<PixelT>::type> {
    typedef typename MaskedPixelType<PixelT>::type return_type;
  public:
    template <class MaskedPixelT>
    inline return_type operator()(PixelT const& value, MaskedPixelT const& mask) const {
      return_type result = value;
      if (is_valid(value) || is_valid(mask))
        validate(result);
      else
        invalidate(result);
      return result;
    }
  };

  template <class ViewT, class MaskViewT>
  BinaryPerPixelView<ViewT,MaskViewT,UnionPixelMask<typename ViewT::pixel_type> >
  union_mask(ImageViewBase<    ViewT> const&      view,
              ImageViewBase<MaskViewT> const& mask_view) {
    typedef BinaryPerPixelView<ViewT,MaskViewT,UnionPixelMask<typename ViewT::pixel_type> > view_type;
    return view_type(view.impl(), mask_view.impl(), UnionPixelMask<typename ViewT::pixel_type>());
  }

  // Is reasonably fast
  template <class ViewT, class MaskViewT>
  struct IsMultiplyAccessible<BinaryPerPixelView<ViewT,MaskViewT,UnionPixelMask<typename ViewT::pixel_type> > >: public boost::mpl::and_<IsMultiplyAccessible<ViewT>,IsMultiplyAccessible<MaskViewT> >::type {};

  //******************************************************************
  /// intersect_mask(view, mask)
  ///
  /// Intersects 'mask' w/ view. View's data is returned
  ///
  template <class PixelT>
  class IntersectPixelMask: public ReturnFixedType<typename MaskedPixelType<PixelT>::type> {
    typedef typename MaskedPixelType<PixelT>::type return_type;
  public:
    template <class MaskedPixelT>
    inline return_type operator()(PixelT const& value, MaskedPixelT const& mask) const {
      return_type result = value;
      if (is_valid(value) && is_valid(mask))
        validate(result);
      else
        invalidate(result);
      return result;
    }
  };

  template <class ViewT, class MaskViewT>
  BinaryPerPixelView<ViewT,MaskViewT,IntersectPixelMask<typename ViewT::pixel_type> >
  intersect_mask(ImageViewBase<    ViewT> const&      view,
                  ImageViewBase<MaskViewT> const& mask_view) {
    typedef BinaryPerPixelView<ViewT,MaskViewT,IntersectPixelMask<typename ViewT::pixel_type> > view_type;
    return view_type(view.impl(), mask_view.impl(),
                     IntersectPixelMask<typename ViewT::pixel_type>());
  }

  // Is reasonable fast
  template <class ViewT, class MaskViewT>
  struct IsMultiplyAccessible<BinaryPerPixelView<ViewT,MaskViewT,IntersectPixelMask<typename ViewT::pixel_type> > >: public boost::mpl::and_<IsMultiplyAccessible<ViewT>,IsMultiplyAccessible<MaskViewT> >::type {};

} // namespace vw

#endif // __VW_IMAGE_MASK_VIEWS_H__
