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
  /// create_mask( view, value )
  ///
  /// Given a view with pixels of type PixelT and a pixel value to
  /// consider as the "no data" or masked value, returns a view with
  /// pixels that are of the PixelMask<PixelT>, with the appropriate
  /// pixels masked.
  /// - Should safely accept inputs which already contain a mask, in 
  ///   which case the input mask is ignored.

  /// Values are valid if they are different than nodata_val
  template <class PixelT>
  class CreatePixelMask : public ReturnFixedType<typename MaskedPixelType<PixelT>::type > {
    PixelT m_nodata_value;
  public:
    CreatePixelMask( PixelT const& nodata_value ) : m_nodata_value(nodata_value) {}
    
    inline typename MaskedPixelType<PixelT>::type operator()( PixelT const& value ) const {
      typedef typename MaskedPixelType<PixelT>::type MPixelT;
      if ( value != m_nodata_value && value == value ) // need the latter for NaNs 
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
  template <class PixelT>
  class CreatePixelMaskLE;

  /// Values are valid if they min_val <= val <= max_val
  template <class PixelT>
  class CreatePixelRangeMask;

  /// Values are valid if nodata_val < val <= max_val
  template <class PixelT>
  class CreatePixelRangeMask2;

  /// Masks out pixels which are equal to NaN
  /// - Only use this with floats and doubles!
  template <class PixelT>
  class CreatePixelMaskNan;


  /// Simple single value nodata
  template <class ViewT>
  UnaryPerPixelView<ViewT,CreatePixelMask<typename ViewT::pixel_type> >
  create_mask( ImageViewBase<ViewT> const& view, typename ViewT::pixel_type const& value ) {
    typedef UnaryPerPixelView<ViewT,CreatePixelMask<typename ViewT::pixel_type> > view_type;
    return view_type( view.impl(), CreatePixelMask<typename ViewT::pixel_type>(value) );
  }

  /// Valid if data falls within a range
  template <class ViewT>
  UnaryPerPixelView<ViewT,CreatePixelRangeMask<typename ViewT::pixel_type> >
  create_mask( ImageViewBase<ViewT> const& view,
               typename ViewT::pixel_type const& valid_min,
               typename ViewT::pixel_type const& valid_max ) {
    typedef UnaryPerPixelView<ViewT,CreatePixelRangeMask<typename ViewT::pixel_type> > view_type;
    return view_type( view.impl(), CreatePixelRangeMask<typename ViewT::pixel_type>( valid_min, valid_max ));
  }

  /// Default mask zero
  template <class ViewT>
  UnaryPerPixelView<ViewT,CreatePixelMask<typename ViewT::pixel_type> >
  create_mask( ImageViewBase<ViewT> const& view ) {
    return create_mask( view.impl(), typename ViewT::pixel_type() );
  }

  /// Mask values unless nodata_val < val <= max_val
  template <class ViewT>
  UnaryPerPixelView<ViewT,CreatePixelRangeMask2<typename ViewT::pixel_type> >
  create_pixel_range_mask2(ImageViewBase<ViewT> const& view,
			   typename ViewT::pixel_type const& nodata_val,
			   typename ViewT::pixel_type const& max_val
			   ) {
    typedef UnaryPerPixelView<ViewT,CreatePixelRangeMask2<typename ViewT::pixel_type> > view_type;
    return view_type( view.impl(),
		      CreatePixelRangeMask2<typename ViewT::pixel_type>(nodata_val, max_val) );
  }

  /// Mask values less than or equal to the nodata value.
  template <class ViewT>
  UnaryPerPixelView<ViewT,CreatePixelMaskLE<typename ViewT::pixel_type> >
  create_mask_less_or_equal( ImageViewBase<ViewT> const& view, typename ViewT::pixel_type const& value ) {
    typedef UnaryPerPixelView<ViewT,CreatePixelMaskLE<typename ViewT::pixel_type> > view_type;
    return view_type( view.impl(), CreatePixelMaskLE<typename ViewT::pixel_type>(value) );
  }
  
  /// Mask out values which are NaN
  template <class ViewT>
  UnaryPerPixelView<ViewT,CreatePixelMaskNan<typename ViewT::pixel_type> >
  create_mask_nan( ImageViewBase<ViewT> const& view ) {
    typedef UnaryPerPixelView<ViewT,CreatePixelMaskNan<typename ViewT::pixel_type> > view_type;
    return view_type( view.impl(), CreatePixelMaskNan<typename ViewT::pixel_type>() );
  }  
  
  

  // Indicate that create_mask is "reasonably fast" and should never
  // induce an extra rasterization step during prerasterization.
  template <class ViewT>
  struct IsMultiplyAccessible<UnaryPerPixelView<ViewT,CreatePixelMask<typename ViewT::pixel_type> > >
    : public IsMultiplyAccessible<ViewT> {};
  template <class ViewT>
  struct IsMultiplyAccessible<UnaryPerPixelView<ViewT,CreatePixelMaskLE<typename ViewT::pixel_type> > >
    : public IsMultiplyAccessible<ViewT> {};
  template <class ViewT>
  struct IsMultiplyAccessible<UnaryPerPixelView<ViewT,CreatePixelRangeMask<typename ViewT::pixel_type> > >
    : public IsMultiplyAccessible<ViewT> {};

  // *******************************************************************
  /// apply_mask( view, value )
  ///
  /// Given a view with pixels of the type PixelMask<T>, this view
  /// returns an image with pixels of type T where any pixel that was
  /// marked as "invalid" in the mask is replaced with the constant
  /// pixel value passed in as value.  The value is T() by default.
  ///
  template <class PixelT>
  class ApplyPixelMask;

  template <class ViewT>
  UnaryPerPixelView<ViewT,ApplyPixelMask<typename UnmaskedPixelType<typename ViewT::pixel_type>::type> >
  apply_mask( ImageViewBase<ViewT> const& view,
              typename UnmaskedPixelType<typename ViewT::pixel_type>::type const& value ) {
    typedef UnaryPerPixelView<ViewT,ApplyPixelMask<typename UnmaskedPixelType<typename ViewT::pixel_type>::type> > view_type;
    return view_type( view.impl(), ApplyPixelMask<typename UnmaskedPixelType<typename ViewT::pixel_type>::type>(value) );
  }

  // We overload the function rather than defaulting the value
  // argument to work around a compiler issue in MSVC 2005.
  template <class ViewT>
  UnaryPerPixelView<ViewT,ApplyPixelMask<typename UnmaskedPixelType<typename ViewT::pixel_type>::type> >
  apply_mask( ImageViewBase<ViewT> const& view ) {
    return apply_mask( view.impl(), typename UnmaskedPixelType<typename ViewT::pixel_type>::type() );
  }

  // Indicate that apply_mask is "reasonably fast" and should never
  // induce an extra rasterization step during prerasterization.
  template <class ViewT>
  struct IsMultiplyAccessible<UnaryPerPixelView<ViewT,ApplyPixelMask<typename UnmaskedPixelType<typename ViewT::pixel_type>::type> > >
    : public IsMultiplyAccessible<ViewT> {};

  // *******************************************************************
  /// copy_mask(view, mask)
  ///
  /// Copies a mask from one image to another.
  ///
  template <class PixelT>
  class CopyPixelMask;

  /// Return a copy of the first argument with a mask copied from the second argument.
  template <class ViewT, class MaskViewT>
  BinaryPerPixelView<ViewT,MaskViewT,CopyPixelMask<typename ViewT::pixel_type> >
  copy_mask( ImageViewBase<    ViewT> const&      view,
             ImageViewBase<MaskViewT> const& mask_view ) {
    typedef BinaryPerPixelView<ViewT,MaskViewT,CopyPixelMask<typename ViewT::pixel_type> > view_type;
    return view_type( view.impl(), mask_view.impl(), CopyPixelMask<typename ViewT::pixel_type>() );
  }

  // Indicate that copy_mask is "reasonably fast" and should never
  // induce an extra rasterization step during prerasterization.
  template <class ViewT, class MaskViewT>
  struct IsMultiplyAccessible<BinaryPerPixelView<ViewT,MaskViewT,CopyPixelMask<typename ViewT::pixel_type> > >
    : public boost::mpl::and_<IsMultiplyAccessible<ViewT>,IsMultiplyAccessible<MaskViewT> >::type {};

  // *******************************************************************
  /// mask_to_alpha(view)
  ///
  /// Converts a mask channel to an alpha channel, generating an image that
  /// is transparent wherever the data is masked.
  ///
  template <class PixelT>
  class MaskToAlpha;

  template <class ViewT>
  UnaryPerPixelView<ViewT,MaskToAlpha<typename ViewT::pixel_type> >
  mask_to_alpha( ImageViewBase<ViewT> const& view ) {
    typedef UnaryPerPixelView<ViewT,MaskToAlpha<typename ViewT::pixel_type> > view_type;
    return view_type( view.impl(), MaskToAlpha<typename ViewT::pixel_type>() );
  }

  // Indicate that mask_to_alpha is "reasonably fast" and should never
  // induce an extra rasterization step during prerasterization.
  template <class ViewT>
  struct IsMultiplyAccessible<UnaryPerPixelView<ViewT,MaskToAlpha<typename ViewT::pixel_type> > >
    : public IsMultiplyAccessible<ViewT> {};


  // *******************************************************************
  /// alpha_to_mask(view)
  ///
  /// Converts an channel to a mask channel, generating an image that
  /// is transparent wherever the data has an alpha value of 0.
  ///
  template <class PixelT>
  class AlphaToMask;

  template <class ViewT>
  UnaryPerPixelView<ViewT,AlphaToMask<typename ViewT::pixel_type> >
  alpha_to_mask( ImageViewBase<ViewT> const& view ) {
    typedef UnaryPerPixelView<ViewT,AlphaToMask<typename ViewT::pixel_type> > view_type;
    return view_type( view.impl(), AlphaToMask<typename ViewT::pixel_type>() );
  }

  // Indicate that alpha_to_mask is "reasonably fast" and should never
  // induce an extra rasterization step during prerasterization.
  template <class ViewT>
  struct IsMultiplyAccessible<UnaryPerPixelView<ViewT,AlphaToMask<typename ViewT::pixel_type> > >
    : public IsMultiplyAccessible<ViewT> {};


  // *******************************************************************
  /// EdgeMaskView
  ///
  /// Create an image with zero-valued (i.e. default contructor)
  /// pixels around the edges masked out.
  template <class ViewT>
  class EdgeMaskView;

  /// \cond INTERNAL
  template <class ViewT>
  struct IsMultiplyAccessible<EdgeMaskView<ViewT> > : public true_type {};
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
  EdgeMaskView<ViewT> edge_mask( ImageViewBase<ViewT> const& v,
                                 const ProgressCallback &progress_callback = ProgressCallback::dummy_instance() ) {
    return EdgeMaskView<ViewT>( v.impl(), typename ViewT::pixel_type(), 0, progress_callback );
  }

  template <class ViewT>
  EdgeMaskView<ViewT> edge_mask( ImageViewBase<ViewT> const& v,
                                 typename ViewT::pixel_type value,
                                 int32 buffer = 0,
                                 const ProgressCallback &progress_callback = ProgressCallback::dummy_instance() ) {
    return EdgeMaskView<ViewT>( v.impl(), value, buffer, progress_callback );
  }

  //******************************************************************
  /// invert_mask(view)
  ///
  /// Given a view with pixels of type PixelMask<T>, this will toggle
  /// all valids to invalid, and invalid to valids.
  template <class PixelT>
  class InvertPixelMask;

  template <class ViewT>
  UnaryPerPixelView<ViewT,InvertPixelMask<typename ViewT::pixel_type> >
  invert_mask( ImageViewBase<ViewT> const& view ) {
    typedef UnaryPerPixelView<ViewT,InvertPixelMask<typename ViewT::pixel_type> > view_type;
    return view_type( view.impl(), InvertPixelMask<typename ViewT::pixel_type>());
  }

  // Invert Pixel Mask is "reasonably fast"
  template <class ViewT>
  struct IsMultiplyAccessible<UnaryPerPixelView<ViewT,InvertPixelMask<typename ViewT::pixel_type> > >
    : public IsMultiplyAccessible<ViewT> {};

  //*****************************************************************
  /// validate_mask(view)
  ///
  /// Given an image of PixelMasks, this will make all pixels valid
  template <class PixelT>
  class ValidatePixelMask;

  template <class ViewT>
  UnaryPerPixelView<ViewT,ValidatePixelMask<typename ViewT::pixel_type> >
  validate_mask( ImageViewBase<ViewT> const& view ) {
    typedef UnaryPerPixelView<ViewT,ValidatePixelMask<typename ViewT::pixel_type> > view_type;
    return view_type( view.impl(), ValidatePixelMask<typename ViewT::pixel_type>());
  }

  // Validate Pixel Mask is "reasonably fast"
  template <class ViewT>
  struct IsMultiplyAccessible<UnaryPerPixelView<ViewT,ValidatePixelMask<typename ViewT::pixel_type> > >
    : public IsMultiplyAccessible<ViewT> {};

  //*****************************************************************
  /// invalidate_mask(view)
  ///
  /// Given an image of PixelMasks, this will make all pixels invalid
  template <class PixelT>
  class InvalidatePixelMask;

  template <class ViewT>
  UnaryPerPixelView<ViewT,InvalidatePixelMask<typename ViewT::pixel_type> >
  invalidate_mask( ImageViewBase<ViewT> const& view ) {
    typedef UnaryPerPixelView<ViewT,InvalidatePixelMask<typename ViewT::pixel_type> > view_type;
    return view_type( view.impl(), InvalidatePixelMask<typename ViewT::pixel_type>());
  }

  // Invalidate Pixel Mask is "reasonably fast"
  template <class ViewT>
  struct IsMultiplyAccessible<UnaryPerPixelView<ViewT,InvalidatePixelMask<typename ViewT::pixel_type> > >
    : public IsMultiplyAccessible<ViewT> {};

  //******************************************************************
  /// union_mask(view, mask)
  ///
  /// Unions 'mask' w/ view. View's data is returned
  ///
  template <class PixelT>
  class UnionPixelMask;

  template <class ViewT, class MaskViewT>
  BinaryPerPixelView<ViewT,MaskViewT,UnionPixelMask<typename ViewT::pixel_type> >
  union_mask( ImageViewBase<    ViewT> const&      view,
              ImageViewBase<MaskViewT> const& mask_view ) {
    typedef BinaryPerPixelView<ViewT,MaskViewT,UnionPixelMask<typename ViewT::pixel_type> > view_type;
    return view_type( view.impl(), mask_view.impl(), UnionPixelMask<typename ViewT::pixel_type>() );
  }

  // Is reasonably fast
  template <class ViewT, class MaskViewT>
  struct IsMultiplyAccessible<BinaryPerPixelView<ViewT,MaskViewT,UnionPixelMask<typename ViewT::pixel_type> > >
    : public boost::mpl::and_<IsMultiplyAccessible<ViewT>,IsMultiplyAccessible<MaskViewT> >::type {};

  //******************************************************************
  /// intersect_mask(view, mask)
  ///
  /// Intersects 'mask' w/ view. View's data is returned
  ///
  template <class PixelT>
  class IntersectPixelMask;

  template <class ViewT, class MaskViewT>
  BinaryPerPixelView<ViewT,MaskViewT,IntersectPixelMask<typename ViewT::pixel_type> >
  intersect_mask( ImageViewBase<    ViewT> const&      view,
                  ImageViewBase<MaskViewT> const& mask_view ) {
    typedef BinaryPerPixelView<ViewT,MaskViewT,IntersectPixelMask<typename ViewT::pixel_type> > view_type;
    return view_type( view.impl(), mask_view.impl(), IntersectPixelMask<typename ViewT::pixel_type>() );
  }

  // Is reasonable fast
  template <class ViewT, class MaskViewT>
  struct IsMultiplyAccessible<BinaryPerPixelView<ViewT,MaskViewT,IntersectPixelMask<typename ViewT::pixel_type> > >
    : public boost::mpl::and_<IsMultiplyAccessible<ViewT>,IsMultiplyAccessible<MaskViewT> >::type {};

} // namespace vw

#include "MaskViews.tcc"

#endif // __VW_IMAGE_MASK_VIEWS_H__
