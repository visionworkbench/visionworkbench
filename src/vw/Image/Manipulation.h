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

/// \file Manipulation.h
///
/// Simple image manipulation functions, such as flipping and cropping.  All of
/// the functions in this file except copy() return shallow views of the
/// ImageView.  That is, they do not copy the data underneath but instead they
/// refer to the same data, indexing and accessing it in a different way.
///
/// The first collection of functions in this file perform basic
/// transformations to the domain of the image, such as transposition,
/// rotation by 90 degree increments, flipping, and cropping.
///
/// This file also provides views and functions that take simple
/// "slices" of images, returning a new view composed from individual
/// channels or planes of the source image. These include:
///
/// - select_col          () : takes a single-column slice of an image
/// - select_row          () : takes a single-row slice of an image
/// - select_plane        () : takes a single-plane slice of an image
/// - select_channel      () : takes a single-channel slice of an image
/// - channels_to_planes  () : reinterprets a multi-channel image as a multi-plane image
/// - planes_to_channels  () : reinterprets a multi-plane image as a multi-channel image
/// - pixel_cast          () : casts the pixels of an image to a new pixel type
/// - pixel_cast_rescale  () : same as above, but rescales pixel values appropriately
/// - channel_cast        () : casts the channels of an image while retaining the pixel format
/// - channel_cast_rescale() : same as above, but rescales pixel values appropriately
/// - channel_cast_round_and_clamp(): round, clamp, then cast. 
/// See PixelMath.h for round(), Algorithms.h for clamp(), etc.

#ifndef __VW_IMAGE_MANIPULATION_H__
#define __VW_IMAGE_MANIPULATION_H__

#include <boost/mpl/logical.hpp>

#include <vw/Image/ImageView.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Image/PerPixelViews.h>
#include <vw/Image/AlgorithmFunctions.h>

namespace vw {

  // *******************************************************************
  // copy()
  // *******************************************************************

  // Class definition
  template <class ImageT>
  class CopyView;

  template <class ImageT>
  struct IsMultiplyAccessible<CopyView<ImageT> > : public true_type {};

  /// Make a (deep) copy of an image.
  template <class ImageT>
  CopyView<ImageT> copy( ImageViewBase<ImageT> const& v ) {
    return CopyView<ImageT>( v.impl() );
  }


  // *******************************************************************
  // Transpose
  // *******************************************************************

  // Specialized pixel accessor
  template <class ChildT>
  class TransposePixelAccessor;

  // Class definition
  template <class ImageT>
  class TransposeView;

  /// \cond INTERNAL
  template <class ImageT>
  struct IsMultiplyAccessible<TransposeView<ImageT> > : public IsMultiplyAccessible<ImageT> {};
  /// \endcond

  /// Transpose an image.
  template <class ImageT>
  TransposeView<ImageT> transpose( ImageViewBase<ImageT> const& v ) {
    return TransposeView<ImageT>( v.impl() );
  }


  // *******************************************************************
  // Rotate180
  // *******************************************************************

  // Specialized pixel accessor
  template <class ChildT>
  class Rotate180PixelAccessor;

  // Image View Class
  template <class ImageT>
  class Rotate180View;

  /// \cond INTERNAL
  template <class ImageT>
  struct IsMultiplyAccessible<Rotate180View<ImageT> > : public IsMultiplyAccessible<ImageT> {};
  /// \endcond

  /// Rotate an image 180 degrees.
  template <class ImageT>
  Rotate180View<ImageT> rotate_180( ImageViewBase<ImageT> const& v ) {
    return Rotate180View<ImageT>( v.impl() );
  }


  // *******************************************************************
  // Rotate90CW
  // *******************************************************************

  // Specialized pixel accessor
  template <class ChildT>
  class Rotate90CWPixelAccessor;

  // Class definition
  template <class ImageT>
  class Rotate90CWView;

  /// \cond INTERNAL
  template <class ImageT>
  struct IsMultiplyAccessible<Rotate90CWView<ImageT> > : public IsMultiplyAccessible<ImageT> {};
  /// \endcond

  /// Rotate an image 90 degrees clockwise.
  template <class ImageT>
  Rotate90CWView<ImageT> rotate_90_cw( ImageViewBase<ImageT> const& v ) {
    return Rotate90CWView<ImageT>( v.impl() );
  }


  // *******************************************************************
  // Rotate90CCW
  // *******************************************************************

  // Specialized pixel accessor
  template <class ChildT>
  class Rotate90CCWPixelAccessor;

  // Class definition
  template <class ImageT>
  class Rotate90CCWView;

  /// \cond INTERNAL
  template <class ImageT>
  struct IsMultiplyAccessible<Rotate90CCWView<ImageT> > : public IsMultiplyAccessible<ImageT> {};
  /// \endcond

  /// Rotate an image 90 degrees counter-clockwise.
  template <class ImageT>
  Rotate90CCWView<ImageT> rotate_90_ccw( ImageViewBase<ImageT> const& v ) {
    return Rotate90CCWView<ImageT>( v.impl() );
  }


  // *******************************************************************
  // FlipVertical
  // *******************************************************************

  // Specialized pixel accessor
  template <class ChildT>
  class FlipVerticalPixelAccessor;

  // Class definition
  template <class ImageT>
  class FlipVerticalView;

  /// \cond INTERNAL
  template <class ImageT>
  struct IsMultiplyAccessible<FlipVerticalView<ImageT> > : public IsMultiplyAccessible<ImageT> {};
  /// \endcond

  /// Flip an image vertically.
  template <class ImageT>
  FlipVerticalView<ImageT> flip_vertical( ImageViewBase<ImageT> const& v ) {
    return FlipVerticalView<ImageT>( v.impl() );
  }


  // *******************************************************************
  // FlipHorizontal
  // *******************************************************************

  // Specialized pixel accessor
  template <class ChildT>
  class FlipHorizontalPixelAccessor;

  // Class definition
  template <class ImageT>
  class FlipHorizontalView;

  /// \cond INTERNAL
  template <class ImageT>
  struct IsMultiplyAccessible<FlipHorizontalView<ImageT> > : public IsMultiplyAccessible<ImageT> {};
  /// \endcond

  /// Flip an image horizontally.
  template <class ImageT>
  FlipHorizontalView<ImageT> flip_horizontal( ImageViewBase<ImageT> const& v ) {
    return FlipHorizontalView<ImageT>( v.impl() );
  }


  // *******************************************************************
  // crop()
  // *******************************************************************

  // Class definition
  template <class ImageT>
  class CropView;

  /// \cond INTERNAL
  // Type traits
  template <class ImageT>
  struct IsFloatingPointIndexable<CropView<ImageT> >  : public IsFloatingPointIndexable<ImageT> {};

  template <class ImageT>
  struct IsMultiplyAccessible<CropView<ImageT> > : public IsMultiplyAccessible<ImageT> {};
  /// \endcond

  /// Crop an image.
  template <class ImageT>
  inline CropView<ImageT> crop( ImageViewBase<ImageT> const& v, int32 upper_left_x, int32 upper_left_y, int32 width, int32 height ) {
    return CropView<ImageT>( v.impl(), upper_left_x, upper_left_y, width, height );
  }

  /// Crop an image.
  template <class ImageT, class BBoxRealT>
  inline CropView<ImageT> crop( ImageViewBase<ImageT> const& v, BBox<BBoxRealT,2> const& bbox ) {
    return CropView<ImageT>( v.impl(), bbox );
  }


  // *******************************************************************
  // subsample()
  // *******************************************************************

  // Specialized image accessor
  template <class ChildT>
  class SubsamplePixelAccessor;

  // Class definition
  template <class ImageT>
  class SubsampleView;

  /// \cond INTERNAL
  template <class ImageT>
  struct IsMultiplyAccessible<SubsampleView<ImageT> > : public IsMultiplyAccessible<ImageT> {};
  /// \endcond

  /// Subsample an image by an integer factor.  Note that this
  /// function does not pre-smooth the image prior to subsampling: it
  /// simply selects every Nth pixel.  You will typically want to
  /// apply some sort of anti-aliasing filter prior to calling this
  /// function.
  template <class ImageT>
  inline SubsampleView<ImageT> subsample( ImageT const& v, int32 subsampling_factor ) {
    return SubsampleView<ImageT>( v, subsampling_factor );
  }

  /// Subsample an image by integer factors in x and y.  Note that
  /// this function does not pre-smooth the image prior to
  /// subsampling: it simply selects every Nth pixel.  You will
  /// typically want to apply some sort of anti-aliasing filter prior
  /// to calling this function.
  template <class ImageT>
  inline SubsampleView<ImageT> subsample( ImageT const& v, int32 xfactor, int32 yfactor ) {
    return SubsampleView<ImageT>( v, xfactor, yfactor );
  }


  // *******************************************************************
  // select_col()
  // *******************************************************************

  /// Return a single column from an image
  /// \see vw::select_col
  template <class ImageT>
  class SelectColView;

  /// \cond INTERNAL
  // View type Traits
  template <class ImageT>
  struct IsMultiplyAccessible<SelectColView<ImageT> > : public IsMultiplyAccessible<ImageT> {};
  /// \endcond

  /// Extracts a single column of an image.  This function returns a
  /// writeable view of a single column of a multi-column image.
  /// \see vw::SelectColView
  template <class ImageT>
  SelectColView<ImageT> select_col( ImageViewBase<ImageT> const& v, int32 col ) {
    return SelectColView<ImageT>( v.impl(), col );
  }


  // *******************************************************************
  // select_row()
  // *******************************************************************

  /// Return a single row from an image
  /// \see vw::select_row
  template <class ImageT>
  class SelectRowView;

  // View type Traits
  template <class ImageT>
  struct IsMultiplyAccessible<SelectRowView<ImageT> > : public IsMultiplyAccessible<ImageT> {};

  /// Extracts a single row of an image.  This function returns a
  /// writeable view of a single row of a multi-row image.
  /// \see vw::SelectRowView
  template <class ImageT>
  SelectRowView<ImageT> select_row( ImageViewBase<ImageT> const& v, int32 row ) {
    return SelectRowView<ImageT>( v.impl(), row );
  }


  // *******************************************************************
  // select_plane()
  // *******************************************************************

  /// Return a single plane from a multi-plane image
  /// \see vw::select_plane
  template <class ImageT>
  class SelectPlaneView;

  /// \cond INTERNAL
  // View type Traits
  template <class ImageT>
  struct IsMultiplyAccessible<SelectPlaneView<ImageT> > : public IsMultiplyAccessible<ImageT> {};
  /// \endcond

  /// Extracts a single plane of a multi-plane image.  This function
  /// returns a writeable view of a single plane of a multi-plane
  /// image.  \see vw::SelectPlaneView
  template <class ImageT>
  SelectPlaneView<ImageT> select_plane( ImageViewBase<ImageT> const& v, int32 plane ) {
    return SelectPlaneView<ImageT>( v.impl(), plane );
  }


  // *******************************************************************
  // select_channel()
  // *******************************************************************

  /// A channel selecting functor, used by \ref select_channel().
  template <class ImageT>
  struct SelectChannelFunctor;

  /// Extracts a single channel of a multi-channel image.  This function
  /// returns a writeable view of a single channel of a multi-channel image.
  template <class ImageT>
  UnaryPerPixelView<ImageT,SelectChannelFunctor<ImageT> >
  inline select_channel( ImageViewBase<ImageT>& image, int32 channel ) {
    return UnaryPerPixelView<ImageT,SelectChannelFunctor<ImageT> >( image.impl(), SelectChannelFunctor<ImageT>(channel) );
  }

  /// Extracts a single channel of a multi-channel image (const overload).
  /// This function returns a writeable view of a single channel of a
  /// multi-channel image.
  template <class ImageT>
  UnaryPerPixelView<ImageT,SelectChannelFunctor<const ImageT> >
  inline select_channel( ImageViewBase<ImageT> const& image, int32 channel ) {
    return UnaryPerPixelView<ImageT,SelectChannelFunctor<const ImageT> >( image.impl(), SelectChannelFunctor<const ImageT>(channel) );
  }

  /// A convenience function to select the alpha channel of an image.
  template <class ImageT>
  UnaryPerPixelView<ImageT,SelectChannelFunctor<ImageT> >
  inline select_alpha_channel( ImageViewBase<ImageT>& image );

  /// A convenience function to select the alpha channel of an image
  /// (const overload).
  template <class ImageT>
  UnaryPerPixelView<ImageT,SelectChannelFunctor<const ImageT> >
  inline select_alpha_channel( ImageViewBase<ImageT> const& image );


  // *******************************************************************
  // channels_to_planes()
  // *******************************************************************

  /// A channels-to-planes pixel accessor adaptor.
  ///
  /// This is a special wrapper pixel accessor type, used by
  /// \ref vw::ChannelsToPlanesView, that treats the channels in a
  /// multi-channel image as planes.
  template <class ChildT>
  class ChannelsToPlanesAccessor;

  /// A view that turns a one plane, multi-channel view into a mulit-plane, one channel view.
  /// \see vw::channels_to_planes
  template <class ImageT>
  class ChannelsToPlanesView;

  /// \cond INTERNAL
  // View type Traits
  template <class ImageT>
  struct IsMultiplyAccessible<ChannelsToPlanesView<ImageT> > : public IsMultiplyAccessible<ImageT> {};
  /// \endcond

  /// Adapts a multi-channel image view so the channels are treated as
  /// planes.  This function returns a writeable view of of a
  /// single-plane, multi-channel image that treats the channels as
  /// planes.  This is primarily intended to simplify interfacing with
  /// legacy non-Vision-Workbench code that can only support
  /// fundamental pixel types.  If you are thinking about using this
  /// function and you are not trying to interface to legacy code then
  /// you are almost certainly doing something wrong.
  /// \see vw::ChannelsToPlanesView
  template <class ImageT>
  ChannelsToPlanesView<ImageT> channels_to_planes( ImageViewBase<ImageT> const& v ) {
    return ChannelsToPlanesView<ImageT>( v.impl() );
  }


  // *******************************************************************
  // planes_to_channels()
  // *******************************************************************

  /// A view that turns a multi-plane, single-channel view into a
  /// one-plane, multi-channel view.
  /// \see vw::planes_to_channels
  template <class PixelT, class ImageT>
  class PlanesToChannelsView;

  /// Adapts a multi-plane image view so the planes are treated as
  /// channels of the given pixel type.  This is primarily intended to
  /// simplify interfacing with legacy non-Vision-Workbench code that
  /// can only support fundamental pixel types.  If you are thinking
  /// about using this function and you are not trying to interface to
  /// legacy code then you are almost certainly doing something wrong.
  /// \see vw::PlanesToChannelsView
  /// \see vw::channels_to_planes
  template <class PixelT, class ImageT>
  PlanesToChannelsView<PixelT,ImageT> planes_to_channels( ImageViewBase<ImageT> const& v ) {
    return PlanesToChannelsView<PixelT,ImageT>( v.impl() );
  }


  // *******************************************************************
  // pixel_cast()
  // *******************************************************************

  /// A pixel casting functor, used by \ref pixel_cast().
  template <class PixelT>
  struct PixelCastFunctor;

  /// Create a new image view by statically casting the pixels to a
  /// new type.
  template <class PixelT, class ImageT>
  inline UnaryPerPixelView<ImageT,PixelCastFunctor<PixelT> > pixel_cast( ImageViewBase<ImageT> const& image ) {
    return UnaryPerPixelView<ImageT,PixelCastFunctor<PixelT> >( image.impl() );
  }


  // *******************************************************************
  // pixel_cast_rescale()
  // *******************************************************************

  /// A pixel casting functor, used by \ref pixel_cast_rescale().
  template <class PixelT>
  struct PixelCastRescaleFunctor;

  /// Create a new image view by casting and rescaling the pixels to a
  /// new type.
  template <class PixelT, class ImageT>
  inline UnaryPerPixelView<ImageT,PixelCastRescaleFunctor<PixelT> > pixel_cast_rescale( ImageViewBase<ImageT> const& image ) {
    return UnaryPerPixelView<ImageT,PixelCastRescaleFunctor<PixelT> >( image.impl() );
  }


  // *******************************************************************
  // channel_cast()
  // *******************************************************************

  /// A pixel channel casting functor, used by \ref channel_cast().
  template <class ChannelT>
  struct PixelChannelCastFunctor;

  /// Create a new image view by statically casting the channels of the pixels to a new type.
  template <class ChannelT, class ImageT>
  inline UnaryPerPixelView<ImageT,PixelChannelCastFunctor<ChannelT> > channel_cast( ImageViewBase<ImageT> const& image ) {
    return UnaryPerPixelView<ImageT,PixelChannelCastFunctor<ChannelT> >( image.impl() );
  }


  // *******************************************************************
  // channel_cast_round_and_clamp()
  // *******************************************************************

  /// A pixel channel casting, rounding and clamping functor, used by
  /// \ref channel_cast_round_and_clamp().
  template <class ChannelT>
  struct PixelChannelCastRoundClampFunctor;

  /// Create a new image view by casting, rounding, and clamping the channels of
  /// the pixels to a new type.
  template <class ChannelT, class ImageT>
  inline UnaryPerPixelView<ImageT,PixelChannelCastRoundClampFunctor<ChannelT> > channel_cast_round_and_clamp( ImageViewBase<ImageT> const& image ) {
    return UnaryPerPixelView<ImageT,PixelChannelCastRoundClampFunctor<ChannelT> >( image.impl() );
  }


  // *******************************************************************
  // channel_cast_rescale()
  // *******************************************************************

  /// A pixel channel casting and rescaling functor, used by
  /// \ref channel_cast_rescale().
  template <class ChannelT>
  struct PixelChannelCastRescaleFunctor;

  /// Create a new image view by casting and rescaling the channels of
  /// the pixels to a new type.
  template <class ChannelT, class ImageT>
  inline UnaryPerPixelView<ImageT,PixelChannelCastRescaleFunctor<ChannelT> > channel_cast_rescale( ImageViewBase<ImageT> const& image ) {
    return UnaryPerPixelView<ImageT,PixelChannelCastRescaleFunctor<ChannelT> >( image.impl() );
  }


  // *******************************************************************
  // weighted_rgb_to_gray()
  // *******************************************************************

  /// A weighted rgb-to-gray pixel conversion functor.
  //class WeightedRGBToGrayFunctor; // TODO: Why does forward declaration not work?

  /// A weighted rgb-to-gray pixel conversion functor.
  class WeightedRGBToGrayFunctor {
    double m_rw, m_gw, m_bw;
  public:
    template <class ArgsT> struct result {};
    template <class FuncT, class ChannelT> struct result<FuncT(PixelRGB<ChannelT>)> { typedef PixelGray<ChannelT> type; };
    template <class FuncT, class ChannelT> struct result<FuncT(PixelRGBA<ChannelT>)> { typedef PixelGrayA<ChannelT> type; };
    WeightedRGBToGrayFunctor( double rw, double gw, double bw ) : m_rw(rw), m_gw(gw), m_bw(bw) {}
    template <class ChannelT> inline PixelGrayA<ChannelT> operator()( PixelRGBA<ChannelT> const& rgb ) const {
      return weighted_rgb_to_gray( rgb, m_rw, m_gw, m_bw );
    }
    template <class ChannelT> inline PixelGray<ChannelT> operator()( PixelRGB<ChannelT> const& rgb ) const {
      return weighted_rgb_to_gray( rgb, m_rw, m_gw, m_bw );
    }
  };


  /// Weighted conversion from PixelRGBA to PixelGrayA using user-specified weights.
  template <class ImageT>
  inline UnaryPerPixelView<ImageT,WeightedRGBToGrayFunctor> weighted_rgb_to_gray( ImageViewBase<ImageT> const& image, double rw, double gw, double bw ) {
    return UnaryPerPixelView<ImageT,WeightedRGBToGrayFunctor>( image.impl(), WeightedRGBToGrayFunctor(rw,gw,bw) );
  }

  /// Weighted conversion from PixelRGBA to PixelGrayA using the default weights.
  template <class ImageT>
  inline UnaryPerPixelView<ImageT,WeightedRGBToGrayFunctor> weighted_rgb_to_gray( ImageViewBase<ImageT> const& image ) {
    WeightedRGBToGrayFunctor func( VW_RGB_TO_GRAY_R_WEIGHT, VW_RGB_TO_GRAY_G_WEIGHT, VW_RGB_TO_GRAY_B_WEIGHT );
    return UnaryPerPixelView<ImageT,WeightedRGBToGrayFunctor>( image.impl(), func );
  }

} // namespace vw

#include "Manipulation.tcc"

#endif // __VW_IMAGE_MANIPULATION_H__
