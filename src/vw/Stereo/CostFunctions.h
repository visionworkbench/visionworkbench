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


#ifndef __VW_STEREO_COSTFUNCTIONS_H__
#define __VW_STEREO_COSTFUNCTIONS_H__

#include <vw/Image/Manipulation.h>
#include <vw/Image/ImageMath.h>
#include <vw/Stereo/Algorithms.h>

namespace vw {
namespace stereo {

  template <template<class> class ResultT>
  struct BinaryReturnUnaryTemplateBind1st {
    template <class Args> struct result { typedef void type; };
    template <class FuncT, class Arg1T, class Arg2T>
    struct result<FuncT(Arg1T,Arg2T)> {
      typedef typename ResultT<Arg1T>::type type;
    };
  };

  /// Given a type, these traits classes help to determine a suitable
  /// working type for accumulation operations or other intermediate
  /// results that require computational headroom.
  //
  // The base type attempts to extract accumulator type from ImageView. For bool and int
  // types use an int64 accumulator to avoid overflow. All this logic is usually
  // used with float images.
  // TODO(oalexan1): The accumulators are best to be always float64.
  template <class ImageT>
  struct AbsAccumulatorType {
    typedef typename AbsAccumulatorType<typename ImageChannelType<ImageT>::type>::type type;
  };
  template <class PixelT> struct AbsAccumulatorType<PixelMathBase<PixelT> > {
    typedef typename AbsAccumulatorType<typename CompoundChannelType<PixelT>::type>::type type;
  };
  template <> struct AbsAccumulatorType<bool>        { typedef vw::int64   type; };
  template <> struct AbsAccumulatorType<vw::uint8>   { typedef vw::int64   type; };
  template <> struct AbsAccumulatorType<vw::int8>    { typedef vw::int64   type; };
  template <> struct AbsAccumulatorType<vw::uint16>  { typedef vw::int64   type; };
  template <> struct AbsAccumulatorType<vw::int16>   { typedef vw::int64   type; };
  template <> struct AbsAccumulatorType<vw::uint32>  { typedef vw::int64   type; };
  template <> struct AbsAccumulatorType<vw::int32>   { typedef vw::int64   type; };
  template <> struct AbsAccumulatorType<vw::uint64>  { typedef vw::int64   type; };
  template <> struct AbsAccumulatorType<vw::int64>   { typedef vw::int64   type; };
  template <> struct AbsAccumulatorType<vw::float32> { typedef vw::float32 type; };
  template <> struct AbsAccumulatorType<vw::float64> { typedef vw::float64 type; };

  // Squared difference accumulator type (needs to be higher to avoid overflow)
  template <class ImageT>
  struct SqrDiffAccumulatorType {
    typedef typename SqrDiffAccumulatorType<typename CompoundChannelType<typename ImageT::pixel_type>::type>::type type;
  };
  template <class PixelT> struct SqrDiffAccumulatorType<PixelMathBase<PixelT> > {
    typedef typename SqrDiffAccumulatorType<typename CompoundChannelType<PixelT>::type>::type type;
  };
  template <> struct SqrDiffAccumulatorType<bool>        { typedef vw::int64   type; };
  template <> struct SqrDiffAccumulatorType<vw::uint8>   { typedef vw::int64   type; };
  template <> struct SqrDiffAccumulatorType<vw::int8>    { typedef vw::int64   type; };
  template <> struct SqrDiffAccumulatorType<vw::uint16>  { typedef vw::int64   type; };
  template <> struct SqrDiffAccumulatorType<vw::int16>   { typedef vw::int64   type; };
  template <> struct SqrDiffAccumulatorType<vw::uint32>  { typedef vw::float32 type; };
  template <> struct SqrDiffAccumulatorType<vw::int32>   { typedef vw::float32 type; };
  template <> struct SqrDiffAccumulatorType<vw::uint64>  { typedef vw::float32 type; };
  template <> struct SqrDiffAccumulatorType<vw::int64>   { typedef vw::float32 type; };
  template <> struct SqrDiffAccumulatorType<vw::float32> { typedef vw::float64 type; };
  template <> struct SqrDiffAccumulatorType<vw::float64> { typedef vw::float64 type; };

  // Simplest Functors
  struct AbsDifferenceFunctor : BinaryReturn1stType {
    struct Helper : BinaryReturn1stType {
      template <class Arg1T, class Arg2T>
      inline Arg1T operator()(Arg1T arg1, Arg2T arg2) const {
        return (arg1 < arg2) ? arg2 - arg1 : arg1 - arg2;
      }

      inline float operator()(float arg1, float arg2) const {
        return fabs(arg1 - arg2);
      }

      inline double operator()(double arg1, double arg2) const {
        return fabs(arg1 - arg2);
      }
    };

    template <class Pixel1T, class Pixel2T>
    Pixel1T operator()(Pixel1T const& arg1, Pixel2T const& arg2) const {
      return BinaryCompoundFunctor<Helper,Pixel1T,Pixel2T>()(arg1, arg2);
    }
  };

  struct SquaredDifferenceFunctor {
    struct Helper : BinaryReturnUnaryTemplateBind1st<SqrDiffAccumulatorType> {
      template <class ArgT>
      inline typename SqrDiffAccumulatorType<ArgT>::type
      operator()(ArgT const& arg1, ArgT const& arg2) const {
        return (arg1 - arg2)*(arg1 - arg2);
      }

      inline SqrDiffAccumulatorType<uint8>::type
      operator()(uint8 const& arg1, uint8 const& arg2) const {
        return (arg1 < arg2) ? (arg2 - arg1)*(arg2 - arg1) :
          (arg1 - arg2)*(arg1 - arg2);
      }

      inline SqrDiffAccumulatorType<uint16>::type
      operator()(uint16 const& arg1, uint16 const& arg2) const {
        return (arg1 < arg2) ? (arg2 - arg1)*(arg2 - arg1) :
          (arg1 - arg2)*(arg1 - arg2);
      }

      inline SqrDiffAccumulatorType<uint32>::type
      operator()(uint32 const& arg1, uint32 const& arg2) const {
        return (arg1 < arg2) ? (arg2 - arg1)*(arg2 - arg1) :
          (arg1 - arg2)*(arg1 - arg2);
      }

      inline SqrDiffAccumulatorType<uint64>::type
      operator()(uint64 const& arg1, uint64 const& arg2) const {
        return (arg1 < arg2) ? (arg2 - arg1)*(arg2 - arg1) :
          (arg1 - arg2)*(arg1 - arg2);
      }
    };

    template <class Pixel1T, class Pixel2T>
    typename PixelChannelCast<Pixel1T, typename SqrDiffAccumulatorType
                              <typename PixelChannelType<Pixel1T>::type>::type>::type
    operator()(Pixel1T const& arg1, Pixel2T const& arg2) const {
      return BinaryCompoundFunctor<Helper,Pixel1T,Pixel2T>()(arg1, arg2);
    }

    // So Boost can determine return type
    template <class Args> struct result { typedef void type; };
    template <class FuncT, class Arg1T, class Arg2T>
    struct result<FuncT(Arg1T,Arg2T)> {
      typedef typename PixelChannelCast<Arg1T,
        typename SqrDiffAccumulatorType<typename PixelChannelType<Arg1T>::type>::type>::type type;
    };
  };

  struct CrossCorrelationFunctor {
    struct Helper : BinaryReturnUnaryTemplateBind1st<SqrDiffAccumulatorType> {
      template <class ArgT>
      inline typename boost::enable_if<boost::is_float<typename SqrDiffAccumulatorType<ArgT>::type>,typename SqrDiffAccumulatorType<ArgT>::type>::type
      operator()(ArgT const& arg1, ArgT const& arg2) const {
        return arg1 * arg2;
      }

      // This is a custom version for integers that saves room in the
      // integer for summing. Realize that we are summing squares! It
      // is very hard to avoid overflowing.
      template <class ArgT>
      inline typename boost::disable_if<boost::is_float<typename SqrDiffAccumulatorType<ArgT>::type>,typename SqrDiffAccumulatorType<ArgT>::type>::type
      operator()(ArgT const& arg1, ArgT const& arg2) const {
        return typename SqrDiffAccumulatorType<ArgT>::type(arg1 * arg2) / 256;
      }
    };

    template <class Pixel1T, class Pixel2T>
    typename PixelChannelCast<Pixel1T,typename SqrDiffAccumulatorType<typename PixelChannelType<Pixel1T>::type>::type>::type
    operator()(Pixel1T const& arg1, Pixel2T const& arg2) const {
      return BinaryCompoundFunctor<Helper,Pixel1T,Pixel2T>()(arg1, arg2);
    }

    // For Boost can determine return type
    template <class Args> struct result { typedef void type; };
    template <class FuncT, class Arg1T, class Arg2T>
    struct result<FuncT(Arg1T,Arg2T)> {
      typedef typename PixelChannelCast<Arg1T,typename SqrDiffAccumulatorType<typename PixelChannelType<Arg1T>::type>::type>::type type;
    };
  };

  // ImageView Functors
  enum CostFunctionType {
    ABSOLUTE_DIFFERENCE,
    SQUARED_DIFFERENCE,
    CROSS_CORRELATION,
    CENSUS_TRANSFORM,
    TERNARY_CENSUS_TRANSFORM
  };

  template <class ImageT>
  struct AbsoluteCost {
    typedef typename AbsAccumulatorType<ImageT>::type accumulator_type;
    typedef typename PixelChannelCast<typename ImageT::pixel_type, accumulator_type>::type pixel_accumulator_type;

    // Does nothing
    AbsoluteCost(ImageViewBase<ImageT> const& /*left*/,
                 ImageViewBase<ImageT> const& /*right*/,
                 Vector2i const& /*kernel_size*/ ) {}

    BinaryPerPixelView<ImageT, ImageT, AbsDifferenceFunctor>
    operator()(ImageViewBase<ImageT> const& left, ImageViewBase<ImageT> const& right) const {
      typedef BinaryPerPixelView<ImageT, ImageT, AbsDifferenceFunctor> result_type;
      return result_type(left.impl(),right.impl());
    }

    // Does nothing
    inline void cost_modification(ImageView<pixel_accumulator_type>& /*cost_metric*/,
                                   Vector2i const& /*disparity*/) const {}

    inline bool quality_comparison(accumulator_type cost,
                                    accumulator_type quality) const {
      return cost < quality;
    }
  };

  template <class ImageT>
  struct SquaredCost {
    typedef typename SqrDiffAccumulatorType<ImageT>::type accumulator_type;
    typedef typename PixelChannelCast<typename ImageT::pixel_type,
                                      accumulator_type>::type pixel_accumulator_type;

    // Does nothing
    SquaredCost(ImageT const& /*left*/, ImageT const& /*right*/,
                            Vector2i const& /*kernel_size*/ ) {}

    BinaryPerPixelView<ImageT, ImageT, SquaredDifferenceFunctor>
    operator()(ImageT const& left, ImageT const& right ) const {
      typedef BinaryPerPixelView<ImageT, ImageT, SquaredDifferenceFunctor> result_type;
      return result_type(left.impl(),right.impl());
    }

    // Does nothing
    inline void cost_modification(ImageView<pixel_accumulator_type>& /*cost_metric*/,
                                   Vector2i const& /*disparity*/) const {}

    inline bool quality_comparison(accumulator_type cost, accumulator_type quality) const {
      return cost < quality;
    }
  };

  // Cross-correlation. Only float32 or float64 images are expected.
  template <class ImageT>
  struct NCCCost {
    typedef typename SqrDiffAccumulatorType<ImageT>::type accumulator_type;
    typedef typename PixelChannelCast<typename ImageT::pixel_type,
                                      accumulator_type>::type pixel_accumulator_type;
    ImageView<pixel_accumulator_type> left_precision, right_precision;

    NCCCost(ImageT const& left, ImageT const& right, Vector2i const& kernel_size ) {
      left_precision = pixel_accumulator_type(1) /
        fast_box_sum<accumulator_type>(square(left.impl()), kernel_size);
      right_precision = pixel_accumulator_type(1) /
        fast_box_sum<accumulator_type>(square(right.impl()), kernel_size);
    }
    
    BinaryPerPixelView<ImageT, ImageT, CrossCorrelationFunctor>
    operator()( ImageT const& left, ImageT const& right ) const {
      typedef BinaryPerPixelView<ImageT, ImageT, CrossCorrelationFunctor> result_type;
      return result_type(left.impl(),right.impl());
    }

    inline void cost_modification(ImageView<pixel_accumulator_type>& cost_metric,
                                   Vector2i const& disparity) const {
      cost_metric *= sqrt(left_precision * crop(right_precision,
                                                bounding_box(left_precision)+disparity));
    }

    inline bool quality_comparison(accumulator_type cost, accumulator_type quality) const {
      return cost > quality;
    }
  };

}}

#endif//__VW_STEREO_COSTFUNCTIONS_H__
