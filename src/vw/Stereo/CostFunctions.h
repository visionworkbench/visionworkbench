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
  // The base type attempts to extract accumulator type from ImageView
  template <class ImageT>
  struct AbsAccumulatorType {
    typedef typename AbsAccumulatorType<typename ImageChannelType<ImageT>::type>::type type;
  };
  template <class PixelT> struct AbsAccumulatorType<PixelMathBase<PixelT> > { typedef typename AbsAccumulatorType<typename CompoundChannelType<PixelT>::type>::type type; };
  template <> struct AbsAccumulatorType<bool>        { typedef vw::int16   type; };
  template <> struct AbsAccumulatorType<vw::uint8>   { typedef vw::int16   type; };
  template <> struct AbsAccumulatorType<vw::int8>    { typedef vw::int16   type; };
  template <> struct AbsAccumulatorType<vw::uint16>  { typedef vw::int32   type; };
  template <> struct AbsAccumulatorType<vw::int16>   { typedef vw::int32   type; };
  template <> struct AbsAccumulatorType<vw::uint32>  { typedef vw::int32   type; };
  template <> struct AbsAccumulatorType<vw::int32>   { typedef vw::int32   type; };
  template <> struct AbsAccumulatorType<vw::uint64>  { typedef vw::int64   type; };
  template <> struct AbsAccumulatorType<vw::int64>   { typedef vw::int64   type; };
  template <> struct AbsAccumulatorType<vw::float32> { typedef vw::float32 type; };
  template <> struct AbsAccumulatorType<vw::float64> { typedef vw::float64 type; };

  // Squared difference accumulator type (needs to be higher to avoid overflow)
  template <class ImageT>
  struct SqrDiffAccumulatorType {
    typedef typename SqrDiffAccumulatorType<typename CompoundChannelType<typename ImageT::pixel_type>::type>::type type;
  };
  template <class PixelT> struct SqrDiffAccumulatorType<PixelMathBase<PixelT> > { typedef typename SqrDiffAccumulatorType<typename CompoundChannelType<PixelT>::type>::type type; };
  template <> struct SqrDiffAccumulatorType<bool>        { typedef vw::int16   type; };
  template <> struct SqrDiffAccumulatorType<vw::uint8>   { typedef vw::int32   type; };
  template <> struct SqrDiffAccumulatorType<vw::int8>    { typedef vw::int32   type; };
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
      inline Arg1T operator()( Arg1T arg1, Arg2T arg2 ) const {
        return ( arg1 < arg2 ) ? arg2 - arg1 : arg1 - arg2;
      }

      inline float operator()( float arg1, float arg2 ) const {
        return fabs( arg1 - arg2 );
      }

      inline double operator()( double arg1, double arg2 ) const {
        return fabs( arg1 - arg2 );
      }
    };

    template <class Pixel1T, class Pixel2T>
    Pixel1T operator()( Pixel1T const& arg1, Pixel2T const& arg2 ) const {
      return BinaryCompoundFunctor<Helper,Pixel1T,Pixel2T>()(arg1, arg2);
    }
  };

  struct SquaredDifferenceFunctor {
    struct Helper : BinaryReturnUnaryTemplateBind1st<SqrDiffAccumulatorType> {
      template <class ArgT>
      inline typename SqrDiffAccumulatorType<ArgT>::type
      operator()( ArgT const& arg1, ArgT const& arg2 ) const {
        return (arg1 - arg2)*(arg1 - arg2);
      }

      inline SqrDiffAccumulatorType<uint8>::type
      operator()( uint8 const& arg1, uint8 const& arg2 ) const {
        return (arg1 < arg2 ) ? (arg2 - arg1)*(arg2 - arg1) :
          (arg1 - arg2)*(arg1 - arg2);
      }

      inline SqrDiffAccumulatorType<uint16>::type
      operator()( uint16 const& arg1, uint16 const& arg2 ) const {
        return (arg1 < arg2 ) ? (arg2 - arg1)*(arg2 - arg1) :
          (arg1 - arg2)*(arg1 - arg2);
      }

      inline SqrDiffAccumulatorType<uint32>::type
      operator()( uint32 const& arg1, uint32 const& arg2 ) const {
        return (arg1 < arg2 ) ? (arg2 - arg1)*(arg2 - arg1) :
          (arg1 - arg2)*(arg1 - arg2);
      }

      inline SqrDiffAccumulatorType<uint64>::type
      operator()( uint64 const& arg1, uint64 const& arg2 ) const {
        return (arg1 < arg2 ) ? (arg2 - arg1)*(arg2 - arg1) :
          (arg1 - arg2)*(arg1 - arg2);
      }
    };

    template <class Pixel1T, class Pixel2T>
    typename PixelChannelCast<Pixel1T,typename SqrDiffAccumulatorType<typename PixelChannelType<Pixel1T>::type>::type>::type
    operator()( Pixel1T const& arg1, Pixel2T const& arg2 ) const {
      return BinaryCompoundFunctor<Helper,Pixel1T,Pixel2T>()(arg1, arg2);
    }

    // So Boost can determine return type
    template <class Args> struct result { typedef void type; };
    template <class FuncT, class Arg1T, class Arg2T>
    struct result<FuncT(Arg1T,Arg2T)> {
      typedef typename PixelChannelCast<Arg1T,typename SqrDiffAccumulatorType<typename PixelChannelType<Arg1T>::type>::type>::type type;
    };
  };

  struct CrossCorrelationFunctor {
    struct Helper : BinaryReturnUnaryTemplateBind1st<SqrDiffAccumulatorType> {
      template <class ArgT>
      inline typename boost::enable_if<boost::is_float<typename SqrDiffAccumulatorType<ArgT>::type>,typename SqrDiffAccumulatorType<ArgT>::type>::type
      operator()( ArgT const& arg1, ArgT const& arg2 ) const {
        return arg1 * arg2;
      }

      // This is a custom version for integers that saves room in the
      // integer for summing. Realize that we are summing squares! It
      // is very hard to avoid overflowing.
      template <class ArgT>
      inline typename boost::disable_if<boost::is_float<typename SqrDiffAccumulatorType<ArgT>::type>,typename SqrDiffAccumulatorType<ArgT>::type>::type
      operator()( ArgT const& arg1, ArgT const& arg2 ) const {
        return typename SqrDiffAccumulatorType<ArgT>::type(arg1 * arg2) / 256;
      }
    };

    template <class Pixel1T, class Pixel2T>
    typename PixelChannelCast<Pixel1T,typename SqrDiffAccumulatorType<typename PixelChannelType<Pixel1T>::type>::type>::type
    operator()( Pixel1T const& arg1, Pixel2T const& arg2 ) const {
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

  template <class ImageT, bool IsInteger>
  struct AbsoluteCost {
    typedef typename AbsAccumulatorType<ImageT>::type accumulator_type;
    typedef typename PixelChannelCast<typename ImageT::pixel_type, accumulator_type>::type pixel_accumulator_type;

    // Does nothing
    template <class ImageT1, class ImageT2>
    AbsoluteCost( ImageViewBase<ImageT1> const& /*left*/,
                        ImageViewBase<ImageT2> const& /*right*/,
                        Vector2i const& /*kernel_size*/ ) {}

    template <class ImageT1, class ImageT2>
    BinaryPerPixelView<ImageT1,ImageT2,AbsDifferenceFunctor>
    operator()( ImageViewBase<ImageT1> const& left,
                ImageViewBase<ImageT2> const& right ) const {
      typedef BinaryPerPixelView<ImageT1,ImageT2,AbsDifferenceFunctor> result_type;
      return result_type(left.impl(),right.impl());
    }

    // Does nothing
    inline void cost_modification( ImageView<pixel_accumulator_type>& /*cost_metric*/,
                                   Vector2i const& /*disparity*/ ) const {}

    inline bool quality_comparison( accumulator_type cost,
                                    accumulator_type quality ) const {
      return cost < quality;
    }
  };

  template <class ImageT, bool IsInteger>
  struct SquaredCost {
    typedef typename SqrDiffAccumulatorType<ImageT>::type accumulator_type;
    typedef typename PixelChannelCast<typename ImageT::pixel_type, accumulator_type>::type pixel_accumulator_type;

    // Does nothing
    template <class ImageT1, class ImageT2>
    SquaredCost( ImageViewBase<ImageT1> const& /*left*/,
                            ImageViewBase<ImageT2> const& /*right*/,
                            Vector2i const& /*kernel_size*/ ) {}

    template <class ImageT1, class ImageT2>
    BinaryPerPixelView<ImageT1,ImageT2,SquaredDifferenceFunctor>
    operator()( ImageViewBase<ImageT1> const& left,
                ImageViewBase<ImageT2> const& right ) const {
      typedef BinaryPerPixelView<ImageT1,ImageT2,SquaredDifferenceFunctor> result_type;
      return result_type(left.impl(),right.impl());
    }

    // Does nothing
    inline void cost_modification( ImageView<pixel_accumulator_type>& /*cost_metric*/,
                                   Vector2i const& /*disparity*/ ) const {}

    inline bool quality_comparison( accumulator_type cost,
                                    accumulator_type quality ) const {
      return cost < quality;
    }
  };

  // The float version of Cross Correlation
  template <class ImageT, bool IsInteger>
  struct NCCCost {
    typedef typename SqrDiffAccumulatorType<ImageT>::type accumulator_type;
    typedef typename PixelChannelCast<typename ImageT::pixel_type, accumulator_type>::type pixel_accumulator_type;
    ImageView<pixel_accumulator_type> left_precision, right_precision;

    template <class ImageT1, class ImageT2>
    NCCCost( ImageViewBase<ImageT1> const& left,
                                 ImageViewBase<ImageT2> const& right,
                                 Vector2i const& kernel_size ) {
      left_precision = pixel_accumulator_type(1) /
        fast_box_sum<accumulator_type>(square(left.impl()), kernel_size);
      right_precision = pixel_accumulator_type(1) /
        fast_box_sum<accumulator_type>(square(right.impl()), kernel_size);
    }

    template <class ImageT1, class ImageT2>
    BinaryPerPixelView<ImageT1,ImageT2,CrossCorrelationFunctor>
    operator()( ImageViewBase<ImageT1> const& left,
                ImageViewBase<ImageT2> const& right ) const {
      typedef BinaryPerPixelView<ImageT1,ImageT2,CrossCorrelationFunctor> result_type;
      return result_type(left.impl(),right.impl());
    }

    inline void cost_modification( ImageView<pixel_accumulator_type>& cost_metric,
                                   Vector2i const& disparity ) const {
      cost_metric *= sqrt( left_precision * crop(right_precision,
                                                 bounding_box(left_precision)+disparity) );
    }

    inline bool quality_comparison( accumulator_type cost,
                                    accumulator_type quality ) const {
      return cost > quality;
    }
  };

  // Specialization of the CrossCorrelation Cost Functor so that it
  // avoids integer overflow.
  template <class ImageT>
  struct NCCCost<ImageT, true> {
    typedef typename SqrDiffAccumulatorType<ImageT>::type accumulator_type;
    typedef typename PixelChannelCast<typename ImageT::pixel_type, accumulator_type>::type pixel_accumulator_type;
    ImageView<pixel_accumulator_type> left_variance, right_variance;

    template <class ImageT1, class ImageT2>
    NCCCost( ImageViewBase<ImageT1> const& left,
                                 ImageViewBase<ImageT2> const& right,
                                 Vector2i const& kernel_size ) {
      left_variance =
        fast_box_sum<accumulator_type>(square(pixel_cast<pixel_accumulator_type>(left.impl()))/256, kernel_size);
      right_variance =
        fast_box_sum<accumulator_type>(square(pixel_cast<pixel_accumulator_type>(right.impl()))/256, kernel_size);
    }

    template <class ImageT1, class ImageT2>
    BinaryPerPixelView<ImageT1,ImageT2,CrossCorrelationFunctor>
    operator()( ImageViewBase<ImageT1> const& left,
                ImageViewBase<ImageT2> const& right ) const {
      typedef BinaryPerPixelView<ImageT1,ImageT2,CrossCorrelationFunctor> result_type;
      return result_type(left.impl(),right.impl());
    }

    inline void cost_modification( ImageView<pixel_accumulator_type>& cost_metric,
                                   Vector2i const& disparity ) const {
      cost_metric =
        (64 * cost_metric) / ( sqrt( left_variance * crop(right_variance,
                                                          bounding_box(left_variance)+disparity) ) / 64 );
    }

    inline bool quality_comparison( accumulator_type cost,
                                    accumulator_type quality ) const {
      return cost > quality;
    }
  };

}}

#endif//__VW_STEREO_COSTFUNCTIONS_H__
