// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#ifndef __VW_STEREO_REWRITE_ALGORITHMS_H__
#define __VW_STEREO_REWRITE_ALGORITHMS_H__

#include <vw/Image/ImageView.h>
#include <numeric>
#include <valarray>

namespace vw {
namespace stereo {
namespace rewrite {

  // This quickly sums pixels in a box. This does not COPY the
  // image. It is okay to crop an image before feeding to this
  // function as crops of an ImageView is still memory striding.
  //
  // This expects the input to already be over cropped.
  template <class AccumulatorType, class ViewT>
  ImageView<typename PixelChannelCast<typename ViewT::pixel_type, AccumulatorType>::type>
  fast_box_sum( ImageViewBase<ViewT> const& image, Vector2i const& kernel ) {
    // Sanity check, constants, and types
    VW_ASSERT( kernel[0] % 2 == 1 && kernel[1] % 2 == 1,
               ArgumentErr() << "fast_box_sum: Kernel input not sized with odd values." );

    typedef typename ViewT::pixel_accessor PAccT;
    typedef typename ViewT::pixel_type PixelT;
    typedef typename PixelChannelCast<PixelT,AccumulatorType>::type AccumT;

    ViewT const& input( image.impl() );
    VW_DEBUG_ASSERT( input.cols() >= kernel[0] && input.rows() >= kernel[1],
                     ArgumentErr() << "fast_box_sum: Image is not big enough for kernel." );

    // Allocating output
    ImageView<AccumT> output( input.cols()-kernel[0]+1,
                              input.rows()-kernel[1]+1 );
    typedef typename ImageView<AccumT>::pixel_accessor OAccT;

    // Start column sum
    std::valarray<AccumT> col_sum(  input.cols() );
    const AccumT* col_sum_end = &col_sum[input.cols()];
    {
      PAccT row_input = input.origin();
      for ( int32 ky = kernel[1]; ky; --ky ) {
        PAccT col_input = row_input;
        AccumT *col = &col_sum[0];
        while ( col != col_sum_end ) {
          *col++ += *col_input;
          col_input.next_col();
        }
        row_input.next_row();
      }
    }

    // start rasterizing rows
    OAccT dst = output.origin();
    PAccT src_row_back = input.origin();
    PAccT src_row_front = input.origin().advance(0,kernel[1]);
    for ( int32 y = output.rows() - 1; y; --y ) {
      // Seed row sum
      AccumT row_sum(0);
      row_sum = std::accumulate(&col_sum[0],&col_sum[kernel[0]],
                                row_sum);

      // Sum down the row line
      AccumT const *cback = &col_sum[0], *cfront = &col_sum[kernel[0]];
      while( cfront != col_sum_end ) {
        *dst = row_sum;
        dst.next_col();
        row_sum += *cfront++ - *cback++;
      }
      *dst = row_sum;
      dst.next_col();

      // Update column sums
      PAccT src_col_back = src_row_back;
      PAccT src_col_front = src_row_front;
      for ( AccumT* col = &col_sum[0]; col != col_sum_end; ++col ) {
        *col += *src_col_front; // We do this in 2 lines to avoid casting.
        *col -= *src_col_back;  // I'm unsure if the assembly is still doing that.
        src_col_back.next_col();
        src_col_front.next_col();
      }

      // Update row iterators
      src_row_back.next_row();
      src_row_front.next_row();
    }

    { // Perform last sum down the line
      // Seed row sum
      AccumT row_sum(0);
      row_sum = std::accumulate(&col_sum[0],&col_sum[kernel[0]],
                                row_sum);

      // Sum down the row line
      AccumT const *cback = &col_sum[0], *cfront = &col_sum[kernel[0]];
      while( cfront != col_sum_end ) {
        *dst = row_sum;
        dst.next_col();
        row_sum += *cfront++ - *cback++;
      }
      *dst = row_sum;
    }

    return output;
  }

}}}

#endif//__VW_STEREO_REWRITE_ALGORITHMS_H__
