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


#ifndef __VW_STEREO_ALGORITHMS_H__
#define __VW_STEREO_ALGORITHMS_H__

#include <vw/Image/ImageView.h>
#include <vw/Image/Filter.h>
#include <vw/Image/EdgeExtension.h>
#include <numeric>
#include <valarray>

#include <vw/Image/PixelTypes.h>
#include <vw/Image/PixelMask.h>

#include <vw/Math/Functors.h>

namespace vw {
namespace stereo {

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

    ViewT const& input( image.impl() ); // This just helps keep the code cleaner
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


/// Apply a median filter to a disparity image
template <typename T>
inline void disparity_median_filter(ImageView<PixelMask<Vector<T, 2> > > const& disparity_in,
                                    ImageView<PixelMask<Vector<T, 2> > >      & disparity_out,
                                    int kernel_size) {

  int half_kernel = (kernel_size - 1) / 2;
  int num_vals = kernel_size*kernel_size;
  disparity_out = disparity_in;

  if (kernel_size < 3) // No smoothing called for
    return;
    
  // Output pixel loop
  for (int row=half_kernel; row<disparity_in.rows()-half_kernel; ++row) {
    for (int col=half_kernel; col<disparity_in.cols()-half_kernel; ++col) {

      if (!is_valid(disparity_in(col,row)))
        continue;

     // Loop through the kernel
     std::vector<T> dx(num_vals), dy(num_vals);     
     int index = 0;
     for (int r=row-half_kernel; r<=row+half_kernel; ++r) {
       for (int c=col-half_kernel; c<=col+half_kernel; ++c) {
         if (is_valid(disparity_in(c,r))) {
           dx[index] = disparity_in(c,r)[0];
           dy[index] = disparity_in(c,r)[1];
           ++index;
         }
       }
     }
     if (index == 0)
       continue;
       
     dx.resize(index);
     dy.resize(index);
     T median_x = math::destructive_median(dx);
     T median_y = math::destructive_median(dy);
      
     disparity_out(col, row) = PixelMask<Vector<T, 2> >(median_x, median_y);
        
    } 
  } // End loop through pixels

} // End disparity_median_filter

/// Replace isolated disparity values with majority surrounding disparity values.
inline void disparity_neighbor_filter(ImageView<PixelMask<Vector2i> > const& disparity_in,
                                    ImageView<PixelMask<Vector2i> >      & disparity_out) {

  const int COPY_COUNT = 5;

  disparity_out = disparity_in;

  std::vector<int> counts(8);
  std::vector<PixelMask<Vector2i> > vals(8);
  for (int row=1; row<disparity_in.rows()-1; ++row) {
    for (int col=1; col<disparity_in.cols()-1; ++col) {
     
      vals[0] = disparity_in(col-1, row-1);
      vals[1] = disparity_in(col,   row-1);
      vals[2] = disparity_in(col+1, row-1);
      vals[3] = disparity_in(col-1, row  );
      vals[4] = disparity_in(col+1, row  );
      vals[5] = disparity_in(col-1, row+1);
      vals[6] = disparity_in(col,   row+1);
      vals[7] = disparity_in(col+1, row+1);
           
      int max_count = 0;
      int max_index = 0;
      for (int i=0; i<8; ++i) {
        counts[i] = 0;
        if (!is_valid(vals[i]))
          continue;
        for (int j=0; j<8; ++j) {
          if (vals[i] == vals[j])
            ++counts[i];
        }
        if (counts[i] > max_count) {
          max_count = counts[i];
          max_index = i;
        }
      } // End max finder
      
      if (max_count >= COPY_COUNT)
        disparity_out(col, row) = vals[max_index];
        
    } 
  } // End loop through pixels

} // End disparity_neighbor_filter


// TODO: Move to Image/AlgorithmFunctions.h
// - Also could speed this up!
template <class ImageT> 
inline void
texture_measure(ImageT          const& input_image, 
                ImageView<float>     & output_image,
                int kernel_size        = 9,
                double gradient_weight = 0.5, 
                double stddev_weight   = 0.5) {

  // Init the output image
  output_image.set_size(input_image.cols(), input_image.rows());
  //set_all(output_image, typename ImageT::pixel_type(0));

  int half_kernel = (kernel_size - 1) / 2;

  //std::cout << "Generating dx/dy images...\n";

  // Derivative images are precomputed here
  EdgeExtensionView<ImageT, ConstantEdgeExtension> input_wrap(input_image);
  EdgeExtensionView<ImageView<float>, ConstantEdgeExtension> dx_image(derivative_filter(input_image, 1, 0));
  EdgeExtensionView<ImageView<float>, ConstantEdgeExtension> dy_image(derivative_filter(input_image, 0, 1));

  //std::cout << "Performing dynamic smoothing...\n";

  for (int row=0; row<output_image.rows(); ++row) {
    for (int col=0; col<output_image.cols(); ++col) {

      // Check for input pixel validity
      if (!is_valid(input_image(col,row)))
        continue; // TODO

      // Compute the mean intensity value within the region
      // - TODO: Speed this computation up
      double mean_value = 0, count=0;
      for (int r=row-half_kernel; r<=row+half_kernel; ++r) {
        for (int c=col-half_kernel; c<=col+half_kernel; ++c) {
          if (is_valid(input_wrap(c,r))) {
            mean_value += input_wrap(c,r);
            count += 1.0;
          }
        }
      }
      if (count < 1.0)
        continue; // TODO
      mean_value /= count;

      double gradient_total = 0, stddev_total = 0;
      for (int r=row-half_kernel; r<=row+half_kernel; ++r) {
        for (int c=col-half_kernel; c<=col+half_kernel; ++c) {
          if (!is_valid(input_wrap(c,r)))
            continue;
          
          gradient_total += std::abs(dx_image(c,r)) + std::abs(dy_image(c,r));
          stddev_total   += pow(input_wrap(c,r) - mean_value, 2);
        }
      }
     
      gradient_total = gradient_total / (2.0*count);  
      stddev_total   = sqrt(stddev_total / count);
     
      float score = gradient_total*gradient_weight + stddev_total*stddev_weight;
      output_image(col, row) = score; 
        
    } 
  } // End loop through pixels
  
} // end function texture_measure()



/// Smooth out a disparity result with intensity inversely proportional to the
///  amount of texture present in the input image.
template <typename T>
inline void texture_preserving_disparity_filter(ImageView<PixelMask<Vector<T, 2> > > const& disparity_in,
                                                ImageView<PixelMask<Vector<T, 2> > >      & disparity_out,
                                                ImageView<float                    > const& texture_image,
                                                float texture_max     = 0.15,
                                                int   max_kernel_size = 11) {

  float texture_scale = max_kernel_size / texture_max;
  
  EdgeExtensionView<ImageView<PixelMask<Vector<T, 2> > >, ConstantEdgeExtension> input_wrap(disparity_in);
  
  // Pixels are unchanged by default
  disparity_out = disparity_in;
  if ((max_kernel_size < 3) || (texture_max <= 0))
    return; // Quit here if no smoothing is needed.

  //std::cout << "Filter output size: " << bounding_box(disparity_in) << std::endl;
  //std::cout << "Texture_scale: " << texture_scale << std::endl;

  // Output pixel loop
  for (int row=0; row<disparity_in.rows(); ++row) {
    for (int col=0; col<disparity_in.cols(); ++col) {

      // TODO: Texture image should always be positive!
      // Skip invalid pixels
      if (!is_valid(disparity_in(col,row)) || (texture_image(col, row)<0) )
        continue;

      float adjusted_texture = (texture_max-texture_image(col, row));
      if (adjusted_texture < 0)
        adjusted_texture = 0;
      int kernel_size = floor(adjusted_texture * texture_scale);
      if (kernel_size % 2 == 0) // Make odd
        kernel_size += 1;
      if (kernel_size < 3) // No smoothing case
        continue;

      if (kernel_size > max_kernel_size) {
        std::cout << "texture_preserving_disparity_filter: Oversize kernel error! = " << kernel_size << std::endl;
        std::cout << texture_image(col, row) << std::endl;
        std::cout << adjusted_texture << std::endl;
        continue;
      }

      // Compute the kernel size to use at this location
      int half_kernel = (kernel_size - 1) / 2;

      // Loop through the kernel
      double dx_total=0, dy_total=0;
      double count=0;
      for (int r=row-half_kernel; r<=row+half_kernel; ++r) {
        for (int c=col-half_kernel; c<=col+half_kernel; ++c) {
          if (is_valid(input_wrap(c,r))) {
            dx_total += input_wrap(c,r)[0];
            dy_total += input_wrap(c,r)[1];
            count += 1.0;
          }
        }
      }
      if (count < 1.0)
        continue;
     
      disparity_out(col, row) = PixelMask<Vector<T, 2> >(dx_total/count, dy_total/count);
        
    } 
  } // End loop through pixels                       
} // end texture_preserving_disparity_filter function


}} // end vw::stereo

#endif//__VW_STEREO_ALGORITHMS_H__
