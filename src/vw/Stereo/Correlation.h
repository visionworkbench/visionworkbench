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


#ifndef __VW_STEREO_CORRELATION_H__
#define __VW_STEREO_CORRELATION_H__

#include <vw/Core/FundamentalTypes.h>
#include <vw/Image/PixelMask.h>
#include <vw/Stereo/CostFunctions.h>
#include <vw/Image/ImageViewRef.h>

#include <vector>
#include <algorithm>
#include <utility>

namespace vw {
namespace stereo {

  /// Lower level implementation function for calc_disparity.
  /// - The inputs must already be rasterized to safe sizes!
  /// - Since the inputs are rasterized, the input images must not be too big.
  template <class CostFuncT, class PixelT>
  ImageView<PixelMask<Vector2i>>
  best_of_search_convolution(ImageView<PixelT> const& left_raster,
                             ImageView<PixelT> const& right_raster,
                             BBox2i            const& left_region,
                             Vector2i          const& search_volume,
                             Vector2i          const& kernel_size) {
    
    typedef ImageView<PixelT> ImageType;
    typedef typename CostFuncT::accumulator_type AccumChannelT;
    typedef typename PixelChannelCast<PixelT,AccumChannelT>::type AccumT;
    typedef typename std::pair<AccumT,AccumT> QualT;

    // Build cost function which sometimes has side car data
    CostFuncT cost_function(left_raster, right_raster, kernel_size);

    // Result buffers
    Vector2i result_size = bounding_box(left_raster).size() - kernel_size + Vector2i(1,1);
    ImageView<PixelMask<Vector2i>> disparity_map(result_size[0], result_size[1]);
    std::fill(disparity_map.data(), disparity_map.data() + prod(result_size),
              PixelMask<Vector2i>(Vector2i()));
    // First channel is best, second is worst
    ImageView<QualT> quality_map(result_size[0], result_size[1]);
    
    // Storage buffers
    ImageView<AccumT> cost_metric      (result_size[0], result_size[1]);
    ImageView<AccumT> cost_applied     (left_raster.cols(), left_raster.rows());
    ImageView<PixelT> right_raster_crop(left_raster.cols(), left_raster.rows());

    // Loop across the disparity range we are searching over.
    Vector2i disparity (0, 0);
    for (disparity.y() = 0; disparity.y() != search_volume[1]; ++disparity.y()) {
      for (disparity.x() = 0; disparity.x() != search_volume[0]; ++disparity.x()) {
      
        // Compute correlations quickly by shifting the right image by the
        //  current disparity, computing the pixel difference at each location,
        //  and using fast_box_sum/cost_function to get the final convolution
        //  value at each location in "cost_metric"
      
        // There's only one raster here. Fast box sum calls each pixel
        // individually by pixel accessor. It only calls each pixel
        // once so there's no reason to copy/rasterize the cost result before hand.
        //
        // The cost function should also not be applying an edge
        // extension as we've already over cropped the input.
        
        right_raster_crop = crop(right_raster, bounding_box(left_raster)+disparity);
        cost_applied      = cost_function(left_raster, right_raster_crop);
        cost_metric       = fast_box_sum<AccumChannelT>(cost_applied, kernel_size);
        cost_function.cost_modification(cost_metric, disparity);

        // Loop across the region we want to compute disparities for.
        // - The correlation score for each pixel is located in "cost_metric"
        // - We update the best and worst disparity for each pixel in "quality_map"

        // These conditionals might be served outside of the iteration
        // of dx and dy. It would make the code slightly longer but
        // would avoid a conditional inside a double loop.
        const AccumT* cost_ptr     = cost_metric.data();
        const AccumT* cost_ptr_end = cost_metric.data() + prod(result_size);
        QualT* quality_ptr         = quality_map.data();
        PixelMask<Vector2i>* disparity_ptr = disparity_map.data();
        if (disparity != Vector2i(0,0)) {
          // Normal comparison operations
          while (cost_ptr != cost_ptr_end) {
            if (cost_function.quality_comparison(*cost_ptr, quality_ptr->first)) {
              // Better than best?
              quality_ptr->first = *cost_ptr;
              disparity_ptr->child() = disparity;
            } else if (!cost_function.quality_comparison(*cost_ptr, quality_ptr->second)) {
              // Worse than worse
              quality_ptr->second = *cost_ptr;
            }
            ++cost_ptr;
            ++quality_ptr;
            ++disparity_ptr;
          }
        } else {
          // Initializing quality_map and disparity_map with first result
          while (cost_ptr != cost_ptr_end) {
            quality_ptr->first = quality_ptr->second = *cost_ptr;
            ++cost_ptr;
            ++quality_ptr;
          }
        }
      } // End x loop
    } // End y loop

    // Determine validity of result (detects rare invalid cases)
    size_t invalid_count = 0;
    const QualT* quality_ptr      = quality_map.data();
    const QualT* quality_ptr_end  = quality_map.data() + prod(result_size);
    PixelMask<Vector2i>* disp_ptr = disparity_map.data();
    while (quality_ptr != quality_ptr_end) {
      if (quality_ptr->first == quality_ptr->second) {
        invalidate(*disp_ptr);
        ++invalid_count;
      }
      ++quality_ptr;
      ++disp_ptr;
    }
    //std::cout << "Invalidated " << invalid_count << " pixels in best_of_search_convolution2\n";

    return disparity_map;
  } // End function best_of_search_convolution


  /// This actually RASTERIZES/COPY the input images. It then makes an
  /// allocation to store current costs.
  ///
  /// Users pass us the active region of the left image. Hopefully this
  /// allows them to consider if they need edge extension.
  ///
  /// The return size of this function will be:
  ///     return size = left_region_size - kernel_size + 1.
  ///
  /// This means the user must take in account the kernel size for
  /// deciding the region size.
  ///
  /// The size of the area we're going to access in the right is
  /// calculated as follows:
  ///     right_region = left_region + search_volume - 1.
  ///
  ImageView<PixelMask<Vector2i>> calc_disparity(CostFunctionType cost_type,
                                                ImageViewRef<PixelGray<float>> const& left_in,
                                                ImageViewRef<PixelGray<float>> const& right_in,
                                                // Valid region in the left image
                                                BBox2i                 const& left_region,
                                                // Max disparity to search in right image
                                                Vector2i               const& search_volume,
                                                Vector2i               const& kernel_size);

  // TODO: Add some named accessors and functions to clean up our code!
  // These are useful for pyramid stereo correlation and several other
  // algorithms. The prime result is "SearchParam" where the first
  // bbox in the area of the result. The second bbox is the disparity
  // in the region defined by the first bbox.
  /// The first  BBox specifies a location in an image.
  /// The second BBox describes a 2D disparity range using the min and max.
  struct SearchParam : public std::pair<BBox2i,BBox2i> {
  
    SearchParam(BBox2i const& image_region, BBox2i const& disparity_range)
      : std::pair<BBox2i,BBox2i>(image_region, disparity_range){}
  
    // Simple wrapper functions to make the code easier to read
    BBox2i      & image_region   ()       {return this->first;}
    BBox2i const& image_region   () const {return this->first;}
    BBox2i      & disparity_range()       {return this->second;}
    BBox2i const& disparity_range() const {return this->second;}

    /// The amount of computation for correlation for a given
    /// SearchParam instance is proportional to the product of the
    /// dimensions of the image area and the disparity range.
    double search_volume() const {
      return (double)this->first.width ()*(double)this->first.height()*
        (double)this->second.width()*(double)this->second.height();
    }
  };
  
  /// Functor for finding faster SearchParam objects
  struct SearchParamLessThan {
    bool operator()(SearchParam const& A, SearchParam const& B) const {
      return A.search_volume() < B.search_volume();
    }
  };
  
  inline std::ostream& operator<<(std::ostream& os, SearchParam const& p) {
    os << "SearchParam: " << p.image_region() << ", " << p.disparity_range() << std::endl;
    return os;
  }
  

  /// Create fake left and right images and search volume.  Do a fake
  /// disparity calculation. Divide the run-time of this calculation
  /// by left region size times search box size. This will enable us
  /// to estimate how long disparity calculation takes for given cost
  /// function and kernel size.
  double calc_seconds_per_op(CostFunctionType cost_type, Vector2i const& kernel_size);
  
  /// Given one large region of an image to search and disparity ranges
  ///  for each pixel, try to break the region up into smaller regions
  ///  which contain narrower disparity ranges.
  /// - The goal is to isolate high-disparity range portions of the image
  ///   into small regions so that only those regions get searched over
  ///   the large disparity ranges while at the same time keeping
  ///   the number of seperate regions (list.size()) small.
  /// - This is a recursive function. It might be ideal to make this a
  ///   template. However in all cases so far, I've only applied to
  ///   PixelMask<Vector2i>. 
  /// - This should stay an image view as we'll be
  ///   accessing the image alot and randomly.
  bool subdivide_regions(ImageView<PixelMask<Vector2i>> const& disparity,
                         BBox2i const& current_bbox,
                         std::vector<SearchParam>& list, // Output goes here
                         Vector2i const& kernel_size,
                         int32 fail_count = 0);

}}

#endif//__VW_STEREO_CORRELATION_H__
