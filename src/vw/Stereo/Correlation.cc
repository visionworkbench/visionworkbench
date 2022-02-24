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

#include <vw/Math/BBox.h>
#include <vw/Image/Statistics.h>
#include <vw/Core/Exception.h>
#include <vw/Stereo/Correlation.h>
#include <vw/Stereo/Algorithms.h>
#include <vw/Core/Stopwatch.h>
#include <vw/Image/AlgorithmFunctions.h>

namespace vw {
namespace stereo {

 
bool subdivide_regions(ImageView<PixelMask<Vector2i> > const& disparity,
                       BBox2i const& current_bbox,
                       std::vector<SearchParam>& list,
                       Vector2i const& kernel_size,
                       int32 fail_count) {

    // Looking at the 2d disparity vectors inside current_bbox

    const int MIN_REGION_SIZE = 16;

    // 1.) Is this region too small? Must we stop?
    if ( prod(current_bbox.size()) <= 200 ||
         current_bbox.width() < MIN_REGION_SIZE || current_bbox.height() < MIN_REGION_SIZE ){
      BBox2i expanded = current_bbox;
      expanded.expand(1);
      expanded.crop( bounding_box( disparity ) );
      PixelAccumulator<EWMinMaxAccumulator<Vector2i> > accumulator;
      for_each_pixel( crop(disparity, expanded), accumulator );
      if ( !accumulator.is_valid() ) return true;

      list.push_back( SearchParam( current_bbox,
                                   BBox2i(accumulator.minimum(),
                                          accumulator.maximum() + Vector2i(1,1) ) ) );
      return true;
    }

    // 2) Divide the current_bbox into 4 quadrants, does it reduce total search?
    Vector2i split_pt = current_bbox.size()/2;
    BBox2i q1( current_bbox.min(), current_bbox.min()+split_pt );
    BBox2i q4( current_bbox.min()+split_pt, current_bbox.max() );
    BBox2i q2( current_bbox.min() + Vector2i(split_pt[0],0),
               Vector2i(current_bbox.max()[0],current_bbox.min()[1]+split_pt[1]) );
    BBox2i q3( current_bbox.min() + Vector2i(0,split_pt[1]),
               Vector2i(current_bbox.min()[0]+split_pt[0],current_bbox.max()[1]) );
    BBox2i q1_search, q2_search, q3_search, q4_search;

    // Inside each of the four quadrants, find the min and max disparity.
    // - Masked out pixels are ignored
    // - Accumulate product of disparity search region + pixel area
    // - TODO: Should get some of this logic into class functions.
    int32 split_search = 0;
    { // Q1
      PixelAccumulator<EWMinMaxAccumulator<Vector2i> > accumulator;
      for_each_pixel( crop(disparity,q1), accumulator );
      if ( accumulator.is_valid() ) {
        q1_search = BBox2i(accumulator.minimum(),
                           accumulator.maximum()+Vector2i(1,1));
        split_search += q1_search.area() * prod(q1.size()+kernel_size);
      }
    }
    { // Q2
      PixelAccumulator<EWMinMaxAccumulator<Vector2i> > accumulator;
      for_each_pixel( crop(disparity,q2), accumulator );
      if ( accumulator.is_valid() ) {
        q2_search = BBox2i(accumulator.minimum(),
                           accumulator.maximum()+Vector2i(1,1));
        split_search += q2_search.area() * prod(q2.size()+kernel_size);
      }
    }
    { // Q3
      PixelAccumulator<EWMinMaxAccumulator<Vector2i> > accumulator;
      for_each_pixel( crop(disparity,q3), accumulator );
      if ( accumulator.is_valid() ) {
        q3_search = BBox2i(accumulator.minimum(),
                           accumulator.maximum()+Vector2i(1,1));
        split_search += q3_search.area() * prod(q3.size()+kernel_size);
      }
    }
    { // Q4
      PixelAccumulator<EWMinMaxAccumulator<Vector2i> > accumulator;
      for_each_pixel( crop(disparity,q4), accumulator );
      if ( accumulator.is_valid() ) {
        q4_search = BBox2i(accumulator.minimum(),
                           accumulator.maximum()+Vector2i(1,1));
        split_search += q4_search.area() * prod(q4.size()+kernel_size);
      }
    }
    // Now we have an estimate of the cost of processing these four
    // quadrants seperately

    // 3) Find current search v2
    //    - Get the min and max disparity search range that we just calculated
    //      for the four quadrants.  This is faster than recomputing the min/max.
    BBox2i current_search_region;
    if ( q1_search != BBox2i() )
      current_search_region = q1_search;
    if ( q2_search != BBox2i() && current_search_region == BBox2i() )
      current_search_region = q2_search;
    else
      current_search_region.grow(q2_search);
    if ( q3_search != BBox2i() && current_search_region == BBox2i() )
      current_search_region = q3_search;
    else
      current_search_region.grow(q3_search);
    if ( q4_search != BBox2i() && current_search_region == BBox2i() )
      current_search_region = q4_search;
    else
      current_search_region.grow(q4_search);
    
    int32 current_search = current_search_region.area() * prod(current_bbox.size()+kernel_size);

    const double IMPROVEMENT_RATIO = 0.8;

    if ( split_search > current_search*IMPROVEMENT_RATIO && fail_count == 0 ) {
      // Splitting up the disparity region did not reduce our workload.
      // This is our first failure, so see if we can still improve by
      //  subdividing the quadrants one more time.
      std::vector<SearchParam> failed;
      if (!subdivide_regions( disparity, q1, list, kernel_size, fail_count + 1 ) )
        failed.push_back(SearchParam(q1,q1_search));
      if (!subdivide_regions( disparity, q2, list, kernel_size, fail_count + 1 ) )
        failed.push_back(SearchParam(q2,q2_search));
      if (!subdivide_regions( disparity, q3, list, kernel_size, fail_count + 1 ) )
        failed.push_back(SearchParam(q3,q3_search));
      if (!subdivide_regions( disparity, q4, list, kernel_size, fail_count + 1 ) )
        failed.push_back(SearchParam(q4,q4_search));
              
      if ( failed.size() == 4 ) {
        // All failed, push back this region as a whole (what we started with)
        list.push_back( SearchParam( current_bbox,
                                     current_search_region ) );
        return true;
      } else if ( failed.size() == 3 ) {
        // 3 failed to split can I merge ?
        // - See the failed==2 case for description!
        std::vector<SearchParam>::const_iterator it1 = failed.begin(), it2 = failed.begin();
        ++it2;
        if ( ( it1->first.min().x() == it2->first.min().x() ||
               it1->first.min().y() == it2->first.min().y() ) &&
             it1->second == it2->second ) {
          BBox2i merge = it1->first;
          merge.grow(it2->first);
          list.push_back( SearchParam( merge, it1->second ) );
          list.push_back( *++it2 );
          return true;
        }
        ++it1; ++it2;
        if ( ( it1->first.min().x() == it2->first.min().x() ||
               it1->first.min().y() == it2->first.min().y() ) &&
             it1->second == it2->second ) {
          BBox2i merge = it1->first;
          merge.grow(it2->first);
          list.push_back( SearchParam( merge, it1->second ) );
          list.push_back( failed.front() );
          return true;
        }
        it1 = failed.begin();
        if ( ( it1->first.min().x() == it2->first.min().x() ||
               it1->first.min().y() == it2->first.min().y() ) &&
             it1->second == it2->second ) {
          BBox2i merge = it1->first;
          merge.grow(it2->first);
          list.push_back( SearchParam( merge, it1->second ) );
          list.push_back( *++it1 );
          return true;
        }
        // Push only the bombed regions, possibly a merge step could go here
        list.insert( list.end(), failed.begin(), failed.end() );
      } else if ( failed.size() == 2 ) {
        // 2 failed to split.
        // If the quadrants are adjacent and have the same disparity range,
        //  merge them into a single search region.
        // - TODO: How often does this actually work?
        if ( ( failed.front().first.min().x() == failed.back().first.min().x() ||
               failed.front().first.min().y() == failed.back().first.min().y() ) &&
             failed.front().second == failed.back().second ) {
          BBox2i merge = failed.front().first;
          merge.grow(failed.back().first);
          list.push_back( SearchParam( merge, failed.front().second ) );
          return true;
        }
        list.insert( list.end(), failed.begin(), failed.end() );
      } else if ( failed.size() == 1 ) {
        // Only 1 failed to split, use it in its entirety, allowing
        // us to take advantage of the other regions which split well.
        list.push_back( failed.front() );
      }
      return true;
    } else if ( split_search > current_search*IMPROVEMENT_RATIO && fail_count > 0 ) {
      // Second failure trying to split this region, give up!
      return false;
    } else {
      // Good split, Try to keep splitting each of the four quadrants further.
      subdivide_regions( disparity, q1, list, kernel_size );
      subdivide_regions( disparity, q2, list, kernel_size );
      subdivide_regions( disparity, q3, list, kernel_size );
      subdivide_regions( disparity, q4, list, kernel_size );
    }
    return true;
  }
  
  ImageView<PixelMask<Vector2i>>
  calc_disparity(CostFunctionType cost_type,
                 ImageViewRef<PixelGray<float>> const& left_in,
                 ImageViewRef<PixelGray<float>> const& right_in,
                 // Valid region in the left image
                 BBox2i                 const& left_region,
                 // Max disparity to search in right image
                 Vector2i               const& search_volume,
                 Vector2i               const& kernel_size) {

    // Sanity check the input:
    VW_DEBUG_ASSERT(kernel_size[0] % 2 == 1 && kernel_size[1] % 2 == 1,
                    ArgumentErr() << "calc_disparity: Kernel input not sized with odd values.");
    VW_DEBUG_ASSERT(kernel_size[0] <= left_region.width() &&
                    kernel_size[1] <= left_region.height(),
                    ArgumentErr() << "calc_disparity: Kernel size too large of active region.");
    VW_DEBUG_ASSERT(search_volume[0] > 0 && search_volume[1] > 0,
                    ArgumentErr() << "calc_disparity: Search volume must be greater than 0.");
    VW_DEBUG_ASSERT(left_region.min().x() >= 0 &&  left_region.min().y() >= 0 &&
                    left_region.max().x() <= left_in.impl().cols() &&
                    left_region.max().y() <= left_in.impl().rows(),
                    ArgumentErr() << "calc_disparity: Region not inside left image.");

    typedef PixelGray<float> pix_type; // to save some typing

    // Rasterize input so that we can do a lot of processing on it.
    BBox2i right_region = left_region;
    right_region.max() += search_volume - Vector2i(1,1);
    ImageView<pix_type> left (crop(left_in.impl(),  left_region));
    ImageView<pix_type> right(crop(right_in.impl(), right_region));
    
    // Call the lower level function with the appropriate cost function type
    switch (cost_type) {
    case CROSS_CORRELATION:
      return best_of_search_convolution<NCCCost<ImageView<pix_type>>, pix_type>
        (left, right, left_region, search_volume, kernel_size);
    case SQUARED_DIFFERENCE:
      return best_of_search_convolution<SquaredCost<ImageView<pix_type>>, pix_type>
        (left, right, left_region, search_volume, kernel_size);
    default: // case ABSOLUTE_DIFFERENCE:
      return best_of_search_convolution<AbsoluteCost<ImageView<pix_type>>, pix_type>
        (left, right, left_region, search_volume, kernel_size);
    }

    return ImageView<PixelMask<Vector2i>>(); // will not be reached
    
  } // End function calc_disparity

  /// Create fake left and right images and search volume.  Do a fake
  /// disparity calculation. Divide the run-time of this calculation
  /// by left region size times search box size. This will enable us
  /// to estimate how long disparity calculation takes for given cost
  /// function and kernel size.
  double calc_seconds_per_op(CostFunctionType cost_type, Vector2i const& kernel_size){
     
    double elapsed = -1.0;
    double seconds_per_op = -1.0;

    // We don't know what sizes to use to get a reliable time estimate.
    // So increase the size until the time estimate is a second.
    int lsize = 100;
    while (elapsed < 1.0){

      // Below we add kernel_size to ensure the image exceeds the
      // kernel size, for correlation to perform properly.
      lsize = (int)ceil(lsize*1.2) + std::max(kernel_size[0], kernel_size[1]);

      ImageView<PixelGray<float>> fake_left(lsize, lsize);
      for (int col = 0; col < fake_left.cols(); col++){
        for (int row = 0; row < fake_left.rows(); row++){
          fake_left(col, row) = col%2 + 2*(row%5); // some values
        }
      }

      ImageView<PixelGray<float>> fake_right(4*lsize, 4*lsize);
      for (int col = 0; col < fake_right.cols(); col++){
        for (int row = 0; row < fake_right.rows(); row++){
          fake_right(col, row) = 3*(col%7) + row%3; // some values
        }
      }

      BBox2i search_region(0, 0, lsize/5, lsize/5);
      BBox2i left_region = bounding_box(fake_left);

      Stopwatch watch;
      watch.start();
      ImageView<PixelMask<Vector2i>> disparity =
        calc_disparity(cost_type, fake_left, fake_right,
                       left_region, search_region.size(), kernel_size);
      watch.stop();

      // Note: We add an infinitesimal contribution of disparity, lest
      // the compiler tries to optimize away the above calculation due
      // to its result being unused.
      elapsed = watch.elapsed_seconds() + 1e-40*disparity(0, 0).child().x();
      SearchParam params(left_region, search_region);
      seconds_per_op = elapsed/params.search_volume();
    }

    return seconds_per_op;
  }
  
}} // end namespace vw::stereo
