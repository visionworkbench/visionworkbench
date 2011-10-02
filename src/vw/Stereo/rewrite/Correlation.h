// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#ifndef __VW_STEREO_REWRITE_CORRELATION_H__
#define __VW_STEREO_REWRITE_CORRELATION_H__

#include <vw/Image/PixelMask.h>
#include <vw/Image/ImageMath.h>
#include <vw/Stereo/Correlate.h>
#include <vw/Stereo/rewrite/Algorithms.h>
#include <vw/Stereo/rewrite/CostFunctions.h>

namespace vw {
namespace stereo {
namespace rewrite {

  // This actually RASTERIZES/COPY the input images. It then makes an
  // allocation to store current costs.
  //
  // Users pass us the active region of the left image. Hopefully this
  // allows them to consider if they need edge extension.
  //
  // The return size of this function will be:
  //     return size = left_region_size - kernel_size + 1.
  //
  // This means the user must take in account the kernel size for
  // deciding the region size.
  //
  // The size of the area we're going to access in the right is
  // calculated as follows:
  //     right_region = left_region + search_volume - 1.
  //
  template <template<class,bool> class CostFuncT, class ImageT1, class ImageT2>
  ImageView<PixelMask<Vector2i> >
  best_of_search_convolution( ImageViewBase<ImageT1> const& left,
                              ImageViewBase<ImageT2> const& right,
                              BBox2i const& left_region,
                              Vector2i const& search_volume,
                              Vector2i const& kernel_size ) {
    // Sanity check the input:
    VW_DEBUG_ASSERT( kernel_size[0] % 2 == 1 && kernel_size[1] % 2 == 1,
                     ArgumentErr() << "best_of_search_convolution: Kernel input not sized with odd values." );
    VW_DEBUG_ASSERT( kernel_size[0] <= left_region.width() &&
                     kernel_size[1] <= left_region.height(),
                     ArgumentErr() << "best_of_search_convolution: Kernel size too large of active region." );
    VW_DEBUG_ASSERT( search_volume[0] > 0 && search_volume[1] > 0,
                     ArgumentErr() << "best_of_search_convolution: Search volume must be greater than 0" );
    VW_DEBUG_ASSERT( left_region.min().x() >= 0 &&  left_region.min().y() >= 0 &&
                     left_region.max().x() <= left.impl().cols() &&
                     left_region.max().y() <= left.impl().rows(),
                     ArgumentErr() << "best_of_search_convolution: Region not inside left image." );

    typedef typename ImageT1::pixel_type PixelT1;
    typedef typename ImageT2::pixel_type PixelT2;
    typedef typename CostFuncT<ImageT1,boost::is_integral<typename PixelChannelType<PixelT1>::type>::value>::accumulator_type AccumChannelT;
    typedef typename PixelChannelCast<PixelT1,AccumChannelT>::type AccumT;
    typedef typename std::pair<AccumT,AccumT> QualT;

    // Now then .. actually rasterize input so that we are not
    // repeatedly applying preprocess filters.
    BBox2i right_region = left_region;
    right_region.max() += search_volume - Vector2i(1,1);
    ImageView<PixelT1> left_raster( crop(left.impl(), left_region) );
    ImageView<PixelT2> right_raster( crop(right.impl(), right_region) );

    // Build cost function which sometimes has side car data
    CostFuncT<ImageT1,boost::is_integral<typename PixelChannelType<PixelT1>::type>::value> cost_function( left_raster, right_raster, kernel_size);

    // Result buffers
    Vector2i result_size = left_region.size() - kernel_size + Vector2i(1,1);
    ImageView<PixelMask<Vector2i> > disparity_map( result_size[0],
                                                   result_size[1] );
    // First channel is best, second is worst.
    ImageView<QualT> quality_map( result_size[0], result_size[1] );
    ImageView<AccumT> cost_metric( result_size[0], result_size[1] );

    // Convolve across search volume
    for ( int32 dx = 0; dx < search_volume[0]; ++dx ) {
      for ( int32 dy = 0; dy < search_volume[1]; ++dy ) {
        Vector2i disparity(dx,dy);

        // There's only one raster here. Fast box sum calls each pixel
        // individually by pixel accessor. It only calls each pixel
        // once so there's no reason to copy/rasterize the cost result
        // before hand.
        //
        // The cost function should also no be applying an edge
        // extension as we've already over cropped the input.
        cost_metric =
          fast_box_sum<AccumChannelT>(cost_function( left_raster,
                                                     crop(right_raster,
                                                          bounding_box(left_raster)+disparity) ),
                                      kernel_size );
        cost_function.cost_modification( cost_metric, disparity );

        // These conditionals might be served outside of the iteration
        // of dx and dy. It would make the code slightly longer but
        // would avoid a conditional inside a double loop.
        AccumT* cost_ptr     = &cost_metric(0,0);
        AccumT* cost_ptr_end = &cost_metric(cost_metric.cols()-1,
                                            cost_metric.rows()-1)+1;
        QualT* quality_ptr   = &quality_map(0,0);
        PixelMask<Vector2i>* disparity_ptr = &disparity_map(0,0);
        if ( disparity != Vector2i() ) {
          // Normal comparison operations
          while ( cost_ptr != cost_ptr_end ) {
            if ( cost_function.quality_comparison( *cost_ptr, quality_ptr->first ) ) {
              // Better than best?
              quality_ptr->first = *cost_ptr;
              *disparity_ptr = PixelMask<Vector2i>(disparity);
            } else if ( !cost_function.quality_comparison( *cost_ptr, quality_ptr->second ) ) {
              // Worse than worse
              quality_ptr->second = *cost_ptr;
            }
            cost_ptr++;
            quality_ptr++;
            disparity_ptr++;
          }
        } else {
          // Initializing quality_map and disparity_map with first result
          while ( cost_ptr != cost_ptr_end ) {
            *disparity_ptr = PixelMask<Vector2i>(disparity);
            quality_ptr->first = quality_ptr->second = *cost_ptr;
            cost_ptr++;
            quality_ptr++;
            disparity_ptr++;
          }
        }
      }
    }

    // Determine validity of result
    QualT* quality_ptr            = &quality_map(0,0);
    QualT* quality_ptr_end        = &quality_map(quality_map.cols()-1,
                                                 quality_map.rows()-1)+1;
    PixelMask<Vector2i>* disp_ptr = &disparity_map(0,0);
    while ( quality_ptr != quality_ptr_end ) {
      if ( quality_ptr->first == quality_ptr->second )
        invalidate( *disp_ptr );
      quality_ptr++;
      disp_ptr++;
    }

    return disparity_map;
  }

  typedef std::pair<BBox2i,BBox2i> SearchParam;

  inline int32 area( BBox2i const& a ) {
    int32 width = a.width();
    int32 heigh = a.height();
    if ( width < 0 || heigh < 0 )
      return 0;
    return width * heigh;
  }

  inline void expand_bbox( BBox2i& a, BBox2i const& b ) {
    a.min().x() = std::min(a.min().x(),b.min().x());
    a.min().y() = std::min(a.min().y(),b.min().y());
    a.max().x() = std::max(a.max().x(),b.max().x());
    a.max().y() = std::max(a.max().y(),b.max().y());
  }

  bool subdivide_regions( ImageView<PixelMask<Vector2i> > const& disparity,
                          BBox2i const& current_bbox,
                          std::list<SearchParam>& list,
                          Vector2i const& kernel_size,
                          int32 fail_count = 0 ) {

    // 1.) Is this region too small? Must we stop?
    if ( prod(current_bbox.size()) <= 200 ||
         current_bbox.width() < 16 || current_bbox.height() < 16 ){
      PixelAccumulator<EWMinMaxAccumulator<Vector2i> > accumulator;
      for_each_pixel( crop(disparity,current_bbox), accumulator );
      if ( !accumulator.is_valid() ) return true;

      list.push_back( SearchParam( current_bbox,
                                   BBox2i(accumulator.minimum(),
                                          accumulator.maximum() + Vector2i(1,1) ) ) );
      return true;
    }

    // 2) Divide the section into 4 quadrants, does it reduce total search?
    Vector2i split_pt = current_bbox.size()/2;
    BBox2i q1( current_bbox.min(), current_bbox.min()+split_pt );
    BBox2i q4( current_bbox.min()+split_pt, current_bbox.max() );
    BBox2i q2( current_bbox.min() + Vector2i(split_pt[0],0),
               Vector2i(current_bbox.max()[0],current_bbox.min()[1]+split_pt[1]) );
    BBox2i q3( current_bbox.min() + Vector2i(0,split_pt[1]),
               Vector2i(current_bbox.min()[0]+split_pt[0],current_bbox.max()[1]) );
    BBox2i q1_search, q2_search, q3_search, q4_search;

    int32 split_search = 0;
    { // Q1
      PixelAccumulator<EWMinMaxAccumulator<Vector2i> > accumulator;
      for_each_pixel( crop(disparity,q1), accumulator );
      if ( accumulator.is_valid() ) {
        q1_search = BBox2i(accumulator.minimum(),
                           accumulator.maximum()+Vector2i(1,1));
        split_search += area(q1_search) * prod(q1.size()+kernel_size);
      }
    }
    { // Q2
      PixelAccumulator<EWMinMaxAccumulator<Vector2i> > accumulator;
      for_each_pixel( crop(disparity,q2), accumulator );
      if ( accumulator.is_valid() ) {
        q2_search = BBox2i(accumulator.minimum(),
                           accumulator.maximum()+Vector2i(1,1));
        split_search += area(q2_search) * prod(q2.size()+kernel_size);
      }
    }
    { // Q3
      PixelAccumulator<EWMinMaxAccumulator<Vector2i> > accumulator;
      for_each_pixel( crop(disparity,q3), accumulator );
      if ( accumulator.is_valid() ) {
        q3_search = BBox2i(accumulator.minimum(),
                           accumulator.maximum()+Vector2i(1,1));
        split_search += area(q3_search) * prod(q3.size()+kernel_size);
      }
    }
    { // Q4
      PixelAccumulator<EWMinMaxAccumulator<Vector2i> > accumulator;
      for_each_pixel( crop(disparity,q4), accumulator );
      if ( accumulator.is_valid() ) {
        q4_search = BBox2i(accumulator.minimum(),
                           accumulator.maximum()+Vector2i(1,1));
        split_search += area(q4_search) * prod(q4.size()+kernel_size);
      }
    }

    // 3) Find current search v2
    BBox2i current_search_region;
    if ( q1_search != BBox2i() )
      current_search_region = q1_search;
    if ( q2_search != BBox2i() && current_search_region == BBox2i() )
      current_search_region = q2_search;
    else
      expand_bbox( current_search_region, q2_search );
    if ( q3_search != BBox2i() && current_search_region == BBox2i() )
      current_search_region = q3_search;
    else
      expand_bbox( current_search_region, q3_search );
    if ( q4_search != BBox2i() && current_search_region == BBox2i() )
      current_search_region = q4_search;
    else
      expand_bbox( current_search_region, q4_search );
    int32 current_search = area(current_search_region) * prod(current_bbox.size()+kernel_size);

    if ( split_search > current_search*0.9 && fail_count == 0 ) {
      // Did bad .. maybe next level will have better luck?
      std::list<SearchParam> failed;
      if (!subdivide_regions( disparity, q1, list, kernel_size, fail_count + 1 ) )
        failed.push_back(SearchParam(q1,q1_search));
      if (!subdivide_regions( disparity, q2, list, kernel_size, fail_count + 1 ) )
        failed.push_back(SearchParam(q2,q2_search));
      if (!subdivide_regions( disparity, q3, list, kernel_size, fail_count + 1 ) )
        failed.push_back(SearchParam(q3,q3_search));
      if (!subdivide_regions( disparity, q4, list, kernel_size, fail_count + 1 ) )
        failed.push_back(SearchParam(q4,q4_search));
      if ( failed.size() == 4 ) {
        // all failed, push back this region as a whole
        list.push_back( SearchParam( current_bbox,
                                     current_search_region ) );
        return true;
      } else if ( failed.size() == 3 ) {
        // 3 failed to split can I merge ?
        std::list<SearchParam>::const_iterator it1 = failed.begin(), it2 = failed.begin();
        it2++;
        if ( ( it1->first.min().x() == it2->first.min().x() ||
               it1->first.min().y() == it2->first.min().y() ) &&
             it1->second == it2->second ) {
          BBox2i merge = it1->first;
          expand_bbox( merge, it2->first );
          list.push_back( SearchParam( merge, it1->second ) );
          it2++;
          list.push_back( *it2 );
          return true;
        }
        it1++; it2++;
        if ( ( it1->first.min().x() == it2->first.min().x() ||
               it1->first.min().y() == it2->first.min().y() ) &&
             it1->second == it2->second ) {
          BBox2i merge = it1->first;
          expand_bbox( merge, it2->first );
          list.push_back( SearchParam( merge, it1->second ) );
          list.push_back( failed.front() );
          return true;
        }
        it1 = failed.begin();
        if ( ( it1->first.min().x() == it2->first.min().x() ||
               it1->first.min().y() == it2->first.min().y() ) &&
             it1->second == it2->second ) {
          BBox2i merge = it1->first;
          expand_bbox( merge, it2->first );
          list.push_back( SearchParam( merge, it1->second ) );
          it1++;
          list.push_back( *it1 );
          return true;
        }
        // Push only the bombed regions, possibly a merge step could go here
        list.insert( list.end(), failed.begin(), failed.end() );
      } else if ( failed.size() == 2 ) {
        // 2 failed to split .. can I merge?
        if ( ( failed.front().first.min().x() == failed.back().first.min().x() ||
               failed.front().first.min().y() == failed.back().first.min().y() ) &&
             failed.front().second == failed.back().second ) {
          BBox2i merge = failed.front().first;
          expand_bbox( merge, failed.back().first );
          list.push_back( SearchParam( merge, failed.front().second ) );
          return true;
        }
        list.insert( list.end(), failed.begin(), failed.end() );
      } else if ( failed.size() == 1 ) {
        // Only 1 failed to split .. push it back
        list.push_back( failed.front() );
      }
      return true;
    } else if ( split_search > current_search*0.9 && fail_count > 0 ) {
      // Bad again .. back up
      return false;
    } else {
      // Good split
      subdivide_regions( disparity, q1, list, kernel_size );
      subdivide_regions( disparity, q2, list, kernel_size );
      subdivide_regions( disparity, q3, list, kernel_size );
      subdivide_regions( disparity, q4, list, kernel_size );
    }
    return true;
  }

}}}

#endif//__VW_STEREO_REWRITE_CORRELATION_H__
