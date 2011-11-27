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
    std::fill( disparity_map.data(), disparity_map.data() + prod(result_size),
               PixelMask<Vector2i>(Vector2i()) );
    // First channel is best, second is worst.
    ImageView<QualT> quality_map( result_size[0], result_size[1] );
    ImageView<AccumT> cost_metric( result_size[0], result_size[1] );

    // Convolve across search volume
    Vector2i disparity(0,0);
    for ( ; disparity.y() != search_volume[1]; ++disparity.y() ) {
      for ( disparity.x() = 0; disparity.x() != search_volume[0]; ++disparity.x() ) {
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
        const AccumT* cost_ptr     = cost_metric.data();
        const AccumT* cost_ptr_end = cost_metric.data() + prod(result_size);
        QualT* quality_ptr         = quality_map.data();
        PixelMask<Vector2i>* disparity_ptr = disparity_map.data();
        if ( disparity != Vector2i(0,0) ) {
          // Normal comparison operations
          while ( cost_ptr != cost_ptr_end ) {
            if ( cost_function.quality_comparison( *cost_ptr, quality_ptr->first ) ) {
              // Better than best?
              quality_ptr->first = *cost_ptr;
              disparity_ptr->child() = disparity;
            } else if ( !cost_function.quality_comparison( *cost_ptr, quality_ptr->second ) ) {
              // Worse than worse
              quality_ptr->second = *cost_ptr;
            }
            ++cost_ptr;
            ++quality_ptr;
            ++disparity_ptr;
          }
        } else {
          // Initializing quality_map and disparity_map with first result
          while ( cost_ptr != cost_ptr_end ) {
            quality_ptr->first = quality_ptr->second = *cost_ptr;
            ++cost_ptr;
            ++quality_ptr;
          }
        }
      }
    }

    // Determine validity of result
    const QualT* quality_ptr      = quality_map.data();
    const QualT* quality_ptr_end  = quality_map.data() + prod(result_size);
    PixelMask<Vector2i>* disp_ptr = disparity_map.data();
    while ( quality_ptr != quality_ptr_end ) {
      if ( quality_ptr->first == quality_ptr->second )
        invalidate( *disp_ptr );
      ++quality_ptr;
      ++disp_ptr;
    }

    return disparity_map;
  }

  // The following are tools for subdividing a disparity map into
  // regions that contain near equal disparity.
  //
  // These are useful for pyramid stereo correlation and several other
  // algorithms. The prime result is "SearchParam" where the first
  // bbox in the area of the result. The second bbox is the disparity
  // in the region defined by the first bbox.
  typedef std::pair<BBox2i,BBox2i> SearchParam;

  // Developer tools for modifying bboxes.
  inline int32 area( BBox2i const& a );
  inline void expand_bbox( BBox2i& a, BBox2i const& b );

  // This is a recursive function. It might be ideal to make this a
  // template. However in all cases so far, I've only applied to
  // PixelMask<Vector2i>. This should stay an image view as we'll be
  // accessing the image alot and randomly.
  bool subdivide_regions( ImageView<PixelMask<Vector2i> > const& disparity,
                          BBox2i const& current_bbox,
                          std::list<SearchParam>& list,
                          Vector2i const& kernel_size,
                          int32 fail_count = 0 );

}}}

#endif//__VW_STEREO_REWRITE_CORRELATION_H__
