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

#include <vw/Image/ImageView.h>
#include <vw/Stereo/DisparityMap.h>
#include <vw/Stereo/Correlate.h>
#include <vw/Stereo/PreFilter.h>
#include <vw/Stereo/CostFunctions.h>
#include <vw/Stereo/Correlation.h>
#include <vw/Stereo/ParabolaSubpixelView.h>

#include <boost/foreach.hpp>

// TODO(oalexan1): The integer disparity is in fact float. Convert all to float.
namespace vw { namespace stereo {

  ImageView<PixelMask<Vector2f>> ParabolaSubpixelView::evaluate
  (ImageView<PixelMask<Vector2i>> const& integer_disparity,    // Cropped input disp
   ImageViewRef<PixelGray<float>> const& left_filtered_image,  // uncropped left
   ImageViewRef<PixelGray<float>> const& right_filtered_image, // uncropped right
   BBox2i const& left_region,      ///< The ROI of the left and right images
   BBox2i const& right_region,     ///  to use in order to do our computations.
   BBox2i const& disparity_region, ///< The ROI in the entire output image
   BBox2i const& search_range ) const { /// Range of input disparity values in region
    
    // Rasterize the ROIs of the left and right input images.
    ImageView<PixelGray<float>> left_raster  = crop(left_filtered_image, left_region );
    ImageView<PixelGray<float>> right_raster = crop(right_filtered_image,right_region);

    // Allocate a buffer to store the costs of the 9 nearest disparities for each 
    //  pixel in the integer disparity image.
    // - This will use 2.25 MB for a 256^2 pixel region
    ImageView<Vector<float,9>> cost_patch(integer_disparity.cols(), integer_disparity.rows());
      
    // TODO(oalexan1): Why is this hard-coded to a cost function that we did not use
    //       when we originally computed the integer correlation????
    typedef AbsoluteCost<ImageView<PixelGray<float>>> CostType;
        
    typedef typename CostType::accumulator_type                      AccumChannelT;
    typedef typename PixelChannelCast<PixelGray<float>, AccumChannelT>::type AccumT;
    typedef typename ImageView<Vector<float,9>>::pixel_accessor     PatchAcc;
    typedef typename ImageView<AccumT>::pixel_accessor               MetricAcc;
    typedef typename ImageView<PixelMask<Vector2i> >::pixel_accessor IDispAcc;

    // Subdivide the input disparity into smaller regions that have
    // the same disparity. The goal here is that in each region the disparity
    // varies only little, then the cost functions at a lot of points in the
    // region will have identical values, so we will pre-compute each
    // distinct cost function and replicate it across all points in the regions.
    // If however in each region the disparity varies a lot, this approach
    // will slow things down rather than speed things up, then we
    // do the cost function computation for each individual pixel.

    // If the ratio of range of disparities in each region to the area
    // of the region is bigger than this ratio, handle each point
    // in the region individually, per above.

    double ratio = 1.0; // The exact value here does not affect things much
    std::vector<SearchParam> big_zones, zones;
    subdivide_regions( integer_disparity,
                       bounding_box(integer_disparity),
                       big_zones, m_kernel_size );
    BOOST_FOREACH( SearchParam& zone, big_zones ) {

      double len1 = zone.image_region().area();
      double len2 = zone.disparity_range().area();
      if (len2/len1 < ratio){
        zones.push_back(zone); 
      }else{

        for (int32 dx = zone.image_region().min().x(); dx < zone.image_region().max().x(); ++dx){
          for (int32 dy = zone.image_region().min().y(); dy < zone.image_region().max().y();
               ++dy){
            BBox2i box(dx, dy, 1, 1);
            PixelAccumulator<EWMinMaxAccumulator<Vector2i> > accumulator;
            for_each_pixel( crop(integer_disparity, box), accumulator );
            if ( !accumulator.is_valid() ) 
              continue;
            zones.push_back
              ( SearchParam( box,
                             BBox2i(accumulator.minimum(),
                                    accumulator.maximum() + Vector2i(1,1) ) ) );
              
          }
        } // End loop through pixels in the zone
      }
    } // End loop through initial zones
      
      // Loop through the final zones
    BOOST_FOREACH( SearchParam& zone, zones ) {
        
      zone.disparity_range().expand(1);
        
      ImageView<AccumT> cost_metric( zone.image_region().width(), zone.image_region().height() );

      // Note: The zone locations are relative to the cropped output disparity region
      BBox2i left_zone = zone.image_region();
      left_zone.max() += m_kernel_size - Vector2i(1,1);

      // Loop through all the disparity values for this zone
      for ( int32 dx = 0; dx < zone.disparity_range().width(); ++dx ) {
        for ( int32 dy = 0; dy < zone.disparity_range().height(); ++dy ) {

          Vector2i disparity( dx, dy ); // Relative disparity number.

          Vector2i disparity_abs = disparity + zone.disparity_range().min(); // Absolute disparity

          // Compute the score for this disparity at each pixel in the current zone
          ImageView<PixelGray<float>> left_raster_crop = 
            crop(left_raster,left_zone);
          ImageView<PixelGray<float>> right_raster_crop = 
            crop(right_raster,left_zone+disparity_abs-search_range.min());
        
          CostType cost_function(left_raster_crop, right_raster_crop, m_kernel_size);    
          cost_metric = fast_box_sum<AccumChannelT>(cost_function(left_raster_crop,
                                                                  right_raster_crop),
                                                    m_kernel_size );
          //cost_function.cost_modification( cost_metric, disparity );
          cost_function.cost_modification( cost_metric, Vector2i(0,0) );

          VW_DEBUG_ASSERT( (zone.image_region().width()  == cost_metric.cols()) &&
                           (zone.image_region().height() == cost_metric.rows()),
                           MathErr() << "Cost Metric seems to have wrong size" );

          // Now take the fast evaluate cost function and see where
          // the cost patch they apply
          // ???
            
          // What we are doing is building up a list of 9 costs for each pixel representing
          // the 3x3 matrix of disparity score surrounding the best disparity score.
            
          PatchAcc  patch_row  = cost_patch.origin();
          MetricAcc metric_row = cost_metric.origin();
          IDispAcc  idisp_row  = integer_disparity.origin();
          patch_row.advance( zone.image_region().min().x(), zone.image_region().min().y() );
          idisp_row.advance( zone.image_region().min().x(), zone.image_region().min().y() );
            
          // Loop through rows of the zone
          for ( int32 j = zone.image_region().height(); j; --j ) {
            PatchAcc  patch_col  = patch_row;
            MetricAcc metric_col = metric_row;
            IDispAcc  idisp_col  = idisp_row;
            // Loop through columns of the zone
            for ( int32 i = zone.image_region().width(); i; --i ) {
              Vector2i delta = disparity_abs - (*idisp_col).child();
                
              // Given the current disparity and location within the zone, assign the 
              // correlation score to the appropriate location within our storage buffer.
              // - The storage indices are:  [0] [1] [2]
              //                             [3] [4] [5]
              //                             [6] [7] [8]
              /*
                if ( delta == Vector2i(-1,-1) ) {
                (*patch_col)[0] = *metric_col;
                } else if ( delta == Vector2i( 0, -1 ) ) {
                (*patch_col)[3] = *metric_col;
                } else if ( delta == Vector2i( 1, -1 ) ) {
                (*patch_col)[6] = *metric_col;
                } else if ( delta == Vector2i( -1, 0 ) ) {
                (*patch_col)[1] = *metric_col;
                } else if ( delta == Vector2i( 0, 0 ) ) {
                (*patch_col)[4] = *metric_col;
                } else if ( delta == Vector2i( 1, 0 ) ) {
                (*patch_col)[7] = *metric_col;
                } else if ( delta == Vector2i( -1, 1 ) ) {
                (*patch_col)[2] = *metric_col;
                } else if ( delta == Vector2i( 0, 1 ) ) {
                (*patch_col)[5] = *metric_col;
                } else if ( delta == Vector2i( 1, 1 ) ) {
                (*patch_col)[8] = *metric_col;
                }
              */                
              if ( delta == Vector2i(-1,-1) ) {
                (*patch_col)[0] = *metric_col;
              } else if ( delta == Vector2i( 0, -1 ) ) {
                (*patch_col)[1] = *metric_col;
              } else if ( delta == Vector2i( 1, -1 ) ) {
                (*patch_col)[2] = *metric_col;
              } else if ( delta == Vector2i( -1, 0 ) ) {
                (*patch_col)[3] = *metric_col;
              } else if ( delta == Vector2i( 0, 0 ) ) {
                (*patch_col)[4] = *metric_col;
              } else if ( delta == Vector2i( 1, 0 ) ) {
                (*patch_col)[5] = *metric_col;
              } else if ( delta == Vector2i( -1, 1 ) ) {
                (*patch_col)[6] = *metric_col;
              } else if ( delta == Vector2i( 0, 1 ) ) {
                (*patch_col)[7] = *metric_col;
              } else if ( delta == Vector2i( 1, 1 ) ) {
                (*patch_col)[8] = *metric_col;
              }
              patch_col.next_col();
              metric_col.next_col();
              idisp_col.next_col();
            }
            patch_row.next_row();
            metric_row.next_row();
            idisp_row.next_row();
          } // End loop through zone pixels
          
        }
      } // end loop through disparities
        
    } // end loop through zones

      // Calculate the new floating point location
    ImageView<PixelMask<Vector2f> > result( disparity_region.width(),
                                            disparity_region.height() );
    typedef typename ImageView<PixelMask<Vector2f> >::pixel_accessor ResultAcc;
    ResultAcc result_row = result.origin();
    IDispAcc  idisp_row  = integer_disparity.origin();
    PatchAcc  patch_row  = cost_patch.origin();
    for ( int32 j = disparity_region.height(); j; --j ) {
      ResultAcc result_col = result_row;
      IDispAcc  idisp_col  = idisp_row;
      PatchAcc  patch_col  = patch_row;
      for ( int32 i = disparity_region.width(); i; --i ) {
        // Given our 9 points of cost around our disparity pixel
        // (patch_col), find the 2D minimum and that will be our
        // floating point update.
        if ( is_valid( *idisp_col ) ) {
          if ( !std::equal( (*patch_col).begin()+1, (*patch_col).end(),
                            (*patch_col).begin() ) ) {
            // First, compute the parameters of the hyperbolic surface by
            // fitting the nine points in 'points' using a linear least
            // squares fit.  This process is fairly fast, since we have
            // already pre-computed the inverse of the A matrix in Ax = b.
            Vector<float,6> x( m_p_A_matrix * (*patch_col) );
            // With these parameters, we have a closed form expression for
            // the surface.  We compute the derivative, and find the point
            // where the slope is zero.  This is our maximum.
            //
            //  Max is at [x,y] where:
            //   dz/dx = 2ax + cy + d = 0
            //   dz/dy = 2by + cx + e = 0
            float denom = 4 * x[0] * x[1] - ( x[2] * x[2] ); // = 4ab - c^2
            Vector2f offset( ( x[2] * x[4] - 2 * x[1] * x[3] ) / denom,
                             ( x[2] * x[3] - 2 * x[0] * x[4] ) / denom );
            const float MAX_SUBPIXEL_SHIFT = 5.0; // TODO: Experiment with this!
            if ( norm_2(offset) < MAX_SUBPIXEL_SHIFT )
              *result_col = PixelMask<Vector2f>( remove_mask(*idisp_col) + offset );
            else
              *result_col = PixelMask<Vector2f>(*idisp_col);
          } else {
            *result_col = PixelMask<Vector2f>(*idisp_col);
          }
        } else {
          *result_col = PixelMask<Vector2f>();
        }
        result_col.next_col();
        idisp_col.next_col();
        patch_col.next_col();
      }
      result_row.next_row();
      idisp_row.next_row();
      patch_row.next_row();
    }

    return result;
  } // End of the evaluate() function

    // Block rasterization section does the actual work
  ParabolaSubpixelView::prerasterize_type
  ParabolaSubpixelView::prerasterize(BBox2i const& bbox) const {

    // Prerasterize the disparity for this crop and work out it's
    // range. This could be faster if we subdivided this down by the
    // areas with equal values.
    ImageView<PixelMask<Vector2i>> disparity_subsection = crop(m_disparity,bbox);

    // We calculate the entire search range only so that we can
    // prerasterize the sections we'll need for left and right images.
    BBox2i entire_search_range = stereo::get_disparity_range(disparity_subsection);
    entire_search_range.max() += Vector2i(1,1);
    entire_search_range.expand(1); // Because we are going to check
    // the neighboring pixels

    // Calculate our left an right regions that we need
    Vector2i half_kernel = m_kernel_size/2;
    BBox2i left_region = bbox;
    left_region.min() -= half_kernel;
    left_region.max() += half_kernel;
    BBox2i right_region = left_region + entire_search_range.min();
    right_region.max() += entire_search_range.size();

    // Unfortunately using the more convenient prefilter_image
    // function yields slightly different results!
    if (m_prefilter_mode == PREFILTER_LOG) {
      stereo::LaplacianOfGaussian prefilter(m_prefilter_width);
      return prerasterize_type(evaluate(disparity_subsection,
                                        prefilter.filter(m_left_image ),
                                        prefilter.filter(m_right_image),
                                        left_region, right_region, bbox,
                                        entire_search_range ),
                               -bbox.min().x(), -bbox.min().y(), cols(), rows());
    } else {
      if (m_prefilter_mode == PREFILTER_MEANSUB) {
        stereo::SubtractedMean prefilter(m_prefilter_width);
        return prerasterize_type(evaluate(disparity_subsection,
                                          prefilter.filter(m_left_image ),
                                          prefilter.filter(m_right_image),
                                          left_region, right_region, bbox,
                                          entire_search_range ),
                                 -bbox.min().x(), -bbox.min().y(), cols(), rows());
      } else { // PREFILTER_NONE
        stereo::NullOperation prefilter;
        return prerasterize_type(evaluate(disparity_subsection,
                                          prefilter.filter(m_left_image ),
                                          prefilter.filter(m_right_image),
                                          left_region, right_region, bbox,
                                          entire_search_range ),
                                 -bbox.min().x(), -bbox.min().y(), cols(), rows());
      }
    }

  } // end function prerasterize

  ParabolaSubpixelView
  parabola_subpixel(ImageViewRef<PixelMask<Vector2f>> const& disparity,
                    ImageViewRef<PixelGray<float>> const& left_image,
                    ImageViewRef<PixelGray<float>> const& right_image,
                    PrefilterModeType prefilter_mode, float prefilter_width,
                    Vector2i const& kernel_size ) {
    typedef ParabolaSubpixelView result_type;
    return result_type( disparity.impl(), left_image.impl(), right_image.impl(),
                        prefilter_mode, prefilter_width, kernel_size );
  }

}} // namespace vw::stereo

