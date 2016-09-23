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


#ifndef __VW_STEREO_SUBPIXEL_VIEW__
#define __VW_STEREO_SUBPIXEL_VIEW__

#include <vw/Image.h>
#include <vw/Stereo/DisparityMap.h>
#include <vw/Stereo/Correlate.h>
#include <vw/Stereo/PreFilter.h>
#include <vw/Stereo/CostFunctions.h>
#include <vw/Stereo/Correlation.h>

#include <boost/foreach.hpp>

namespace vw {
namespace stereo {

  template <class DImageT, class Image1T, class Image2T>
  class ParabolaSubpixelView : public ImageViewBase<ParabolaSubpixelView<DImageT,Image1T,Image2T> > {
    DImageT    m_disparity;
    Image1T    m_left_image;
    Image2T    m_right_image;
    Vector2i   m_kernel_size;

    // These two variables pick a prefilter which is applied to each pyramid level
    PrefilterModeType m_prefilter_mode; ///< See Prefilter.h for the types
    float m_prefilter_width;     ///< Preprocessing filter width

    Matrix<float,6,9> m_p_A_matrix;

    /// Compute the subpixel disparity for each input integer disparity
    /// - This function is written in a roundabout way to maximize the benefit from our
    ///   fast_box_sum() function.  Of course, we already performed all of these computations
    ///   back when we generated the integer disparity!!!
    template <class FImage1T, class FImage2T>
    ImageView<PixelMask<Vector2f> >
    evaluate( ImageView<PixelMask<Vector2i> > const& integer_disparity, ///< Input disparity, cropped to disparity_region
              ImageViewBase<FImage1T>         const& left_filtered_image,
              ImageViewBase<FImage2T>         const& right_filtered_image,
              BBox2i const& left_region,      ///< The ROI of the left and right images that we will need
              BBox2i const& right_region,     ///  to use in order to do our computations.
              BBox2i const& disparity_region, ///< The ROI in the entire output image that we are computing
              BBox2i const& search_range ) const { ///< The range of input disparity values in disparity_region
      
      // TODO: Try moving this step outside the function to eliminate the template parameters!        
      // Rasterize the ROIs of the left and right input images.
      ImageView<typename FImage1T::pixel_type> left_raster  = crop(left_filtered_image, left_region );
      ImageView<typename FImage2T::pixel_type> right_raster = crop(right_filtered_image,right_region);

      // Allocate a buffer to store the costs of the 9 nearest disparities for each 
      //  pixel in the integer disparity image.
      // - This will use 2.25 MB for a 256^2 pixel region
      ImageView<Vector<float,9> > cost_patch( integer_disparity.cols(),
                                              integer_disparity.rows() );
      
      // TODO: Why is this hard-coded to a cost function that we did not use
      //       when we originally computed the integer correlation????
      typedef AbsoluteCost<FImage1T,boost::is_integral<typename PixelChannelType<typename FImage1T::pixel_type>::type>::value> CostType;
      //typedef NCCCost<FImage1T,boost::is_integral<typename PixelChannelType<typename FImage1T::pixel_type>::type>::value> CostType;
      
      typedef typename CostType::accumulator_type                      AccumChannelT;
      typedef typename PixelChannelCast<typename FImage1T::pixel_type,AccumChannelT>::type AccumT;
      typedef typename ImageView<Vector<float,9> >::pixel_accessor     PatchAcc;
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
            for (int32 dy = zone.image_region().min().y(); dy < zone.image_region().max().y(); ++dy){
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
            ImageView<typename FImage2T::pixel_type> left_raster_crop = 
              crop(left_raster,left_zone);
            ImageView<typename FImage2T::pixel_type> right_raster_crop = 
              crop(right_raster,left_zone+disparity_abs-search_range.min());
        
            CostType cost_function( left_raster_crop, right_raster_crop, m_kernel_size );    
            cost_metric = fast_box_sum<AccumChannelT>(cost_function( left_raster_crop,
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

  public:
    typedef PixelMask<Vector2f> pixel_type;
    typedef pixel_type result_type;
    typedef ProceduralPixelAccessor<ParabolaSubpixelView> pixel_accessor;

    ParabolaSubpixelView(ImageViewBase<DImageT>    const& disparity,
                         ImageViewBase<Image1T>    const& left_image,
                         ImageViewBase<Image2T>    const& right_image,
                         PrefilterModeType prefilter_mode, float prefilter_width,
                         Vector2i const& kernel_size ) :
      m_disparity( disparity.impl() ), m_left_image( left_image.impl() ),
      m_right_image( right_image.impl() ), 
      m_kernel_size( kernel_size ),
      m_prefilter_mode(prefilter_mode), m_prefilter_width(prefilter_width)
      {
      VW_ASSERT( m_disparity.cols() == m_left_image.cols() &&
                 m_disparity.rows() == m_left_image.rows(),
                 ArgumentErr() << "SubpixelView: Disparity image must match left image." );

      // We get a considerable speedup in our 2d subpixel correlation if
      // we go ahead and compute the pseudoinverse of the A matrix (where
      // each row in A is [ x^2 y^2 xy x y 1] (our 2d parabolic surface)
      // for the range of x = [-1:1] and y = [-1:1].
      static float pinvA_data[] =/*
        { 1.0/6,  1.0/6,  1.0/6, -1.0/3, -1.0/3, -1.0/3,  1.0/6,  1.0/6,  1.0/6,
          1.0/6, -1.0/3,  1.0/6,  1.0/6, -1.0/3,  1.0/6,  1.0/6, -1.0/3,  1.0/6,
          1.0/4,    0.0, -1.0/4,    0.0,    0.0,    0.0, -1.0/4,    0.0,  1.0/4,
          -1.0/6, -1.0/6, -1.0/6,    0.0,    0.0,   0.0,  1.0/6,  1.0/6,  1.0/6,
          -1.0/6,    0.0,  1.0/6, -1.0/6,    0.0, 1.0/6, -1.0/6,    0.0,  1.0/6,
          -1.0/9,  2.0/9, -1.0/9,  2.0/9,  5.0/9, 2.0/9, -1.0/9,  2.0/9, -1.0/9 };*/
        {  1.0/6, -1.0/3,  1.0/6,  1.0/6, -1.0/3,  1.0/6,   1.0/6, -1.0/3,  1.0/6,  // = a
         1.0/6,  1.0/6,  1.0/6, -1.0/3, -1.0/3, -1.0/3,   1.0/6,  1.0/6,  1.0/6,  // = b
         1.0/4,    0.0, -1.0/4,    0.0,    0.0,    0.0,  -1.0/4,    0.0,  1.0/4,  // = c
        -1.0/6,    0.0,  1.0/6, -1.0/6,    0.0,  1.0/6,  -1.0/6,    0.0,  1.0/6,  // = d
        -1.0/6, -1.0/6, -1.0/6,    0.0,    0.0,    0.0,   1.0/6,  1.0/6,  1.0/6,  // = e
        -1.0/9,  2.0/9, -1.0/9,  2.0/9,   5.0/9, 2.0/9,  -1.0/9,  2.0/9, -1.0/9 };// = f
      m_p_A_matrix = Matrix<float,6,9>( pinvA_data );
    }

    inline int32 cols  () const { return m_disparity.cols(); }
    inline int32 rows  () const { return m_disparity.rows(); }
    inline int32 planes() const { return 1; }

    inline pixel_accessor origin() const { return pixel_accessor( *this, 0, 0 ); }
    inline pixel_type operator() ( int32 /*i*/, int32 /*j*/, int32 /*p*/ = 0 ) const {
      vw_throw( NoImplErr() << "SubpixelView:operator() has not been implemented." );
      return pixel_type();
    }

    // Block rasterization section does the actual work
    typedef CropView<ImageView<pixel_type> > prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i const& bbox ) const {

      // Prerasterize the disparity for this crop and work out it's
      // range. This could be faster if we subdivided this down by the
      // areas with equal values.
      ImageView<PixelMask<Vector2i> > disparity_subsection = crop(m_disparity,bbox);

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

      // Unfortunately using the more convenient prefilter_image function yields slightly different results!
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

    template <class DestT>
    inline void rasterize( DestT const& dest, BBox2i const& bbox ) const {
      vw::rasterize(prerasterize(bbox), dest, bbox );
    }
  }; // End class ParabolaSubpixelView

  template <class DImageT, class Image1T, class Image2T>
  ParabolaSubpixelView<DImageT, Image1T, Image2T>
  parabola_subpixel( ImageViewBase<DImageT>    const& disparity,
                     ImageViewBase<Image1T>    const& left_image,
                     ImageViewBase<Image2T>    const& right_image,
                     PrefilterModeType prefilter_mode, float prefilter_width,
                     Vector2i const& kernel_size ) {
    typedef ParabolaSubpixelView<DImageT, Image1T, Image2T> result_type;
    return result_type( disparity.impl(), left_image.impl(), right_image.impl(),
                        prefilter_mode, prefilter_width, kernel_size );
  }

//----------------------------------------------------------------

  enum PyramidSubpixelView_Algorithm {
    SUBPIXEL_LUCAS_KANADE = 0,
    SUBPIXEL_FAST_AFFINE  = 1,
    SUBPIXEL_BAYES_EM     = 2
  };

  /// An image view for performing image correlation using affine sub-pixel correlation.
  template <class ImageT1, class ImageT2, class ImageTD>
  class PyramidSubpixelView : public ImageViewBase<PyramidSubpixelView<ImageT1, ImageT2, ImageTD> > {

    ImageTD m_disparity_map;
    ImageT1 m_left_image;
    ImageT2 m_right_image;

    // General Settings - Could any of these be shared with the other classes?
    Vector2i       m_kernel_size;
    int32          m_max_pyramid_levels;
    PyramidSubpixelView_Algorithm m_algorithm;

    // These two variables pick a prefilter which is applied to each pyramid level
    PrefilterModeType m_prefilter_mode; ///< See Prefilter.h for the types
    float m_prefilter_width;     ///< Preprocessing filter width


  public:
    typedef PixelMask<Vector2f> pixel_type;
    typedef pixel_type          result_type;
    typedef ProceduralPixelAccessor<PyramidSubpixelView> pixel_accessor;

    PyramidSubpixelView(ImageViewBase<ImageTD> const& disparity_map,
                        ImageViewBase<ImageT1> const& left_image,
                        ImageViewBase<ImageT2> const& right_image,
                        PrefilterModeType prefilter_mode, float prefilter_width,
                        Vector2i const& kernel_size,
                        int32 max_pyramid_levels,
                        PyramidSubpixelView_Algorithm algorithm) :
      m_disparity_map(disparity_map.impl()),
      m_left_image(left_image.impl()), m_right_image(right_image.impl()),
      m_kernel_size(kernel_size), 
      m_max_pyramid_levels(max_pyramid_levels), m_algorithm(algorithm),
      m_prefilter_mode(prefilter_mode), m_prefilter_width(prefilter_width) {

      // Basic assertions
      VW_ASSERT( m_disparity_map.cols() == m_left_image.cols() &&
                 m_disparity_map.rows() == m_left_image.rows(),
                 ArgumentErr() << "PyramidSubpixelView::PyramidSubpixelView():  Disparity image must match left image.\n");

      VW_ASSERT((m_left_image.channels() == 1) && (m_left_image.planes() == 1) &&
                (m_right_image.channels() == 1) && (m_right_image.planes() == 1),
                ArgumentErr() << "PyramidSubpixelView::PyramidSubpixelView(): multi-channel, multi-plane images are not supported.\n");

      // Max pyramid can't go below 0 ... or bayes em won't process anything
      if ( m_max_pyramid_levels < 0 ) m_max_pyramid_levels = 0;
    }

    // Standard ImageView interface methods
    inline int32 cols  () const { return m_left_image.cols(); }
    inline int32 rows  () const { return m_left_image.rows(); }
    inline int32 planes() const { return 1; }

    inline pixel_accessor origin() const { return pixel_accessor( *this, 0, 0 ); }

    inline pixel_type operator()(float /*x*/, float /*y*/, int32 /*p*/ = 0) const {
      vw_throw(NoImplErr() << "PyramidSubpixelView::operator() is not yet implemented.");
      return PixelMask<Vector2f>(); // Never reached
    }


    /// \cond INTERNAL
    typedef CropView<ImageView<pixel_type> > prerasterize_type;
    inline prerasterize_type prerasterize(BBox2i const& bbox) const {

#if VW_DEBUG_LEVEL > 0
      Stopwatch watch;
      watch.start();
#endif
      // Find the range of disparity values for this patch.
      ImageView<pixel_type > disparity_map_patch =
        crop(m_disparity_map, bbox);
      BBox2i search_range = get_disparity_range(disparity_map_patch);

      // The area in the right image that we'll be searching is
      // determined by the bbox of the left image plus the search range.
      BBox2i left_crop_bbox(bbox);
      BBox2i right_crop_bbox(bbox.min() + search_range.min(),
                             bbox.max() + search_range.max());

      // The correlator requires the images to be the same size. The
      // search bbox will always be larger than the given left image
      // bbox, so we just make the left bbox the same size as the right bbox.
      left_crop_bbox.max() = left_crop_bbox.min() + right_crop_bbox.size();

      // Finally, we must adjust both bounding boxes to account for
      // the size of the kernel itself.
      right_crop_bbox.min() -= m_kernel_size;
      right_crop_bbox.max() += m_kernel_size;
      left_crop_bbox.min() -=  m_kernel_size;
      left_crop_bbox.max() +=  m_kernel_size;

      // We crop the images to the expanded bounding box and edge
      // extend in case the new bbox extends past the image bounds.
      ImageView<float> left_image_patch, right_image_patch;
      // parabola subpixel does the same preprocessing as the pyramid correlator
      // Using the convenient prefilter function yields slightly different results here
      //  for some reason so we use this method to avoid having to redo all the tests.
      if (m_prefilter_mode == PREFILTER_LOG) {
        stereo::LaplacianOfGaussian prefilter(m_prefilter_width);
        left_image_patch  = crop(prefilter.filter(m_left_image ), left_crop_bbox );
        right_image_patch = crop(prefilter.filter(m_right_image), right_crop_bbox);
      } else {
        if (m_prefilter_mode == PREFILTER_MEANSUB) {
          stereo::SubtractedMean prefilter(m_prefilter_width);
          left_image_patch  = crop(prefilter.filter(m_left_image ), left_crop_bbox );
          right_image_patch = crop(prefilter.filter(m_right_image), right_crop_bbox);
        } else { // PREFILTER_NONE
          stereo::NullOperation prefilter;
          left_image_patch  = crop(prefilter.filter(m_left_image ), left_crop_bbox );
          right_image_patch = crop(prefilter.filter(m_right_image), right_crop_bbox);
        }
      }
      disparity_map_patch = crop(edge_extend(m_disparity_map, ZeroEdgeExtension()),
                                 left_crop_bbox);

      // Adjust the disparities to be relative to the cropped
      // image pixel locations
      PixelMask<Vector2f> disparity_patch_translation( search_range.min() );
      disparity_map_patch -= disparity_patch_translation;

      std::vector< ImageView<float> > l_patches, r_patches;
      std::vector<BBox2i> rois;
      ImageView<pixel_type > d_subpatch;

      // This initialization makes a difference only if the number of
      // levels is 0, which is a valid situation.
      d_subpatch = disparity_map_patch;

      // I'd like for image subsampling to use a gaussian when
      // downsampling however it was introducing some edge effects
      // that I couldn't figure out within a reasonable time frame.
      for ( int32 i = 0; i < m_max_pyramid_levels; i++ ) {
        if ( i > 0 ) {
          // Building all other levels
          l_patches.push_back( subsample( l_patches.back(), 2 ) );
          r_patches.push_back( subsample( r_patches.back(), 2 ) );
          ImageView<PixelMask<Vector2f> > d_subpatch_buf = disparity_subsample( d_subpatch );
          d_subpatch = d_subpatch_buf;
          rois.push_back( rois.back()/2 );
        } else {
          // First step down from native resolution
          l_patches.push_back( subsample( left_image_patch, 2 ) );
          r_patches.push_back( subsample( right_image_patch, 2 ) );
          d_subpatch = disparity_subsample( disparity_map_patch );
          rois.push_back( BBox2i( m_kernel_size[0], m_kernel_size[1],
                                  bbox.width(), bbox.height() ) );
        }
      }

      // Loop through all but final pyramid levels.
      for ( int32 i = m_max_pyramid_levels-1; i >= 0; i-- ) {

        switch(m_algorithm) {
        case SUBPIXEL_LUCAS_KANADE:
          subpixel_optimized_LK_2d(d_subpatch,
                                   l_patches[i], r_patches[i],
                                   m_kernel_size[0], m_kernel_size[1],
                                   rois[i], true, true, false );
          break;
        case SUBPIXEL_FAST_AFFINE:
          subpixel_optimized_affine_2d(d_subpatch,
                                       l_patches[i], r_patches[i],
                                       m_kernel_size[0], m_kernel_size[1],
                                       rois[i], true, true, false );
          break;
        case SUBPIXEL_BAYES_EM:
          subpixel_optimized_affine_2d_EM(d_subpatch,
                                          l_patches[i], r_patches[i],
                                          m_kernel_size[0], m_kernel_size[1],
                                          rois[i], true, true, false );
          break;
        default:
          vw_throw(NoImplErr() << "Invalid algorithm selection passed to PyramidSubpixelView.");
        }

        BBox2i crop_bbox;
        if ( i > 0 )
          crop_bbox = BBox2i(0,0,l_patches[i-1].cols(), l_patches[i-1].rows());
        else
          crop_bbox = BBox2i(0,0,left_image_patch.cols(), left_image_patch.rows());
        ImageView<pixel_type > d_subpatch_buf =
          crop(disparity_upsample(edge_extend(d_subpatch)), crop_bbox);
        d_subpatch = d_subpatch_buf;
      } // End loop through pyramid levels

      disparity_map_patch = d_subpatch;

      // Perfrom final pass at native resolution
      switch(m_algorithm) {
      case SUBPIXEL_LUCAS_KANADE:
        subpixel_optimized_LK_2d(disparity_map_patch,
                                 left_image_patch, right_image_patch,
                                 m_kernel_size[0], m_kernel_size[1],
                                 BBox2i(m_kernel_size[0],m_kernel_size[1],
                                 bbox.width(), bbox.height()),
                                 true, true, false );
        break;
      case SUBPIXEL_FAST_AFFINE:
        subpixel_optimized_affine_2d(disparity_map_patch,
                                     left_image_patch, right_image_patch,
                                     m_kernel_size[0], m_kernel_size[1],
                                     BBox2i(m_kernel_size[0],m_kernel_size[1],
                                            bbox.width(), bbox.height()),
                                     true, true, false );
        break;
      case SUBPIXEL_BAYES_EM:
        subpixel_optimized_affine_2d_EM(disparity_map_patch,
                                        left_image_patch, right_image_patch,
                                        m_kernel_size[0], m_kernel_size[1],
                                        BBox2i(m_kernel_size[0],m_kernel_size[1],
                                               bbox.width(), bbox.height()),
                                        true, true, false );
        break;
      default:
        vw_throw(NoImplErr() << "Invalid algorithm selection passed to PyramidSubpixelView.");
      }

      // Undo the above adjustment
      disparity_map_patch += disparity_patch_translation;

#if VW_DEBUG_LEVEL > 0
      watch.stop();
      vw_out(DebugMessage,"stereo") << "Tile " << bbox << " processed in " << watch.elapsed_seconds() << " s\n";
#endif

      // This may seem confusing, but we must crop here so that the
      // good pixel data is placed into the coordinates specified by
      // the bbox.  This allows rasterize to touch those pixels
      // using the coordinates inside the bbox.  The pixels outside
      // those coordinates are invalid, and they never get accessed.
      return crop(disparity_map_patch, BBox2i(m_kernel_size[0]-bbox.min()[0],
                                              m_kernel_size[1]-bbox.min()[1],
                                              m_left_image.cols(),
                                              m_left_image.rows() ));
    }

    template <class DestT> inline void rasterize(DestT const& dest, BBox2i const& bbox) const {
      vw::rasterize(prerasterize(bbox), dest, bbox);
    }
    /// \endcond
  };

  // Set of wrapper functions to help use PyramidSubpixelView
  template <class ImageT1, class ImageT2, class DisparityT>
  PyramidSubpixelView<ImageT1, ImageT2, DisparityT>
  lk_subpixel( ImageViewBase<DisparityT> const& disparity_map,
               ImageViewBase<ImageT1   > const& left_image,
               ImageViewBase<ImageT2   > const& right_image,
               PrefilterModeType prefilter_mode, float prefilter_width,
               Vector2i const& kernel_size,
               int max_pyramid_levels = 2 ) {
    typedef PyramidSubpixelView<ImageT1,ImageT2,DisparityT> result_type;

    return result_type( disparity_map.impl(), left_image.impl(),
                        right_image.impl(), 
                        prefilter_mode, prefilter_width,
                        kernel_size,
                        max_pyramid_levels,  SUBPIXEL_LUCAS_KANADE);
  }

  template <class ImageT1, class ImageT2, class DisparityT>
  PyramidSubpixelView<ImageT1, ImageT2, DisparityT>
  affine_subpixel( ImageViewBase<DisparityT> const& disparity_map,
                   ImageViewBase<ImageT1   > const& left_image,
                   ImageViewBase<ImageT2   > const& right_image,
                   PrefilterModeType prefilter_mode, float prefilter_width,
                   Vector2i const& kernel_size,
                   int max_pyramid_levels = 2 ) {
    typedef PyramidSubpixelView<ImageT1,ImageT2,DisparityT> result_type;

    return result_type( disparity_map.impl(), left_image.impl(),
                        right_image.impl(), 
                        prefilter_mode, prefilter_width,
                        kernel_size,
                        max_pyramid_levels,  SUBPIXEL_FAST_AFFINE);
  }

  template <class ImageT1, class ImageT2, class DisparityT>
  PyramidSubpixelView<ImageT1, ImageT2, DisparityT>
  bayes_em_subpixel( ImageViewBase<DisparityT> const& disparity_map,
                     ImageViewBase<ImageT1   > const& left_image,
                     ImageViewBase<ImageT2   > const& right_image,
                     PrefilterModeType prefilter_mode, float prefilter_width,
                     Vector2i const& kernel_size,
                     int max_pyramid_levels = 2 ) {
    typedef PyramidSubpixelView<ImageT1,ImageT2,DisparityT> result_type;

    return result_type( disparity_map.impl(), left_image.impl(),
                        right_image.impl(),
                        prefilter_mode, prefilter_width,
                        kernel_size,
                        max_pyramid_levels,  SUBPIXEL_BAYES_EM);
  }

  // End of components for Pyramid subpixel view









}} // namespace vw::stereo

#endif // __VW_STEREO_CORRELATOR_VIEW__
