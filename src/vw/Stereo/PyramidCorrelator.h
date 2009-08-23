// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_STEREO_CORRELATOR_H__
#define __VW_STEREO_CORRELATOR_H__

#include <vw/Image/ImageView.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Stereo/DisparityMap.h>
#include <vw/Stereo/OptimizedCorrelator.h>

namespace vw {
namespace stereo {
  class PyramidCorrelator {

    BBox2 m_initial_search_range;
    Vector2i m_kernel_size;
    float m_cross_correlation_threshold;
    float m_corrscore_rejection_threshold;
    int m_cost_blur;
    stereo::CorrelatorType m_correlator_type;    
    int m_pyramid_min_image_dimension;

    std::string m_debug_prefix;

    /*
    // Reduce the disparity size by a factor of two by averaging the pixels
    template <class PixelDisparity>
    ImageView<PixelDisparity> subsample_by_two(ImageView<PixelDisparity> &img) {

      //ImageView<PixelDisparity<float> >

      ImageView<PixelDisparity> outImg(img.cols()/2, img.rows()/2,img.planes());		
      int32 i, j, p;
      
      for (p = 0; p < outImg.planes() ; p++) {
        for (i = 0; i < outImg.cols(); i++) {
          for (j = 0; j < outImg.rows(); j++) {  
            outImg(i,j,p) = 0.0f;
            outImg(i,j,p) += img(2*i     , 2*j    ,p);
            outImg(i,j,p) += img(2*i + 1 , 2*j    ,p);
            outImg(i,j,p) += img(2*i     , 2*j + 1,p);
            outImg(i,j,p) += img(2*i + 1 , 2*j + 1,p);
            outImg(i,j,p) /= 4;
          }
        }
      }
    
      return outImg;
    }
    */

    // Reduce the image size by a factor of two by averaging the pixels
    template <class PixelT>
    ImageView<PixelT> subsample_by_two(ImageView<PixelT> &img) {
      
      ImageView<PixelT> outImg(img.cols()/2, img.rows()/2,img.planes());		
      int32 i, j, p;
    
      for (p = 0; p < outImg.planes() ; p++) {
        for (i = 0; i < outImg.cols(); i++) {
          for (j = 0; j < outImg.rows(); j++) {  
            outImg(i,j,p) = 0.0f;
            outImg(i,j,p) += img(2*i     , 2*j    ,p);
            outImg(i,j,p) += img(2*i + 1 , 2*j    ,p);
            outImg(i,j,p) += img(2*i     , 2*j + 1,p);
            outImg(i,j,p) += img(2*i + 1 , 2*j + 1,p);
            outImg(i,j,p) /= 4;
          }
        }
      }
      return outImg;
    }

    // Reduce the image size by a factor of two by averaging the pixels
    template <class MaskPixelT>
    ImageView<MaskPixelT> subsample_mask_by_two(ImageView<MaskPixelT> &img) {
      
      ImageView<MaskPixelT> outImg(img.cols()/2, img.rows()/2,img.planes());		
      int32 i, j, p;
    
      for (p = 0; p < outImg.planes() ; p++) {
        for (i = 0; i < outImg.cols(); i++) {
          for (j = 0; j < outImg.rows(); j++) {  
            if (!img(2*i     , 2*j    ,p) ||
                !img(2*i + 1 , 2*j    ,p) ||
                !img(2*i     , 2*j + 1,p) || 
                !img(2*i + 1 , 2*j + 1,p) ) 
              outImg(i,j,p) =  MaskPixelT();
            else 
              outImg(i,j,p) = ScalarTypeLimits<MaskPixelT>::highest();
          }
        }
      }
      return outImg;
    }

    // Iterate over the nominal blocks, creating output blocks for correlation
    BBox2 compute_matching_blocks(BBox2i const& nominal_block, BBox2 search_range, 
                                  BBox2i &left_block, BBox2i &right_block);

    std::vector<BBox2> compute_search_ranges(ImageView<PixelDisparity<float> > const& prev_disparity_map, 
                                             std::vector<BBox2i> nominal_blocks);
    
    void write_debug_images(int n, ImageViewRef<PixelDisparity<float> > const& disparity_map, std::vector<BBox2i> nominal_blocks);
    std::vector<BBox2i> subdivide_bboxes(ImageView<PixelDisparity<float> > const& disparity_map, BBox2i const& box);

    // do_correlation()
    //
    // Takes an image pyramid of images and conducts dense stereo
    // matching using a pyramid based approach.
    template <class ChannelT, class PreProcFilterT>
    vw::ImageView<vw::PixelDisparity<float> > do_correlation(std::vector<ImageView<ChannelT> > left_pyramid, 
                                                             std::vector<ImageView<ChannelT> > right_pyramid,
                                                             std::vector<ImageView<uint8> > left_masks, 
                                                             std::vector<ImageView<uint8> > right_masks,
                                                             PreProcFilterT const& preproc_filter) {

      int pyramid_levels = left_pyramid.size();
      BBox2 initial_search_range = m_initial_search_range / pow(2, pyramid_levels-1);
      ImageView<PixelDisparity<float> > disparity_map;
    
      // Refined the disparity map by searching in the local region
      // where the last good disparity value was found.
      for (int n = pyramid_levels - 1; n >=0; --n) {
        std::ostringstream current_level;
        current_level << n;
        TerminalProgressCallback prog(InfoMessage,"\tLevel " + current_level.str() );

        ImageView<PixelDisparity<float> > new_disparity_map(left_pyramid[n].cols(), left_pyramid[n].rows());
    
        // 1. Subdivide disparity map into subregions.  We build up
        //    the disparity map for the level, one subregion at a
        //    time.  For now, subregions that are 512x512 pixels seems
        //    to be an efficient size.
        //
        //    We also build a list of search ranges from the previous
        //    level's disparity map.  If this is the first level of the
        //    pyramid, we go with the full search range.
        std::vector<BBox2> search_ranges;
        std::vector<BBox2i> nominal_blocks;
        if (n == (pyramid_levels-1) ) {
          nominal_blocks.push_back(BBox2i(0,0,left_pyramid[n].cols(), left_pyramid[n].rows()));
          search_ranges.push_back(initial_search_range);
        } else { 
          //nominal_blocks = image_blocks(left_pyramid[n], 512, 512);
          nominal_blocks = subdivide_bboxes(disparity_map, BBox2i(0,0,left_pyramid[n].cols(), left_pyramid[n].rows()));
          search_ranges = compute_search_ranges(disparity_map, nominal_blocks);
        }

        for (unsigned r = 0; r < nominal_blocks.size(); ++r) {
          prog.report_progress((float)r/nominal_blocks.size());
      
          // Given a block from the left image, compute the bounding
          // box of pixels we will be searching in the right image
          // given the disparity range for the current left image
          // bbox.
          //
          // There's no point in correlating in areas where the second
          // image has no data, so we adjust the block sizes here to avoid
          // doing unnecessary work.
          BBox2i left_block, right_block;
          BBox2i right_image_workarea = BBox2i(Vector2i(nominal_blocks[r].min().x()+int(floor(search_ranges[r].min().x())),
                                                        nominal_blocks[r].min().y()+int(floor(search_ranges[r].min().y()))),
                                               Vector2i(nominal_blocks[r].max().x()+int(ceil(search_ranges[r].max().x())),
                                                        nominal_blocks[r].max().y()+int(ceil(search_ranges[r].max().y()))));
          BBox2i right_image_bounds = BBox2i(0,0,
                                             right_pyramid[n].cols(),
                                             right_pyramid[n].rows());
          right_image_workarea.crop(right_image_bounds);
          if (right_image_workarea.width() == 0 || right_image_workarea.height() == 0) { continue; }
          //           BBox2i left_image_workarea = BBox2i(Vector2i(right_image_workarea.min().x()-int(floor(search_ranges[r].min().x()))),
          //                                                        rights_image_workarea.min().y()-int(floor(search_ranges[r].min().y()))),
          //                                               Vector2i(right_image_workarea.max().x()-int(ceil(search_ranges[r].max().x())),
          //                                                        right_image_workarea.max().y()-int(ceil(search_ranges[r].max().y()))));
          //           if (left_image_workarea.width() <= 0 || left_image_workarea.height() <= 0) { continue; }
          //           nominal_blocks[r] = left_image_workarea;
          BBox2 adjusted_search_range = compute_matching_blocks(nominal_blocks[r], search_ranges[r], left_block, right_block);
      
          //   2. Run the correlation for this level.  We pass in the
          //      offset (difference) between the adjusted_search_range
          //      and original search_ranges[r] so that this can be added
          //      back in when setting the final disparity.
          float h_disp_offset = search_ranges[r].min().x() - adjusted_search_range.min().x();
          float v_disp_offset = search_ranges[r].min().y() - adjusted_search_range.min().y();
      
          // Place this block in the proper place in the complete
          // disparity map.
          //
          // FIXME: this is an unecessary extra copy that is done as a
          // workaround for a constness bug. -mbroxton
          ImageView<ChannelT> block1 = crop(edge_extend(left_pyramid[n],ReflectEdgeExtension()),left_block);
          ImageView<ChannelT> block2 = crop(edge_extend(right_pyramid[n],ReflectEdgeExtension()),right_block);
          ImageView<PixelDisparity<float> > disparity_block;

          disparity_block = this->correlate( block1, block2,
                                             adjusted_search_range, 
                                             Vector2(h_disp_offset, v_disp_offset), 
                                             preproc_filter );

          crop(new_disparity_map, nominal_blocks[r]) = crop(disparity_block, 
                                                            m_kernel_size[0], 
                                                            m_kernel_size[1],
                                                            nominal_blocks[r].width(),
                                                            nominal_blocks[r].height());
        }
        prog.report_finished();
        
    
        // Clean up the disparity map by rejecting outliers in the lower
        // resolution levels of the pyramid.  These are some settings that
        // seem to work well in practice.
        int32 rm_half_kernel = 5;
        double rm_min_matches_percent = 0.5;
        double rm_threshold = 0.0;

        // At the lower levels, we want to do a little bit of
        // disparity map clean-up to prevent spurious matches from
        // impacting the search range estimates, but at the top level
        // of the pyramid, we would rather return the "raw" disparity
        // map straight from the correlator.
        if (n != 0)
          disparity_map = disparity::mask(disparity::clean_up(new_disparity_map,
                                                              rm_half_kernel, 
                                                              rm_half_kernel,
                                                              rm_threshold,
                                                              rm_min_matches_percent),
                                          left_masks[n], right_masks[n]); 
        else
          disparity_map = disparity::mask(new_disparity_map, left_masks[n], right_masks[n]);
        
        // Debugging output at each level
        if (m_debug_prefix.size() > 0)
          write_debug_images(n, disparity_map, nominal_blocks);
      }

      return disparity_map;
    }

    template <class ViewT, class PreProcFilterT>
    vw::ImageView<vw::PixelDisparity<float> > correlate(ImageViewBase<ViewT> const& left_image, 
                                                        ImageViewBase<ViewT> const& right_image, 
                                                        BBox2 search_range, Vector2 offset, 
                                                        PreProcFilterT const& preproc_filter) {
      
      vw::stereo::OptimizedCorrelator correlator( BBox2i(Vector2i(int(floor(search_range.min().x())), int(ceil(search_range.min().y()))),
                                                         Vector2i(int(floor(search_range.max().x())), int(ceil(search_range.max().y()))) ),
                                                  m_kernel_size[0],
                                                  m_cross_correlation_threshold,
                                                  m_corrscore_rejection_threshold,
                                                  m_cost_blur, m_correlator_type);
      ImageView<PixelDisparity<float> > result = correlator( left_image.impl(), right_image.impl(), preproc_filter);
      
      for (int j = 0; j < result.rows(); ++j) 
        for (int i = 0; i < result.cols(); ++i) 
          if (!result(i,j).missing()) {
            result(i,j).h() += offset[0];
            result(i,j).v() += offset[1];
          }
      return result;
    }
    
  public:

    /// Correlator Constructor
    ///
    /// Set pyramid_levels to 0 to force the use of a single pyramid
    /// level (essentially disabling pyramid correlation).
    PyramidCorrelator(BBox2 initial_search_range, 
                      Vector2i kernel_size, 
                      float cross_correlation_threshold = 1, 
                      float corrscore_rejection_threshold = 1.0, 
                      int cost_blur = 1,
                      stereo::CorrelatorType correlator_type = ABS_DIFF_CORRELATOR,
                      int pyramid_min_image_dimension = 256) :
      m_initial_search_range(initial_search_range), 
      m_kernel_size(kernel_size), 
      m_cross_correlation_threshold(cross_correlation_threshold), 
      m_corrscore_rejection_threshold(corrscore_rejection_threshold), 
      m_cost_blur(cost_blur),
      m_correlator_type(correlator_type),
      m_pyramid_min_image_dimension(pyramid_min_image_dimension) {
      m_debug_prefix = "";
    }
    
    /// Turn on debugging output.  The debug_file_prefix string is
    /// used as a prefix for all debug image files.
    void set_debug_mode(std::string const& debug_file_prefix) { m_debug_prefix = debug_file_prefix; }
    
    template <class ViewT, class MaskViewT, class PreProcFilterT>
    ImageView<PixelDisparity<float> > operator() (ImageViewBase<ViewT> const& left_image,
                                                  ImageViewBase<ViewT> const& right_image,
                                                  ImageViewBase<MaskViewT> const& left_mask,
                                                  ImageViewBase<MaskViewT> const& right_mask,
                                                  PreProcFilterT const& preproc_filter) {

      typedef typename ViewT::pixel_type pixel_type;
      typedef typename PixelChannelType<pixel_type>::type channel_type;
      
      VW_ASSERT(left_image.impl().cols() == right_image.impl().cols() && 
                left_image.impl().rows() == right_image.impl().rows(),
                ArgumentErr() << "Correlator(): input image dimensions do not match.");

      VW_ASSERT(left_image.impl().cols() == left_mask.impl().cols() && 
                left_image.impl().rows() == left_mask.impl().rows(),
                ArgumentErr() << "Correlator(): input image and mask dimensions do not match.");

      VW_ASSERT(left_image.impl().cols() == right_mask.impl().cols() && 
                left_image.impl().rows() == right_mask.impl().rows(),
                ArgumentErr() << "Correlator(): input image and mask dimensions do not match.");

      // Compute the number of pyramid levels needed to reach the
      // minimum image resolution set by the user.  
      // 
      // Level 0 : finest (high resolution) image
      // Level n : coarsest (low resolution) image ( n == pyramid_levels - 1 )
      int pyramid_levels = 1;
      if (m_pyramid_min_image_dimension != 0) {
        int dimension = std::min(left_image.impl().cols(), left_image.impl().rows());
        while (dimension > m_pyramid_min_image_dimension) { pyramid_levels++; dimension /= 2;}
      } 
      vw_out(InfoMessage, "stereo") << "Initializing pyramid correlator with " << pyramid_levels << " levels.\n";

      // Build the image pyramid
      std::vector<ImageView<channel_type> > left_pyramid(pyramid_levels), right_pyramid(pyramid_levels);
      std::vector<ImageView<uint8> > left_masks(pyramid_levels), right_masks(pyramid_levels);

      left_pyramid[0] = channels_to_planes(left_image);
      right_pyramid[0] = channels_to_planes(right_image);
      left_masks[0] = channels_to_planes(left_mask);
      right_masks[0] = channels_to_planes(right_mask);
      
      // Produce the image pyramid
      for (int n = 1; n < pyramid_levels; ++n) {
        std::ostringstream ostr;
        ostr << n;

        left_pyramid[n] = subsample_by_two(left_pyramid[n-1]);
        right_pyramid[n] = subsample_by_two(right_pyramid[n-1]);
        left_masks[n] = subsample_mask_by_two(left_masks[n-1]);
        right_masks[n] = subsample_mask_by_two(right_masks[n-1]);
      }

      return do_correlation(left_pyramid, right_pyramid, left_masks, right_masks, preproc_filter);
    }
      
  };
  
}} // namespace vw::stereo

#endif // __VW_STEREO_CORRELATOR_H__
