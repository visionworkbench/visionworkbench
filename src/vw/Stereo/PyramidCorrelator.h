// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_STEREO_CORRELATOR_H__
#define __VW_STEREO_CORRELATOR_H__

#include <vw/Image/ImageView.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Image/MaskViews.h>
#include <vw/Image/Transform.h>
#include <vw/Image/Filter.h>
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
    int m_pyramid_levels;
    int m_min_subregion_dim;
    typedef PixelMask<Vector2f> PixelDisp;

    std::string m_debug_prefix;

    // Reduce the image size by a factor of two by averaging the pixels
    template <class PixelT>
    ImageView<PixelT> subsample_by_two(ImageView<PixelT> &img) {

      ImageView<PixelT> outImg(img.cols()/2, img.rows()/2,img.planes());
      int32 i, j, p;

      for (p = 0; p < outImg.planes() ; p++) {
        for (i = 0; i < outImg.cols(); i++) {
          for (j = 0; j < outImg.rows(); j++) {
            outImg(i,j,p) = PixelT(0);
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

    template <class PixelT>
    ImageView<PixelT> upsample_by_two(ImageView<PixelT> &img) {
      ImageView<PixelT> outImg(img.cols()*2, img.rows()*2,img.planes());
      int32 i, j, p;

      for (p = 0; p < outImg.planes() ; p++) {
        for (i = 0; i < outImg.cols(); i++) {
          for (j = 0; j < outImg.rows(); j++) {
            outImg(i,j,p) = img(i/2, j/2, p);
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
    BBox2 compute_matching_blocks(BBox2i const& nominal_block,
                                  BBox2 search_range,
                                  BBox2i &left_block, BBox2i &right_block);

    std::vector<BBox2>
    compute_search_ranges(ImageView<PixelDisp> const& prev_disparity_map,
                          std::vector<BBox2i> const& nominal_blocks);

    void write_debug_images(int n, ImageViewRef<PixelDisp> const& disparity_map,
                            std::vector<BBox2i> nominal_blocks);
    std::vector<BBox2i>
    subdivide_bboxes(ImageView<PixelDisp> const& disparity_map,
                     ImageView<PixelMask<uint32> > const& valid_pad,
                     BBox2i const& box);

    template <class ViewT>
    int count_valid_pixels(ImageViewBase<ViewT> const& img) {
      typedef typename ViewT::const_iterator view_iter;

      int count = 0;
      for (view_iter i = img.impl().begin(); i != img.impl().end(); i++) {
        if (i->valid())
          count++;
      }

      return count;
    }

    // do_correlation()
    //
    // Takes an image pyramid of images and conducts dense stereo
    // matching using a pyramid based approach.
    template <class ChannelT, class PreProcFilterT>
    ImageView<PixelDisp>
    do_correlation(std::vector<ImageView<ChannelT> > left_pyramid,
                   std::vector<ImageView<ChannelT> > right_pyramid,
                   std::vector<ImageView<uint8> > left_masks,
                   std::vector<ImageView<uint8> > right_masks,
                   PreProcFilterT const& preproc_filter) {

      BBox2 initial_search_range =
        m_initial_search_range / pow(2.0, m_pyramid_levels-1);
      ImageView<PixelMask<Vector2f> > disparity_map;

      // Overall Progress Bar
      TerminalProgressCallback prog( "stereo", "Pyr Search:", DebugMessage);

      // Refined the disparity map by searching in the local region
      // where the last good disparity value was found.
      for (int n = m_pyramid_levels - 1; n >=0; --n) {
        std::ostringstream current_level;
        current_level << n;

        SubProgressCallback subbar(prog,float(m_pyramid_levels-1-n)/float(m_pyramid_levels), float(m_pyramid_levels-n)/float(m_pyramid_levels));

        ImageView<PixelDisp> new_disparity_map(left_pyramid[n].cols(),
                                               left_pyramid[n].rows());

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
        if (n == (m_pyramid_levels-1) ) {
          nominal_blocks.push_back(BBox2i(0,0,left_pyramid[n].cols(),
                                          left_pyramid[n].rows()));
          search_ranges.push_back(initial_search_range);
        } else {
          std::vector<vw::uint32> x_kern(m_kernel_size.x()),
            y_kern(m_kernel_size.y());
          std::fill(x_kern.begin(), x_kern.end(), 1);
          std::fill(y_kern.begin(), y_kern.end(), 1);

          // valid_pad masks all the pixels already masked by
          // disparity_map, with the addition of a m_kernel_size/2 pad
          // around each pixel.  This is used to prevent
          // subdivide_bboxes from rejecting subregions that may
          // actually later get filled with valid pixels at a higher
          // scale (which helps prevent 'cutting' into the disparity map)
          ImageViewRef<uint32> valid =
            apply_mask(copy_mask(constant_view<uint32>(1, disparity_map.cols(),
                                                       disparity_map.rows()),
                                 disparity_map));
          ImageView<PixelMask<uint32> > valid_pad =
            create_mask(separable_convolution_filter(valid, x_kern, y_kern));

          nominal_blocks = subdivide_bboxes(disparity_map, valid_pad,
                                            BBox2i(0,0,left_pyramid[n].cols(),
                                                   left_pyramid[n].rows()));
          search_ranges = compute_search_ranges(disparity_map, nominal_blocks);
        }

        for (unsigned r = 0; r < nominal_blocks.size(); ++r) {
          subbar.report_progress((float)r/nominal_blocks.size());

          // Given a block from the left image, compute the bounding
          // box of pixels we will be searching in the right image
          // given the disparity range for the current left image
          // bbox.
          //
          // There's no point in correlating in areas where the second
          // image has no data, so we adjust the block sizes here to avoid
          // doing unnecessary work.
          BBox2i left_block, right_block;
          BBox2i right_image_workarea =
            BBox2i(Vector2i(nominal_blocks[r].min().x()+int(floor(search_ranges[r].min().x())),
                            nominal_blocks[r].min().y()+int(floor(search_ranges[r].min().y()))),
                   Vector2i(nominal_blocks[r].max().x()+int(ceil(search_ranges[r].max().x())),
                            nominal_blocks[r].max().y()+int(ceil(search_ranges[r].max().y()))));
          BBox2i right_image_bounds =
            BBox2i(0,0, right_pyramid[n].cols(), right_pyramid[n].rows());
          right_image_workarea.crop(right_image_bounds);
          if (right_image_workarea.width() == 0 ||
              right_image_workarea.height() == 0) { continue; }
          BBox2 adjusted_search_range =
            compute_matching_blocks(nominal_blocks[r],
                                    search_ranges[r], left_block, right_block);

          //   2. Run the correlation for this level.  We pass in the
          //      offset (difference) between the adjusted_search_range
          //      and original search_ranges[r] so that this can be added
          //      back in when setting the final disparity.
          float h_disp_offset =
            search_ranges[r].min().x() - adjusted_search_range.min().x();
          float v_disp_offset =
            search_ranges[r].min().y() - adjusted_search_range.min().y();

          // Place this block in the proper place in the complete
          // disparity map.
          ImageViewRef<ChannelT> block1 =
            crop(edge_extend(left_pyramid[n],ReflectEdgeExtension()),left_block);
          ImageViewRef<ChannelT> block2 =
            crop(edge_extend(right_pyramid[n],ReflectEdgeExtension()),right_block);
          ImageView<PixelDisp> disparity_block;

          disparity_block =
            this->correlate( block1, block2, adjusted_search_range,
                             Vector2(h_disp_offset, v_disp_offset),
                             preproc_filter );

          crop(new_disparity_map, nominal_blocks[r]) =
            crop(disparity_block, m_kernel_size[0], m_kernel_size[1],
                 nominal_blocks[r].width(), nominal_blocks[r].height());
        }
        subbar.report_finished();


        // Clean up the disparity map by rejecting outliers in the lower
        // resolution levels of the pyramid.  These are some settings that
        // seem to work well in practice.
        int32 rm_half_kernel = 5;
        double rm_min_matches_percent = 0.5;
        double rm_threshold = 3.0;


        ImageView<PixelDisp> disparity_map_clean;
        disparity_map_clean =
          disparity_mask(disparity_clean_up(new_disparity_map,
                                            rm_half_kernel, rm_half_kernel,
                                            rm_threshold,
                                            rm_min_matches_percent),
                         left_masks[n], right_masks[n]);

        if (n == m_pyramid_levels - 1) {
          // At the highest level of the pyramid, use the cleaned version
          // of the disparity map just obtained (since there are no
          // previous results to learn from)
          disparity_map = disparity_map_clean;
        } else if (n == 0) {
          // At the last level, return the raw results from the correlator
          disparity_map = disparity_mask(new_disparity_map,
                                         left_masks[n], right_masks[n]);
        } else {
          // If we have a missing pixel that correlated properly in
          // the previous pyramid level, use the disparity found at
          // the previous pyramid level
          ImageView<PixelDisp > disparity_map_old;
          disparity_map_old =
            2*crop(edge_extend(upsample_by_two(disparity_map),
                               ZeroEdgeExtension()),
                   BBox2i(0, 0, disparity_map_clean.cols(),
                          disparity_map_clean.rows()));
          ImageView<PixelDisp > disparity_map_old_diff =
            invert_mask(intersect_mask(disparity_map_old, disparity_map_clean));
          disparity_map =
            disparity_mask(create_mask(apply_mask(disparity_map_clean) +
                                       apply_mask(disparity_map_old_diff)),
                           left_masks[n], right_masks[n]);
        }

        // Debugging output at each level
        if (m_debug_prefix.size() > 0)
          write_debug_images(n, disparity_map, nominal_blocks);
      }
      prog.report_finished();

      return disparity_map;
    }

    template <class ViewT, class PreProcFilterT>
    ImageView<PixelDisp > correlate(ImageViewBase<ViewT> const& left_image,
                                    ImageViewBase<ViewT> const& right_image,
                                    BBox2 search_range, Vector2 offset,
                                    PreProcFilterT const& preproc_filter) {

      stereo::OptimizedCorrelator correlator( BBox2i(Vector2i(int(floor(search_range.min().x())), int(ceil(search_range.min().y()))),
                                                     Vector2i(int(floor(search_range.max().x())), int(ceil(search_range.max().y()))) ),
                                              m_kernel_size[0],
                                              m_cross_correlation_threshold,
                                              m_corrscore_rejection_threshold,
                                              m_cost_blur, m_correlator_type);
      ImageView<PixelDisp > result = correlator( left_image.impl(),
                                                 right_image.impl(),
                                                 preproc_filter);

      for (int j = 0; j < result.rows(); ++j)
        for (int i = 0; i < result.cols(); ++i)
          if ( is_valid(result(i,j)) ) {
            result(i,j)[0] += offset[0];
            result(i,j)[1] += offset[1];
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
                      int pyramid_levels = 4) :
      m_initial_search_range(initial_search_range),
      m_kernel_size(kernel_size),
      m_cross_correlation_threshold(cross_correlation_threshold),
      m_corrscore_rejection_threshold(corrscore_rejection_threshold),
      m_cost_blur(cost_blur),
      m_correlator_type(correlator_type),
      m_pyramid_levels(pyramid_levels) {
      m_debug_prefix = "";
      m_min_subregion_dim = 128;
    }

    /// Turn on debugging output.  The debug_file_prefix string is
    /// used as a prefix for all debug image files.
    void set_debug_mode(std::string const& debug_file_prefix) { m_debug_prefix = debug_file_prefix; }

    template <class ViewT, class MaskViewT, class PreProcFilterT>
    ImageView<PixelDisp > operator() (ImageViewBase<ViewT> const& left_image,
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
      vw_out(DebugMessage, "stereo") << "Initializing pyramid correlator with "
                                     << m_pyramid_levels << " levels.\n";

      // Build the image pyramid
      std::vector<ImageView<channel_type> > left_pyramid(m_pyramid_levels), right_pyramid(m_pyramid_levels);
      std::vector<ImageView<uint8> > left_masks(m_pyramid_levels), right_masks(m_pyramid_levels);

      left_pyramid[0] = channels_to_planes(left_image);    // Is this really what we want
      right_pyramid[0] = channels_to_planes(right_image);  // shouldn't we channel cast to
      left_masks[0] = channels_to_planes(left_mask);       // to scalar?
      right_masks[0] = channels_to_planes(right_mask);

      // Produce the image pyramid
      for (int n = 1; n < m_pyramid_levels; ++n) {
        left_pyramid[n] = subsample_by_two(left_pyramid[n-1]);
        right_pyramid[n] = subsample_by_two(right_pyramid[n-1]);
        left_masks[n] = subsample_mask_by_two(left_masks[n-1]);
        right_masks[n] = subsample_mask_by_two(right_masks[n-1]);
      }

      int mask_padding = std::max(m_kernel_size[0], m_kernel_size[1])/2;
      for (int n = 0; n < m_pyramid_levels; ++n) {
        left_masks[n] = apply_mask(edge_mask(left_masks[n], 0, mask_padding),0);
        right_masks[n] = apply_mask(edge_mask(right_masks[n], 0, mask_padding),0);
      }

      return do_correlation(left_pyramid, right_pyramid,
                            left_masks, right_masks, preproc_filter);
    }

  };

}} // namespace vw::stereo

#endif // __VW_STEREO_CORRELATOR_H__
