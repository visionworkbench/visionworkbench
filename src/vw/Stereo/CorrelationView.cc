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

#include <vw/Core/Stopwatch.h>
#include <vw/Stereo/CorrelationView.h>
#include <vw/Stereo/PreFilter.h>
#include <vw/Stereo/Correlate.h>
#include <vw/Stereo/DisparityMap.h>
#include <vw/Core/Thread.h>
#include <vw/Image/Algorithms.h>
#include <vw/Image/ErodeView.h>
#include <vw/FileIO/DiskImageView.h>
#include <vw/Stereo/Correlation.h>

#include <boost/foreach.hpp>
#include <ctime>
#include <vw/Image/Manipulation.h>

namespace vw { namespace stereo {

  /// Downsample a mask by two. If at least two mask pixels in a 2x2
  /// region are on, the output pixel is on.
  struct SubsampleMaskByTwoFunc: public ReturnFixedType<uint8> {
    BBox2i work_area() const { return BBox2i(0,0,2,2); }
    
    template <class PixelAccessorT>
    typename boost::remove_reference<typename PixelAccessorT::pixel_type>::type
    operator()(PixelAccessorT acc) const {
      
      typedef typename PixelAccessorT::pixel_type PixelT;
      
      uint8 count = 0;
      if (*acc) count++;
      acc.next_col();
      if (*acc) count++;
      acc.advance(-1,1);
      if (*acc) count++;
      acc.next_col();
      if (*acc) count++;
      if (count > 1)
        return PixelT(ScalarTypeLimits<PixelT>::highest());
      return PixelT();
    }
  }; // End struct SubsampleMaskByTwoFunc
  
  ImageView<uint8> subsample_mask_by_two(ImageView<uint8> const& input) {
    return subsample(per_pixel_accessor_filter(input.impl(), SubsampleMaskByTwoFunc()), 2);
  }
  
  //=========================================================================
  // Correlation with several pyramid levels
  bool PyramidCorrelationView::
  build_image_pyramids(BBox2i const& bbox, int32 const max_pyramid_levels,
                       std::vector<ImageView<PixelGray<float>>> & left_pyramid,
                       std::vector<ImageView<PixelGray<float>>> & right_pyramid,
                       std::vector<ImageView<uint8>>            & left_mask_pyramid,
                       std::vector<ImageView<uint8>>            & right_mask_pyramid) const {
    
    Vector2i half_kernel = m_kernel_size/2;

    // Init the pyramids: Highest resolution image is stored at index zero.
    left_pyramid.resize      (max_pyramid_levels + 1);
    right_pyramid.resize     (max_pyramid_levels + 1);
    left_mask_pyramid.resize (max_pyramid_levels + 1);
    right_mask_pyramid.resize(max_pyramid_levels + 1);
  
    int32 max_upscaling = 1 << max_pyramid_levels;
    vw_out(VerboseDebugMessage, "stereo") << "max_upscaling = " << max_upscaling << std::endl;
    BBox2i left_global_region, right_global_region;
    // Region in the left image is the input bbox, expanded by half the kernel size
    //  to make sure that we have full base of support for the stereo correlation.
    // - Make sure that the kernel buffer size is large enough to support the kernel
    //   buffer at the lower resolution levels!
    Vector2i region_offset = half_kernel * max_upscaling;
    vw_out(VerboseDebugMessage, "stereo") << "pyramid region offset = "
                                          << region_offset << std::endl;
    left_global_region = bbox;
    left_global_region.expand(region_offset);
    // Region in the right image is the left region plus search range offsets
    // - What is the last upscaling term for?
    right_global_region = left_global_region + m_search_region.min();
    right_global_region.max() += m_search_region.size();// + Vector2i(max_upscaling,max_upscaling); 
    // Total increase in right size = 2*region_offset + search_region_size + max_upscaling

    // TODO(oalexan1): left and right global regions must be multiples of 2^num_levels.
    // Be careful with m_search_region.
    
    vw_out(VerboseDebugMessage, "stereo") << "Left pyramid base bbox:  "
                                          << left_global_region  << std::endl;
    vw_out(VerboseDebugMessage, "stereo") << "Right pyramid base bbox: "
                                          << right_global_region << std::endl;
  
    // Extract the lowest resolution layer
    // - Constant extension is used here to help the correlator make matches near the image edge.
    // TODO(oalexan1): This must be extension with no-data values as otherwise we get junk.
    left_pyramid      [0] = crop(edge_extend(m_left_image,  ConstantEdgeExtension()),
                                 left_global_region);
    right_pyramid     [0] = crop(edge_extend(m_right_image, ConstantEdgeExtension()),
                                 right_global_region);
    // TODO(oalexan1): Put below ZeroEdgeExtension!
    left_mask_pyramid [0] = crop(edge_extend(m_left_mask,   ConstantEdgeExtension()),
                                 left_global_region);
    right_mask_pyramid[0] = crop(edge_extend(m_right_mask,  ConstantEdgeExtension()),
                                 right_global_region);

#if VW_DEBUG_LEVEL > 0
    VW_OUT(DebugMessage,"stereo") << " > Left ROI: "    << left_global_region
                                  << "\n > Right ROI: " << right_global_region << "\n";
#endif

    // Fill in the nodata of the left and right images with a mean
    // pixel value. This helps with the edge quality of a DEM.
    // - Note that this will not fill in the edge-extended values.
    // TODO(oalexan1): This is not accurate! Fill instead with avg from neighbors
    PixelGray<float> left_mean;
    PixelGray<float> right_mean;
    try {
      left_mean  = mean_pixel_value(subsample(copy_mask(left_pyramid [0],
                                                        create_mask(left_mask_pyramid [0], 0)), 2));
      right_mean = mean_pixel_value(subsample(copy_mask(right_pyramid[0],
                                                        create_mask(right_mask_pyramid[0], 0)), 2));
    } catch (const ArgumentErr& err) {
      // Mean pixel value will throw an argument error if there
      // are no valid pixels. If that happens, it means either the
      // left or the right image is full masked.
      return false;
    }
    // Now paste the mean value into the masked pixels
    left_pyramid [0] = apply_mask(copy_mask(left_pyramid[0],
                                            create_mask(left_mask_pyramid[0], 0)),
                                  left_mean);
    right_pyramid[0] = apply_mask(copy_mask(right_pyramid[0],
                                            create_mask(right_mask_pyramid[0], 0)),
                                  right_mean);

    vw_out(DebugMessage, "stereo") << "Left  pyramid base size = "
                                   << bounding_box(left_pyramid[0]) << std::endl;
    vw_out(DebugMessage, "stereo") << "Right pyramid base size = "
                                   << bounding_box(right_pyramid[0]) << std::endl;

#if 0
    // This is useful for debugging
    std::ostringstream ostr;
    ostr << "tile_" << bbox.min().x() << "_" << bbox.min().y() << "_"
         << bbox.width() << "_" << bbox.height();
    std::string suff = ostr.str();
    
    // Left image and mask
    std::string left_image_file = "left_image_" + suff + ".tif";
    std::cout << std::endl; // To deal with the progress dialog
    std::cout << "Writing: " << left_image_file << std::endl;
    write_image(left_image_file, left_pyramid[0]);
    std::string left_mask_file = "left_mask_" + suff + ".tif";
    std::cout << "Writing: " << left_mask_file << std::endl;
    write_image(left_mask_file, left_mask_pyramid[0]);

    // Right image and mask
    std::string right_image_file = "right_image_" + suff + ".tif";
    std::cout << "Writing: " << right_image_file << std::endl;
    write_image(right_image_file, right_pyramid[0]);
    std::string right_mask_file = "right_mask_" + suff + ".tif";
    std::cout << "Writing: " << right_mask_file << std::endl;
    std::cout << std::endl; // To deal with the progress dialog
    write_image(right_mask_file, right_mask_pyramid[0]);
#endif
    
    // Reduce the mask images from the expanded-size region to the
    //   actual sized region.
    // - The actual sized region is just the input bbox plus the
    //   search range, no expanded base of support or anything.
    // - The mask is not used in the expanded base of support region
    //   of the image, only the main region.
    // - The larger mask size before is used to compute the mean color
    //   values.
    // - Zero edge extension is used here so we don't treat the edge
    //   extended pixels as valid.
    //   Currently this mainly affects the SGM algorithm which needs
    //   to know not to solve for those pixels.
    BBox2i right_mask = bbox + m_search_region.min();
    right_mask.max() += m_search_region.size();
    left_mask_pyramid [0] = crop(edge_extend(m_left_mask, ZeroEdgeExtension()), bbox);
    right_mask_pyramid[0] = crop(edge_extend(m_right_mask,ZeroEdgeExtension()), right_mask);

    vw_out(DebugMessage, "stereo") << "Left  pyramid mask base size = " << bbox       << std::endl;
    vw_out(DebugMessage, "stereo") << "Right pyramid mask base size = " << right_mask << std::endl;

    // Set up smoothing kernel used before downsampling.
    std::vector<typename DefaultKernelT<PixelGray<float>>::type > kernel = 
      generate_pyramid_smoothing_kernel();
    std::vector<uint8> mask_kern(max(m_kernel_size));
    std::fill(mask_kern.begin(), mask_kern.end(), 1);

    // Smooth and downsample to build the pyramid (don't smooth the masks)
    // TODO(oalexan1): Use a PixelMask rather handling the image and its mask separately!
    for (int32 i = 1; i <= max_pyramid_levels; i++) {
      left_pyramid      [i] = subsample(separable_convolution_filter(left_pyramid [i-1],
                                                                     kernel, kernel),2);
      right_pyramid     [i] = subsample(separable_convolution_filter(right_pyramid[i-1],
                                                                      kernel, kernel),2);
      left_mask_pyramid [i] = subsample_mask_by_two(left_mask_pyramid [i-1]);
      right_mask_pyramid[i] = subsample_mask_by_two(right_mask_pyramid[i-1]);
    
      vw_out(DebugMessage, "stereo") << "--- Created pyramid level " << i << std::endl;    
      vw_out(DebugMessage, "stereo") << "Left  pyramid size = "
                                     << bounding_box(left_pyramid[i]     ) << std::endl;
      vw_out(DebugMessage, "stereo") << "Right pyramid size = "
                                     << bounding_box(right_pyramid[i]    ) << std::endl;
      vw_out(DebugMessage, "stereo") << "Left  pyramid mask size = "
                                     << bounding_box(left_mask_pyramid[i]) << std::endl;
      vw_out(DebugMessage, "stereo") << "Right pyramid mask size = "
                                     << bounding_box(right_mask_pyramid[i]) << std::endl;
      vw_out(DebugMessage, "stereo") << "Level search size = "
                                     << (m_search_region.size() / (1 << i)) << std::endl;
    }

    // Apply the prefilter to each pyramid level
    // TODO(oalexan1): Use a PixelMask rather handling the image and its mask separately!
    for (int32 i = 0; i <= max_pyramid_levels; i++) {
      left_pyramid [i] = prefilter_image(left_pyramid [i], m_prefilter_mode, m_prefilter_width);
      right_pyramid[i] = prefilter_image(right_pyramid[i], m_prefilter_mode, m_prefilter_width);
    }
    
    return true;
  }

  /// Filter out small blobs of valid pixels (they are usually bad)
  void PyramidCorrelationView::
  disparity_blob_filter(ImageView<pixel_typeI> &disparity, int level,
                        int max_blob_area) const {

    // Throw out blobs with this many pixels or fewer
    int scaling = 1 << level;
    int area    = max_blob_area / scaling;
    if (area < 1)
      return; // Skip if erode turned off
    //vw_out() << "Removing blobs smaller than: " << area << std::endl;

    if (0) { // DEBUG
      vw_out() << "Writing pre-blob image...\n";
      std::ostringstream ostr;
      ostr << "disparity_preblob_" << level;
      write_image(ostr.str() + ".tif", pixel_cast<PixelMask<Vector2f>>(disparity));
      vw_out() << "Finished writing DEBUG data...\n";
    } // End DEBUG

      // Do the entire image at once!
    BBox2i tile_size = bounding_box(disparity);
    int big_size = tile_size.width();
    if (tile_size.height() > big_size) 
      big_size = tile_size.height();
    
    BlobIndexThreaded smallBlobIndex(disparity, area, big_size);
    ImageView<pixel_typeI> filtered_image = applyErodeView(disparity, smallBlobIndex);

    disparity = filtered_image;
  }

  PyramidCorrelationView::prerasterize_type
  PyramidCorrelationView::prerasterize (BBox2i const& bbox) const {

    // Sanity check, the currently processed bbox must fit in the lr_diff buffer.
    if (m_lr_disp_diff != NULL) {
      BBox2i lr_diff_box(0, 0, m_lr_disp_diff->cols(), m_lr_disp_diff->rows());
      lr_diff_box += m_region_ul;
      if (!lr_diff_box.contains(bbox)) 
        vw_throw(ArgumentErr() << "The L-R to R-L difference image domain "
                 << "does not contain the current tile.");
    }
    
    time_t start, end;
    if (m_corr_timeout){
      std::time (&start);
    }

#if VW_DEBUG_LEVEL > 0
    Stopwatch watch;
    watch.start();
#endif

    // 1.0) Determining the number of levels to process
    //      There's a maximum based on kernel size. There's also
    //      maximum defined by the search range. Here we determine
    //      the maximum based on kernel size and current bbox.
    // - max_pyramid_levels is the number of levels not including the original resolution level.
    // - Each pyramid level is shrunk by the (reduced size) kernel size on that level.
    int32 smallest_bbox      = math::min(bbox.size()); // Get smallest/largest of height/width
    int32 largest_kernel     = math::max(m_kernel_size);
    int32 max_pyramid_levels = std::floor(log(smallest_bbox)/log(2.0f)
                                          - log(largest_kernel)/log(2.0f));
    if (m_max_level_by_search < max_pyramid_levels)
      max_pyramid_levels = m_max_level_by_search;
    if (max_pyramid_levels < 1)
      max_pyramid_levels = 0;
    Vector2i half_kernel = m_kernel_size/2;
    int32 max_upscaling = 1 << max_pyramid_levels;

    // 2.0) Build the pyramids
    //      - Highest resolution image is stored at index zero.
    //      - Remember that these pyramid images are larger than just the input bbox!
    std::vector<ImageView<PixelGray<float>>> left_pyramid;
    std::vector<ImageView<PixelGray<float>>> right_pyramid;
    std::vector<ImageView<uint8>> left_mask_pyramid;
    std::vector<ImageView<uint8>> right_mask_pyramid;

    if (!build_image_pyramids(bbox, max_pyramid_levels, left_pyramid, right_pyramid, 
                              left_mask_pyramid, right_mask_pyramid)){
#if VW_DEBUG_LEVEL > 0
      watch.stop();
      double elapsed = watch.elapsed_seconds();
      vw_out(DebugMessage,"stereo") << "Tile " << bbox << " has no data. Processed in "
                                    << elapsed << " s\n";
#endif
      return prerasterize_type(ImageView<result_type>(bbox.width(), bbox.height()),
                               -bbox.min().x(), -bbox.min().y(),
                               cols(), rows());
    }
    
    // 3.0) Actually perform correlation now
    ImageView<pixel_typeI> disparity, prev_disparity, disparity_rl, prev_disparity_rl;
    std::vector<stereo::SearchParam> zones; 
    // Start off the search at the lowest resolution pyramid level.  This zone covers
    // the entire image and uses the disparity range that was loaded into the class.
    BBox2i initial_disparity_range = BBox2i(0, 0, 
                                            m_search_region.width()  / max_upscaling + 1,
                                            m_search_region.height() / max_upscaling + 1);
    zones.push_back(SearchParam(bounding_box(left_mask_pyramid[max_pyramid_levels]),
                                initial_disparity_range));
    vw_out(DebugMessage,"stereo") << "initial_disparity_range = " 
      << initial_disparity_range << "\n";

    // Perform correlation. Keep track of how much time elapsed
    // since we started and stop if we estimate that doing one more
    // image chunk will bring us over time.

    // To not slow us down with timing, we use some heuristics to
    // estimate how much time elapsed, as time to do an image chunk
    // is proportional with image area times search range area. This
    // is not completely accurate, so every now and then do actual
    // timing, no more often than once in measure_spacing seconds.
    double estim_elapsed   = 0.0;
    int    measure_spacing = 2; // seconds
    double prev_estim      = estim_elapsed;

    ImageView<result_type> subpixel_disparity;
    const bool use_sgm = (m_algorithm != VW_CORRELATION_BM); // Anything but block matching

    // Loop down through all of the pyramid levels, low res to high res.
    for (int32 level = max_pyramid_levels; level >= 0; --level) {

      const bool use_mgm = ((m_algorithm == VW_CORRELATION_MGM) || 
                            ((m_algorithm == VW_CORRELATION_FINAL_MGM) && (level == 0)));

      const bool on_last_level = (level == 0);    
      bool check_rl = false;

      int32 scaling = 1 << level; // scaling = 2^level
      if (use_sgm) {
        prev_disparity    = disparity;
        prev_disparity_rl = disparity_rl;
      }

      disparity.set_size(left_mask_pyramid[level]); // Note, no kernel padding here.
      
      // This is the number of padding pixels added to expand the base of support for the
      //  correlation kernel.  It matches the amount computed in build_image_pyramid().
      Vector2i region_offset = max_upscaling*half_kernel/scaling;
      
      vw_out(DebugMessage,"stereo") << "\n\nProcessing level: " << level 
                                    << " with size " << disparity.get_size() << "\n\n";
      vw_out(DebugMessage,"stereo") << "region_offset = " << region_offset << std::endl;
      vw_out(DebugMessage,"stereo") << "Number of zones = " << zones.size() << std::endl;

      ImageView<uint8> right_rl_mask, left_rl_mask;

      // SGM method
      if (use_sgm) {

        // Mimic processing in normal case with a single zone
        BBox2i disparity_range = BBox2i(0,0,m_search_region.width()/scaling,
                                        m_search_region.height()/scaling);
        SearchParam zone(bounding_box(left_mask_pyramid[level]), // Non-padded size
                         disparity_range);
        
        // The input zone is in the normal pixel coordinates for this level.
        // We need to convert it to a bbox in the expanded base of support image at this level.
        // - For right region of support, remember that the right image was already shifted
        //   to account for the disparity locations so only size is accounted for here.
        BBox2i left_region = zone.image_region() + region_offset;
        left_region.expand(half_kernel);
        BBox2i right_region = left_region;
        right_region.max() += zone.disparity_range().size();
        
        ImageView<pixel_typeI> *prev_disp_ptr=0; // Pass in upper level disparity
        if (level < max_pyramid_levels) {
          prev_disp_ptr = &prev_disparity;
          vw_out(VerboseDebugMessage, "stereo") << "Disparity size      = "
                                                << bounding_box(disparity    ) << std::endl;
          vw_out(VerboseDebugMessage, "stereo") << "Prev Disparity size = "
                                                << bounding_box(prev_disparity) << std::endl;
        }
        
        // Note: The masks contain exactly the region of interests from the input masks,
        //       with the right mask containing the region offset by the search range.
        //       The left mask size should exactly equal the output size here.
        // - To be fully accurate, should crop the right mask slightly
        //   but SGM does not require this.
        
        boost::shared_ptr<SemiGlobalMatcher> sgm_matcher_ptr;
        crop(disparity, zone.image_region()) // This crop is not needed in SGM case!
          = calc_disparity_sgm(m_cost_type,
                               crop(left_pyramid [level], left_region), 
                               crop(right_pyramid[level], right_region),
                               // Specify that the whole cropped region is valid
                               left_region - left_region.min(), 
                               zone.disparity_range().size(), 
                               m_kernel_size, use_mgm, m_sgm_subpixel_mode,
                               m_sgm_search_buffer, m_memory_limit_mb,
                               sgm_matcher_ptr,
                               &(left_mask_pyramid[level]), &(right_mask_pyramid[level]),
                               prev_disp_ptr);
        // Delete the matcher pointer right after we use it to free up its large buffers.
        // - On the last level we need to generate the subpixel view before we delete it.
        // - Note that the subpixel image is created BEFORE filtering out bad pixels at the
        //   integer level.  This is ok, we just apply the integer filter results before 
        //   returning the subpixel disparity.  Doing things in this order avoids having
        //   to keep both the LR and the RL large accumulation buffers in memory at the 
        //   same time but it does mean we waste time computing subpixel values for pixels
        //   that will get invalidated later.
        if (level == 0)
          subpixel_disparity = sgm_matcher_ptr->create_disparity_view_subpixel(disparity);
        sgm_matcher_ptr.reset();

        // If the user requested a left<->right consistency check at this level,
        //   compute right to left disparity.
        if (m_consistency_threshold >= 0.0 && level >= m_min_consistency_level) {

          check_rl = true;

          // Update m_search_region for this level
          BBox2i search_region_level = m_search_region;
          search_region_level /= scaling;
          
          // To properly perform the reverse correlation, we need to
          // fix the ROIs to account for the different sizes of the
          // left and right images and make sure they line up with the
          // previous disparity image and the input masks.

          BBox2i right_reverse_region = right_region;
          // Shift to right
          BBox2i left_reverse_region = left_region - zone.disparity_range().size();
          // Enlarge to fit the search range
          left_reverse_region.max() += 2*zone.disparity_range().size();

          if (m_write_debug_images) { // DEBUG
            std::cout << "\n====== RL CHECK =====\n";
            std::cout << "left region = " << left_region << std::endl;
            std::cout << "right region = " << right_region << std::endl;
            std::cout << "scaling       = " << scaling << std::endl;
            std::cout << "half_kernel   = " << half_kernel << std::endl;
            std::cout << "region_offset = " << region_offset << std::endl;
            std::cout << "right_reverse_region  = " << right_reverse_region << std::endl;
            std::cout << "left_reverse_region   = " << left_reverse_region << std::endl;
            std::cout << "m_search_region     = " << m_search_region << std::endl;
            std::cout << "search_region_level = " << search_region_level << std::endl;
            std::cout << "search_region_level.max() = " << search_region_level.max() << std::endl;
            std::cout << "zone.image_region() = " << zone.image_region() << std::endl;
            std::cout << "zone.disparity_range() = " << zone.disparity_range() << std::endl;
            std::cout << "zone.disparity_range().max() = "
                      << zone.disparity_range().max() << std::endl;
            std::cout << "right_pyramid[level].size()   = "
                      << bounding_box(right_pyramid[level]) << std::endl;
            std::cout << "left_pyramid[level].size()   = "
                      << bounding_box(left_pyramid[level]) << std::endl;
          }

          ImageView<pixel_typeI> *prev_disp_ptr_rl=0; // Pass in upper level disparity
          if (level < max_pyramid_levels) {
            prev_disp_ptr_rl = &prev_disparity_rl;
            vw_out(VerboseDebugMessage, "stereo") << "Prev Disparity size RL = "
                                                  << bounding_box(prev_disparity_rl) << std::endl;
          }

          
          // Set the masks to the exact size needed and adjust the
          // position of the left one to match the image shift.
          BBox2i right_mask_bbox = BBox2i(0,0, right_reverse_region.width ()-2*half_kernel[0],
                                          right_reverse_region.height()-2*half_kernel[1]);
          BBox2i left_mask_bbox  = BBox2i(0,0, left_reverse_region.width  ()-2*half_kernel[0],
                                          left_reverse_region.height ()-2*half_kernel[1])
            - zone.disparity_range().size();
          
          // Rasterization needed because we pass these by address below.
          // - Could be avoided with some refactoring.
          right_rl_mask = crop(edge_extend(right_mask_pyramid[level], ZeroEdgeExtension()),
                               right_mask_bbox); 
          left_rl_mask  = crop(edge_extend(left_mask_pyramid [level], ZeroEdgeExtension()),
                               left_mask_bbox);

          if (m_write_debug_images) { // DEBUG
            std::cout << "left mask input: "
                      << bounding_box(left_mask_pyramid [level]) << std::endl;
            std::cout << "right mask input: "
                      << bounding_box(right_mask_pyramid[level]) << std::endl;
            std::cout << "left mask: "  << left_mask_bbox << std::endl;
            std::cout << "right mask: " << right_mask_bbox << std::endl;
          }

          //write_image("lr_cropl.tif", crop(left_pyramid [level], left_region));
          //write_image("lr_cropr.tif", crop(right_pyramid[level], right_region));
          //write_image("rl_cropl.tif", crop(edge_extend(left_pyramid[level]),
          // left_reverse_region));
          //write_image("rl_cropr.tif", crop(right_pyramid[level], right_reverse_region));


          // Write out masks with borders inserted - should line up with images!
          //right_mask_bbox.expand(half_kernel);
          //write_image("rl_cropRmaskPad.tif", crop(edge_extend(right_rl_mask,
          // ZeroEdgeExtension()), right_mask_bbox));
          //BBox2i temp = bounding_box(left_rl_mask);
          //temp.expand(half_kernel);
          //write_image("rl_cropLmaskPad.tif", crop(edge_extend(left_rl_mask,
          // ZeroEdgeExtension()), temp));
          //write_image("lr_result.tif", crop(disparity, zone.image_region()));

          boost::shared_ptr<SemiGlobalMatcher> sgm_right_matcher_ptr;
          disparity_rl = calc_disparity_sgm(m_cost_type,
                                            crop(right_pyramid[level], right_reverse_region),
                                            crop(edge_extend(left_pyramid[level]),
                                                 left_reverse_region),
                                            // Full RR region
                                            right_reverse_region - right_reverse_region.min(),
                                            zone.disparity_range().size(), 
                                            m_kernel_size, use_mgm, m_sgm_subpixel_mode,
                                            m_sgm_search_buffer, m_memory_limit_mb,
                                            sgm_right_matcher_ptr,
                                            &(right_rl_mask), 
                                            &(left_rl_mask),
                                            prev_disp_ptr_rl);
          sgm_right_matcher_ptr.reset(); // Immediately delete this to clear memory.

          //write_image("rl_result.tif", disparity_rl);

          // Convert from RL to negative LR values
          pixel_typeI offset(zone.disparity_range().size());
          disparity_rl -= offset;

          //write_image("rl_result2.tif", disparity_rl);

          // Prepare to save the L-R to R-L discrepancy
          // discrepancy. Do it only at level 0. Find the upper-left
          // corner offset. Take into account that m_lr_disp_diff does
          // not span the full image, but only the portion needed by
          // the parent.
          Vector2i ul_corner_offset(0, 0);
          ImageView<PixelMask<float>> * lr_disp_diff = NULL;
          if (level == 0 && m_lr_disp_diff != NULL) {
            ul_corner_offset = zone.image_region().min() + bbox.min() - m_region_ul;
            lr_disp_diff = m_lr_disp_diff;
          }
          
          // Find pixels where the disparity distance is greater than m_consistency_threshold
          //  and flag those pixels as invalid.
          // Crop not needed for SGM!
          const bool verbose = true;
          stereo::cross_corr_consistency_check(crop(disparity, zone.image_region()), 
                                               disparity_rl, m_consistency_threshold,
                                               lr_disp_diff, ul_corner_offset,
                                               verbose);
          
          //write_image("lr_result_cons.tif", crop(disparity, zone.image_region()));
                                               
          //disparity_rl *= -1;
          //write_image("rl_result3.tif", disparity_rl);
          //disparity_rl *= -1;
          
          // This offset was used for the cross-corr check but must be
          // reset for the previous disparity image
          disparity_rl += offset;                                               
                                               
        } // End of last level right to left disparity check


      } else { // Normal block matching method

        // 3.1) Process each zone with their refined search estimates
        // - The zones are subregions of the image with similar disparities
        //   that we identified in previous iterations.
        // - Prioritize the zones which take less time so we don't miss
        //   a bunch of tiles because we spent all our time on a slow one.
        // Sort the zones, smallest to largest.
        std::sort(zones.begin(), zones.end(), SearchParamLessThan()); 
        BOOST_FOREACH(SearchParam const& zone, zones) {

          // The input zone is in the normal pixel coordinates for this level.
          // We need to convert it to a bbox in the expanded base of support image at this level.
          BBox2i left_region = zone.image_region() + region_offset; // Kernel width offset
          left_region.expand(half_kernel);

          // Make right region contain all of the needed match area.
          BBox2i right_region = left_region + zone.disparity_range().min();
          right_region.max() += zone.disparity_range().size();
          
          // Setting up the ROIs in this way means that the range of
          // disparities calculated is always >=0

          // Check timing estimate to see if we should go ahead with this zone or quit.
          SearchParam params(left_region, zone.disparity_range());
          double next_elapsed = m_seconds_per_op * params.search_volume();
          if (m_corr_timeout > 0.0 && estim_elapsed + next_elapsed > m_corr_timeout){
            vw_out() << "Tile: " << bbox << " reached timeout: "
                     << m_corr_timeout << " s" << std::endl;
            break;
          }else
            estim_elapsed += next_elapsed;

          // See if it is time to actually accurately compute the time
          if (m_corr_timeout > 0.0 && estim_elapsed - prev_estim > measure_spacing){
            std::time (&end);
            double diff = std::difftime(end, start);
            estim_elapsed = diff;
            prev_estim = estim_elapsed;
          }

          // Compute left to right disparity vectors in this zone.
          // - The cropped regions we pass in have padding for the kernel.
          crop(disparity, zone.image_region())
            = calc_disparity(m_cost_type,
                             crop(left_pyramid [level], left_region), 
                             crop(right_pyramid[level], right_region),
                             // Specify that the whole cropped region is valid
                             left_region - left_region.min(), 
                             zone.disparity_range().size(), 
                             m_kernel_size);

          // TODO(oalexan1):  Support checks at higher levels like with SGM!
          // If at the last level and the user requested a left<->right consistency check,
          //   compute right to left disparity.
          if (m_consistency_threshold >= 0 && level == 0) {

            check_rl = true;

            // Check the time again before moving on with this
            SearchParam params2(right_region, zone.disparity_range());
            double next_elapsed = m_seconds_per_op * params2.search_volume();
            if (m_corr_timeout > 0.0 && estim_elapsed + next_elapsed > m_corr_timeout){
              vw_out() << "Tile: " << bbox << " reached timeout: "
                       << m_corr_timeout << " s" << std::endl;
              break;
            }else{
              estim_elapsed += next_elapsed;
            }
            // Compute right to left disparity in this zone       
            // TODO(oalexan1): Below use masks, and edge extend with zero edge extension!
            disparity_rl = calc_disparity(m_cost_type,
                                          crop(edge_extend(right_pyramid[level]), right_region),
                                          crop(edge_extend(left_pyramid [level]),
                                               left_region - zone.disparity_range().size()),
                                          right_region - right_region.min(),
                                          zone.disparity_range().size(), m_kernel_size)
              - pixel_typeI(zone.disparity_range().size());

            // Prepare to save the L-R to R-L disparity
            // discrepancy. Do it only at level 0. Find the upper-left
            // corner offset. Take into account that m_lr_disp_diff
            // does not span the full image, but only the portion
            // needed by the parent.
            Vector2i ul_corner_offset(0, 0);
            ImageView<PixelMask<float>> * lr_disp_diff = NULL;
            if (level == 0 && m_lr_disp_diff != NULL) {
              ul_corner_offset = zone.image_region().min() + bbox.min() - m_region_ul;
              lr_disp_diff = m_lr_disp_diff;
            }

            // Find pixels where the disparity distance is greater than m_consistency_threshold
            const bool verbose = true;
            stereo::cross_corr_consistency_check(crop(disparity, zone.image_region()),
                                                 disparity_rl, m_consistency_threshold,
                                                 lr_disp_diff, ul_corner_offset,
                                                 verbose);
          } // End of last level right to left disparity check

            // Fix the offsets to account for cropping.
          crop(disparity, zone.image_region()) += pixel_typeI(zone.disparity_range().min());
        } // End of zone loop
      } // End non-SGM case

        // 3.2a) Filter the disparity so we are not processing more
        // than we need to. Inner function filtering is only to
        // catch "speckle" type noise of individual outliers.  Outer
        // function just merges the masks over the filtered disparity
        // image.  //const int32 rm_half_kernel = 5;
      const float rm_min_matches_percent = 0.5;
      const float rm_threshold = 3.0;

      // R-L image only needs to be filtered when using SGM because it is used in that case
      //  to initialize the search range of the following pyramid level.

      if (m_filter_half_kernel > 0) { // Skip filtering if zero radius passed in
        if (!on_last_level) {
          disparity = disparity_mask(disparity_cleanup_using_thresh
                                     (disparity,
                                      m_filter_half_kernel, m_filter_half_kernel,
                                      rm_threshold,
                                      rm_min_matches_percent),
                                     left_mask_pyramid [level],
                                     right_mask_pyramid[level]);
          if (check_rl && use_sgm) {
            disparity_rl = disparity_mask(disparity_cleanup_using_thresh
                                          (disparity_rl,
                                           m_filter_half_kernel, m_filter_half_kernel,
                                           rm_threshold,
                                           rm_min_matches_percent),
                                          right_rl_mask, 
                                          left_rl_mask);
          }
        } else { // On the last level
        
          // We don't do a single hot pixel check on the final level as it leaves a border.
          disparity = disparity_mask(rm_outliers_using_thresh
                                     (disparity,
                                      m_filter_half_kernel, m_filter_half_kernel,
                                      rm_threshold,
                                      rm_min_matches_percent),
                                     left_mask_pyramid [level],
                                     right_mask_pyramid[level]);

          // No need to filter R-L disparity with SGM on the last level.
        }
      } // End of 
      
      // The kernel based filtering tends to leave isolated blobs behind.
      disparity_blob_filter(disparity, level, m_blob_filter_area);
      if (check_rl && !on_last_level)
        disparity_blob_filter(disparity_rl, level, m_blob_filter_area);

      // 3.2b) Refine search estimates but never let them go beyond
      // the search region defined by the user
      // - SGM method does not use zones.
      if (!on_last_level && !use_sgm) {
        const size_t next_level = level-1;
        zones.clear();
        
        // On the next resolution level, break up the image area into multiple
        // smaller zones with similar disparities.  This helps minimize
        // the total amount of searching done on the image.
        subdivide_regions(disparity, bounding_box(disparity),
                          zones, m_kernel_size);
      
        scaling >>= 1;
        
        // Scale search range defines the maximum search range that
        // is possible in the next step. This (at lower levels) will
        // actually be larger than the search range that the user
        // specified. We are able to do this because we are taking
        // advantage of the half kernel padding needed at the hight
        // level of the pyramid.
        BBox2i scale_search_region(0,0,
                                   right_pyramid[next_level].cols()
                                   - left_pyramid[next_level].cols(),
                                   right_pyramid[next_level].rows()
                                   - left_pyramid[next_level].rows());
        BBox2i next_zone_size = bounding_box(left_mask_pyramid[level-1]);
        
        BBox2i default_disparity_range = BBox2i(0,0,m_search_region.width(),
                                                m_search_region.height());
        
        BOOST_FOREACH(SearchParam& zone, zones) {
        
          zone.image_region() *= 2;
          zone.image_region().crop(next_zone_size);
          zone.disparity_range() *= 2;
          zone.disparity_range().expand(2); // This is practically required. Our
          // correlation will fail if the search has only one solution.
          // - Increasing this expansion number improves results slightly but
          //   significantly increases the processing times.
          
          zone.disparity_range().crop(scale_search_region);
          
          if (zone.disparity_range().empty()) {
            zone.disparity_range() = default_disparity_range; // Reset invalid disparity!
          }
        } // End zone update loop

      } // End zone handling
      
      if (m_write_debug_images) { // DEBUG
        vw_out() << "Writing DEBUG data...\n";
        BBox2i scaled = bbox/1; // Why divide by 2??
        std::ostringstream ostr;
        ostr << "disparity_" << scaled.min()[0] << "_"
             << scaled.min()[1] << "_" << scaled.max()[0] << "_"
             << scaled.max()[1] << "_" << level;
        write_image(ostr.str() + ".tif", pixel_cast<PixelMask<Vector2f>>(disparity));

        if (use_sgm && check_rl)
          write_image(ostr.str() + "_rl.tif", pixel_cast<PixelMask<Vector2f>>(disparity_rl));
          
        if (!use_sgm) { // SGM does not use zones
          std::ofstream f((ostr.str() + "_zone.txt").c_str());
          BOOST_FOREACH(SearchParam& zone, zones) {
            f << zone.image_region() << " " << zone.disparity_range() << "\n";
          }
          f.close();
        }
        write_image(ostr.str() + "left.tif",  left_pyramid [level]);
        write_image(ostr.str() + "right.tif", right_pyramid[level]);
        write_image(ostr.str() + "lmask.tif", left_mask_pyramid [level]);
        write_image(ostr.str() + "rmask.tif", right_mask_pyramid[level]);
        vw_out() << "Finished writing DEBUG data...\n";
      } // End DEBUG
      
        //if (level == 2)
        //  vw_throw(NoImplErr() << "DEBUG");
      
    } // End of the level loop

    VW_ASSERT(bbox.size() == bounding_box(disparity).size(),
              MathErr() << "PyramidCorrelation: Solved disparity "
              << "doesn't match requested bbox size.");

#if VW_DEBUG_LEVEL > 0
    watch.stop();
    double elapsed = watch.elapsed_seconds();
    vw_out(DebugMessage,"stereo") << "Tile " << bbox << " processed in "
                                  << elapsed << " s\n";
    if (m_corr_timeout > 0.0){
      vw_out(DebugMessage,"stereo")
        << "Elapsed (actual/estimated/ratio): " << elapsed << ' '
        << estim_elapsed << ' ' << elapsed/estim_elapsed << std::endl;
    }
#endif

    // If filtering removed disparities, also invalidate m_lr_disp_diff at the same pixels.
    if (m_lr_disp_diff != NULL) {
      Vector2i ul_corner_offset = bbox.min() - m_region_ul;
      for (int r = 0; r < disparity.rows(); r++){
        for (int c = 0; c < disparity.cols(); c++){
          if (!is_valid(disparity(c, r))) 
            (*m_lr_disp_diff)(c + ul_corner_offset[0], r + ul_corner_offset[1]).invalidate();
        }
      }
    }

    // 5.0) Reposition our result back into the global solution. Also
    // we need to correct for the offset we applied to the search
    // region. At this point we either cast to floating point or run a
    // subpixel refinement algorithm.

    if (m_algorithm != VW_CORRELATION_BM) {
    
      // Copy the filtered out pixels to the subpixel view.
      for (int r = 0; r < disparity.rows(); r++){
        for (int c = 0; c < disparity.cols(); c++){
          if (!is_valid(disparity(c,r)))
            invalidate(subpixel_disparity(c,r));
        }
      }
    
      // For SGM, subpixel correlation is performed here, not in stereo_rfne.     
      return prerasterize_type(subpixel_disparity + result_type(m_search_region.min()),
                               -bbox.min().x(), -bbox.min().y(),
                               cols(), rows());      
    } else {
      // TODO CLEANUP
      ImageView<pixel_typeI> temp = disparity + pixel_typeI(m_search_region.min());
      ImageView<result_type> float_type = pixel_cast<result_type, ImageView<pixel_typeI>>(temp);
      return prerasterize_type(float_type,
                               -bbox.min().x(), -bbox.min().y(),
                               cols(), rows());
    }
  } // End function prerasterize

}} // namespace stereo

