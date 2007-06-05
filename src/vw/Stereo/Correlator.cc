
#include <vw/Stereo/Correlator.h>
#include <vw/Stereo/OptimizedCorrelator.h>
#include <vw/Core/FundamentalTypes.h>
#include <vw/Core/ProgressCallback.h>
#include <vw/Image/Transform.h>

/// A progress monitor that prints a progress bar on STDOUT.
namespace vw {
  class PyramidCorrelatorProgressCallback : public vw::ProgressCallback {
    MessageLevel m_level;
    std::string m_prefix;
  public:
    PyramidCorrelatorProgressCallback( std::string prefix, MessageLevel level = InfoMessage ) : m_level(level), m_prefix(prefix) {}
    virtual ~PyramidCorrelatorProgressCallback() {}
    
    virtual void report_progress(double progress) const {
      int pi = static_cast<int>(progress * 60);
      vw_out(m_level) << "\r" << m_prefix << " [";
      for( int i=0; i<pi; ++i ) vw_out(m_level) << "*";
      for( int i=60; i>pi; --i ) vw_out(m_level) << ".";
      vw_out(m_level) << "] " << (int)(progress*100) << "%" << std::flush;
    }
    
    virtual void report_aborted(std::string why="") const {
      vw_out(m_level) << " Aborted: " << why << std::endl;
    }
    
    virtual void report_finished() const {
      vw_out(m_level) << "\r" << m_prefix << " [************************************************************] Complete!" << std::endl;
    }
  };
} // namespace vw

// Iterate over the nominal blocks, creating output blocks for correlation
//
// To compute the block size, we must 1) Adjust the and size of
// the blocks for the right image to account for the disparity
// search range, 2) adjust the size of the left image blocks to
// be equal to the size of the right image blocks, and 3) buffer
// both blocks by the kernel size.
//
// Returns an updated search range that has been recentered for use in
// comparing the left_blocks and right_blocks.
vw::BBox2 vw::stereo::Correlator::compute_matching_blocks(BBox2i const& nominal_block, BBox2 search_range,
                                                           BBox2i &left_block, BBox2i &right_block) {
  
  left_block = nominal_block;
    
  // The bounds of the right box depend on the size of the requested search range.
  right_block = BBox2i(Vector2i(nominal_block.min().x()+int(floor(search_range.min().x())),
                                nominal_block.min().y()+int(floor(search_range.min().y()))),
                       Vector2i(nominal_block.max().x()+int(ceil(search_range.max().x())),
                                nominal_block.max().y()+int(ceil(search_range.max().y()))));
  
  // Ensure that the left and right blocks are the same size.
  left_block.max() = Vector2i(nominal_block.min().x() + right_block.width(),
                              nominal_block.min().y() + right_block.height());
  
  // Pad all blocks by the kernel size
  right_block.min() -= Vector2i(m_kernel_size[0], m_kernel_size[1]);
  right_block.max() += Vector2i(m_kernel_size[0], m_kernel_size[1]);
  left_block.min() -= Vector2i(m_kernel_size[0], m_kernel_size[1]);
  left_block.max() += Vector2i(m_kernel_size[0], m_kernel_size[1]);
  
  return BBox2(0,0,int(ceil(search_range.width())),int(ceil(search_range.height())));
}

//  Use the previous level's disparity map to narrow the search range
//  for each nominal block. At each new level, we buffer the search
//  range slightly so that we cover all feasible disparity levels at
//  the higher level of resolution.
std::vector<vw::BBox2> vw::stereo::Correlator::compute_search_ranges(ImageView<PixelDisparity<float> > const& prev_disparity_map, 
                                                                      std::vector<BBox2i> nominal_blocks) {
  std::vector<BBox2> search_ranges(nominal_blocks.size());
  std::vector<int> good_pixel_vec(nominal_blocks.size());
  std::vector<int> fixed_blocks;

  // Step two: compute the search ranges from the disparity map.
  for (int r = 0; r < nominal_blocks.size(); ++r)
    search_ranges[r] = disparity::get_disparity_range(crop(prev_disparity_map,(nominal_blocks[r])/2), good_pixel_vec[r], false);
      

  // Step two: adjust the search range or fix it if it came from a
  // block weith zero valid pixels.
  for (int r = 0; r < nominal_blocks.size(); ++r) {

    // Pick a reasonable search range based on neighboring tiles rather
    // than skipping out entirely here.
    if (good_pixel_vec[r] == 0) {
      // Find adjancent tiles and use the search ranges of adjacent
      // tiles to seed us with a reasonable search range of our own.
      //
      // To do this test, we expand this block by 1 pixel in every
      // direction and then test for intersections with other blocks.
      BBox2i this_bbox = nominal_blocks[r];
      this_bbox.expand(1);
      bool found_one_match = false;
      for (int b = 0; b < nominal_blocks.size(); ++b) 
        if (good_pixel_vec[b] != 0 && this_bbox.intersects(nominal_blocks[b])) {
          if (found_one_match) {
            search_ranges[r].grow(search_ranges[b]);
          } else {
            found_one_match = true;
            fixed_blocks.push_back(r);
            search_ranges[r] = search_ranges[b];
          }
        }
    }
  } 

  // Set the good_pixel_vec[r] to a positive value for the fixed
  // blocks so that the search range is adjusted correctly be the code
  // below.
  for (int i = 0; i < fixed_blocks.size(); ++i) 
    good_pixel_vec[fixed_blocks[i]] = 1;

  // Step three: scale up the search range for the next pyramid
  // level and pad it here.
  for (int r = 0; r < nominal_blocks.size(); ++r) {
    if (good_pixel_vec[r] != 0) {
      search_ranges[r] *= 2;
      search_ranges[r].min().x() -= 2;
      search_ranges[r].max().x() += 2;
      search_ranges[r].min().y() -= 2;
      search_ranges[r].max().y() += 2;
    }
  } 
  return search_ranges;
}

vw::ImageView<vw::PixelDisparity<float> > vw::stereo::Correlator::correlate(ImageView<uint8> left_image, ImageView<uint8> right_image, 
                                                                            BBox2 search_range, Vector2 offset) {
  vw::stereo::OptimizedCorrelator correlator( int(floor(search_range.min().x())), int(ceil(search_range.max().x())),
                                              int(floor(search_range.min().y())), int(ceil(search_range.max().y())),
                                              m_kernel_size[0], m_kernel_size[1],
                                              false, m_cross_correlation_threshold,
                                              false, false );
  ImageView<PixelDisparity<float> > result = correlator( left_image, right_image, true);
  
  for (int j = 0; j < result.rows(); ++j) 
    for (int i = 0; i < result.cols(); ++i) 
      if (!result(i,j).missing()) {
        result(i,j).h() += offset[0];
        result(i,j).v() += offset[1];
      }

  return result;
}



// do_correlation()
//
// Takes an image pyramid of SLOG images and conducts dense stereo
// matching using a pyramid based approach.
vw::ImageView<vw::PixelDisparity<float> > vw::stereo::Correlator::do_correlation(std::vector<ImageView<uint8> > left_slog_pyramid, 
                                                                                 std::vector<ImageView<uint8> > right_slog_pyramid,
                                                                                 std::vector<ImageView<bool> > left_masks, 
                                                                                 std::vector<ImageView<bool> > right_masks) {

  int pyramid_levels = left_slog_pyramid.size();
  
  // Clean up the disparity map by rejecting outliers in the lower
  // resolution levels of the pyramid.  These are some settings that
  // seem to work well in practice.
  int32 rm_half_kernel = 5;
  double rm_min_matches_percent = 0.7;
  double rm_threshold = 0;
  
  // Run the full correlation complete with left-right/right-left
  // checks.  This produces a rough, low res version of the disparity
  // map.
  BBox2 initial_search_range = m_initial_search_range / pow(2, pyramid_levels-1);
  std::ostringstream current_level0;
  current_level0 << (pyramid_levels-1);
  PyramidCorrelatorProgressCallback prog0("\tLevel " + current_level0.str());
  prog0.report_finished();
  ImageView<PixelDisparity<float> > disparity_map = disparity::clean_up(this->correlate(left_slog_pyramid[pyramid_levels-1], right_slog_pyramid[pyramid_levels-1], initial_search_range, Vector2(0,0)),
                                                                        rm_half_kernel, 
                                                                        rm_half_kernel,
                                                                        rm_threshold,
                                                                        rm_min_matches_percent); 
  disparity::mask(disparity_map, left_masks[pyramid_levels-1], right_masks[pyramid_levels-1]);

  // Debugging: Print out the disparity map at the lowest resolution
  if (m_debug_prefix.size() > 0) {
    std::ostringstream current_level;
    current_level << pyramid_levels-1;
    BBox2 disp_range = disparity::get_disparity_range(disparity_map);
    write_image( m_debug_prefix + "-H-" + current_level.str() + ".jpg", normalize(clamp(select_channel(disparity_map,0), disp_range.min().x(), disp_range.max().x() )));
    write_image( m_debug_prefix + "-V-" + current_level.str() + ".jpg", normalize(clamp(select_channel(disparity_map,1), disp_range.min().y(), disp_range.max().y() )));
  }
  
  // Refined the disparity map by searching in the local region where the last good disparity value was found
  for (int n = pyramid_levels - 2; n >=0; --n) {
    std::ostringstream current_level;
    current_level << n;
    PyramidCorrelatorProgressCallback prog("\tLevel " + current_level.str() );
    
    // 1. Subdivide disparity map into subregions.  We build up the
    //    disparity map for the level, one subregion at a time.
    std::vector<BBox2i> nominal_blocks = image_blocks(left_slog_pyramid[n], 512, 512);
    ImageView<PixelDisparity<float> > new_disparity_map(left_slog_pyramid[n].cols(), left_slog_pyramid[n].rows());
    
    // First we build a list of search ranges from the previous
    // level's disparity map.
    std::vector<BBox2> search_ranges = compute_search_ranges(disparity_map, nominal_blocks);

    for (int r = 0; r < nominal_blocks.size(); ++r) {
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
      BBox2i right_image_workarea = BBox2i(-int(floor(search_ranges[r].min().x())),
                                           -int(floor(search_ranges[r].min().y())),
                                           right_slog_pyramid[n].cols() - int(ceil(search_ranges[r].max().x())),
                                           right_slog_pyramid[n].rows() - int(ceil(search_ranges[r].max().y())));
      BBox2i right_image_bounds = BBox2i(0,0,
                                         right_slog_pyramid[n].cols(),
                                         right_slog_pyramid[n].rows());
      right_image_workarea.crop(right_image_bounds);
      nominal_blocks[r].crop(right_image_workarea);
      if (nominal_blocks[r].width() == 0 || nominal_blocks[r].height() == 0) { continue; }
      BBox2 adjusted_search_range = compute_matching_blocks(nominal_blocks[r], search_ranges[r], left_block, right_block);
      
      //   2. Run the correlation for this level.  We pass in the
      //      offset (difference) between the adjusted_search_range
      //      and original search_ranges[r] so that this can be added
      //      back in when setting the final disparity.
      float h_disp_offset = search_ranges[r].min().x() - adjusted_search_range.min().x();
      float v_disp_offset = search_ranges[r].min().y() - adjusted_search_range.min().y();
      
      // Place this block in the proper place in the complete
      // disparity map.
      ImageView<PixelDisparity<float> > disparity_block = this->correlate( crop(edge_extend(left_slog_pyramid[n],ReflectEdgeExtension()),left_block),
                                                                           crop(edge_extend(right_slog_pyramid[n],ReflectEdgeExtension()),right_block),
                                                                           adjusted_search_range, 
                                                                           Vector2(h_disp_offset, v_disp_offset) );
      crop(new_disparity_map, nominal_blocks[r]) = crop(disparity_block, 
                                                        m_kernel_size[0], 
                                                        m_kernel_size[1],
                                                        nominal_blocks[r].width(),
                                                        nominal_blocks[r].height());
    }
    prog.report_finished();

    
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
    if (m_debug_prefix.size() > 0) {
      std::ostringstream current_level;
      current_level << n;
      double min_h_disp, min_v_disp, max_h_disp, max_v_disp;
      BBox2 disp_range = disparity::get_disparity_range(disparity_map);
      write_image( m_debug_prefix+"-H-" + current_level.str() + ".jpg", normalize(clamp(select_channel(disparity_map,0), disp_range.min().x(), disp_range.max().x() )));
      write_image( m_debug_prefix+"-V-" + current_level.str() + ".jpg", normalize(clamp(select_channel(disparity_map,1), disp_range.min().y(), disp_range.max().y() )));
    }
    
  }
  // Do subpixel correlation
  subpixel_correlation(disparity_map, left_slog_pyramid[0], right_slog_pyramid[0], m_kernel_size[0], m_kernel_size[1], m_do_h_subpixel, m_do_v_subpixel);
  return disparity_map;

}

