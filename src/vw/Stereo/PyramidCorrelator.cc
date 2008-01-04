#include <vw/Stereo/PyramidCorrelator.h>
#include <vw/Stereo/OptimizedCorrelator.h>
#include <vw/Core/FundamentalTypes.h>
#include <vw/Core/ProgressCallback.h>
#include <vw/Image/Transform.h>

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
vw::BBox2 vw::stereo::PyramidCorrelator::compute_matching_blocks(BBox2i const& nominal_block, BBox2 search_range,
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
std::vector<vw::BBox2> vw::stereo::PyramidCorrelator::compute_search_ranges(ImageView<PixelDisparity<float> > const& prev_disparity_map, 
                                                                      std::vector<BBox2i> nominal_blocks) {
  std::vector<BBox2> search_ranges(nominal_blocks.size());
  std::vector<int> good_pixel_vec(nominal_blocks.size());
  std::vector<int> fixed_blocks;

  // Step two: compute the search ranges from the disparity map.
  for (unsigned r = 0; r < nominal_blocks.size(); ++r)
    search_ranges[r] = disparity::get_disparity_range(crop(prev_disparity_map,(nominal_blocks[r])/2), good_pixel_vec[r], false);
      

  // Step two: adjust the search range or fix it if it came from a
  // block weith zero valid pixels.
  for (unsigned r = 0; r < nominal_blocks.size(); ++r) {

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
      for (unsigned b = 0; b < nominal_blocks.size(); ++b) 
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
  for (unsigned i = 0; i < fixed_blocks.size(); ++i) 
    good_pixel_vec[fixed_blocks[i]] = 1;

  // Step three: scale up the search range for the next pyramid
  // level and pad it here.
  for (unsigned r = 0; r < nominal_blocks.size(); ++r) {
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





