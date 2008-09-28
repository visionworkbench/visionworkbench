#include <vw/Stereo/PyramidCorrelator.h>
#include <vw/Stereo/OptimizedCorrelator.h>
#include <vw/Core/FundamentalTypes.h>
#include <vw/Core/ProgressCallback.h>


std::vector<vw::BBox2i> vw::stereo::PyramidCorrelator::subdivide_bboxes(vw::ImageView<vw::PixelDisparity<float> > const& disparity_map, vw::BBox2i const& box) {
  //  std::cout << "SUBDIVIDING BBOXES: " << disparity_map.cols() << " " << disparity_map.rows() << "    " << box <<"\n";
  std::vector<BBox2i> result;
  BBox2i box_div_2 = box;
  box_div_2.min() = box.min()/2;  box_div_2.max() = box.max()/2;
  BBox2 disp_range = stereo::disparity::get_disparity_range(crop(disparity_map, box_div_2));
  double cost  = box.width()*box.height()*disp_range.width()*disp_range.height();

  if (cost == 0 && 
      disp_range.max().x() == 0 && disp_range.max().y() == 0 &&
      disp_range.min().x() == 0 && disp_range.min().y() == 0) {
    // Search range is zero.  Reject this bounding box.
    //    std::cout << "\t Rejecting: " << box << " " << disp_range << " " << cost << "  " << disp_range << "\n";
    return result;
  } else if (cost < 4e5 || box.width() < 128 || box.height() < 128) {
    // The bounding box is small enough.  
    //    std::cout << "\t      Leaf: " << box << " " << disp_range << " " << cost << "\n";
    result.push_back(box);
    return result;
  } else {
    BBox2i subbox1 = box, subbox2 = box;
    if (box.width() > box.height()) {
      subbox1.max().x() = box.min().x() + box.width()/2;
      subbox2.min().x() = box.min().x() + box.width()/2;
    } else {
      subbox1.max().y() = box.min().y() + box.height()/2;
      subbox2.min().y() = box.min().y() + box.height()/2;
    }
    //    std::cout << "\tSubdividing: " << box << " " << subbox1 << " " << subbox2 << "  " << cost << "\n";
    result = subdivide_bboxes(disparity_map, subbox1);
    std::vector<BBox2i> l2 = subdivide_bboxes(disparity_map, subbox2);
    result.insert(result.end(), l2.begin(), l2.end());
    return result;
  }
}

void draw_bbox(vw::ImageView<vw::PixelRGB<float> > &view, vw::BBox2i const& bbox) {
  int u,v; 
  // Top
  v = bbox.min().y();
  for (u = bbox.min().x(); u < bbox.max().x(); ++u) 
    view(u,v) = vw::PixelRGB<float>(1.0,0.0,0.0);
  
  // Bottom
  v = bbox.max().y()-1;
  for (u = bbox.min().x(); u < bbox.max().x(); ++u) 
    view(u,v) = vw::PixelRGB<float>(1.0,0.0,0.0);
  
  // Left
  u = bbox.min().x();
  for (v = bbox.min().y(); v < bbox.max().y(); ++v) 
    view(u,v) = vw::PixelRGB<float>(1.0,0.0,0.0);
  
  // Left
  u = bbox.max().x()-1;
  for (v = bbox.min().y(); v < bbox.max().y(); ++v) 
    view(u,v) = vw::PixelRGB<float>(1.0,0.0,0.0);
}

void vw::stereo::PyramidCorrelator::write_debug_images(int n, ImageViewRef<PixelDisparity<float> > const& disparity_map, std::vector<BBox2i> nominal_blocks) {
  std::ostringstream current_level;
  current_level << n;
  BBox2 disp_range = disparity::get_disparity_range(disparity_map);
  ImageView<PixelRGB<float> > horz = normalize(clamp(select_channel(disparity_map,0), disp_range.min().x(), disp_range.max().x() ));
  ImageView<PixelRGB<float> > vert = normalize(clamp(select_channel(disparity_map,1), disp_range.min().y(), disp_range.max().y() ));
  
  for (unsigned i = 0; i < nominal_blocks.size(); ++i) {
    draw_bbox(horz, nominal_blocks[i]);
    draw_bbox(vert, nominal_blocks[i]);
  }
  
  write_image( m_debug_prefix+current_level.str()+"-H.jpg", horz);
  write_image( m_debug_prefix+current_level.str()+"-V.jpg", vert);
  
  // For debugging:
  //           write_image( m_debug_prefix+"-L-" + current_level.str() + "-mask.jpg", normalize(channel_cast<uint8>(left_masks[n])));
  //           write_image( m_debug_prefix+"-R-" + current_level.str() + "-mask.jpg", normalize(channel_cast<uint8>(right_masks[n])));
}

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





