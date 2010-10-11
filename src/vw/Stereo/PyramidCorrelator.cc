// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Stereo/PyramidCorrelator.h>
#include <vw/Stereo/OptimizedCorrelator.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/FileIO/DiskImageResource.h>

using namespace vw;
using namespace stereo;

std::vector<vw::BBox2i>
PyramidCorrelator::subdivide_bboxes(ImageView<PixelMask<Vector2f> > const& disparity_map,
                                    ImageView<PixelMask<uint32> > const& valid_pad,
                                    BBox2i const& box) {
  std::vector<BBox2i> result;
  BBox2i box_div_2 = box;
  box_div_2.min() = box.min()/2;  box_div_2.max() = box.max()/2;
  BBox2 disp_range;
  try {
    disp_range = get_disparity_range(crop(disparity_map, box_div_2));
  } catch ( std::exception &/*e*/ ) {
    // If there are no good pixels, don't add this box
    if (count_valid_pixels(crop(valid_pad, box_div_2)) == 0)
      return result;

    // Return a large bbox so that we keep dividing
    disp_range = BBox2();
  }

  if (disp_range.width()*disp_range.height() <= 4 ||
      (box.width() < m_min_subregion_dim && box.height() < m_min_subregion_dim)) {
    // The bounding box is small enough.
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

    result = subdivide_bboxes(disparity_map, valid_pad, subbox1);
    std::vector<BBox2i> l2 = subdivide_bboxes(disparity_map, valid_pad, subbox2);
    result.insert(result.end(), l2.begin(), l2.end());
    return result;
  }
}

void draw_bbox(ImageView<PixelRGB<float> > &view, BBox2i const& bbox) {
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

void PyramidCorrelator::write_debug_images(int n, ImageViewRef<PixelMask<Vector2f> > const& disparity_map,
                                           std::vector<BBox2i> nominal_blocks) {
  std::ostringstream current_level;
  current_level << n;
  BBox2 disp_range;
  try {
    disp_range = get_disparity_range(disparity_map);
  } catch ( std::exception & /*e*/ ) {
    // There was no good pixels
    disp_range = BBox2();
  }
  ImageView<PixelRGB<float> > horz = normalize(clamp(select_channel(disparity_map,0), disp_range.min().x(), disp_range.max().x() ));
  ImageView<PixelRGB<float> > vert = normalize(clamp(select_channel(disparity_map,1), disp_range.min().y(), disp_range.max().y() ));

  for (unsigned i = 0; i < nominal_blocks.size(); ++i) {
    draw_bbox(horz, nominal_blocks[i]);
    draw_bbox(vert, nominal_blocks[i]);
  }

  write_image( m_debug_prefix+current_level.str()+"-H.jpg", horz);
  write_image( m_debug_prefix+current_level.str()+"-V.jpg", vert);

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
BBox2 PyramidCorrelator::compute_matching_blocks(BBox2i const& nominal_block, BBox2 search_range,
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
std::vector<BBox2>
PyramidCorrelator::compute_search_ranges(ImageView<PixelMask<Vector2f> > const& prev_disparity_map,
                                         std::vector<BBox2i> const& nominal_blocks) {
  std::vector<BBox2> search_ranges(nominal_blocks.size());
  std::vector<bool> is_good(nominal_blocks.size());

  // Step 1: compute the search ranges from the disparity map.
  for (unsigned i = 0; i < nominal_blocks.size(); ++i) {
    ImageViewRef<PixelMask<Vector2f> > crop_disparity = crop(prev_disparity_map,(nominal_blocks[i])/2);
    if (count_valid_pixels(crop_disparity) > 20 * 20) {
      search_ranges[i] = get_disparity_range(crop_disparity);
      is_good[i]=true;
    } else {
      // Not enough good pixels available
      search_ranges[i] = BBox2i(0, 0, 0, 0);
      is_good[i]=false;
    }
  }

  // Step 2: adjust the search range or fix it if it came from a
  // block with zero valid pixels.
  std::list<unsigned> was_corrected;
  for (unsigned r = 0; r < nominal_blocks.size(); ++r) {

    // Pick a reasonable search range based on neighboring tiles rather
    // than skipping out entirely here.
    if ( !is_good[r] ) {
      // Find adjancent tiles and use the search ranges of adjacent
      // tiles to seed us with a reasonable search range of our own.
      //
      // To do this test, we expand this block by 1 pixel in every
      // direction and then test for intersections with other blocks.
      BBox2i this_bbox = nominal_blocks[r];
      this_bbox.expand(1);
      bool found_one_match = false;
      for (unsigned b = 0; b < nominal_blocks.size(); ++b)
        if ( is_good[b] && this_bbox.intersects(nominal_blocks[b])) {
          if (found_one_match) {
            search_ranges[r].grow(search_ranges[b]);
          } else {
            found_one_match = true;
            was_corrected.push_back( r );
            search_ranges[r] = search_ranges[b];
          }
        }
    }
  }

  // Step 2 1/2: Mark good the search ranges that have been corrected
  for ( std::list<unsigned>::const_iterator iter = was_corrected.begin();
        iter != was_corrected.end(); iter++ )
    is_good[*iter] = true;

  // Step 3: scale up the search range for the next pyramid
  // level and pad it here.
  for (unsigned r = 0; r < nominal_blocks.size(); ++r)
    if ( is_good[r] ) {
      search_ranges[r] *= 2;
      search_ranges[r].expand(2);
    }

  return search_ranges;
}
