// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
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
                                    ImageView<uint8> const& dmask,
                                    BBox2i const& box) {
  if (count_valid_pixels(crop(create_mask(dmask), box)) == 0) {
    return std::vector<BBox2i>();
  }

  // No need for try catch on disparity range, because subdivide_bboxes 
  // should only be called for valid subregions
  BBox2 disp_range = get_disparity_range(crop(disparity_map, box));
  bool bad_range = (disp_range.width()+1) * (disp_range.height()+1) > 4;
  
  BBox2i subbox1 = box, subbox2 = box;
  if (box.width() > box.height()) {
    subbox1.max().x() = box.min().x() + box.width()/2;
    subbox2.min().x() = box.min().x() + box.width()/2;
  } else {
    subbox1.max().y() = box.min().y() + box.height()/2;
    subbox2.min().y() = box.min().y() + box.height()/2;
  }

  bool good_size = box.width() > m_min_subregion_dim || box.height() > m_min_subregion_dim;
  bool good_subregion = is_good_subregion(disparity_map, dmask, subbox1) && is_good_subregion(disparity_map, dmask, subbox2);
  bool kill_invalid = count_valid_pixels(crop(create_mask(dmask), subbox1)) == 0 || count_valid_pixels(crop(create_mask(dmask), subbox2)) == 0;
  
  if (kill_invalid || (bad_range && good_size && good_subregion)) {
    std::vector<BBox2i> r = subdivide_bboxes(disparity_map, dmask, subbox1);
    std::vector<BBox2i> l = subdivide_bboxes(disparity_map, dmask, subbox2);
    r.insert(r.end(), l.begin(), l.end());
    return r;
  }
  
  std::vector<vw::BBox2i> result;
  result.push_back(box);
  return result;
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
  } catch ( std::exception & e ) {
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
                                         ImageView<uint8> const& dmask,
                                         std::vector<BBox2i> const& nominal_blocks,
                                         BBox2i const& default_search_range) {
  std::vector<BBox2> search_ranges(nominal_blocks.size());

  for (unsigned i = 0; i < nominal_blocks.size(); i++) {
    if (is_good_subregion(prev_disparity_map, dmask, nominal_blocks[i])) {
      // get disparity range should never throw when is_good_subregion is true
      search_ranges[i] = get_disparity_range(crop(prev_disparity_map, nominal_blocks[i]));
      search_ranges[i].expand(2);
    } else {
      search_ranges[i] = default_search_range;
    }
  }

  return search_ranges;
}
