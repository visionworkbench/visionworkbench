// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Plate/TileManipulation.h>

// Given a bbox, returns a list of smaller bboxes that perfectly
// tile the space of the larger bbox.
std::list<vw::BBox2i> vw::platefile::bbox_tiles(vw::BBox2i const& bbox, int width, int height) {
  std::list<vw::BBox2i> bboxes;

  vw::int32 j_offset = bbox.min().y();
  while ( j_offset < bbox.max().y() ) {
    vw::int32 j_dim = (bbox.max().y() - j_offset) < height ? (bbox.max().y() - j_offset) : height;
    vw::int32 i_offset = bbox.min().x();
    while ( i_offset < bbox.max().x() ) {
      vw::int32 i_dim = (bbox.max().x() - i_offset) < width ? (bbox.max().x() - i_offset) : width;
      bboxes.push_back(vw::BBox2i(i_offset,j_offset,i_dim,j_dim));
      i_offset += i_dim;
    }
    j_offset += j_dim;
  }
  return bboxes;
}
