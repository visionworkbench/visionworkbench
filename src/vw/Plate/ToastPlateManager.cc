// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Plate/ToastPlateManager.h>

using namespace vw::platefile;

std::vector<ToastPlateManager::TileInfo> 
ToastPlateManager::wwt_image_tiles( BBox2i const& image_bbox, 
                                     int32 const resolution,
                                     int32 const tile_size) {
  std::vector<TileInfo> result;

  // There's no point in starting the search before there is good
  // image data, so we adjust the start point here.
  int32 minx = int(floor(image_bbox.min().x() / (tile_size-1)) * (tile_size-1));
  int32 miny = int(floor(image_bbox.min().y() / (tile_size-1)) * (tile_size-1));
  int x = minx / (tile_size-1);
  int y = miny / (tile_size-1);

  // Iterate over the bounding boxes in the entire TOAST space...
  int curx = minx;
  int cury = miny;
  while (cury < image_bbox.max().y() - 1) {
    while (curx < image_bbox.max().x() - 1) {
      
      TileInfo be(x, y, BBox2i(curx, cury, tile_size, tile_size));
      
      // ...but only add bounding boxes that overlap with the image.
      if (image_bbox.intersects(be.bbox))
        result.push_back(be);
      
      curx += (tile_size-1);
      ++x;
    }
    curx = minx;
    x = minx / (tile_size-1);
    cury += (tile_size-1);
    ++y;
  }
  return result;
}



