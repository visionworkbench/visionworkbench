// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/Plate/PolarStereoPlateManager.h>
using namespace vw::platefile;
using namespace vw;

void
vw::platefile::stereo_image_tiles( BBox2i const& image_bbox,
                                   cartography::GeoTransform const& geotx,
                                   int32 tile_size,
                                   std::list<TileInfo> & tiles) {
  BBox2i pyramid_px_bbox = geotx.forward_bbox(image_bbox);
  tiles.clear();

  int32 min_tile_x =
    boost::numeric_cast<int32>(floor(pyramid_px_bbox.min().x() / tile_size ));
  int32 min_tile_y =
    boost::numeric_cast<int32>(floor(pyramid_px_bbox.min().y() / tile_size ));
  int32 max_tile_x =
    boost::numeric_cast<int32>(ceil(pyramid_px_bbox.max().x()  / tile_size ));
  int32 max_tile_y =
    boost::numeric_cast<int32>(ceil(pyramid_px_bbox.max().y()  / tile_size ));

  for ( int32 tile_x = min_tile_x; tile_x <= max_tile_x; tile_x++ ) {
    for ( int32 tile_y = min_tile_y; tile_y <= max_tile_y; tile_y++ ) {
      TileInfo tile( tile_x, tile_y,
                     BBox2i(tile_x*tile_size,tile_y*tile_size,
                            tile_size, tile_size) );

      // See if it intersects
      bool intersects = false;

      // Check top boundry of bbox
      for ( int32 px_x = tile.bbox.min()[0];
            px_x < tile.bbox.max()[0]-1 && !intersects; px_x++ )
        if ( image_bbox.contains( geotx.reverse( Vector2(px_x,tile.bbox.min()[1]))))
          intersects = true;

      // Check right boundry of bbox
      for ( int32 px_y = tile.bbox.min()[1];
            px_y < tile.bbox.max()[1]-1 && !intersects; px_y++ )
        if ( image_bbox.contains( geotx.reverse( Vector2(tile.bbox.max()[0],px_y))))
          intersects = true;

      // Check bottom boundry of bbox
      for ( int32 px_x = tile.bbox.max()[0]-1;
            px_x > tile.bbox.min()[0] && !intersects; px_x-- )
        if ( image_bbox.contains( geotx.reverse( Vector2(px_x,tile.bbox.max()[1]))))
          intersects = true;

      // Check left boundry of bbox
      for ( int32 px_y = tile.bbox.max()[1]-1;
            px_y > tile.bbox.min()[1] && !intersects; px_y-- )
        if ( image_bbox.contains( geotx.reverse( Vector2(tile.bbox.min()[0],px_y))))
          intersects = true;

      // If it intersects, its worth rendering
      if ( intersects )
        tiles.push_back( tile );
    }
  }
}

