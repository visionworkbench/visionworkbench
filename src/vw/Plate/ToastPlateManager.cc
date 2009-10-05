// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
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
  int32 minx = floor(image_bbox.min().x() / (tile_size-1)) * (tile_size-1);
  int32 miny = floor(image_bbox.min().y() / (tile_size-1)) * (tile_size-1);
  int x = minx / (tile_size-1);
  int y = miny / (tile_size-1);

  // Similarly, there's no point in searching past the end.
  int32 maxx = (floor(image_bbox.max().x() / (tile_size-1)) + 1) * (tile_size-1);
  int32 maxy = (floor(image_bbox.max().y() / (tile_size-1)) + 1) * (tile_size-1);
  
  // Iterate over the bounding boxes in the entire TOAST space...
  int curx = minx;
  int cury = miny;
  while (cury <= maxy) {
    while (curx <= maxx) {
      
      TileInfo be(x, y, BBox2i(curx, cury, tile_size, tile_size));
      
      // ...but only add bounding boxes that overlap with the image.
      if (image_bbox.contains(be.bbox))
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



vw::ImageView<vw::PixelRGBA<vw::uint8> > 
ToastPlateManager::load_tile( int32 level, int32 x, int32 y ) {
  int32 num_tiles = 1 << level;
  if( x==-1 ) {
    if( y==-1 ) {
      return load_tile(level, num_tiles-1, num_tiles-1);
    }
    if( y==num_tiles ) {
      return load_tile(level, num_tiles-1, 0);
    }
    ImageView<PixelRGBA<uint8> > tile = load_tile(level, 0, num_tiles-1-y);
    if( tile ) return rotate_180(tile);
    else return tile;
  }
  if( x==num_tiles ) {
    if( y==-1 ) {
      return load_tile(level, 0, num_tiles-1);
    }
    if( y==num_tiles ) {
      return load_tile(level, 0, 0);
    }
    ImageView<PixelRGBA<uint8> > tile = load_tile(level, num_tiles-1, num_tiles-1-y);
    if( tile ) return rotate_180(tile);
    else return tile;
  }
  if( y==-1 ) {
    ImageView<PixelRGBA<uint8> > tile = load_tile(level, num_tiles-1-x, 0);
    if( tile ) return rotate_180(tile);
    else return tile;
  }
  if( y==num_tiles ) {
    ImageView<PixelRGBA<uint8> > tile = load_tile(level, num_tiles-1-x, num_tiles-1);
    if( tile ) return rotate_180(tile);
    else return tile;
  }
    
  // TODO: Reenable cache
  //
  // // Check the cache
  // for( typename cache_t::iterator i=m_cache.begin(); i!=m_cache.end(); ++i ) {
  //   if( i->level==level && i->x==x && i->y==y ) {
  //     CacheEntry e = *i;
  //     m_cache.erase(i);
  //     m_cache.push_front(e);
  //     return e.tile;
  //   }
  // }


  // First we try to access the indexrecord for this tile.  If that
  // fails, then we must be trying to access a node in the tree that
  // simply doesn't exist.  In this case, we create a totally empty
  // tile and return it.
  ImageView<PixelRGBA<uint8> > tile;
  IndexRecord rec;
  try {
    rec = m_platefile->read_record(x, y, level);
  } catch (TileNotFoundErr &e) {
    return tile;
  }

  // If the record lookup succeded, we look at the current status of
  // the tile to decide what to do next.


  if (rec.status() == INDEX_RECORD_VALID) {

    // CASE 1 : Valid tiles can be returned without any further processing.
    m_platefile->read(tile, x, y, level);
    return tile;

  } else if (rec.status() == INDEX_RECORD_EMPTY || 
             rec.status() == INDEX_RECORD_STALE) {
    
    // CASE 2 : Empty tiles need to be regenerated from scratch.

    // Create an image large enough to store all of the child nodes
    int tile_size = m_platefile->default_tile_size();
    ImageView<PixelRGBA<uint8> > super(4*tile_size-3, 4*tile_size-3);
        
    // Iterate over the children, gathering them and (recursively)
    // regenerating them if necessary.
    for( int j=-1; j<3; ++j ) {
      for( int i=-1; i<3; ++i ) {
        ImageView<PixelRGBA<uint8> > child = load_tile(level+1,2*x+i,2*y+j);
        if( child ) crop(super,(tile_size-1)*(i+1),(tile_size-1)*(j+1),tile_size,tile_size) = child;	    
      }
    }
        
    // In the WWT implementation of TOAST the pixel centers
    // (rather than the than pixel corners) are grid-aligned, so
    // we need to use an odd-sized antialiasing kernel instead of
    // the usual 2x2 box filter.  The following 5-pixel kernel was
    // optimized to avoid the extra blurring associated with using
    // a kernel wider than 2 pixels.  Math was involved.
    std::vector<float> kernel(5);
    kernel[0] = kernel[4] = -0.0344;
    kernel[1] = kernel[3] = 0.2135;
    kernel[2] = 0.6418;

    tile = subsample( crop( separable_convolution_filter( super, 
                                                          kernel, 
                                                          kernel, 
                                                          NoEdgeExtension() ),
                            tile_size-1, tile_size-1, 2*tile_size, 2*tile_size ), 2 );

    if (rec.status() == INDEX_RECORD_STALE) {
      std::cout << "\t    [ " << x << " " << y << " @ " << level << " ] -- Regenerating tile.\n";

      if( ! is_transparent(tile) ) {
        ImageView<PixelRGBA<uint8> > old_data(tile.cols(), tile.rows());
        try {
          m_platefile->read(old_data, x, y, level);
        } catch (TileNotFoundErr &e) { 
          // Do nothing... we already have a default constructed empty image above! 
        }

        VW_ASSERT(old_data.cols() == tile.cols() && old_data.rows() == tile.rows(),
                  LogicErr() << "WritePlateFileTask::operator() -- new tile dimensions do not " 
                  << "match old tile dimensions.");
        
        vw::platefile::CompositeView<PixelRGBA<uint8> > composite;
        composite.insert(old_data, 0, 0);
        composite.insert(tile, 0, 0);
        composite.prepare();

        ImageView<PixelRGBA<uint8> > composite_tile = composite;
        if( ! is_transparent(composite_tile) ) 
          m_platefile->write(composite_tile, x, y, level);
      }

    } else {
      std::cout << "\t    [ " << x << " " << y << " @ " << level << " ] -- Creating tile.\n";
      if( ! is_transparent(tile) ) 
        m_platefile->write(tile, x, y, level);
    }
  }

  // TODO: Reenable cache
  //
  // // Save it in the cache.  The cache size of 1024 tiles was chosen
  // // somewhat arbitrarily.
  // if( m_cache.size() >= 1024 )
  //   m_cache.pop_back();
  // CacheEntry e;
  // e.level = level;
  // e.x = x;
  // e.y = y;
  // e.tile = tile;
  // m_cache.push_front(e);
  return tile;
}
