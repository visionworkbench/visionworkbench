// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/Plate/ToastPlateManager.h>
using namespace vw::platefile;
using namespace vw;


template <class PixelT>
ImageView<PixelT> vw::platefile::ToastPlateManager<PixelT>::fetch_child_tile(int x, int y, int level, 
                                                                             int transaction_id) const {

  //  std::cout << "Fetching child tile " << x << " " << y << "\n";

  int32 num_tiles = 1 << level;
  if( x==-1 ) {
    if( y==-1 ) {
      return fetch_child_tile(num_tiles-1, num_tiles-1, level, transaction_id);
    }
    if( y==num_tiles ) {
      return fetch_child_tile(num_tiles-1, 0, level, transaction_id);
    }
    ImageView<PixelT> tile = fetch_child_tile(0, num_tiles-1-y, level, transaction_id);
    if( tile ) return rotate_180(tile);
    else return tile;
  }
  if( x==num_tiles ) {
    if( y==-1 ) {
      return fetch_child_tile(0, num_tiles-1, level, transaction_id);
    }
    if( y==num_tiles ) {
      return fetch_child_tile(0, 0, level, transaction_id);
    }
    ImageView<PixelT> tile = fetch_child_tile(num_tiles-1, num_tiles-1-y, level, transaction_id);
    if( tile ) return rotate_180(tile);
    else return tile;
  }
  if( y==-1 ) {
    ImageView<PixelT> tile = fetch_child_tile(num_tiles-1-x, 0, level, transaction_id);
    if( tile ) return rotate_180(tile);
    else return tile;
  }
  if( y==num_tiles ) {
    ImageView<PixelT> tile = fetch_child_tile(num_tiles-1-x, num_tiles-1, level, transaction_id);
    if( tile ) return rotate_180(tile);
    else return tile;
  }
   
  // Tile accesses during mipmapping of a TOAST mosaic are highly
  // localized, so it speeds things up considerably to cache the tiles
  // as we access them. We check the cache for this tile here.
  //
  for( typename cache_t::iterator i=m_cache.begin(); i!=m_cache.end(); ++i ) {
    if( i->level==level && i->x==x && i->y==y ) {
      CacheEntry e = *i;
      m_cache.erase(i);
      m_cache.push_front(e);
      vw_out(VerboseDebugMessage, "platefile") << "Found cached tile at " << x << " " << y << " " << level << "\n";
      return e.tile;
    }
  }

  ImageView<PixelT> tile;
  try {
    
    // If the tile is not in the cache, we attempt to access it in the index.
    vw_out(VerboseDebugMessage, "platefile") << "Reading tile at " << x << " " << y << " " << level << "\n";
    m_platefile->read(tile, x, y, level, transaction_id, true); // exact_transaction_match == true

  } catch (TileNotFoundErr &e) {
    // If that fails, then there is no tile.  We return an empty image.
  }

  // Save the tile in the cache.  The cache size of 1024 tiles was chosen
  // somewhat arbitrarily.
  //
  if( m_cache.size() >= 1024 )
    m_cache.pop_back();
  CacheEntry e;
  e.level = level;
  e.x = x;
  e.y = y;
  e.transaction_id = transaction_id;
  e.tile = tile;
  m_cache.push_front(e);

  // Regardless of what happened above, we return the tile here.  If
  // the read failed, the tile will be an empty image.
  return tile;
}


template <class PixelT>
void vw::platefile::ToastPlateManager<PixelT>::generate_mipmap_tile(int col, int row, 
                                                                    int level, int transaction_id) const {

  // Create an image large enough to store all of the child nodes
  int tile_size = m_platefile->default_tile_size();
  ImageView<PixelT> super(4*tile_size-3, 4*tile_size-3);
        
  // Iterate over the children, gathering them and (recursively)
  // regenerating them if necessary.
  for( int j=-1; j<3; ++j ) {
    for( int i=-1; i<3; ++i ) {
      ImageView<PixelT> child = fetch_child_tile(2*col+i, 2*row+j, level+1, transaction_id);
      if(child) {
        crop(super,(tile_size-1)*(i+1),(tile_size-1)*(j+1),tile_size,tile_size) = child;
      }
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
  
  ImageView<PixelT> new_tile = subsample( crop( separable_convolution_filter( super, 
                                                                              kernel, 
                                                                              kernel, 
                                                                              NoEdgeExtension() ),
                                                tile_size-1, tile_size-1, 2*tile_size, 2*tile_size ), 2 );

  if (!is_transparent(new_tile)) {
    //    std::cout << "Writing " << col << " " << row << " @ " << level << "\n";
    //m_platefile->write_request();
    m_platefile->write_update(new_tile, col, row, level, transaction_id);
    //m_platefile->write_complete();
  }
}


// Explicit template instatiation
namespace vw { 
namespace platefile {

  template 
  void vw::platefile::ToastPlateManager<PixelGrayA<uint8> >::generate_mipmap_tile(int col, int row, 
                                                                                  int level, 
                                                                                  int transaction_id) const;
  template 
  void vw::platefile::ToastPlateManager<PixelGrayA<int16> >::generate_mipmap_tile(int col, int row, 
                                                                                 int level, 
                                                                                 int transaction_id) const;
  template 
  void vw::platefile::ToastPlateManager<PixelRGBA<uint8> >::generate_mipmap_tile(int col, int row, 
                                                                                 int level, 
                                                                                 int transaction_id) const;

  template 
  ImageView<PixelGrayA<uint8> > vw::platefile::ToastPlateManager<PixelGrayA<uint8> >::fetch_child_tile(int col, int row, 
                                                                                                       int level, 
                                                                                                       int transaction_id) const;
  template 
  ImageView<PixelGrayA<uint8> > vw::platefile::ToastPlateManager<PixelGrayA<int16> >::fetch_child_tile(int col, int row, 
                                                                                                       int level, 
                                                                                                       int transaction_id) const;
  template 
  ImageView<PixelRGBA<uint8> > vw::platefile::ToastPlateManager<PixelRGBA<uint8> >::fetch_child_tile(int col, int row, 
                                                                                                     int level, 
                                                                                                     int transaction_id) const;
  
}}
