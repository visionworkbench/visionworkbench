// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/Plate/ToastPlateManager.h>
using namespace vw::platefile;
using namespace vw;

// // Ok. This is one of those really annoying and esoteric c++
// // template problems: we can't call load_tile_impl<> directly from
// // the CRTP superclass because the template appears in the return
// // type of this method.  Instead, we add the extra layer of
// // indirection (load_tile<>, above), which has the return value in
// // the function arguments.  
// template <class PixelT>
// ImageView<PixelT> ToastPlateManager<PixelT>::load_tile_impl( int32 level, int32 x, int32 y, 
//                                                              int transaction_id, 
//                                                              int max_level ) {
//   int32 num_tiles = 1 << level;
//   if( x==-1 ) {
//     if( y==-1 ) {
//       return load_tile_impl(level, num_tiles-1, num_tiles-1, 
//                                     transaction_id, max_level);
//     }
//     if( y==num_tiles ) {
//       return load_tile_impl(level, num_tiles-1, 0, 
//                                     transaction_id, max_level);
//     }
//     ImageView<PixelT> tile = load_tile_impl(level, 0, num_tiles-1-y, 
//                                                     transaction_id, max_level);
//     if( tile ) return rotate_180(tile);
//     else return tile;
//   }
//   if( x==num_tiles ) {
//     if( y==-1 ) {
//       return load_tile_impl(level, 0, num_tiles-1, 
//                                     transaction_id, max_level);
//     }
//     if( y==num_tiles ) {
//       return load_tile_impl(level, 0, 0, transaction_id, max_level);
//     }
//     ImageView<PixelT> tile = load_tile_impl(level, num_tiles-1, num_tiles-1-y, 
//                                                     transaction_id, max_level);
//     if( tile ) return rotate_180(tile);
//     else return tile;
//   }
//   if( y==-1 ) {
//     ImageView<PixelT> tile = load_tile_impl(level, num_tiles-1-x, 0, 
//                                                     transaction_id, max_level);
//     if( tile ) return rotate_180(tile);
//     else return tile;
//   }
//   if( y==num_tiles ) {
//     ImageView<PixelT> tile = load_tile_impl(level, num_tiles-1-x, num_tiles-1, 
//                                                     transaction_id, max_level);
//     if( tile ) return rotate_180(tile);
//     else return tile;
//   }
   
//   // Check the cache
//   for( typename cache_t::iterator i=m_cache.begin(); i!=m_cache.end(); ++i ) {
//     if( i->level==level && i->x==x && i->y==y && i->transaction_id == transaction_id ) {
//       CacheEntry e = *i;
//       m_cache.erase(i);
//       m_cache.push_front(e);
//       return e.tile;
//     }
//   }

//   // First we try to access the indexrecord for this tile.  If
//   // that fails, then we must be trying to access a node in the
//   // tree that simply doesn't exist.  In this case, we create a
//   // totally empty tile and return it.  Note that we are looking
//   // for the tile using the transaction_id here because we
//   // want to read tiles that have just been inserted in this round
//   // of processing.
//   ImageView<PixelT> tile;
//   IndexRecord rec;
//   bool found_record = true;
//   try {
//     rec = m_platefile->read_record(x, y, level, transaction_id, true);  // exact_transaction_match = true
//   } catch (TileNotFoundErr &e) {
//     found_record = false;
//   }

//   // STEP 1 : Combining and subsampling to create a new tile
//   // 

//   // If the record lookup succeded, we look at the current status of
//   // the tile to decide what to do next.
//   if (found_record && rec.status() == INDEX_RECORD_VALID) {

//     // CASE 1 : Valid tiles can be returned without any further processing.
//     m_platefile->read(tile, x, y, level, transaction_id, true);          // exact_transaction_match = true

//   } else if (found_record && (rec.status() == INDEX_RECORD_EMPTY || 
//                               rec.status() == INDEX_RECORD_STALE)) {
    
//     // CASE 2 : Empty tiles need to be regenerated from scratch.

//     // Create an image large enough to store all of the child nodes
//     int tile_size = m_platefile->default_tile_size();
//     ImageView<PixelT> super(4*tile_size-3, 4*tile_size-3);
        
//     // Iterate over the children, gathering them and (recursively)
//     // regenerating them if necessary.
//     for( int j=-1; j<3; ++j ) {
//       for( int i=-1; i<3; ++i ) {
//         ImageView<PixelT> child = load_tile_impl(level+1,2*x+i,2*y+j,transaction_id, max_level);
//         if( child ) crop(super,(tile_size-1)*(i+1),(tile_size-1)*(j+1),tile_size,tile_size) = child;	    
//       }
//     }
        
//     // In the WWT implementation of TOAST the pixel centers
//     // (rather than the than pixel corners) are grid-aligned, so
//     // we need to use an odd-sized antialiasing kernel instead of
//     // the usual 2x2 box filter.  The following 5-pixel kernel was
//     // optimized to avoid the extra blurring associated with using
//     // a kernel wider than 2 pixels.  Math was involved.
//     std::vector<float> kernel(5);
//     kernel[0] = kernel[4] = -0.0344;
//     kernel[1] = kernel[3] = 0.2135;
//     kernel[2] = 0.6418;

//     ImageView<PixelT> new_tile = subsample( crop( separable_convolution_filter( super, 
//                                                           kernel, 
//                                                           kernel, 
//                                                           NoEdgeExtension() ),
//                             tile_size-1, tile_size-1, 2*tile_size, 2*tile_size ), 2 );

//     // STEP 2 : MIPMAPPING
//     // 
//     // The composite_mosaic_tile() function looks for any tiles at
//     // equal or lower resolution in the mosaic, and composites this
//     // tile on top of those tiles, supersampling the low-res tile if
//     // necessary.
//     vw_out(0) << "\t    [ " << x << " " << y << " @ " << level << " ] -- Mipmapping tile.\n";
//     tile = m_compositor.composite_mosaic_tile(m_platefile, new_tile, x, y, 
//                                               level, max_level, transaction_id);
//   }

//   // Save the tile in the cache.  The cache size of 1024 tiles was chosen
//   // somewhat arbitrarily.
//   if( m_cache.size() >= 1024 )
//     m_cache.pop_back();
//   CacheEntry e;
//   e.level = level;
//   e.x = x;
//   e.y = y;
//   e.transaction_id = transaction_id;
//   e.tile = tile;
//   m_cache.push_front(e);

//   return tile;
// }


template <class PixelT>
ImageView<PixelT> vw::platefile::ToastPlateManager<PixelT>::composite_child_tile(int &num_composited, 
                                                                                 int x, int y, 
                                                                                 int level, 
                                                                                 int read_transaction_id,
                                                                                 int write_transaction_id) const {

  int32 num_tiles = 1 << level;
  if( x==-1 ) {
    if( y==-1 ) {
      return composite_child_tile(num_composited, num_tiles-1, num_tiles-1, level,
                                  read_transaction_id, write_transaction_id);
    }
    if( y==num_tiles ) {
      return composite_child_tile(num_composited, num_tiles-1, 0, level,
                                  read_transaction_id, write_transaction_id);
    }
    ImageView<PixelT> tile = composite_child_tile(num_composited, 0, num_tiles-1-y, level, 
                                                  read_transaction_id, write_transaction_id);
    if( tile ) return rotate_180(tile);
    else return tile;
  }
  if( x==num_tiles ) {
    if( y==-1 ) {
      return composite_child_tile(num_composited, 0, num_tiles-1, level,
                                  read_transaction_id, write_transaction_id);
    }
    if( y==num_tiles ) {
      return composite_child_tile(num_composited, 0, 0, level,
                                  read_transaction_id, write_transaction_id);
    }
    ImageView<PixelT> tile = composite_child_tile(num_composited, num_tiles-1, num_tiles-1-y, level,
                                                  read_transaction_id, write_transaction_id);
    if( tile ) return rotate_180(tile);
    else return tile;
  }
  if( y==-1 ) {
    ImageView<PixelT> tile = composite_child_tile(num_composited, num_tiles-1-x, 0, level,
                                                  read_transaction_id, write_transaction_id);
    if( tile ) return rotate_180(tile);
    else return tile;
  }
  if( y==num_tiles ) {
    ImageView<PixelT> tile = composite_child_tile(num_composited, num_tiles-1-x, num_tiles-1, level,
                                                  read_transaction_id, write_transaction_id);
    if( tile ) return rotate_180(tile);
    else return tile;
  }
   
  // Check the cache
  //
  // XXX : disabled for now 
  //
  //   for( typename cache_t::iterator i=m_cache.begin(); i!=m_cache.end(); ++i ) {
  //     if( i->level==level && i->x==x && i->y==y ) {
  //       CacheEntry e = *i;
  //       m_cache.erase(i);
  //       m_cache.push_front(e);
  //       return e.tile;
  //     }
  //   }

  // STEP 1 : Combining and subsampling to create a new tile
  // 

  // If the record lookup succeded, we look at the current status of
  // the tile to decide what to do next.
  std::cout << "Attempting to access tiles at " << x << " " << y << " " << level << "\n";
  std::list<ImageView<PixelT> > tiles;
  try {
    std::list<TileHeader> headers = m_platefile->multi_read(tiles, 
                                                            x, y, level, 
                                                            read_transaction_id,
                                                            write_transaction_id);
    num_composited = tiles.size();
    std::cout << "There were " << num_composited << " tiles here.\n";

    if (tiles.size() == 1) {

      // First case: there is just one tile.  No need to composite.
      // Just return the tile.
      return *(tiles.begin());

    } else if (tiles.size() > 1) {

      // Second case: multiple tiles.  Composite them together in the
      // order that they were returned.  The _reversed_ list of images
      // should be in order of increasing transaction id.
      tiles.reverse();
      vw::mosaic::ImageComposite<PixelT> composite;
      for (typename std::list<ImageView<PixelT> >::iterator tile_iter = tiles.begin();
           tile_iter != tiles.end(); ++tile_iter) {
        composite.insert(*tile_iter, 0, 0);
      }
      composite.set_draft_mode( true );
      composite.prepare();
      ImageView<PixelT> result_tile = composite;
      return result_tile;

    } else {

      // Third case: no tiles were returned.  Do nothing.

    }

  } catch (TileNotFoundErr &e) {
    
    std::cout << "Missing tile: " << e.what() << "\n";
    // Do nothing here.  There are no tiles to composite.

  }

  // Save the tile in the cache.  The cache size of 1024 tiles was chosen
  // somewhat arbitrarily.
  //
  // XXX : disabled for now
  //   if( m_cache.size() >= 1024 )
  //     m_cache.pop_back();
  //   CacheEntry e;
  //   e.level = level;
  //   e.x = x;
  //   e.y = y;
  //   e.transaction_id = transaction_id;
  //   e.tile = tile;
  //   m_cache.push_front(e);

  // If we reach this point, then there was no data to composite.
  // Return an empty image.
  return ImageView<PixelT>();
}

template <class PixelT>
void vw::platefile::ToastPlateManager<PixelT>::regenerate_tile(int col, int row, 
                                                               int level, 
                                                               int start_transaction_id, 
                                                               int end_transaction_id, 
                                                               int write_transaction_id) const {
  std::cout << "\t\tGenerating " << col << " " << row << "\n";

  // Create an image large enough to store all of the child nodes
  int tile_size = m_platefile->default_tile_size();
  ImageView<PixelT> super(4*tile_size-3, 4*tile_size-3);
        
  // Iterate over the children, gathering them and (recursively)
  // regenerating them if necessary.
  for( int j=-1; j<3; ++j ) {
    for( int i=-1; i<3; ++i ) {
      int num_composited;
      ImageView<PixelT> child = composite_child_tile(num_composited, 2*col+i, 2*row+j, level+1, 
                                                     start_transaction_id, end_transaction_id);
      if( num_composited > 0 && child) {
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
    std::cout << "\twriting tile\n";
    m_platefile->write(new_tile, col, row, level, write_transaction_id);
  }
}


// Explicit template instatiation
namespace vw { 
namespace platefile {

  // template 
  // ImageView<PixelGrayA<uint8> > ToastPlateManager<PixelGrayA<uint8> >::load_tile_impl( int32 level, int32 x, int32 y, 
  //                                                                                      int transaction_id, 
  //                                                                                      int max_level);
  // template
  // ImageView<PixelRGBA<uint8> > ToastPlateManager<PixelRGBA<uint8> >::load_tile_impl( int32 level, int32 x, int32 y, 
  //                                                                                    int transaction_id,
  //                                                                                    int max_level );
  
  template 
  void vw::platefile::ToastPlateManager<PixelGrayA<uint8> >::regenerate_tile(int col, int row, 
                                                                             int level, 
                                                                             int start_transaction_id, 
                                                                             int end_transaction_id,
                                                                             int write_transaction_id) const;
  template 
  void vw::platefile::ToastPlateManager<PixelRGBA<uint8> >::regenerate_tile(int col, int row, 
                                                                            int level, 
                                                                            int start_transaction_id, 
                                                                            int end_transaction_id,
                                                                            int write_transaction_id) const;

  template 
  ImageView<PixelGrayA<uint8> > vw::platefile::ToastPlateManager<PixelGrayA<uint8> >::composite_child_tile(int &num_composited, 
                                                                                                           int col, int row, 
                                                                                                           int level, 
                                                                                                           int start_transaction_id,
                                                                                                           int end_transaction_id) const;
  template 
  ImageView<PixelRGBA<uint8> > vw::platefile::ToastPlateManager<PixelRGBA<uint8> >::composite_child_tile(int &num_composited, 
                                                                                                         int col, int row, 
                                                                                                         int level, 
                                                                                                         int start_transaction_id,
                                                                                                         int end_transaction_id) const;

  
}}
