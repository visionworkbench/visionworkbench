// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Plate/ToastPlateManager.h>
using namespace vw::platefile;
using namespace vw;

std::vector<vw::platefile::TileInfo> 
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


// Ok. This is one of those really annoying and esoteric c++
// template problems: we can't call load_tile_impl<> directly from
// the CRTP superclass because the template appears in the return
// type of this method.  Instead, we add the extra layer of
// indirection (load_tile<>, above), which has the return value in
// the function arguments.  
template <class PixelT>
ImageView<PixelT> ToastPlateManager::load_tile_impl( int32 level, int32 x, int32 y, 
                                                     int write_transaction_id, 
                                                     int max_depth ) {
  int32 num_tiles = 1 << level;
  if( x==-1 ) {
    if( y==-1 ) {
      return load_tile_impl<PixelT>(level, num_tiles-1, num_tiles-1, 
                                    write_transaction_id, max_depth);
    }
    if( y==num_tiles ) {
      return load_tile_impl<PixelT>(level, num_tiles-1, 0, 
                                    write_transaction_id, max_depth);
    }
    ImageView<PixelT> tile = load_tile_impl<PixelT>(level, 0, num_tiles-1-y, 
                                                    write_transaction_id, max_depth);
    if( tile ) return rotate_180(tile);
    else return tile;
  }
  if( x==num_tiles ) {
    if( y==-1 ) {
      return load_tile_impl<PixelT>(level, 0, num_tiles-1, 
                                    write_transaction_id, max_depth);
    }
    if( y==num_tiles ) {
      return load_tile_impl<PixelT>(level, 0, 0, write_transaction_id, max_depth);
    }
    ImageView<PixelT> tile = load_tile_impl<PixelT>(level, num_tiles-1, num_tiles-1-y, 
                                                    write_transaction_id, max_depth);
    if( tile ) return rotate_180(tile);
    else return tile;
  }
  if( y==-1 ) {
    ImageView<PixelT> tile = load_tile_impl<PixelT>(level, num_tiles-1-x, 0, 
                                                    write_transaction_id, max_depth);
    if( tile ) return rotate_180(tile);
    else return tile;
  }
  if( y==num_tiles ) {
    ImageView<PixelT> tile = load_tile_impl<PixelT>(level, num_tiles-1-x, num_tiles-1, 
                                                    write_transaction_id, max_depth);
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


  // First we try to access the indexrecord for this tile.  If
  // that fails, then we must be trying to access a node in the
  // tree that simply doesn't exist.  In this case, we create a
  // totally empty tile and return it.  Note that we are looking
  // for the tile using the write_transaction_id here because we
  // want to read tiles that have just been inserted in this round
  // of processing.
  ImageView<PixelT> tile;
  IndexRecord rec;
  try {
    rec = m_platefile->read_record(x, y, level, write_transaction_id);
    //    std::cout << "XXX : " << x << " " << y << " @ " << level << " \n" << rec.DebugString() << "\n";
  } catch (TileNotFoundErr &e) {
    return tile;
  }

  // STEP 1 : Combining and subsampling to create a new tile
  // 

  // If the record lookup succeded, we look at the current status of
  // the tile to decide what to do next.
  if (rec.status() == INDEX_RECORD_VALID) {

    // CASE 1 : Valid tiles can be returned without any further processing.
    m_platefile->read(tile, x, y, level, write_transaction_id);
    return tile;

  } else if (rec.status() == INDEX_RECORD_EMPTY || 
             rec.status() == INDEX_RECORD_STALE) {
    
    // CASE 2 : Empty tiles need to be regenerated from scratch.

    // Create an image large enough to store all of the child nodes
    int tile_size = m_platefile->default_tile_size();
    ImageView<PixelT> super(4*tile_size-3, 4*tile_size-3);
        
    // Iterate over the children, gathering them and (recursively)
    // regenerating them if necessary.
    for( int j=-1; j<3; ++j ) {
      for( int i=-1; i<3; ++i ) {
        ImageView<PixelT> child = load_tile_impl<PixelT>(level+1,2*x+i,2*y+j,
                                                         write_transaction_id, max_depth);
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


    // STEP 2 : MIPMAPPING
    // 
    // The composite_mosaic_tile() function looks for any tiles at
    // equal or lower resolution in the mosaic, and composites this
    // tile on top of those tiles, supersampling the low-res tile if
    // necessary.
    return composite_mosaic_tile(m_platefile, tile, x, y, level, write_transaction_id);
  }
  //   // Next we need to assess whether there is a tile
  //   try {
  //     IndexRecord old_rec = m_platefile->read_record(x, y, level, write_transaction_id-1);

  //     // Mipmapping must happen strictly in order of transaction
  //     // ID.  Here we check to see if the underlying data tile is
  //     // available yet.  It may be "locked" by another mipmapping
  //     // session running on a different instance of image2plate.
  //     // If that's the case, we sleep for a short time and then
  //     // try again.
  //     while (old_rec.status() == INDEX_RECORD_LOCKED) {
  //       vw_out(0) << "WAITING for tile [ " << x << " " << y << " @ " << level << "]\n";
  //       sleep(5.0);
  //       old_rec = m_platefile->read_record(x, y, level, write_transaction_id-1);
  //     }
  //   } catch (TileNotFoundErr &e) { 
  //     // Do nothing... we already have a default constructed empty image above! 
  //   }

  //   if (old_rec.status() != INDEX_RECORD_VALID) {
  //     // Regenerate a tile by overlaying it on top of the existing tile.
  //     std::cout << "\t    [ " << x << " " << y << " @ " << level << " ] -- Regenerating tile.\n";

  //     ImageView<PixelT> old_data(tile.cols(), tile.rows());
  //     m_platefile->read(old_data, x, y, level, write_transaction_id - 1);

  //     if( ! is_transparent(tile) ) {
  //       try {

  //         // Mipmapping must happen strictly in order of transaction
  //         // ID.  Here we check to see if the underlying data tile is
  //         // available yet.  It may be "locked" by another mipmapping
  //         // session running on a different instance of image2plate.
  //         // If that's the case, we sleep for a short time and then
  //         // try again.
  //         IndexRecord old_rec = m_platefile->read_record(x, y, level, write_transaction_id-1);
  //         while (old_rec.status() == INDEX_RECORD_LOCKED) {
  //           vw_out(0) << "WAITING for tile [ " << x << " " << y << " @ " << level << "]\n";
  //           sleep(5.0);
  //           old_rec = m_platefile->read_record(x, y, level, write_transaction_id-1);
  //         }
  //         m_platefile->read(old_data, x, y, level, write_transaction_id - 1);
  //       } catch (TileNotFoundErr &e) { 
  //         // Do nothing... we already have a default constructed empty image above! 
  //       }

  //       VW_ASSERT(old_data.cols() == tile.cols() && old_data.rows() == tile.rows(),
  //                 LogicErr() << "WritePlateFileTask::operator() -- new tile dimensions do not " 
  //                 << "match old tile dimensions.");
        
  //       vw::mosaic::ImageComposite<PixelT> composite;
  //       composite.insert(old_data, 0, 0);
  //       composite.insert(tile, 0, 0);
  //       composite.set_draft_mode( true );
  //       composite.prepare();

  //       ImageView<PixelT> composite_tile = composite;
  //       if( ! is_transparent(composite_tile) ) 
  //         m_platefile->write(composite_tile, x, y, level, write_transaction_id);
  //     }

  //   } else {
  //     // Create a new tile from scratch.
  //     if( ! is_transparent(tile) ) {
  //       std::cout << "\t    [ " << x << " " << y << " @ " << level << " ] -- Creating tile.\n";
  //       m_platefile->write(tile, x, y, level, write_transaction_id);
  //     }
  //   }
  // }

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

  vw_out(0) << "WARNING: Encountered a LOCKED tile while mipmapping. "
            <<"This isn\'t supposed to happen!\n";
  return tile;
}

// Explicit template instatiation
namespace vw { 
namespace platefile {


template 
ImageView<PixelGrayA<uint8> > ToastPlateManager::load_tile_impl<PixelGrayA<uint8> >( int32 level, int32 x, int32 y, 
                                                                                     int write_transaction_id, 
                                                                                     int max_depth);

template 
ImageView<PixelGrayA<float> > ToastPlateManager::load_tile_impl<PixelGrayA<float> >( int32 level, int32 x, int32 y, 
                                                                                     int write_transaction_id,
                                                                                     int max_depth);

template
ImageView<PixelRGBA<uint8> > ToastPlateManager::load_tile_impl<PixelRGBA<uint8> >( int32 level, int32 x, int32 y, 
                                                                                   int write_transaction_id,
                                                                                   int max_depth );


template 
ImageView<PixelRGBA<uint16> > ToastPlateManager::load_tile_impl<PixelRGBA<uint16> >( int32 level, int32 x, int32 y, 
                                                                                     int write_transaction_id,
                                                                                     int max_depth );


}}
