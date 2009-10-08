// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/Plate/KmlPlateManager.h>

using namespace vw::platefile;

std::vector<PlateManager::TileInfo> 
KmlPlateManager::kml_image_tiles( BBox2i const& image_bbox, 
                                  int32 const resolution,
                                  int32 const tile_size) {
  std::vector<TileInfo> result;

  // There's no point in starting the search before there is good
  // image data, so we adjust the start point here.
  int32 minx = floor(image_bbox.min().x() / tile_size) * tile_size;
  int32 miny = floor(image_bbox.min().y() / tile_size) * tile_size;
  int x = minx / tile_size;
  int y = miny / tile_size;

  // Iterate over the bounding boxes in the entire KML space...
  int curx = minx;
  int cury = miny;
  while (cury < image_bbox.max().y()) {
    while (curx <= image_bbox.max().x()) {
      
      TileInfo be(x, y, BBox2i(curx, cury, tile_size, tile_size));
      
      // ...but only add bounding boxes that overlap with the image.
      if (image_bbox.intersects(be.bbox))
        result.push_back(be);
      
      curx += tile_size;
      ++x;
    }
    curx = minx;
    x = minx / tile_size;
    cury += tile_size;
    ++y;
  }
  return result;
}



vw::ImageView<vw::PixelRGBA<vw::uint8> > 
KmlPlateManager::load_tile( int32 level, int32 x, int32 y ) {
    
  // First we try to access the indexrecord for this tile.  If that
  // fails, then we must be trying to access a node in the tree that
  // simply doesn't exist.  In this case, we create a totally empty
  // tile with zero pixels and return it.
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
    ImageView<PixelRGBA<uint8> > super(2*tile_size, 2*tile_size);
        
    // Iterate over the children, gathering them and (recursively)
    // regenerating them if necessary.
    for( int j=0; j<2; ++j ) {
      for( int i=0; i<2; ++i ) {
        ImageView<PixelRGBA<uint8> > child = load_tile(level+1,2*x+i,2*y+j);
        if( child ) crop(super,tile_size*i,tile_size*j,tile_size,tile_size) = child;	    
      }
    }
        
    // We subsample after blurring with a standard 2x2 box filter.
    std::vector<float> kernel(2);
    kernel[0] = kernel[1] = 0.5;
    
    tile = subsample( separable_convolution_filter( super, 
                                                    kernel, 
                                                    kernel, 
                                                    1, 1,
                                                    ConstantEdgeExtension() ), 2);
    
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

  return tile;
}
