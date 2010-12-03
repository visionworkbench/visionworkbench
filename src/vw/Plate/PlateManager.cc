// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Plate/PlateManager.h>
#include <vw/Plate/TileManipulation.h>
#include <vw/Image/Transform.h>
#include <vw/Image/Filter.h>
#include <vw/Plate/PlateCarreePlateManager.h>
#include <vw/Plate/PolarStereoPlateManager.h>
#include <vw/Plate/ToastPlateManager.h>

using namespace vw;
using namespace vw::platefile;

namespace {
  class BresenhamLine {
    vw::int32 x0, y0, x1, y1;
    vw::int32 x, y;
    bool steep;
    vw::int32 deltax, deltay, error, ystep;
  public:
    BresenhamLine( vw::Vector2i const& start, vw::Vector2i const& stop ) :
    x0(start[0]), y0(start[1]), x1(stop[0]), y1(stop[1]) {
      steep = abs(y1-y0) > abs(x1-x0);
      if (steep) {
        std::swap(x0,y0);
        std::swap(x1,y1);
      }
      if ( x0 > x1 ) {
        std::swap(x0,x1);
        std::swap(y0,y1);
      }
      deltax = x1 - x0;
      deltay = abs(y1-y0);
      error = deltax / 2;
      ystep = y0 < y1 ? 1 : -1;
      x = x0; y = y0;
    }

    vw::Vector2i operator*() const {
      if (steep)
        return vw::Vector2i(y,x);
      else
        return vw::Vector2i(x,y);
    }

    void operator++() {
      x++;
      error -= deltay;
      if ( error < 0 ) {
        y += ystep;
        error += deltax;
      }
    }

    bool is_good() const { return x < x1; }
  };
}

// mipmap() generates mipmapped (i.e. low resolution) tiles in the mosaic.
template <class PixelT>
void PlateManager<PixelT>::mipmap(int starting_level, BBox2i const& bbox,
                                  TransactionOrNeg read_transaction_id,
                                  Transaction output_transaction_id,
                                  bool preblur,
                                  const ProgressCallback &progress_callback,
                                  int stopping_level) const {

  // Adjust the size of the bbox for the first mipmapping level, which
  // is one level up from the starting_level.  We also pad it
  // slightly here, since we want to make sure we don't trucate at
  // tiles when we divide by 2.
  BBox2i level_bbox = bbox;
  level_bbox.min().x() =
    boost::numeric_cast<int32>(floor( float(level_bbox.min().x()) / 2.0 ));
  level_bbox.min().y() =
    boost::numeric_cast<int32>(floor( float(level_bbox.min().y()) / 2.0 ));
  level_bbox.max().x() =
    boost::numeric_cast<int32>(ceil( float(level_bbox.max().x()+1) / 2.0 ));
  level_bbox.max().y() =
    boost::numeric_cast<int32>(ceil( float(level_bbox.max().y()+1) / 2.0 ));

  // Compute the range of progress for this SubProgressCallback.  The
  // geometric sum below probably could be computed more easily, but
  // I'm tired and a little lazy at the moment.
  float total_num_tiles = 0.0;
  float sum_denom = 4.0;
  for ( int level = starting_level-1; level >= 0; --level) {
    total_num_tiles += float(bbox.width() * bbox.height()) / sum_denom;
    sum_denom *= 4.0f;
  }

  float current_num_tiles = 0;
  sum_denom = 4.0;
  float prev_num_tiles = 0;
  for ( int level = starting_level-1;
        level >= (stopping_level >= 0 ? stopping_level : 0); --level) {

    // Do a little progress callback math.
    current_num_tiles += float(bbox.width() * bbox.height()) / sum_denom;
    sum_denom *= 4.0f;
    SubProgressCallback sub_progress(progress_callback,
                                     prev_num_tiles / total_num_tiles,
                                     current_num_tiles / total_num_tiles);
    prev_num_tiles = current_num_tiles;

    // Subdivide the bbox into smaller workunits if necessary.
    // This helps to keep operations efficient.
    std::list<BBox2i> tile_workunits = bbox_tiles(level_bbox, 16, 16);
    int prog_counter = 0;
    for ( std::list<BBox2i>::iterator iter = tile_workunits.begin(); iter != tile_workunits.end(); ++iter) {
      SubProgressCallback sub_sub_progress(sub_progress,
                                           float(prog_counter)/float(tile_workunits.size()),
                                           float(prog_counter+1)/float(tile_workunits.size()));
      prog_counter++;

      // The original bbox passed into the mipmap function only serves
      // as a hint to where we might find tiles.  Some tiles will be
      // missing due to transparency, and it's also possible that a
      // very large number of tiles might be missing if the bounding
      // box happens to cross one edge of the mosaic.  In both cases,
      // the call to valid_tiles() here fetches the tiles that
      // actually contain good data.
      BBox2i parent_region = *iter;
      parent_region.min() *= 2;
      parent_region.max() *= 2;
      std::list<TileHeader> valid_tile_records =
        m_platefile->search_by_region(level+1, parent_region,
                                      read_transaction_id, read_transaction_id, 1);

      // Debugging:
      // std::cout << "Queried for valid_tiles in " << parent_region << " @ " << (level+1)
      //           << "   found " << valid_tile_records.size() << "\n";

      if (!valid_tile_records.empty()) {

        // Once we compute the valid tiles at the parent level, we
        // translate that back down into valid tiles at this level.
        BBox2i trimmed_region;
        for (std::list<TileHeader>::iterator trim_iter = valid_tile_records.begin();
             trim_iter != valid_tile_records.end(); ++trim_iter) {
          trimmed_region.grow(Vector2i(trim_iter->col()/2, trim_iter->row()/2));
        }

        // Pad the bottom right edge of the bbox to make sure we get
        // tiles from the parent that only hald overlapped with the
        // child.
        trimmed_region.max().x() += 1;
        trimmed_region.max().y() += 1;

        ImageView<PixelT> img;

        float inc_amt = 1.0f/float(trimmed_region.width() * trimmed_region.height());
        for (int row = trimmed_region.min().y(); row < trimmed_region.max().y(); ++row) {
          for (int col = trimmed_region.min().x(); col < trimmed_region.max().x(); ++col) {
            this->generate_mipmap_tile(img,col,row,level,read_transaction_id,preblur);

            if (img && !is_transparent(img)) {
              vw_out(VerboseDebugMessage, "platefile") << "Writing " << col << " " << row << " @ " << level << "\n";
              this->m_platefile->write_update(img, col, row, level, output_transaction_id);
            }

            sub_sub_progress.report_incremental_progress(inc_amt);
          }
        }
      }
      sub_sub_progress.report_finished();
    }
    sub_progress.report_finished();

    // Adjust the size of the bbox for this level
    level_bbox.min().x() =
      boost::numeric_cast<int32>(floor( float(level_bbox.min().x()) / 2 ));
    level_bbox.min().y() =
      boost::numeric_cast<int32>(floor( float(level_bbox.min().y()) / 2 ));
    level_bbox.max().x() =
      boost::numeric_cast<int32>(ceil( float(level_bbox.max().x()) / 2 ));
    level_bbox.max().y() =
      boost::numeric_cast<int32>(ceil( float(level_bbox.max().y()) / 2 ));
  }
  progress_callback.report_finished();
}

template <class PixelT>
void PlateManager<PixelT>::generate_mipmap_tile(
      ImageView<PixelT>& dest, int col, int row, int level,
      TransactionOrNeg transaction_id, bool preblur) const
{
  // Create an image large enough to store all of the child nodes
  uint32 tile_size = this->m_platefile->default_tile_size();
  ImageView<PixelT> super(2*tile_size, 2*tile_size);

  bool found = false;

  // Gather the children and stick them into a supertile.
  for( uint32 j=0; j<2; ++j ) {
    for( uint32 i=0; i<2; ++i ) {
      try {
        uint32 child_col = 2*col+i;
        uint32 child_row = 2*row+j;
        this->m_platefile->read(dest, child_col, child_row, level+1, transaction_id, true); // exact_transaction
        crop(super,tile_size*i,tile_size*j,tile_size,tile_size) = dest;
        found = true;
      } catch (TileNotFoundErr &e) { /*Do Nothing*/ }
    }
  }

  if (!found) {
    dest.reset();
    return;
  }

  // We subsample after blurring with a standard 2x2 box filter.
  std::vector<float> kernel(2);
  kernel[0] = kernel[1] = 0.5;

  if (preblur)
    dest = subsample( separable_convolution_filter( super, kernel, kernel, 1, 1, ConstantEdgeExtension() ), 2);
  else
    dest = subsample( super, 2 );
}


// Calculate the TileInfo objects of the tiles affected by
// transforming an image of size 'image_size' with transform 'tx'.
template <class PixelT>
void PlateManager<PixelT>::affected_tiles(BBox2i const& image_size,
                                          TransformRef const& tx, int tile_size,
                                          int /*lvl*/, std::list<TileInfo>& tiles ) const {
  BBox2f pyramid_px_bbox = tx.forward_bbox(image_size);
  tiles.clear();

  int32 min_tile_x =
    boost::numeric_cast<int32>(floor(pyramid_px_bbox.min().x() / tile_size ));
  int32 min_tile_y =
    boost::numeric_cast<int32>(floor(pyramid_px_bbox.min().y() / tile_size ));
  int32 max_tile_x =
    boost::numeric_cast<int32>(ceil(pyramid_px_bbox.max().x()  / tile_size ));
  int32 max_tile_y =
    boost::numeric_cast<int32>(ceil(pyramid_px_bbox.max().y()  / tile_size ));

  for ( int32 tile_x = min_tile_x; tile_x < max_tile_x; tile_x++ ) {
    for ( int32 tile_y = min_tile_y; tile_y < max_tile_y; tile_y++ ) {
      TileInfo tile( tile_x, tile_y,
                     BBox2i(tile_x*tile_size,tile_y*tile_size,
                            tile_size, tile_size) );

      // See if it intersects
      bool intersects = false;

      // We're testing a pattern that looks like
      //      +--+--+
      //      |\ | /|
      //      ---+---
      //      |/ | \|
      //      +--+--+

      // Testing \'
      BresenhamLine line( tile.bbox.min(), tile.bbox.max() );
      while ( line.is_good() && !intersects ) {
        if ( image_size.contains( tx.reverse( *line ) ) )
          intersects = true;
        ++line;
      }

      // Testing /
      line = BresenhamLine( tile.bbox.min()+Vector2i(tile.bbox.width(),0),
                            tile.bbox.max()-Vector2i(tile.bbox.width(),0) );
      while ( line.is_good() && !intersects ) {
        if ( image_size.contains( tx.reverse( *line ) ) )
          intersects = true;
        ++line;
      }

      // Testing |
      line = BresenhamLine( tile.bbox.min()+Vector2i(tile.bbox.width()/2,0),
                            tile.bbox.max()-Vector2i(tile.bbox.width()/2,0) );
      while ( line.is_good() && !intersects ) {
        if ( image_size.contains( tx.reverse( *line ) ) )
          intersects = true;
        ++line;
      }

      // Testing -
      line = BresenhamLine( tile.bbox.min()+Vector2i(0,tile.bbox.height()/2),
                            tile.bbox.max()-Vector2i(0,tile.bbox.height()/2) );
      while ( line.is_good() && !intersects ) {
        if ( image_size.contains( tx.reverse( *line ) ) )
          intersects = true;
        ++line;
      }

      // Testing far left
      line = BresenhamLine( Vector2(0,0), Vector2(0,tile.bbox.height()) );
      while ( line.is_good() && !intersects ) {
        if ( image_size.contains( tx.reverse( *line ) ) )
          intersects = true;
        ++line;
      }

      // Testing far right
      line = BresenhamLine( Vector2(tile.bbox.width()-1,0),
                            Vector2(tile.bbox.width()-1,tile.bbox.height()) );
      while ( line.is_good() && !intersects ) {
        if ( image_size.contains( tx.reverse( *line ) ) )
          intersects = true;
        ++line;
      }

      // Testing top
      line = BresenhamLine( Vector2(0,0), Vector2(tile.bbox.width(),0) );
      while ( line.is_good() && !intersects ) {
        if ( image_size.contains( tx.reverse( *line ) ) )
          intersects = true;
        ++line;
      }

      // Testing bottom
      line = BresenhamLine( Vector2(0,tile.bbox.height()-1),
                            Vector2(tile.bbox.width(),tile.bbox.height()-1) );
      while ( line.is_good() && !intersects ) {
        if ( image_size.contains( tx.reverse( *line ) ) )
          intersects = true;
        ++line;
      }

      // If it intersects, its worth rendering
      if ( intersects )
        tiles.push_back( tile );
    }
  }
}

template <class PixelT>
PlateManager<PixelT>*
PlateManager<PixelT>::make( std::string const& mode,
                            boost::shared_ptr<PlateFile> platefile ) {
  std::string mode_l = boost::to_lower_copy( mode );

  if ( mode_l == "equi" ) {
    return new PlateCarreePlateManager<PixelT>(platefile);
  } else if ( mode_l == "toast" ) {
    return new ToastPlateManager<PixelT>(platefile);
  } else if ( mode_l == "polar" ) {
    return new PolarStereoPlateManager<PixelT>(platefile);
  } else {
    vw_throw( ArgumentErr() << "Unknown option: \"" << mode << "\".\n" );
  }
}

// Explicit template instantiation
namespace vw {
namespace platefile {

#define VW_INSTANTIATE_PLATE_MANAGER_TYPES(PIXELT)                             \
  template void                                                                \
  PlateManager<PIXELT >::mipmap(int starting_level, BBox2i const& bbox,        \
                                TransactionOrNeg transaction_id,               \
                                Transaction output_transaction_id,             \
                                bool preblur,                                  \
                                const ProgressCallback &progress_callback,     \
                                int stopping_level) const;                     \
  template void                                                                \
  PlateManager<PIXELT >::affected_tiles(BBox2i const& image_size,              \
                                        TransformRef const& tx, int tile_size, \
                                        int level, std::list<TileInfo>& tiles ) const; \
  template PlateManager<PIXELT >*                                              \
  PlateManager<PIXELT >::make( std::string const& mode,                        \
                               boost::shared_ptr<PlateFile> platefile );

  VW_INSTANTIATE_PLATE_MANAGER_TYPES(PixelGrayA<uint8>)
  VW_INSTANTIATE_PLATE_MANAGER_TYPES(PixelGrayA<int16>)
  VW_INSTANTIATE_PLATE_MANAGER_TYPES(PixelGrayA<float32>)
  VW_INSTANTIATE_PLATE_MANAGER_TYPES(PixelRGBA<uint8>)
}}
