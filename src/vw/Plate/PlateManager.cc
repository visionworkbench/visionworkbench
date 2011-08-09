// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Plate/PlateManager.h>
#include <vw/Plate/TileManipulation.h>
#include <vw/Plate/PlateCarreePlateManager.h>
#include <vw/Plate/PolarStereoPlateManager.h>
#include <vw/Plate/ToastPlateManager.h>
#include <vw/Plate/detail/MipmapHelpers.h>

using namespace vw;
using namespace vw::platefile;
namespace d = vw::platefile::detail;

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

namespace {
  using namespace vw;

  BBox2i move_up(const BBox2i& input) {
    return BBox2i( input.min().x()/2,   input.min().y()/2,
                  (input.width()+1)/2, (input.height()+1)/2);
  }

  typedef std::map<d::rowcol_t, std::vector<d::rowcol_t> > composite_map_t;

  template <typename PixelT>
  struct tile_cache_t : public std::map<d::rowcol_t, ImageView<PixelT> > {};

  template <typename PixelT>
  struct level_data {
    composite_map_t hdrs;
    tile_cache_t<PixelT> cache;
    void clear() {
      hdrs.clear();
      cache.clear();
    }
  };
}

namespace {

// Given a bottom tile count, approximate the total count of tiles for the pyramid
uint32 approximate_total_tiles(const BBox2i& region, uint32 bottom_tile_count) {
  // total number of tiles that need to be generated from this is a geometic
  // sum with r between 0.5 (for rectangle base) and 0.25 (for square base).
  // (geometic sum: E(n=0..inf)(a*r^n) = a / (1-r))
  float w = region.width(), h = region.height();
  //squareness will range from 1 (square) to asmptotic 0 (rect)
  float squareness = (w < h) ? w/h : h/w;
  // they tend toward rect, so range it from 0.5 to 0.3
  float r = ((1-squareness) * 0.2) + 0.3;
  // take off the top layer, we already have it
  return bottom_tile_count / (1-r);
}

template <typename PixelT>
void load_tiles(PlateFile& plate, Datastore::TileSearch& headers, tile_cache_t<PixelT>& cache, d::RememberCallback& pc) {
  BOOST_FOREACH(const Tile& t, plate.batch_read(headers)) {
    ImageView<PixelT>& image = cache[d::rowcol_t(t.hdr.row(), t.hdr.col())];
    boost::scoped_ptr<SrcImageResource> r(SrcMemoryImageResource::open(t.hdr.filetype(), &t.data->operator[](0), t.data->size()));
    read_image(image, *r);
    pc.tick();
  }
}

template <class PixelT>
void fast_mipmap( PlateFile& plate, uint32 starting_level, level_data<PixelT>* src_data, bool preblur, d::RememberCallback& mosaic_pc)
{
  typedef ImageView<PixelT> image_t;
  level_data<PixelT> scratch;
  level_data<PixelT> *curr = src_data, *prev = &scratch;

  for (int32 level = starting_level-1; level >= 0; --level)
  {
    std::swap(curr, prev); // point prev at 'current'
    curr->clear();         // and dump the old 'prev'

    // assign the previous-level hdrs to their parents
    BOOST_FOREACH(const composite_map_t::value_type& v, prev->hdrs)
      curr->hdrs[d::parent_tile(d::therow(v.first), d::thecol(v.first))].push_back(v.first);

    // build the new tiles
    BOOST_FOREACH(const composite_map_t::value_type& v, curr->hdrs) {
      const d::rowcol_t& parent = v.first;
      const std::vector<d::rowcol_t>& children = v.second;
      VW_ASSERT(children.size() > 0,  LogicErr() << "How can there be zero here?");
      VW_ASSERT(children.size() <= 4, LogicErr() << "How can there be more than four here?");

      std::map<uint32, image_t> c;
      BOOST_FOREACH(const d::rowcol_t& child, children)
        c[d::calc_composite_id(parent, child)] = prev->cache[child];

      image_t& new_image = curr->cache[parent];
      mipmap_one_tile(new_image, plate.default_tile_size(), c[0], c[1], c[2], c[3], preblur);
      plate.write_update(new_image, d::thecol(parent), d::therow(parent), level);
      mosaic_pc.tick();
    }
  }
}

}

template <class PixelT>
void PlateManager<PixelT>::mipmap(uint32 starting_level, BBox2i const& region,
                                  TransactionOrNeg read_transaction_id,
                                  bool preblur,
                                  const ProgressCallback &progress_callback) const
{
  const uint64 TILE_BYTES  = m_platefile->default_tile_size() * m_platefile->default_tile_size() * uint32(PixelNumBytes<PixelT>::value);
  const uint64 CACHE_BYTES = vw_settings().system_cache_size();
  const uint64 CACHE_TILES = CACHE_BYTES / TILE_BYTES;
  VW_ASSERT(CACHE_TILES > 100, LogicErr() << "You will have many problems if you can't cache at least 100 tiles (you can only store " << CACHE_TILES << ")");

  level_data<PixelT> scratch;
  level_data<PixelT> *curr = &scratch;
  uint32 total_tiles;
  {
    Datastore::TileSearch tiles;
    BOOST_FOREACH(const TileHeader& hdr, m_platefile->search_by_region(starting_level, region, read_transaction_id)) {
      tiles.push_back(hdr);
      // prime the entry and just use the default constructor (we don't need
      // the data, just the region, because the data comes from the batch_read)
      curr->hdrs[d::rowcol_t(hdr.row(), hdr.col())] = composite_map_t::mapped_type();
    }

    // Approximate that loading the tiles from disk is 3% of the total work
    d::RememberCallback tile_load_pc(progress_callback, 0.03, tiles.size());
    load_tiles<PixelT>(*m_platefile, tiles, curr->cache, tile_load_pc);
    tile_load_pc.report_finished();

    // take the top layer off the count, we already have it
    total_tiles = approximate_total_tiles(region, tiles.size()) - tiles.size();
  }

  {
    d::RememberCallback mosaic_pc(progress_callback, 0.97, total_tiles);
    fast_mipmap(*m_platefile, starting_level-1, curr, preblur, mosaic_pc);
  }

  progress_callback.report_finished();
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
    vw_throw( ArgumentErr() << "Unknown option: \"" << mode << "\"." );
  }
}

// Explicit template instantiation
namespace vw {
namespace platefile {

#define VW_INSTANTIATE_PLATE_MANAGER_TYPES(PIXELT)                             \
  template void                                                                \
  PlateManager<PIXELT >::mipmap(uint32 starting_level, BBox2i const& bbox,     \
                                TransactionOrNeg transaction_id,               \
                                bool preblur,                                  \
                                const ProgressCallback &progress_callback) const; \
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
