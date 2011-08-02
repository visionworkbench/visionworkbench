// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Plate/SnapshotManager.h>
#include <vw/Plate/detail/MipmapHelpers.h>
#include <vw/Plate/PlateFile.h>
#include <vw/Plate/TileManipulation.h>
#include <vw/Mosaic/ImageComposite.h>
#include <vw/Image/ImageView.h>
#include <boost/foreach.hpp>
#include <boost/lambda/construct.hpp>

namespace d = vw::platefile::detail;

namespace {
  using namespace vw;
  using namespace vw::platefile;

  // row, col at target level, tile headers from next level down (sorted in Tid order)
  typedef std::map<d::rowcol_t, std::vector<TileHeader> > composite_map_t;
  // Cache of parsed images
  template <typename PixelT>
  struct tile_cache_t : public std::map<d::rowcoltid_t, ImageView<PixelT> > {};

  template <typename PixelT>
  void parse_image_and_store(const Tile& t, tile_cache_t<PixelT>& tile_cache) {
    ImageView<PixelT>& image = tile_cache[d::rowcoltid_t(t.hdr.row(), t.hdr.col(), t.hdr.transaction_id())];
    boost::scoped_ptr<SrcImageResource> r(SrcMemoryImageResource::open(t.hdr.filetype(), &t.data->operator[](0), t.data->size()));
    read_image(image, *r);
  }

  template <typename PixelT>
  void mosaic_in_tid_order(ImageView<PixelT>& out, const size_t size, const tile_cache_t<PixelT>& tile_cache, std::vector<TileHeader>& tiles) {
    typedef ImageView<PixelT> image_t;
    typedef tile_cache_t<PixelT> cache_t;

    std::sort(tiles.begin(), tiles.end(), d::SortByTidDesc());
    mosaic::ImageComposite<PixelT> composite;
    composite.set_draft_mode(true);
    // Insert the images into the composite from highest to lowest tid (already sorted due to SortByTidDesc)
    BOOST_FOREACH(const TileHeader& hdr, tiles) {
      typename cache_t::const_iterator i = tile_cache.find(d::rowcoltid_t(d::therow(hdr), d::thecol(hdr), d::thetid(hdr)));
      if (i == tile_cache.end()) {
        vw_out(WarningMessage, "platefile.snapshot") << "Failed to load image for " << hdr << std::endl;
        continue;
      }
      composite.insert(i->second, 0, 0);
      if (is_opaque(i->second))
        break;
    }
    if (composite.cols() == 0 || composite.rows() == 0) {
      out.reset();
      return;
    }
    composite.prepare(BBox2i(0,0,size,size));
    out = composite;
  }
}

namespace vw {
namespace platefile {

template <class PixelT>
void SnapshotManager<PixelT>::snapshot(uint32 level, BBox2i const& tile_region, TransactionRange range, const ProgressCallback &progress_) const {
  typedef ImageView<PixelT> image_t;

  // The byte size of an uncompressed image
  const uint64 TILE_BYTES  = m_write_plate->default_tile_size() * m_write_plate->default_tile_size() * uint32(PixelNumBytes<PixelT>::value);
  const uint64 CACHE_BYTES = vw_settings().system_cache_size();
  const uint64 CACHE_TILES = CACHE_BYTES / TILE_BYTES;
  // This is an arbitrary value, to hopefully catch pathlogically-small cache sizes
  VW_ASSERT(CACHE_TILES > 100, LogicErr() << "You will have many problems if you can't cache at least 100 tiles (you can only store " << CACHE_TILES << ")");

  Datastore::TileSearch tile_lookup;
  tile_lookup.reserve(CACHE_TILES);

  // Divide up the region into moderately-sized chunks
  std::list<BBox2i> regions = bbox_tiles(tile_region, 1024, 1024);
  BOOST_FOREACH(const BBox2i& region, regions) {
    SubProgressCallback progress(progress_, progress_.progress(), progress_.progress() + 1./regions.size());

    tile_lookup.clear();

    // This map holds the headers for the current level
    composite_map_t composite_map;

    // Grab the tiles in the zone
    BOOST_FOREACH(const TileHeader& hdr, m_read_plate->search_by_region(level, region, range))
      composite_map[d::rowcol_t(hdr.row(), hdr.col())].push_back(hdr);

    d::RememberCallback pc(progress, 1, composite_map.size());

    // queue up CACHE_TILES worth of tiles
    BOOST_FOREACH(composite_map_t::value_type& t, composite_map) {
      if (tile_lookup.size() + t.second.size() >= CACHE_TILES) {
        VW_ASSERT(tile_lookup.size() > 0, LogicErr() << "Your cache size in tiles (" << CACHE_TILES << ") is smaller than the required minimum of " << t.second.size());

        tile_cache_t<PixelT> tile_cache;
        BOOST_FOREACH(const Tile& t, m_read_plate->batch_read(tile_lookup))
          parse_image_and_store(t, tile_cache);
        image_t tile;
        mosaic_in_tid_order(tile, m_write_plate->default_tile_size(), tile_cache, t.second);
        if (!tile) {
          vw_out(WarningMessage, "platefile.snapshot") << "Empty tile list, skipping writing tile row=" << d::therow(t.first) << " col=" << d::thecol(t.first) << std::endl;
          continue;
        }
        m_write_plate->write_update(tile, d::thecol(t.first), d::therow(t.first), level);
        pc.tick(tile_lookup.size());
        tile_lookup.clear();
      }
      std::transform(t.second.begin(), t.second.end(), tile_lookup.end(), boost::lambda::constructor<Tile>());
    }
  }
  progress_.report_finished();
}

// Create a full snapshot of every level and every region in the mosaic.
template <class PixelT>
void SnapshotManager<PixelT>::full_snapshot(TransactionRange read_transaction_range) const {
  for (int32 level = m_read_plate->num_levels()-1; level >= 0; --level) {
    TerminalProgressCallback prog("plate.snapshot", std::string("Snapshot level ") + stringify(level));
    snapshot(level, d::move_down(BBox2i(0,0,1,1), level), read_transaction_range, prog);
  }
}

}} // namespace vw::platefile

// Explicit template instatiation
#define _VW_INSTANTIATE(Px)\
  template\
  void vw::platefile::SnapshotManager<Px >::snapshot(uint32 level, BBox2i const& bbox, TransactionRange, const ProgressCallback &progress) const;\
  template\
  void vw::platefile::SnapshotManager<Px >::full_snapshot(TransactionRange) const;

_VW_INSTANTIATE(vw::PixelGrayA<vw::uint8>);
_VW_INSTANTIATE(vw::PixelGrayA<vw::int16>);
_VW_INSTANTIATE(vw::PixelGrayA<vw::float32>);
_VW_INSTANTIATE(vw::PixelRGBA<vw::uint8>);
