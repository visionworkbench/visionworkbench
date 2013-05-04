// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NASA Vision Workbench is licensed under the Apache License,
//  Version 2.0 (the "License"); you may not use this file except in
//  compliance with the License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__


#include <vw/Plate/SnapshotManager.h>
#include <vw/Plate/detail/MipmapHelpers.h>
#include <vw/Plate/PlateFile.h>
#include <vw/Plate/TileManipulation.h>
#include <vw/Mosaic/ImageComposite.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/UtilityViews.h>
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

    // We arrange in descending order, so that when the higher level
    // tiles fill the output completely, we'll stop inserting.
    out = constant_view(PixelT(), size, size);
    std::sort(tiles.begin(), tiles.end(), d::SortByTidDesc());
    // Insert the images into the composite from lowest to highest tid
    BOOST_FOREACH(const TileHeader& hdr, tiles) {
      typename cache_t::const_iterator i = tile_cache.find(d::rowcoltid_t(d::therow(hdr), d::thecol(hdr), d::thetid(hdr)));
      if (i == tile_cache.end()) {
        VW_OUT(WarningMessage, "platefile.snapshot") << "Failed to load image for " << hdr << std::endl;
        continue;
      }

      // Insert the next lower tile and then lay the previous
      // composite on top. If this new combination is opaque, we've
      // finished the snapshot.
      mosaic::ImageComposite<PixelT> composite;
      composite.set_draft_mode(true);
      composite.insert(i->second, 0, 0);
      composite.insert( out, 0, 0 );
      composite.prepare(BBox2i(0,0,size,size));
      out = composite;
      if ( is_opaque( out ) )
        return;
    }
  }
}

namespace vw {
namespace platefile {

template <class PixelT>
void SnapshotManager<PixelT>::snapshot(uint32 level, BBox2i const& tile_region, TransactionRange range, const ProgressCallback &progress) const {
  // The byte size of an uncompressed image
  const uint64 TILE_BYTES  = m_write_plate->default_tile_size() * m_write_plate->default_tile_size() * uint32(PixelNumBytes<PixelT>::value);
  const uint64 CACHE_BYTES = vw_settings().system_cache_size();
  const uint64 CACHE_TILES = CACHE_BYTES / TILE_BYTES;
  VW_ASSERT( CACHE_TILES > 1, LogicErr() << "Cache too small to process any tiles in snapshot." );
  // This is an arbitrary value, to hopefully catch
  // pathlogically-small cache sizes
  if ( CACHE_TILES < 100 )
    VW_OUT(WarningMessage) << "You will lose a lot speed to thrashing if you can't cache at least 100 tiles (you can only store " << CACHE_TILES << ")\n";

  // Divide up the region into moderately-sized chunks
  std::list<BBox2i> regions = bbox_tiles(tile_region, 1024, 1024);
  BOOST_FOREACH(const BBox2i& region, regions) {
    SubProgressCallback region_pc(progress, progress.progress(), progress.progress() + 1./regions.size());

    // This map holds the headers for the current level
    composite_map_t composite_map;

    // Grab the tiles in the zone
    BOOST_FOREACH(const TileHeader& hdr, m_read_plate->search_by_region(level, region, range))
      composite_map[d::rowcol_t(hdr.row(), hdr.col())].push_back(hdr);

    if (composite_map.size() == 0) {
      region_pc.report_finished();
      continue;
    }

    d::RememberCallback pc(region_pc, 1, composite_map.size());

    composite_map_t::const_iterator i = composite_map.begin(), end = composite_map.end();

    do {
      // queue up CACHE_TILES worth of tiles into composite_batch
      size_t size = 0;
      composite_map_t composite_batch;

      composite_map_t::const_iterator i_before = i;
      for (; i != end; ++i) {
        if (size + i->second.size() <= CACHE_TILES) {
          composite_batch.insert(*i);
          size += i->second.size();
        } else {
          VW_ASSERT( i_before != i,
                     LogicErr() << "Cache too small to perform composite. Stuck in infinite loop." );
          break;
        }
      }
      if (composite_batch.size() == 0)
        continue;

      Datastore::TileSearch tile_lookup;
      tile_cache_t<PixelT> tile_cache;

      BOOST_FOREACH(const composite_map_t::value_type& t, composite_batch) {
        //std::cerr << "Batching tiles for " << t.first << ": ";
        //std::transform(t.second.begin(), t.second.end(), std::ostream_iterator<uint32>(std::cerr, ", "), boost::bind(&TileHeader::transaction_id, _1));
        //std::cerr << " (currently have " << tile_lookup.size() << "/" << CACHE_TILES << ")" << std::endl;
        std::copy(t.second.begin(), t.second.end(), std::back_inserter(tile_lookup));
      }

      BOOST_FOREACH(const Tile& t, m_read_plate->batch_read(tile_lookup))
        parse_image_and_store(t, tile_cache);
      tile_lookup.clear();

      BOOST_FOREACH(composite_map_t::value_type& t, composite_batch) {

        // Early exit condition if there is only one tile
        if (t.second.size() == 1 ) {
          typename tile_cache_t<PixelT>::const_iterator i =
            tile_cache.find(d::rowcoltid_t(d::therow(t.second[0]), d::thecol(t.second[0]), d::thetid(t.second[0])));
          if ( i == tile_cache.end() ) {
            VW_OUT(WarningMessage, "platefile.snapshot") << "Failed to load image for " << t.second[0] << std::endl;
            continue;
          }
          m_write_plate->write_update(i->second, d::thecol(t.first), d::therow(t.first), level);
          continue;
        }

        // Perform an actual composite if there are multiple tiles to be inserted
        ImageView<PixelT> tile;
        mosaic_in_tid_order<PixelT>(tile, m_write_plate->default_tile_size(), tile_cache, t.second);
        if (!tile) {
          VW_OUT(WarningMessage, "platefile.snapshot") << "Empty tile list, skipping writing tile row=" << d::therow(t.first) << " col=" << d::thecol(t.first) << std::endl;
          continue;
        }
        m_write_plate->write_update(tile, d::thecol(t.first), d::therow(t.first), level);
        pc.tick();
      }
      composite_batch.clear();
      size = 0;
    } while (i != end);
    region_pc.report_finished();
  }
  progress.report_finished();
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
