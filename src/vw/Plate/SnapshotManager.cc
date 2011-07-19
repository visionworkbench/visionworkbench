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

namespace d = vw::platefile::detail;

namespace {
  using namespace vw;
  using namespace vw::platefile;

  // row, col at target level, tile headers from next level down (sorted in Tid order)
  typedef std::map<d::rowcol_t, std::vector<d::tile_order_t> > composite_map_t;
  // Cache of parsed images
  template <typename PixelT>
  struct tile_cache_t : public std::map<d::rowcoltid_t, ImageView<PixelT> > {};

  // composite_order_override lets you override the composite order
  void schedule_tile(const TileHeader& hdr, std::vector<d::tile_order_t>& stack, Datastore::TileSearch& tile_lookup, uint32 order) {
    stack.push_back(boost::make_tuple(order, d::rowcoltid_t(hdr.row(), hdr.col(), hdr.transaction_id())));
    tile_lookup.push_back(hdr);
  }

  template <typename PixelT>
  void parse_image_and_store(const Tile& t, tile_cache_t<PixelT>& tile_cache) {
    ImageView<PixelT>& image = tile_cache[d::rowcoltid_t(t.hdr.row(), t.hdr.col(), t.hdr.transaction_id())];
    boost::scoped_ptr<SrcImageResource> r(SrcMemoryImageResource::open(t.hdr.filetype(), &t.data->operator[](0), t.data->size()));
    read_image(image, *r);
  }
}

namespace vw {
namespace platefile {

template <class PixelT>
void SnapshotManager<PixelT>::snapshot(uint32 level, BBox2i const& tile_region, TransactionRange range, const ProgressCallback &progress_) const {
  typedef ImageView<PixelT> image_t;

  // Divide up the region into moderately-sized chunks
  std::list<BBox2i> regions = bbox_tiles(tile_region, 1024, 1024);
  BOOST_FOREACH(const BBox2i& region, regions) {
    SubProgressCallback progress(progress_, progress_.progress(), progress_.progress() + 1./regions.size());

    // Declare these inside the loop so we keep the memory use down to a dull roar
    // This map holds the headers for the parent/dest tiles and the tiles used to generate them
    composite_map_t composite_map;
    // This map holds the headers for the parent/dest tiles and the tiles that are a level below that need to be composed and downsampled first
    // The key is parent tile (same as the composite_map) and the value is a
    // 4-elt vector, the index of which identifies where to compose the image
    // [i.e. 0==UL, 1==UR, 2==LL, 3==LR]
    composite_map_t intermediate_map;

    tile_cache_t<PixelT> tile_cache;

    {
      // List of tiles we need to fetch
      Datastore::TileSearch tile_lookup;

      // Grab the tiles in the zone
      BOOST_FOREACH(const TileHeader& hdr, m_platefile->search_by_region(level, region, range))
        schedule_tile(hdr, composite_map[d::rowcol_t(hdr.row(), hdr.col())], tile_lookup, hdr.transaction_id());

      // Grab the result of the previous level snapshot
      BOOST_FOREACH(const TileHeader& hdr, m_platefile->search_by_region(level+1, d::move_down(region, 1), TransactionRange(m_platefile->transaction_id()))) {
        d::rowcol_t parent = d::parent_tile(hdr.row(), hdr.col());
        uint32 composite_id = d::calc_composite_id(parent, hdr);
        schedule_tile(hdr, intermediate_map[parent], tile_lookup, composite_id);
      }

      // Now load the images from disk, decode them, and store them back to the cache
      d::RememberCallback tile_load_pc(progress, 0.2, tile_lookup.size());
      BOOST_FOREACH(const Tile& t, m_platefile->batch_read(tile_lookup)) {
        parse_image_and_store(t, tile_cache);
        tile_load_pc.tick();
      }
    }

    {
      d::RememberCallback subtile_composite_pc(progress, 0.64, intermediate_map.size());
      // Now look up all the child tiles, and compose/downsample them to the target layer
      BOOST_FOREACH(const composite_map_t::value_type& t, intermediate_map) {
        const d::rowcol_t    dest_loc(d::therow(t.first), d::thecol(t.first));
        const d::rowcoltid_t dest_loc_tid(d::therow(t.first), d::thecol(t.first), 0);
        const std::vector<d::tile_order_t>& children = t.second;
        VW_ASSERT(children.size() > 0,  LogicErr() << "How can there be zero here?");
        VW_ASSERT(children.size() <= 4, LogicErr() << "How can there be more than four here?");

        std::map<uint32, image_t> c;
        BOOST_FOREACH(const d::tile_order_t& order, children)
          c[order.get<0>()] = tile_cache[order.get<1>()];

        // use tid 0 as the dest, to make sure it ends up under the other tiles
        mipmap_one_tile(tile_cache[dest_loc_tid], m_platefile->default_tile_size(), c[0], c[1], c[2], c[3]);
        composite_map[dest_loc].push_back(boost::make_tuple(0, dest_loc_tid));
        subtile_composite_pc.tick();
      }
    }

    {
      d::RememberCallback composite_pc(progress, 0.16, composite_map.size());
      // now iterate the output tile set and create the composites
      BOOST_FOREACH(composite_map_t::value_type& t, composite_map) {
        std::sort(t.second.begin(), t.second.end(), d::SortByTidDesc());
        mosaic::ImageComposite<PixelT> composite;
        composite.set_draft_mode(true);
        // Insert the images into the composite from highest to lowest tid (already sorted due to sort_heap)
        BOOST_FOREACH(const d::tile_order_t& order, t.second) {
          const d::rowcoltid_t& hdr = order.get<1>();
          image_t& img = tile_cache[d::rowcoltid_t(d::therow(hdr), d::thecol(hdr), d::thetid(hdr))];
          if (!img) {
            vw_out(WarningMessage, "platefile.snapshot") << "Failed to load image for " << hdr << std::endl;
            continue;
          }
          composite.insert(img, 0, 0);
          if (is_opaque(img))
            break;
        }
        if (composite.cols() == 0 || composite.rows() == 0) {
          vw_out(WarningMessage, "platefile.snapshot") << "Empty tile list, skipping writing tile row=" << d::therow(t.first) << " col=" << d::thecol(t.first) << std::endl;
          continue;
        }
        composite.prepare(BBox2i(0,0,m_platefile->default_tile_size(), m_platefile->default_tile_size()));
        image_t tile = composite;
        m_platefile->write_update(tile, d::thecol(t.first), d::therow(t.first), level);
        composite_pc.tick();
      }
    }
  }
  progress_.report_finished();
}

// Create a full snapshot of every level and every region in the mosaic.
template <class PixelT>
void SnapshotManager<PixelT>::full_snapshot(TransactionRange read_transaction_range) const {
  for (int32 level = m_platefile->num_levels()-1; level >= 0; --level) {
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
