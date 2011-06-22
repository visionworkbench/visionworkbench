// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Plate/SnapshotManager.h>
#include <vw/Plate/PlateFile.h>
#include <vw/Plate/TileManipulation.h>
#include <vw/Mosaic/ImageComposite.h>
#include <vw/Image/ImageView.h>
#include <boost/foreach.hpp>
#include <boost/assign/list_of.hpp>
#include <set>

namespace {

  using namespace vw;
  using namespace vw::platefile;

  typedef boost::tuple<uint32, uint32, uint32> rowcoltid_t;
  typedef boost::tuple<uint32, uint32>         rowcol_t;
  typedef boost::tuple<uint32, rowcoltid_t>    tile_order_t;

  BBox2i move_down(const BBox2i& input, uint32 level_change) {
    return input * (1 << level_change);
  }

  rowcol_t parent_tile(const uint32 row, const uint32 col) {
    return rowcol_t(row/2, col/2);
  }

  struct CmpTuple {
    template <typename T1, typename T2>
    bool operator()(const boost::tuple<T1, T2>& a, const boost::tuple<T1, T2>& b) const {
      if (a.get<0>() != b.get<0>())
        return a.get<0>() < b.get<0>();
      return a.get<1>() < b.get<1>();
    }
    template <typename T1, typename T2, typename T3>
    bool operator()(const boost::tuple<T1, T2, T3>& a, const boost::tuple<T1, T2, T3>& b) const {
      if (a.get<0>() != b.get<0>())
        return a.get<0>() < b.get<0>();
      else if (a.get<1>() != b.get<1>())
        return a.get<1>() < b.get<1>();
      return a.get<2>() < b.get<2>();
    }
  };

  // given a parent tile (1 level up) and a hdr, calculate the composite id;
  // [0==UL, 1==UR, 2==LL, 3==LR]
  uint32 calc_composite_id(const rowcol_t& parent, const TileHeader& hdr) {
    typedef std::map<rowcol_t, uint32, CmpTuple> map_t;
    static const map_t lookup = boost::assign::map_list_of
      (rowcol_t(0,0), 0)
      (rowcol_t(0,1), 1)
      (rowcol_t(1,0), 2)
      (rowcol_t(1,1), 3);

    rowcol_t offset(hdr.row() - parent.get<0>() * 2,
                    hdr.col() - parent.get<1>() * 2);

    map_t::const_iterator i = lookup.find(offset);
    VW_ASSERT(i != lookup.end(), LogicErr() << "Cannot determine composite id for hdr " << hdr << " and parent "
        << parent.get<1>() << "," << parent.get<0>());
    return i->second;
  }

  struct SortByTidDesc {
    bool operator()(const tile_order_t& a, const tile_order_t& b)
    {
      return b.get<0>() < a.get<0>();
    }
  };

  // row, col at target level, tile headers from next level down (sorted in Tid order)
  typedef std::map<rowcol_t, std::vector<tile_order_t>, CmpTuple> composite_map_t;
  // Cache of parsed images
  template <typename PixelT>
  struct tile_cache_t : public std::map<rowcoltid_t, ImageView<PixelT>, CmpTuple> {};

  // composite_order_override lets you override the composite order
  void schedule_tile(const TileHeader& hdr, std::vector<tile_order_t>& stack, Datastore::TileSearch& tile_lookup, uint32 order) {
    stack.push_back(boost::make_tuple(order, rowcoltid_t(hdr.row(), hdr.col(), hdr.transaction_id())));
    tile_lookup.push_back(hdr);
  }

  template <typename PixelT>
  void parse_image_and_store(const Tile& t, tile_cache_t<PixelT>& tile_cache) {
    ImageView<PixelT>& image = tile_cache[rowcoltid_t(t.hdr.row(), t.hdr.col(), t.hdr.transaction_id())];
    boost::scoped_ptr<SrcImageResource> r(SrcMemoryImageResource::open(t.hdr.filetype(), &t.data->operator[](0), t.data->size()));
    read_image(image, *r);
  }
} // anon

namespace vw {
namespace platefile {

std::ostream& operator<<(std::ostream& o, const rowcoltid_t& hdr) {
  return (o << hdr.get<1>() << "," << hdr.get<0>() << " (t_id = " << hdr.get<2>() << ")");
}

template <class PixelT>
void SnapshotManager<PixelT>::snapshot(uint32 level, BBox2i const& tile_region, TransactionRange range) const {
  typedef ImageView<PixelT> image_t;

  // Divide up the region into moderately-sized chunks
  BOOST_FOREACH(const BBox2i& region, bbox_tiles(tile_region, 1024, 1024)) {
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
        schedule_tile(hdr, composite_map[rowcol_t(hdr.row(), hdr.col())], tile_lookup, hdr.transaction_id());

      // Grab the result of the previous level snapshot
      BOOST_FOREACH(const TileHeader& hdr, m_platefile->search_by_region(level+1, move_down(region, 1), TransactionRange(m_platefile->transaction_id()))) {
        rowcol_t parent = parent_tile(hdr.row(), hdr.col());
        uint32 composite_id = calc_composite_id(parent, hdr);
        schedule_tile(hdr, intermediate_map[parent], tile_lookup, composite_id);
      }

      // Now load the images from disk, decode them, and store them back to the cache
      BOOST_FOREACH(const Tile& t, m_platefile->batch_read(tile_lookup))
        parse_image_and_store(t, tile_cache);

      // Now look up all the child tiles, and compose/downsample them to the target layer
      BOOST_FOREACH(const composite_map_t::value_type& t, intermediate_map) {
        const rowcol_t    dest_loc(t.first.get<0>(), t.first.get<1>());
        const rowcoltid_t dest_loc_tid(t.first.get<0>(), t.first.get<1>(), 0);
        const std::vector<tile_order_t>& children = t.second;
        VW_ASSERT(children.size() > 0,  LogicErr() << "How can there be zero here?");
        VW_ASSERT(children.size() <= 4, LogicErr() << "How can there be more than four here?");

        std::map<uint32, image_t> c;
        BOOST_FOREACH(const tile_order_t& order, children)
          c[order.get<0>()] = tile_cache[order.get<1>()];

        // use tid 0 as the dest, to make sure it ends up under the other tiles
        mipmap_one_tile(tile_cache[dest_loc_tid], m_platefile->default_tile_size(), c[0], c[1], c[2], c[3]);
        composite_map[dest_loc].push_back(boost::make_tuple(0, dest_loc_tid));
      }
    }

    // now iterate the output tile set and create the composites
    BOOST_FOREACH(composite_map_t::value_type& t, composite_map) {
      std::sort(t.second.begin(), t.second.end(), SortByTidDesc());
      mosaic::ImageComposite<PixelT> composite;
      composite.set_draft_mode(true);
      // Insert the images into the composite from highest to lowest tid (already sorted due to sort_heap)
      BOOST_FOREACH(const tile_order_t& order, t.second) {
        const rowcoltid_t& hdr = order.get<1>();
        image_t& img = tile_cache[rowcoltid_t(hdr.get<0>(), hdr.get<1>(), hdr.get<2>())];
        if (!img) {
          vw_out(WarningMessage, "platefile.snapshot") << "Failed to load image for " << hdr << std::endl;
          continue;
        }
        composite.insert(img, 0, 0);
        if (is_opaque(img))
          break;
      }
      if (composite.cols() == 0 || composite.rows() == 0) {
        vw_out(WarningMessage, "platefile.snapshot") << "Empty tile list, skipping writing tile row=" << t.first.get<0>() << " col=" << t.first.get<1>() << std::endl;
        continue;
      }
      composite.prepare(BBox2i(0,0,m_platefile->default_tile_size(), m_platefile->default_tile_size()));
      image_t tile = composite;
      m_platefile->write_update(tile, t.first.get<1>(), t.first.get<0>(), level);
    }
  }
}

// Create a full snapshot of every level and every region in the mosaic.
template <class PixelT>
void SnapshotManager<PixelT>::full_snapshot(TransactionRange read_transaction_range) const {

  for (uint32 level = 0; level < m_platefile->num_levels(); ++level) {

    // For debugging:
    // for (int level = 0; level < 10; ++level) {

    // Snapshot the entire region at each level.  These region will be
    // broken down into smaller work units in snapshot().
    uint32 region_size = 1 << level;
    uint32 subdivided_region_size = region_size / 16;
    if (subdivided_region_size < 1024) subdivided_region_size = 1024;
    BBox2i full_region(0,0,region_size,region_size);
    std::list<BBox2i> workunits = bbox_tiles(full_region,
                                             subdivided_region_size,
                                             subdivided_region_size);
    for ( std::list<BBox2i>::iterator region_iter = workunits.begin();
          region_iter != workunits.end(); ++region_iter) {
      snapshot(level, *region_iter, read_transaction_range);
    }
  }
}
}} // namespace vw::platefile

// Explicit template instatiation
#define _VW_INSTANTIATE(Px)\
  template\
  void vw::platefile::SnapshotManager<Px >::snapshot(uint32 level, BBox2i const& bbox, TransactionRange) const;\
  template\
  void vw::platefile::SnapshotManager<Px >::full_snapshot(TransactionRange) const;

_VW_INSTANTIATE(vw::PixelGrayA<vw::uint8>);
_VW_INSTANTIATE(vw::PixelGrayA<vw::int16>);
_VW_INSTANTIATE(vw::PixelGrayA<vw::float32>);
_VW_INSTANTIATE(vw::PixelRGBA<vw::uint8>);
