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
#include <set>

namespace {

  using namespace vw;
  using namespace vw::platefile;

  typedef boost::tuple<uint32, uint32, uint32> rowcoltid_t;
  typedef boost::tuple<uint32, uint32>         rowcol_t;
  typedef boost::tuple<uint32, TileHeader>     tile_order_t;

  BBox2i move_down(const BBox2i& input, uint32 level_change) {
    return input * (1 << level_change);
  }

  rowcol_t parent_tile(const uint32 row, const uint32 col) {
    return rowcol_t(row/2, col/2);
  }

  struct CmpTuple {
    template <typename T1, typename T2>
    bool operator()(const boost::tuple<T1, T2>& a, const boost::tuple<T1, T2>& b) {
      if (a.get<0>() != b.get<0>())
        return a.get<0>() < b.get<0>();
      return a.get<1>() < b.get<1>();
    }
    template <typename T1, typename T2, typename T3>
    bool operator()(const boost::tuple<T1, T2, T3>& a, const boost::tuple<T1, T2, T3>& b) {
      if (a.get<0>() != b.get<0>())
        return a.get<0>() < b.get<0>();
      else if (a.get<1>() != b.get<1>())
        return a.get<1>() < b.get<1>();
      return a.get<2>() < b.get<2>();
    }
  };

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
  void schedule_tile(const TileHeader& hdr, std::vector<tile_order_t>& stack, Datastore::TileSearch& tile_lookup, int32 composite_order_override) {
    stack.push_back(boost::make_tuple(composite_order_override == -1 ? hdr.transaction_id() : boost::numeric_cast<uint32>(composite_order_override), hdr));
    push_heap(stack.begin(), stack.end(), SortByTidDesc());
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
    composite_map_t composite_map;
    tile_cache_t<PixelT> tile_cache;

    {
      // List of tiles we need to fetch
      Datastore::TileSearch tile_lookup;

      // Grab the tiles in the zone
      BOOST_FOREACH(const TileHeader& hdr, m_platefile->search_by_region(level, region, range))
        schedule_tile(hdr, composite_map[rowcol_t(hdr.row(), hdr.col())], tile_lookup, -1);

      // Grab the result of the previous level snapshot
      BOOST_FOREACH(const TileHeader& hdr, m_platefile->search_by_region(level+1, move_down(region, 1), TransactionRange(m_platefile->transaction_id())))
        schedule_tile(hdr, composite_map[parent_tile(hdr.row(), hdr.col())], tile_lookup, 0);

      // Now load the images from disk, decode them, and store them back to the cache
      BOOST_FOREACH(const Tile& t, m_platefile->batch_read(tile_lookup)) {
        if (t.data)
          parse_image_and_store(t, tile_cache);
      }
    }

    // now iterate the output tile set and create the composites
    BOOST_FOREACH(composite_map_t::value_type& t, composite_map) {
      sort_heap(t.second.begin(), t.second.end(), SortByTidDesc());
      mosaic::ImageComposite<PixelT> composite;
      composite.set_draft_mode(true);
      // Insert the images into the composite from highest to lowest tid (already sorted due to sort_heap)
      BOOST_FOREACH(const tile_order_t& order, t.second) {
        const TileHeader& hdr = order.get<1>();
        image_t& img = tile_cache[rowcoltid_t(hdr.row(), hdr.col(), hdr.transaction_id())];
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
