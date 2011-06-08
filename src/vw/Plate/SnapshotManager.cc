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

namespace vw {
namespace platefile {

#if 0

  typedef std::map<uint32, TileHeader> tilemap_t;

  // ----------------------- WeightedAvg blending -------------------------------

  // Performs alpha blending with un-premultiplied alpha
  template <class PixelT>
  struct WeightedAvgBlendFunctor : ReturnFixedType<PixelT> {
    inline PixelT operator()( PixelT const& pixel_a, PixelT const& pixel_b ) const {
      typedef typename PixelChannelType<PixelT>::type channel_type;

      // Extract the alpha channels
      channel_type weight_a = alpha_channel(pixel_a);
      channel_type weight_b = alpha_channel(pixel_b);

      // Compute the new pixel values
      //      PixelT result = pixel_a * weight_a + pixel_b * weight_b * (1-weight_a); // For debugging: alpha blending
      PixelT result = pixel_a * weight_a + pixel_b * (1-weight_a);


      // This check in necessary to make sure that NaNs and Inf, -Inf
      // don't slip through the cracks.  This can happen sometimes
      // when you mask out bad values, but leave the NaNs in the
      // grayscale channel of the image.  The NaN will propagate
      // through the computation regardless of whether they are
      // "masked out", and appear in the output DEM.  This prevents
      // this from ever happening!
      if (result[0] != result[0]) {
        return PixelT();
      }

      // Compute the new weight value as the max of the individual weight values.
      //      result.a() = weight_a + weight_b * (1-weight_a);   // For debugging: alpha blending
      result.a() = std::max(weight_a, weight_b);

      // This check ensures that the data value in the snapshot is
      // zero if the alpha value is zero.
      if (result.a() == 0)
        return PixelT();

      return result;
    }
  };


  /// Create a new view which is the alpha blended version of two
  /// image views.  Image A ends up on top of image B.
  template <class ImageT>
  inline BinaryPerPixelView<ImageT, ImageT, WeightedAvgBlendFunctor<typename ImageT::pixel_type> > weighted_avg_blend( ImageViewBase<ImageT> const& image_a, ImageViewBase<ImageT> const& image_b ) {
    VW_ASSERT( image_a.impl().rows() == image_b.impl().rows() &&
               image_a.impl().cols() == image_b.impl().cols(),
               ArgumentErr()
               << "platefile::weighted_avg_blend() failed.  Image dimensions do not match." );

    return BinaryPerPixelView<ImageT,ImageT,WeightedAvgBlendFunctor<typename ImageT::pixel_type> >( image_a.impl(), image_b.impl() );
  }

  // -----------------------   Utilities   -------------------------------

  bool is_leaf(boost::shared_ptr<PlateFile> platefile, TileHeader const& tile_header) {

    uint32 current_col = tile_header.col();
    uint32 current_row = tile_header.row();

    for (uint32 new_row = current_row*2; new_row < current_row*2+2; ++new_row) {
      for (uint32 new_col = current_col*2; new_col < current_col*2+2; ++new_col) {

        std::list<TileHeader> tile_records =
          platefile->search_by_location(new_col, new_row, tile_header.level() + 1, TransactionRange(tile_header.transaction_id()));

          // If this node has child tiles, then it is not a leaf node.
          if (tile_records.size() > 0)
            return false;
      }
    }

    // If not child tiles were found, then this is a leaf node.
    return true;
  }
#endif

  namespace detail {
    BBox2i move_down(const BBox2i& input, uint32 level_change) {
      return input * (1 << level_change);
    }

#if 0
    void collect(boost::shared_ptr<PlateFile> plate, tilemap_t& all_tiles, uint32 col, uint32 row, uint32 level, vw::BBox2i const& target_region, uint32 target_level, const TransactionRange& range)
    {
      std::list<TileHeader> tiles = plate->search_by_location(col, row, level, range);
      // This branch is empty. nothing to find.
      if (tiles.empty())
        return;

      // Add the current level to the list (we already verified that these are in-region)
      BOOST_FOREACH(const TileHeader& t, tiles)
        all_tiles.insert(std::make_pair(t.transaction_id(), t));

      if (target_level == level)
        return;

      for (uint32 row2 = row*2; row2 < row*2+2; ++row2)
        for (uint32 col2 = col*2; col2 < col*2+2; ++col2)
          if (move_down(BBox2i(col2, row2, 1, 1), target_level - level).intersects(target_region))
            collect(plate, all_tiles, col2, row2, level+1, target_region, target_level, range);
    }
#endif
  }

#if 0
  template <class PixelT>
  int SnapshotManager<PixelT>::snapshot_helper(uint32 current_col,
        uint32 current_row,
        uint32 current_level,
        vw::BBox2i const& target_region,
        uint32 target_level,
        TransactionRange read_transaction_range) const
  {

    tilemap_t composite_tiles;
    detail::collect(m_platefile, composite_tiles, current_col, current_row, current_level, target_region, target_level, read_transaction_range);
    vw_out(DebugMessage, "plate::snapshot") << "Collected tiles for snapshot: " << composite_tiles.size() << std::endl;

    // If we arrive at the target level and find that we have fewer
    // than two tiles to composite, then there's no actual
    // compositing that needs to be done.  We can safely bail on
    // this tile.
    if (composite_tiles.size() < 2)
      return 0;

    // Create a tile into which we will accumulate composited data.
    ImageView<PixelT> composite_tile(m_platefile->default_tile_size(), m_platefile->default_tile_size());

    vw_out(DebugMessage, "plate::snapshot") << "Starting Compositing run...\n";
    int num_composited = 0;

    BOOST_REVERSE_FOREACH(const tilemap_t::value_type& tile, composite_tiles) {
      const TileHeader& hdr = tile.second;
      ImageView<PixelT> new_tile(m_platefile->default_tile_size(), m_platefile->default_tile_size());

      // Tile access is highly localized during snapshotting, so it
      // speeds things up considerably to cache them here.  We check
      // the cache first and then read the tile from the platefile
      // if nothing is found.
      if ( !(this->restore_tile(new_tile, hdr.col(), hdr.row(), hdr.level(), hdr.transaction_id())) ) {
        try {
          // If the tile isn't in our cache, then we take the
          // performance hit and read the tile from the platefile.
          m_platefile->read(new_tile, hdr.col(), hdr.row(), hdr.level(), hdr.transaction_id(), true); // exact_transaction_match
        } catch (const BlobIoErr &e) {
          m_platefile->error_log() << "BlobIoErr, tile[" << hdr << "]: " << e.what() << "\n";
        } catch (const IOErr &e) {
          m_platefile->error_log() << "IOErr, tile[" << hdr << "]: " << e.what() << "\n";
        }

        // Save the tile to the read cache.
        this->save_tile(new_tile, hdr.col(), hdr.row(), hdr.level(), hdr.transaction_id());
      }

      vw_out(DebugMessage, "plate::snapshot") << "\t--> Adding tile: " << hdr << "\n";

      // If this is the first tile in this location, it's already at
      // the target_level, and it's opaque (which is a very common
      // case), then we don't need to save it as part of the
      // snapshot since its already the top tile in the snapshot.
      // This optimization saves us both space and time!
      if ( hdr.level() == target_level && &tile == &(*composite_tiles.rbegin()) && is_opaque(new_tile) ) {
        return 0;
      }

      // If we are compositing a tile from a lower resolution in the
      // pyramid, it needs to be supersampled and cropped before we
      // can use it here.
      if (hdr.level() != target_level) {

        new_tile = resample_img_from_level(
            new_tile, hdr.col(), hdr.row(), hdr.level(), current_col, current_row, target_level);


      }

      // Add the new tile UNDERNEATH the tiles we have already composited.
      //
      // If we are snapshotting terrain, we do this with a weighted
      // average.  Otherwise we use normal alpha blending.
      if (m_tweak_settings_for_terrain) {
        composite_tile = weighted_avg_blend(composite_tile, new_tile);
      } else {
        vw::mosaic::ImageComposite<PixelT> composite;
        composite.insert(new_tile, 0, 0);
        composite.insert(composite_tile, 0, 0);
        composite.set_draft_mode( true );
        composite.prepare();
        composite_tile = composite;
      }
      ++num_composited;

      // Check to see if the composite tile is opaque.  If it is,
      // then there's no point in compositing any more tiles because
      // they will end up hidden anyway.
      if ( is_opaque(composite_tile) ) {
        break;
      }
    }

    if (num_composited > 0) {
      // -- debug
      vw_out(DebugMessage, "plate::snapshot")  << "\t    Compositing " << current_col << " "
        << current_row << " @ " << current_level
        << "  [ " << num_composited << " ]\n";

      m_platefile->write_update(composite_tile, current_col, current_row, current_level);
      return 1;
    } else {
      return 0;
    }
  }
#endif

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
  bool operator()(const TileHeader& a, const TileHeader& b)
  {
    return b.transaction_id() < a.transaction_id();
  }
};

template <class PixelT>
void SnapshotManager<PixelT>::snapshot(uint32 level, BBox2i const& tile_region, TransactionRange range) const {
  typedef boost::tuple<uint32, uint32, uint32> rowcoltid_t;
  typedef boost::tuple<uint32, uint32> rowcol_t;
  typedef ImageView<PixelT> image_t;

  // row, col at target level, tile headers from next level down (sorted in Tid order)
  typedef std::map<rowcol_t, std::vector<TileHeader>, CmpTuple>    composite_map_t;
  // Cache of parsed images
  typedef std::map<rowcoltid_t, image_t, CmpTuple>  tile_cache_t;

  composite_map_t composite_map;
  tile_cache_t tile_cache;

  // XXX: This will make it so the snapshot tiles go on top of the existing tiles at a level
  // This ensures zoom consistency at the expense of possible canonical data.
  // This can be switched by demoting the priority of the output transaction

  // Divide up the region into moderately-sized chunks
  BOOST_FOREACH(const BBox2i& region, bbox_tiles(tile_region, 1024, 1024)) {
    // Grab all tiles that the chunk depends on (one level below) and sort them by tid
    BOOST_FOREACH(const TileHeader& hdr, m_platefile->search_by_region(level, detail::move_down(region, 1), range)) {
      composite_map_t::mapped_type& h = composite_map[rowcol_t(hdr.row(), hdr.col())];
      h.push_back(hdr);
      push_heap(h.begin(), h.end(), SortByTidDesc());
      tile_cache[rowcoltid_t(hdr.row(), hdr.col(), hdr.transaction_id())] = image_t();
    }
    // XXX: Doesn't look up tiles

    // now iterate the output tile set and create the composites
    BOOST_FOREACH(composite_map_t::value_type& t, composite_map) {
      sort_heap(t.second.begin(), t.second.end(), SortByTidDesc());
      mosaic::ImageComposite<PixelT> composite;
      composite.set_draft_mode(true);
      // Insert the images into the composite from highest to lowest tid (already sorted due to sort_heap)
      BOOST_FOREACH(TileHeader& hdr, t.second) {
        image_t& img = tile_cache[rowcoltid_t(hdr.row(), hdr.col(), hdr.transaction_id())];
        if (!img) {
          vw_out(WarningMessage, "platefile.snapshot") << "Failed to load image for " << hdr << std::endl;
          continue;
        }
        composite.insert(img, 0, 0);
        if (is_opaque(img))
          break;
      }
      composite.prepare();
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
