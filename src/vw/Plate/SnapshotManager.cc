// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Plate/SnapshotManager.h>
#include <vw/Plate/PlateFile.h>
#include <vw/Plate/TileManipulation.h>
#include <vw/Mosaic/ImageComposite.h>
#include <vw/Image/ImageView.h>
#include <set>

namespace vw {
namespace platefile {

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

    int current_col = tile_header.col();
    int current_row = tile_header.row();

    for (int new_row = current_row*2; new_row < current_row*2+2; ++new_row) {
      for (int new_col = current_col*2; new_col < current_col*2+2; ++new_col) {

          std::list<TileHeader> tile_records;
          tile_records = platefile->search_by_location(new_col, new_row,
                                                       tile_header.level() + 1,
                                                       tile_header.transaction_id(),
                                                       tile_header.transaction_id(),
                                                       false); // fetch_one_additional_entry

          // If this node has child tiles, then it is not a leaf node.
          if (tile_records.size() > 0)
            return false;
      }
    }

    // If not child tiles were found, then this is a leaf node.
    return true;
  }

  template <class PixelT>
  int SnapshotManager<PixelT>::snapshot_helper(int current_col,
        int current_row,
        int current_level,
        std::map<TransactionOrNeg, TileHeader> composite_tiles,
        vw::BBox2i const& target_region,
        int target_level,
        TransactionOrNeg start_transaction_id,
        TransactionOrNeg end_transaction_id,
        Transaction write_transaction_id,
        bool tweak_settings_for_terrain) const
  {

    typedef std::map<TransactionOrNeg, TileHeader> TileMap;
    typedef std::list<TileHeader> TileList;

    // Check to see if there are any tiles at this level that need to be
    // snapshotted.  We fetch one additional tile outside of the
    // specified range so that we get any tiles that were not actually
    // composited during the last snapshot.
    TileList tile_records;
    tile_records = m_platefile->search_by_location(current_col, current_row, current_level,
                                                   start_transaction_id, end_transaction_id,
                                                   true); // fetch_one_additional_entry


    // If there are no valid tiles at this level, then there is nothing
    // further for us to do here on this branch of the recursion.  We
    // return immediately.
    if (tile_records.empty())
      return 0;

    // Insert the TileHeaders for this level into the composite_tiles
    // map<>.  Elements will be inserted in order, and any duplicate
    // transaction_ids will be overwritten by new tile records.  In
    // this fashion, we build up a list of the highest resolution
    // valid tiles in order of transaction id.
    for (TileList::const_iterator iter = tile_records.begin();
         iter != tile_records.end(); ++iter) {

      // We always add any valid tiles at the target_level, since
      // they should always be included in the composite.
      //
      // When adding tiles that may need to be supersampled, we are
      // only interested in leaf nodes.  All others will have higher
      // resolution tiles available as their children that we should
      // use instead.
      if (current_level == target_level || is_leaf(m_platefile, *iter) ) {

          vw_out(DebugMessage, "plate::snapshot")
            << "Adding tile @ " << iter->transaction_id() << " : "
            << " [ " << iter->transaction_id() << " ]  "
            << iter->col() << " " << iter->row() << " @ " << iter->level() << "\n";
          composite_tiles[iter->transaction_id()] = *iter;

      }

    }

    // If we have reached the target level, then we use the
    // accumulated tiles to generate a snapshot for this tile.
    if (current_level == target_level) {

      // Cull out "extra" extra tiles.  We need to do this here
      // because *each* search_by_location() call above can
      // potentially return it's own additional tile, and we only need
      // the most recent additional tile in the stack to include in
      // the composite.
      TileMap::iterator cull_iter = composite_tiles.begin();
      TransactionOrNeg extra_tid_id = -1;
      while (cull_iter != composite_tiles.end()) {

        // The rest of the tiles are not "extra," so we stop searching
        // for extras here.
        if (cull_iter->first >= start_transaction_id)
          break;

        // Search for the most recent transaction id that isn't
        // greater than the start_transaction_id.  This will be the
        // "extra" transaction id that we may want to include from an
        // earlier snapshot.
        if (cull_iter->first > extra_tid_id) {
          extra_tid_id = cull_iter->first;
        }
        ++cull_iter;
      }

      if (extra_tid_id != -1) {
        cull_iter  = composite_tiles.begin();
        while (cull_iter != composite_tiles.end()) {
          TileMap::iterator current_iter = cull_iter;
          ++cull_iter;

          // Once again, we are only interested in considering "extra"
          // tiles here, and not any that actually fall wihin the
          // valid range.
          if (current_iter->first >= start_transaction_id)
            break;

          // Delete any "extra" extra tiles, keeping only one that is
          // the best match for this snapshot.
          if (current_iter->first != extra_tid_id) {
            vw_out(DebugMessage, "plate::snapshot")
              << "Culling extraneous extra tiles that fall outside of transaction range:"
              << current_iter->first << " @ " << current_iter->second.level()
              << " (extra_tid_id = " << extra_tid_id << ")\n";
            composite_tiles.erase(current_iter);
          }

        }
      }

      // Check to see if the second entry is the start_transaction_id.
      // If this is the case, then we can safely skip the first tile
      // (which is then assumed to be the "additional_entry" that
      // falls earlier than the given transaction range) since it
      // will have already been incorporated into the previous
      // snapshot.
      if (composite_tiles.size() >= 2) {
        TileMap::iterator lowest_tid_iter = composite_tiles.begin();
        TileMap::const_iterator second_tid_iter = lowest_tid_iter;
        ++second_tid_iter;
        if (second_tid_iter->first == start_transaction_id) {
          vw_out(DebugMessage, "plate::snapshot")
            << "Culling tile that was already part of the last snapshot: "
            << lowest_tid_iter->first << " @ " << lowest_tid_iter->second.level() << "\n";
          composite_tiles.erase(lowest_tid_iter);
        }
      }

      // If we arrive at the target level and find that we have fewer
      // than two tiles to composite, then there's no actual
      // compositing that needs to be done.  We can safely bail on
      // this tile.
      if (composite_tiles.size() < 2)
        return 0;

      // Create a tile into which we will accumulate composited data.
      ImageView<PixelT> composite_tile(m_platefile->default_tile_size(),
                                       m_platefile->default_tile_size());

      vw_out(DebugMessage, "plate::snapshot") << "Starting Compositing run...\n";
      int num_composited = 0;
      for (TileMap::reverse_iterator iter = composite_tiles.rbegin();
           iter != composite_tiles.rend(); ++iter) {

        TileHeader &current_hdr = iter->second;
        ImageView<PixelT> new_tile(m_platefile->default_tile_size(),
                                   m_platefile->default_tile_size());

        // Tile access is highly localized during snapshotting, so it
        // speeds things up considerably to cache them here.  We check
        // the cache first and then read the tile from the platefile
        // if nothing is found.
        if ( !(this->restore_tile(new_tile, current_hdr.col(), current_hdr.row(),
                                  current_hdr.level(), current_hdr.transaction_id())) ) {
          try {
            // If the tile isn't in our cache, then we take the
            // performance hit and read the tile from the platefile.
            m_platefile->read(new_tile,
                            current_hdr.col(), current_hdr.row(),
                            current_hdr.level(), current_hdr.transaction_id(),
                            true); // exact_transaction_match
          } catch (BlobIoErr &e) {
            // If we get a BlobIO error, that's bad news, but not worth
            // killing the snapshot for.  Instead we log the error here and move
            // onto the next location.
            std::ostringstream ostr;
            ostr << "WARNING: a BlobIoErr occured while reading tile [" << current_hdr.col()
                 << " " << current_hdr.row() << "] @ " << current_hdr.level()
                 << " (t_id = " << current_hdr.transaction_id() << "): " << e.what();
            m_platefile->log(ostr.str());
            vw_out(ErrorMessage) << ostr.str() << "\n";
          } catch (IOErr &e) {
            // If we get a IoErr error, that's bad news, but not worth
            // killing the snapshot for.  Instead we log the error here and move
            // onto the next location.
            std::ostringstream ostr;
            ostr << "WARNING: an IOErr occurred while reading tile [" << current_hdr.col()
                 << " " << current_hdr.row() << "] @ " << current_hdr.level()
                 << " (t_id = " << current_hdr.transaction_id() << "): " << e.what();
            m_platefile->log(ostr.str());
            vw_out(ErrorMessage) << ostr.str() << "\n";
          }

          // Save the tile to the read cache.
          this->save_tile(new_tile, current_hdr.col(), current_hdr.row(),
                          current_hdr.level(), current_hdr.transaction_id());

        }

        vw_out(DebugMessage, "plate::snapshot")
          << "\t--> Adding tile: " << current_hdr.col() << " " << current_hdr.row()
          << " @ " << current_hdr.level() << "\n";

        // If this is the first tile in this location, it's already at
        // the target_level, and it's opaque (which is a very common
        // case), then we don't need to save it as part of the
        // snapshot since its already the top tile in the snapshot.
        // This optimization saves us both space and time!
        if ( current_hdr.level() == target_level &&
             iter == composite_tiles.rbegin() &&
             is_opaque(new_tile) ) {
          return 0;
        }

        // If we are compositing a tile from a lower resolution in the
        // pyramid, it needs to be supersampled and cropped before we
        // can use it here.
        if (current_hdr.level() != target_level) {

          new_tile = resample_img_from_level(
                      new_tile, current_hdr.col(), current_hdr.row(), current_hdr.level(),
                      current_col, current_row, target_level);


        }

        // Add the new tile UNDERNEATH the tiles we have already composited.
        //
        // If we are snapshotting terrain, we do this with a weighted
        // average.  Otherwise we use normal alpha blending.
        if (tweak_settings_for_terrain) {
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
        //---

        // for testing purposes
        // std::cout << "Writing dummy tiles " << current_col << " " << current_row << " @ " << current_level << " for t_id = " << write_transaction_id << "\n";
        // ImageView<PixelRGBA<float> > test_tile(256,256);
        // for (int j = 0; j < test_tile.rows(); ++j) {
        //   for (int i = 0; i < test_tile.cols(); ++i) {
        //     if (abs(i-j) < 10)
        //       test_tile(i,j) = PixelRGBA<float>(1.0,0,0,1.0);
        //   }
        // }
        //        m_platefile->write_update(test_tile,

        m_platefile->write_update(composite_tile,
                                  current_col, current_row, current_level,
                                  write_transaction_id);
        return 1;
      } else {
        return 0;
      }

    } else {

      // If we haven't reached the target level, then are four
      // possible children that we can explore for more snapshotting
      // information.
      int num_tiles_updated = 0;
      for (int new_row = current_row*2; new_row < current_row*2+2; ++new_row) {
        for (int new_col = current_col*2; new_col < current_col*2+2; ++new_col) {

          // We should only explore those child tiles that overlap
          // with the target_region at the target level.  To test
          // this, we create a bounding box for the child tile we are
          // testing, and scale that bounding box up to see how big it
          // is when projected into the target_level of the pyramid.
          BBox2i parent_region(new_col, new_row, 1, 1);
          parent_region.min() *= pow(2,target_level-(current_level+1));
          parent_region.max() *= pow(2,target_level-(current_level+1));


          // With the scaled up parent region and the target region in
          // hand, we can do a simple intersection test to assess
          // whether this child tile is worth exploring.
          if (parent_region.intersects(target_region)) {
            num_tiles_updated += snapshot_helper(new_col, new_row, current_level+1,
                                                 composite_tiles, target_region, target_level,
                                                 start_transaction_id, end_transaction_id,
                                                 write_transaction_id,
                                                 tweak_settings_for_terrain);
          }
        }
      }
      return num_tiles_updated;
    }
  }
}} // namespace vw::platefile


template <class PixelT>
void vw::platefile::SnapshotManager<PixelT>::snapshot(int level, BBox2i const& tile_region,
                                                      TransactionOrNeg start_transaction_id,
                                                      TransactionOrNeg end_transaction_id,
                                                      Transaction write_transaction_id,
                                                      bool tweak_settings_for_terrain) const {

  // Subdivide the bbox into smaller workunits if necessary.
  // This helps to keep operations efficient.
  std::list<BBox2i> tile_workunits = bbox_tiles(tile_region, 1024, 1024);
  for ( std::list<BBox2i>::iterator region_iter = tile_workunits.begin();
        region_iter != tile_workunits.end(); ++region_iter) {

    // Create an empty list of composite tiles and then kick off the
    // recursive snapshotting process.
    std::map<TransactionOrNeg, TileHeader> composite_tiles;
    int num_tiles_updated = snapshot_helper(0, 0, 0, composite_tiles, *region_iter, level,
                                            start_transaction_id, end_transaction_id,
                                            write_transaction_id,
                                            tweak_settings_for_terrain);

    if (num_tiles_updated > 0)
      vw_out() << "\t--> Snapshot " << *region_iter << " @ level " << level
               << " (" << num_tiles_updated << " tiles updated).\n";

  }
}

// Create a full snapshot of every level and every region in the mosaic.
template <class PixelT>
void vw::platefile::SnapshotManager<PixelT>::full_snapshot(TransactionOrNeg start_transaction_id,
                                                           TransactionOrNeg end_transaction_id,
                                                           Transaction write_transaction_id,
                                                           bool tweak_settings_for_terrain) const {

  for (int level = 0; level < m_platefile->num_levels(); ++level) {

    // For debugging:
    // for (int level = 0; level < 10; ++level) {

    // Snapshot the entire region at each level.  These region will be
    // broken down into smaller work units in snapshot().
    int region_size = 1 << level;
    int subdivided_region_size = region_size / 16;
    if (subdivided_region_size < 1024) subdivided_region_size = 1024;
    BBox2i full_region(0,0,region_size,region_size);
    std::list<BBox2i> workunits = bbox_tiles(full_region,
                                             subdivided_region_size,
                                             subdivided_region_size);
    for ( std::list<BBox2i>::iterator region_iter = workunits.begin();
          region_iter != workunits.end(); ++region_iter) {
      snapshot(level, *region_iter, start_transaction_id,
               end_transaction_id, write_transaction_id,
               tweak_settings_for_terrain);
    }
  }
}

// Explicit template instatiation
namespace vw {
namespace platefile {

  template
  void vw::platefile::SnapshotManager<PixelGrayA<uint8> >::snapshot(int level,
                                                                    BBox2i const& bbox,
                                                                    TransactionOrNeg start_transaction_id,
                                                                    TransactionOrNeg end_transaction_id,
                                                                    Transaction write_transaction_id,
                                                                    bool tweak_settings_for_terrain) const;
  template
  void vw::platefile::SnapshotManager<PixelGrayA<int16> >::snapshot(int level,
                                                                    BBox2i const& bbox,
                                                                    TransactionOrNeg start_transaction_id,
                                                                    TransactionOrNeg end_transaction_id,
                                                                    Transaction write_transaction_id,
                                                                    bool tweak_settings_for_terrain) const;
  template
  void vw::platefile::SnapshotManager<PixelGrayA<float32> >::snapshot(int level,
                                                                      BBox2i const& bbox,
                                                                      TransactionOrNeg start_transaction_id,
                                                                      TransactionOrNeg end_transaction_id,
                                                                      Transaction write_transaction_id,
                                                                      bool tweak_settings_for_terrain) const;
  template
  void vw::platefile::SnapshotManager<PixelRGBA<uint8> >::snapshot(int level,
                                                                   BBox2i const& bbox,
                                                                   TransactionOrNeg start_transaction_id,
                                                                   TransactionOrNeg end_transaction_id,
                                                                   Transaction write_transaction_id,
                                                                   bool tweak_settings_for_terrain) const;

  template
  void vw::platefile::SnapshotManager<PixelGrayA<uint8> >::full_snapshot(TransactionOrNeg start_transaction_id,
                                                                         TransactionOrNeg end_transaction_id,
                                                                         Transaction write_transaction_id,
                                                                         bool tweak_settings_for_terrain) const;
  template
  void vw::platefile::SnapshotManager<PixelGrayA<int16> >::full_snapshot(TransactionOrNeg start_transaction_id,
                                                                         TransactionOrNeg end_transaction_id,
                                                                         Transaction write_transaction_id,
                                                                         bool tweak_settings_for_terrain) const;
  template
  void vw::platefile::SnapshotManager<PixelGrayA<float32> >::full_snapshot(TransactionOrNeg start_transaction_id,
                                                                           TransactionOrNeg end_transaction_id,
                                                                           Transaction write_transaction_id,
                                                                           bool tweak_settings_for_terrain) const;
  template
  void vw::platefile::SnapshotManager<PixelRGBA<uint8> >::full_snapshot(TransactionOrNeg start_transaction_id,
                                                                        TransactionOrNeg end_transaction_id,
                                                                        Transaction write_transaction_id,
                                                                        bool tweak_settings_for_terrain) const;



}}
