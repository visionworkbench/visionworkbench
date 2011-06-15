// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_PLATE_SNAPSHOT_MANAGER_H__
#define __VW_PLATE_SNAPSHOT_MANAGER_H__

#include <vw/Plate/IndexData.pb.h>
#include <vw/Plate/Datastore.h>
#include <vw/Plate/FundamentalTypes.h>
#include <vw/Image/ImageView.h>

namespace vw {
namespace platefile {

  class PlateFile;

  template <class PixelT>
  class SnapshotManager {
    protected:
      boost::shared_ptr<PlateFile> m_platefile;

    public:
      SnapshotManager(boost::shared_ptr<PlateFile> platefile, bool tweak_settings_for_terrain)
        : m_platefile(platefile) {
        VW_ASSERT(!tweak_settings_for_terrain, NoImplErr() << "Terrain snapshotting is currently unimplemented");
      }

    // ---------------------------- SNAPSHOTTING --------------------------------

    // snapshot() creates a complete, composited view of the mosaic.
    //
    //   level -- select the pyramid level on which to carry out mipmapping
    //   bbox -- bounding box (in terms of tiles) containing the tiles that need
    //           to be snapshotted at starting_level.  Use to specify affected tiles.
    //   start_transaction_id -- select a transaction_id to use when accessing tiles.
    //   end_transaction_id -- select a transaction_id to use when accessing tiles.
    //
    void snapshot(uint32 level, BBox2i const& bbox, TransactionRange read_transaction_range) const;

    // Create a full snapshot of every level and every region in the mosaic.
    void full_snapshot(TransactionRange read_transaction_range) const;
  };

}} // namespace vw::plate

#endif // __VW_PLATE_SNAPSHOT_MANAGER_H__
