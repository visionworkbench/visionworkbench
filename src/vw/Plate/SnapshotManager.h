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
  class ReadOnlyPlateFile;

  template <class PixelT>
  class SnapshotManager {
    protected:
      boost::shared_ptr<ReadOnlyPlateFile> m_read_plate;
      boost::shared_ptr<PlateFile> m_write_plate;

    public:
      SnapshotManager(boost::shared_ptr<ReadOnlyPlateFile> read_plate, boost::shared_ptr<PlateFile> write_plate, bool tweak_settings_for_terrain)
        : m_read_plate(read_plate), m_write_plate(write_plate) {
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
    void snapshot(uint32 level, BBox2i const& bbox, TransactionRange read_transaction_range, const ProgressCallback &progress = ProgressCallback::dummy_instance()) const;

    // Create a full snapshot of every level and every region in the mosaic.
    void full_snapshot(TransactionRange read_transaction_range) const;
  };

}} // namespace vw::plate

#endif // __VW_PLATE_SNAPSHOT_MANAGER_H__
