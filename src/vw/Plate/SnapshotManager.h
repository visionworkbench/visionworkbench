// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_PLATE_SNAPSHOT_MANAGER_H__
#define __VW_PLATE_SNAPSHOT_MANAGER_H__

#include <vw/Plate/PlateFile.h>

namespace vw {
namespace platefile {

  // -------------------------------------------------------------------------
  //                              SNAPSHOT MANAGER
  // -------------------------------------------------------------------------

  template <class PixelT>
  class SnapshotManager {

  protected:
    boost::shared_ptr<PlateFile> m_platefile;

  public:

    SnapshotManager(boost::shared_ptr<PlateFile> platefile) : m_platefile(platefile) {}

    // ---------------------------- SNAPSHOTTING --------------------------------

    // snapshot() creates a complete, composited view of the mosaic.
    //
    //   level -- select the pyramid level on which to carry out mipmapping
    //   bbox -- bounding box (in terms of tiles) containing the tiles that need 
    //           to be snapshotted at starting_level.  Use to specify affected tiles.
    //   start_transaction_id -- select a transaction_id to use when accessing tiles.
    //   end_transaction_id -- select a transaction_id to use when accessing tiles.
    //   write_transaction_id -- the t_id to use for writing the snapshotted tiles.
    //
    void snapshot(int level, BBox2i const& bbox, 
                  int start_transaction_id, int end_transaction_id, 
                  int write_transaction_id) const;

    // Create a full snapshot of every level and every region in the mosaic.
    void full_snapshot(int start_transaction_id, 
                       int end_transaction_id, 
                       int write_transaction_id) const;
      
  };    

}} // namespace vw::plate

#endif // __VW_PLATE_SNAPSHOT_MANAGER_H__
