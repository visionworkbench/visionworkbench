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
