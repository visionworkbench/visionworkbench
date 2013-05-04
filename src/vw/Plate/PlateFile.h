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


/// \file PlateFile.h
///
/// Once written, tiles in the plate file are never deleted.  New tiles
/// can supercede old ones, but all of the data remains in the blob files,
/// and the index maintains a sorted list of pointers to previous versions
/// of a tile.  The sorting key for tiles is a field in each index entry
/// called 'transaction_id', which works much like a timestamp to keep
/// track of the order with which tiles were added to the platefile
/// system.
///
/// transaction_id is not simply a timestamp.  Instead, it is a unique
/// integer assigned to each new image that is added to the mosaic.
/// transaction_ids are incremented as each new image is added.
///
/// READING TILES using the transaction_id
/// --------------------------------------
///
/// Read requests can come from the mod_plate apache module (in response
/// to a web request) or from a mosaicking client (in the process of
/// compositing new tiles in the mosaic).  In both cases, it is important
/// for a read request to return tiles from a consistent version of the
/// mosaic.  In other words, we don't want to read tiles from a region
/// where tiles are being actively written -- doing so would result in a
/// patchwork of new and old tiles, and that would look very bad(tm).
///
/// transaction_ids allow you to go back to the latest consistent version
/// (or any previous version) of the mosaic.  Each read request can
/// specify an transaction_id, thereby forcing the server to access the
/// index as it would have appeared ON OR BEFORE the requested
/// transaction_id.  The index server keeps track of the latest "complete"
/// transaction ID (lets call it the transaction_cursor), and will furnish
/// that to clients if they are simply interested in accessing the latest
/// complete version of the mosaic.
///
///
/// WRITING TILES using the transaction_id
/// --------------------------------------
///
/// transaction_id's are assigned on a per-image (i.e. per transaction)
/// basis in the mosaic.  Prior to writing any tiles, the mosaicking
/// client is expected to send a TRANSACTION_REQUEST message to obtain a
/// transaction_id from the index server.  (Note to self: that request
/// should contain some indentifying information that is written to the
/// index log so that we can reconstruct what happened later if the
/// transaction fails.)  The index server doles these out in increasing
/// order so that each client gets a unique transaction id, and that these
/// ids are issued in the order that they were requested.
///
/// Mosaicking clients should also request the transaction_id pointed to
/// by the transaction_cursor (see above) so that it has a consistent
/// version of the mosaic to use for its read requests during compositing.
/// Both the read_transaction_id and write_transaction_id are stored for
/// use in subsequent mosaicking steps.
///
/// When the mosaicking client completes its work for a given transaction,
/// it is expected to send a TRANSACTION_COMPLETE message to inform the
/// index server of its succesful completion.  The index server uses this
/// information to update the transaction_cursor that can be used for read
/// requests.
///
///
/// CAVEATS
/// -------
///
/// - What happens if a mosaicking client fails in mid-mosaic and never
///   issues a TRANSACTION_COMPLETE message?
///
///   This is bad!!  The transaction_cursor will get stuck waiting for the
///   half-finished transaction ID.
///
///   Possbile solutions: keep track of the number of completed
///  transactions, and force the transaction ID to complete if more than
///   a certain number of transaction have elapsed.  (i.e. assume it has
///   failed and note this failure in the index server log.)
///
///   (But then we should make sure that a TRANSACTION_COMPLETE message
///   that comes back from the dead doesn't somehow mess up the
///   transaction_cursor!)
///
/// - What happens to the transaction_cursor if the index server crashes?
///
///   Transactions updates don't happen very often.  We can afford to
///   write the new transaction cursor location to the index file on
///   disk before we actually update it in the live index.
///

#ifndef __VW_PLATE_PLATEFILE_H__
#define __VW_PLATE_PLATEFILE_H__

#include <vw/Plate/FundamentalTypes.h>
#include <vw/Plate/Datastore.h>
#include <vw/Plate/Blob.h>
#include <vw/Plate/Exception.h>
#include <vw/Plate/HTTPUtils.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/FileIO/MemoryImageResource.h>

#include <vw/Image/ImageView.h>
#include <vw/Image/Algorithms.h>

#include <sstream>
#include <boost/tuple/tuple.hpp>

namespace vw {
namespace platefile {

  class ReadOnlyPlateFile {
    protected:
      boost::shared_ptr<Datastore> m_data;
      ReadOnlyPlateFile(const Url& url, std::string type, std::string description, uint32 tile_size, std::string tile_filetype, PixelFormatEnum pixel_format, ChannelTypeEnum channel_type);

    public:
      ReadOnlyPlateFile(const Url& url);

      /// The destructor saves the platefile to disk.
      virtual ~ReadOnlyPlateFile() {}

      IndexHeader index_header() const;
      std::string default_file_type() const;
      uint32 default_tile_size() const;
      PixelFormatEnum pixel_format() const;
      ChannelTypeEnum channel_type() const;
      uint32 num_levels() const;

      boost::tuple<std::string, TileHeader, TileData>
      img_file_name(std::string const& base_name, int col, int row, int level,
                    TransactionOrNeg transaction_id, bool exact_transaction_match = false) const;

      /// Read data directly to a file on disk. You supply a base name
      /// (without the file's image extension).  The image extension
      /// will be appended automatically for you based on the filetype
      /// in the TileHeader.
      std::pair<std::string, TileHeader>
      read_to_file(std::string const& base_name, int col, int row, int level,
                   TransactionOrNeg transaction_id, bool exact_transaction_match = false) const;

      std::pair<TileHeader, TileData>
      read(int col, int row, int level, TransactionOrNeg transaction_id, bool exact_transaction_match = false) const;

      Datastore::TileSearch&
      batch_read(Datastore::TileSearch&) const;

      /// Read an image from the specified tile location in the plate file.
      ///
      /// By default, this call to read will return a tile with the MOST
      /// RECENT transaction_id <= to the transaction_id you specify
      /// here in the function arguments (if a tile exists).  However,
      /// setting exact_transaction_match = true will force the
      /// PlateFile to search for a tile that has the EXACT SAME
      /// transaction_id as the one that you specify.
      ///
      /// A transaction ID of -1 indicates that we should return the most recent
      /// tile, regardless of its transaction id. This trumps everything including
      /// exact_transaction_match.
      template <class ViewT>
      TileHeader read(ViewT &view, int col, int row, int level,
                      TransactionOrNeg transaction_id, bool exact_transaction_match = false) const {

        std::pair<TileHeader, TileData> ret = this->read(col, row, level, transaction_id, exact_transaction_match);
        boost::scoped_ptr<SrcImageResource> r(SrcMemoryImageResource::open(ret.first.filetype(), &ret.second->operator[](0), ret.second->size()));
        read_image(view, *r);
        return ret.first;
      }

      /// Returns a list of valid tiles that match this level, region, and
      /// range of transaction_id's.  Returns a list of TileHeaders with
      /// col/row/level and transaction_id of every tile at each valid location within the range.
      /// Note: the region is EXCLUSIVE: i.e. BBox2i(0,0,1,1) does not include the point (1,1)
      std::list<TileHeader> search_by_region(int level, vw::BBox2i const& region, const TransactionRange& range) const;

      /// Read one ore more images at a specified location in the
      /// platefile by specifying a range of transaction ids of
      /// interest.  This range is inclusive at both ends.
      std::list<TileHeader>
      search_by_location(int col, int row, int level, const TransactionRange& range);
  };

  class PlateFile : public ReadOnlyPlateFile {
      boost::shared_ptr<Transaction> m_transaction;
      boost::shared_ptr<WriteState> m_write_state;
    public:
      PlateFile(const Url& url);

      PlateFile(const Url& url, std::string type, std::string description,
                uint32 tile_size, std::string tile_filetype,
                PixelFormatEnum pixel_format, ChannelTypeEnum channel_type);

      void sync() const;

      std::ostream& audit_log();
      std::ostream& error_log();
      void log(std::string message) VW_DEPRECATED;

      /// Writing, pt. 1: Locks a blob and returns the blob id that can
      /// be used to write tiles.
      void write_request();

      /// Writing, pt. 2: Write an image to the specified tile location
      /// in the plate file.
      template <class ViewT>
      void write_update(ImageViewBase<ViewT> const& view, int col, int row, int level) {

        std::string type = this->default_file_type();
        if (type == "auto") {
          // This specialization saves us TONS of space by storing opaque tiles
          // as jpgs.  However it does come at a small cost of having to conduct
          // this extra check to see if the tile is opaque or not.
          if ( is_opaque(view.impl()) )
            type = "jpg";
          else
            type = "png";
        }
        boost::scoped_ptr<DstMemoryImageResource> r(DstMemoryImageResource::create(type, view.format()));
        write_image(*r, view);
        this->write_update(r->data(), r->size(), col, row, level, type);
      }

      /// Writing, pt. 2, alternate: Write raw data (as a tile) to a specified
      /// tile location. Use the filetype to identify the data later; empty type
      /// means "platefile default".
      void write_update(const uint8* data, uint64 data_size, int col, int row, int level, const std::string& type = "");

      /// Writing, pt. 3: Signal the completion of the write operation.
      void write_complete();

      // --------------------- TRANSACTIONS ------------------------

      // Clients are expected to make a transaction request whenever
      // they start a self-contained chunk of mosaicking work.  .
      virtual const Transaction& transaction_begin(const std::string& description, TransactionOrNeg transaction_id_override);

      // Resume a previous transaction (useful for multi-process runs). You are
      // responsible for making sure transaction_begin was run previously.
      virtual void transaction_resume(const Transaction& tid);

      // Once a chunk of work is complete, clients can "commit" their
      // work to the mosaic by issuding a transaction_complete method.
      virtual void transaction_end(bool update_read_cursor);

      virtual const Transaction& transaction_id() const;
  };

}} // namespace vw::plate

#endif // __VW_PLATE_PLATEFILE_H__
