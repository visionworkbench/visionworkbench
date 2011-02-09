// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
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
#include <vw/Plate/Index.h>
#include <vw/Plate/Blob.h>
#include <vw/Plate/Exception.h>
#include <vw/Plate/HTTPUtils.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/FileIO/MemoryImageResource.h>

#include <vw/Image/ImageView.h>
#include <vw/Image/Algorithms.h>

#include <sstream>

namespace vw {
namespace platefile {

  class PlateFile {
    boost::shared_ptr<Index> m_index;
    boost::shared_ptr<Blob> m_write_blob;
    int m_write_blob_id;

  public:
    PlateFile(const Url& url);

    PlateFile(const Url& url, std::string type, std::string description,
              int tile_size, std::string tile_filetype,
              PixelFormatEnum pixel_format, ChannelTypeEnum channel_type);

    /// The destructor saves the platefile to disk.
    virtual ~PlateFile() {}

    /// Returns the name of the root directory containing the plate file.
    std::string name() const { return m_index->platefile_name(); }

    /// Returns the name of the root directory containing the plate file.
    IndexHeader index_header() const { return m_index->index_header(); }

    /// Returns the file type used to store tiles in this plate file.
    std::string default_file_type() const { return m_index->tile_filetype(); }

    int default_tile_size() const { return m_index->tile_size(); }

    PixelFormatEnum pixel_format() const { return m_index->pixel_format(); }

    ChannelTypeEnum channel_type() const { return m_index->channel_type(); }

    int num_levels() const { return m_index->num_levels(); }

    void sync() const { m_index->sync(); }

    void log(std::string message) { m_index->log(message); }

    /// Read data directly to a file on disk. You supply a base name
    /// (without the file's image extension).  The image extension
    /// will be appended automatically for you based on the filetype
    /// in the TileHeader.
    std::pair<std::string, TileHeader>
    read_to_file(std::string const& base_name, int col, int row, int level,
                 TransactionOrNeg transaction_id, bool exact_transaction_match = false) const;

    std::pair<TileHeader, TileData>
    read(int col, int row, int level, TransactionOrNeg transaction_id, bool exact_transaction_match = false) const;

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

    /// Writing, pt. 1: Locks a blob and returns the blob id that can
    /// be used to write tiles.
    void write_request();

    /// Writing, pt. 2: Write an image to the specified tile location
    /// in the plate file.
    template <class ViewT>
    void write_update(ImageViewBase<ViewT> const& view,
                      int col, int row, int level, Transaction transaction_id) {

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
      this->write_update(r->data(), r->size(), col, row, level, transaction_id, type);
    }

    /// Writing, pt. 2, alternate: Write raw data (as a tile) to a specified
    /// tile location. Use the filetype to identify the data later; empty type
    /// means "platefile default".
    void write_update(const uint8* data, uint64 data_size,
                      int col, int row, int level, Transaction transaction_id, const std::string& type = "");

    /// Writing, pt. 3: Signal the completion of the write operation.
    void write_complete();

    /// Read a record out of the platefile.
    ///
    /// By default, this call to read will return a tile with the MOST
    /// RECENT transaction_id <= to the transaction_id you specify
    /// here in the function arguments (if a tile exists).  However,
    /// setting exact_transaction_match = true will force the
    /// PlateFile to search for a tile that has the EXACT SAME
    /// transaction_id as the one that you specify.
    ///
    /// A transaction ID of -1 indicates that we should return the
    /// most recent tile, regardless of its transaction id.
    IndexRecord read_record(int col, int row, int level,
                            TransactionOrNeg transaction_id, bool exact_transaction_match = false);

    // --------------------- TRANSACTIONS ------------------------

    // Clients are expected to make a transaction request whenever
    // they start a self-contained chunk of mosaicking work.  .
    virtual Transaction transaction_request(std::string transaction_description,
                                            TransactionOrNeg transaction_id_override) {
      return m_index->transaction_request(transaction_description, transaction_id_override);
    }

    // Once a chunk of work is complete, clients can "commit" their
    // work to the mosaic by issuding a transaction_complete method.
    virtual void transaction_complete(Transaction transaction_id, bool update_read_cursor) {
      m_index->transaction_complete(transaction_id, update_read_cursor);
    }

    // If a transaction fails, we may need to clean up the mosaic.
    virtual void transaction_failed(Transaction transaction_id) {
      m_index->transaction_failed(transaction_id);
    }

    virtual Transaction transaction_cursor() {
      return m_index->transaction_cursor();
    }

    // ----------------------- UTILITIES --------------------------

    /// Returns a list of valid tiles that match this level, region, and
    /// range of transaction_id's.  Returns a list of TileHeaders with
    /// col/row/level and transaction_id of the most recent tile at each
    /// valid location.  Note: there may be other tiles in the transaction
    /// range at this col/row/level, but valid_tiles() only returns the
    /// first one.
    std::list<TileHeader> search_by_region(int level, vw::BBox2i const& region,
                                           TransactionOrNeg start_transaction_id,
                                           TransactionOrNeg end_transaction_id,
                                           uint32 min_num_matches,
                                           bool fetch_one_additional_entry = false) const {
      return m_index->search_by_region(level, region,
                                  start_transaction_id, end_transaction_id,
                                  min_num_matches, fetch_one_additional_entry);
    }

    /// Read one ore more images at a specified location in the
    /// platefile by specifying a range of transaction ids of
    /// interest.  This range is inclusive of the first entry, but not
    /// the last entry: [ begin_transaction_id, end_transaction_id )
    ///
    /// This is mostly useful when compositing tiles during mipmapping.
    std::list<TileHeader> search_by_location(int col, int row, int level,
                                             TransactionOrNeg begin_transaction_id, TransactionOrNeg end_transaction_id,
                                             bool fetch_one_additional_entry = false) {
      return m_index->search_by_location(col, row, level,
                                         begin_transaction_id, end_transaction_id,
                                         fetch_one_additional_entry);
    }

  };



}} // namespace vw::plate

#endif // __VW_PLATE_PLATEFILE_H__
