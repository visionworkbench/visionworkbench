// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_PLATEFILE_INDEX_H__
#define __VW_PLATEFILE_INDEX_H__

#include <vw/Plate/FundamentalTypes.h>
#include <vw/Plate/IndexData.pb.h>
#include <vw/Image/PixelTypeInfo.h>
#include <vw/Math/BBox.h>
#include <boost/shared_ptr.hpp>
#include <list>

#define VW_PLATE_INDEX_VERSION 3

namespace vw {
namespace platefile {

  class IndexPage;
  class Url;

  // -------------------------------------------------------------------
  //                          INDEX BASE CLASS
  // -------------------------------------------------------------------
  class Index {
  public:

    /// Static methods for generating an index of the proper type
    /// (i.e. local vs. remote) based on a URL string.  Any URL string
    /// starting with "pf://..." will cause a remote index to be
    /// created.  All other URLs will create a local index.
    ///
    /// This method opens an existing index.
    static boost::shared_ptr<Index> construct_open(const Url& url);

    /// Static methods for generating an index of the proper type
    /// (i.e. local vs. remote) based on a URL string.  Any URL string
    /// starting with "pf://..." will cause a remote index to be
    /// created.  All other URLs will create a local index.
    ///
    /// This method creates a new index.
    static boost::shared_ptr<Index> construct_create(const Url& url, const IndexHeader& new_index_info);

    /// Destructor
    virtual ~Index() {}

    /// Sync any unsaved data in the index to disk.
    virtual void sync() = 0;

    /// Log a message to the platefile log.
    virtual void log(std::string message) = 0;

    // -------------------------- I/O ---------------------------

    /// Grab an IndexPage.  Useful if you want to serialize it by hand
    /// to disk.
    virtual boost::shared_ptr<IndexPage> page_request(int col, int row, int level) const = 0;

    /// Attempt to access a tile in the index.  Throws an
    /// TileNotFoundErr if the tile cannot be found.
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
    virtual IndexRecord read_request(int col, int row, int depth, TransactionOrNeg transaction_id, bool exact_transaction_match = false) = 0;

    /// Writing, pt. 1: Locks a blob and returns the blob id that can
    /// be used to write a tile.
    virtual int write_request(uint64 &size) = 0;

    /// Writing, pt. 2: Supply information to update the index and
    /// unlock the blob id.
    virtual void write_update(TileHeader const& header, IndexRecord const& record) = 0;

    /// Writing, pt. 3: Signal the completion of the write operation.
    virtual void write_complete(int blob_id, uint64 blob_offset) = 0;


    // ----------------------- PROPERTIES  ----------------------

    /// Returns a list of valid tiles that match this level, region, and
    /// range of transaction_id's.  Returns a list of TileHeaders with
    /// col/row/level and transaction_id of the most recent tile at each
    /// valid location.  Note: there may be other tiles in the transaction
    /// range at this col/row/level, but valid_tiles() only returns the
    /// first one.
    virtual std::list<TileHeader> search_by_region(int level, vw::BBox2i const& region,
                                                   TransactionOrNeg start_transaction_id,
                                                   TransactionOrNeg end_transaction_id,
                                                   uint32 min_num_matches,
                                                   bool fetch_one_additional_entry = false) const = 0;

    /// Return multiple tile headers that match the specified
    /// transaction id range.  This range is inclusive of the first
    /// entry, but not the last entry: [ begin_transaction_id, end_transaction_id )
    virtual std::list<TileHeader> search_by_location(int col, int row, int level,
                                                     TransactionOrNeg start_transaction_id,
                                                     TransactionOrNeg end_transaction_id,
                                                     bool fetch_one_additional_entry = false) const = 0;

    virtual IndexHeader index_header() const = 0;

    virtual int32 version() const = 0;
    virtual int32 num_levels() const = 0;

    virtual std::string platefile_name() const = 0;

    virtual int32 tile_size() const = 0;
    virtual std::string tile_filetype() const = 0;

    virtual PixelFormatEnum pixel_format() const = 0;
    virtual ChannelTypeEnum channel_type() const = 0;

    // --------------------- TRANSACTIONS ------------------------

    // Clients are expected to make a transaction request whenever
    // they start a self-contained chunk of mosaicking work.  Use
    // transaction_id_override to force the use of a transaction ID
    // for an upcoming transaction.  Setting transaction_id_override
    // to -1 lets the platefile choose its own transaction_id.
    virtual Transaction transaction_request(std::string transaction_description,
                                            TransactionOrNeg transaction_id_override) = 0;

    /// Once a chunk of work is complete, clients can "commit" their
    /// work to the mosaic by issuding a transaction_complete method.
    virtual void transaction_complete(Transaction transaction_id, bool update_read_cursor) = 0;

    // If a transaction fails, we may need to clean up the mosaic.
    virtual void transaction_failed(Transaction transaction_id) = 0;

    virtual Transaction transaction_cursor() = 0;

  };

}} // namespace vw::plate

#endif // __VW_PLATEFILE_INDEX_H__
