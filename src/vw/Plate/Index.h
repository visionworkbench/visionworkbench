// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_PLATEFILE_INDEX_H__
#define __VW_PLATEFILE_INDEX_H__

#include <vw/Core/FundamentalTypes.h>
#include <vw/Core/Log.h>
#include <vw/Image/PixelTypeInfo.h>
#include <vw/Math/BBox.h>
#include <vw/Plate/ProtoBuffers.pb.h>

#include <list>

#define VW_PLATE_INDEX_VERSION 3

namespace vw {
namespace platefile {

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
    static boost::shared_ptr<Index> construct_open(std::string url);

    /// Static methods for generating an index of the proper type
    /// (i.e. local vs. remote) based on a URL string.  Any URL string
    /// starting with "pf://..." will cause a remote index to be
    /// created.  All other URLs will create a local index.
    ///
    /// This method creates a new index.
    static boost::shared_ptr<Index> construct_create(std::string url, IndexHeader new_index_info);

    /// Destructor
    virtual ~Index() {}

    /// Sync any unsaved data in the index to disk.
    virtual void sync() = 0;

    // -------------------------- I/O ---------------------------

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
    virtual IndexRecord read_request(int col, int row, int depth, int transaction_id, bool exact_transaction_match = false) = 0;

    /// Return multiple index entries that match the specified
    /// transaction id range.  This range is inclusive of the first
    /// entry, but not the last entry: [ begin_transaction_id, end_transaction_id )
    ///
    /// Results are return as a std::pair<int32, IndexRecord>.  The
    /// first value in the pair is the transaction id for that
    /// IndexRecord.
    virtual std::list<std::pair<int32, IndexRecord> > multi_read_request(int col, int row, int level, 
                                                                         int begin_transaction_id, 
                                                                         int end_transaction_id) = 0;
  
    /// Writing, pt. 1: Locks a blob and returns the blob id that can
    /// be used to write a tile.
    virtual int write_request(int size) = 0;

    /// Writing, pt. 2: Supply information to update the index and
    /// unlock the blob id.
    virtual void write_complete(TileHeader const& header, IndexRecord const& record) = 0;

    // ----------------------- PROPERTIES  ----------------------

    /// Returns a list of valid tiles that match this level, region, and
    /// range of transaction_id's.  Returns a list of TileHeaders with
    /// col/row/level and transaction_id of the most recent tile at each
    /// valid location.  Note: there may be other tiles in the transaction
    /// range at this col/row/level, but valid_tiles() only returns the
    /// first one.
    virtual std::list<TileHeader> valid_tiles(int level, vw::BBox2i const& region,
                                              int start_transaction_id, 
                                              int end_transaction_id, 
                                              int min_num_matches) const = 0;

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
    virtual int32 transaction_request(std::string transaction_description,
                                      int transaction_id_override) = 0;

    /// Once a chunk of work is complete, clients can "commit" their
    /// work to the mosaic by issuding a transaction_complete method.
    virtual void transaction_complete(int32 transaction_id, bool update_read_cursor) = 0;

    // If a transaction fails, we may need to clean up the mosaic.  
    virtual void transaction_failed(int32 transaction_id) = 0;

    virtual int32 transaction_cursor() = 0;

  };

}} // namespace vw::plate

#endif // __VW_PLATEFILE_INDEX_H__
