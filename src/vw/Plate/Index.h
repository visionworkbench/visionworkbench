// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_PLATEFILE_INDEX_H__
#define __VW_PLATEFILE_INDEX_H__

#include <vw/Core/FundamentalTypes.h>
#include <vw/Core/Exception.h>
#include <vw/Core/Log.h>

#include <vw/Image/PixelTypeInfo.h>

#include <vw/Plate/Tree.h>
#include <vw/Plate/ProtoBuffers.pb.h>

#define VW_PLATE_INDEX_VERSION 2

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
  
    /// Writing, pt. 1: Locks a blob and returns the blob id that can
    /// be used to write a tile.
    virtual int write_request(int size) = 0;

    /// Writing, pt. 2: Supply information to update the index and
    /// unlock the blob id.
    virtual void write_complete(TileHeader const& header, IndexRecord const& record) = 0;

    // ----------------------- PROPERTIES  ----------------------

    virtual IndexHeader index_header() const = 0;

    virtual int32 version() const = 0;
    virtual int32 num_levels() const = 0;

    virtual std::string platefile_name() const = 0;

    virtual int32 tile_size() const = 0;
    virtual std::string tile_filetype() const = 0;

    virtual PixelFormatEnum pixel_format() const = 0;
    virtual ChannelTypeEnum channel_type() const = 0;

    // --------------------- TRANSACTIONS ------------------------

    /// Clients are expected to make a transaction request whenever
    /// they start a self-contained chunk of mosaicking work.  .
    virtual int32 transaction_request(std::string transaction_description,
                                      std::vector<TileHeader> const& tile_headers) = 0;

    /// Called right before the beginning of the mipmapping pass
    virtual void root_complete(int32 transaction_id,
                               std::vector<TileHeader> const& tile_headers) = 0;

    /// Once a chunk of work is complete, clients can "commit" their
    /// work to the mosaic by issuding a transaction_complete method.
    virtual void transaction_complete(int32 transaction_id) = 0;

    // If a transaction fails, we may need to clean up the mosaic.  
    virtual void transaction_failed(int32 transaction_id) = 0;

    virtual int32 transaction_cursor() = 0;

    // --------------------- UTILITIES ------------------------

    /// Iterate over all nodes in a tree, calling func for each
    /// location.  Note: this will only be implemented for local
    /// indexes.  This function will throw an error if called on a
    /// remote index.
    virtual void map(boost::shared_ptr<TreeMapFunc> func) { 
      vw_throw(NoImplErr() << "Index::map() not implemented for this index type.");
    }
  };

}} // namespace vw::plate

#endif // __VW_PLATEFILE_INDEX_H__
