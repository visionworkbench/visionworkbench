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


#ifndef __VW_PLATEFILE_INDEX_H__
#define __VW_PLATEFILE_INDEX_H__

#include <vw/Plate/FundamentalTypes.h>
#include <vw/Plate/IndexData.pb.h>
#include <vw/Plate/IndexDataPrivate.pb.h>
#include <vw/Image/PixelTypeInfo.h>
#include <vw/Math/BBox.h>
#include <boost/shared_ptr.hpp>
#include <list>

#define VW_PLATE_INDEX_VERSION 3

namespace vw {
namespace platefile {
  class Url;

namespace detail {
  class IndexPage;

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
    virtual std::ostream& log() = 0;

    // -------------------------- I/O ---------------------------

    /// Grab an IndexPage.  Useful if you want to serialize it by hand
    /// to disk.
    virtual boost::shared_ptr<IndexPage> page_request(uint32 col, uint32 row, uint32 level) const = 0;

    virtual uint64 page_id(uint32 col, uint32 row, uint32 level) const = 0;

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
    virtual IndexRecord read_request(uint32 col, uint32 row, uint32 depth, TransactionOrNeg transaction_id, bool exact_transaction_match = false) = 0;

    /// Writing, pt. 1: Locks a blob and returns the blob id that can
    /// be used to write a tile.
    virtual uint32 write_request() = 0;

    /// Writing, pt. 2: Supply information to update the index and
    /// unlock the blob id.
    virtual void write_update(TileHeader const& header, IndexRecord const& record) = 0;

    /// Writing, pt. 3: Signal the completion of the write operation.
    virtual void write_complete(uint32 blob_id) = 0;


    // ----------------------- PROPERTIES  ----------------------

    /// Returns a list of valid tiles that match this level, region, and
    /// range of transaction_id's.  Returns a list of TileHeaders with
    /// col/row/level and transaction_id of the most recent tile at each
    /// valid location.  Note: there may be other tiles in the transaction
    /// range at this col/row/level, but valid_tiles() only returns the
    /// first one.
    /// Note: the region is EXCLUSIVE: i.e. BBox2i(0,0,1,1) does not include the point (1,1)
    virtual std::list<TileHeader> search_by_region(uint32 level, vw::BBox2i const& region,
                                                   TransactionOrNeg start_transaction_id,
                                                   TransactionOrNeg end_transaction_id) const = 0;

    /// Return multiple tile headers that match the specified
    /// transaction id range.  This range is inclusive at both ends.
    virtual std::list<TileHeader> search_by_location(uint32 col, uint32 row, uint32 level,
                                                     TransactionOrNeg start_transaction_id,
                                                     TransactionOrNeg end_transaction_id) const = 0;

    virtual IndexHeader index_header() const = 0;

    virtual uint32 platefile_id() const = 0;
    virtual uint32 version() const = 0;
    virtual uint32 num_levels() const = 0;

    virtual std::string platefile_name() const = 0;

    virtual uint32 tile_size() const = 0;
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

}}}

#endif // __VW_PLATEFILE_INDEX_H__
