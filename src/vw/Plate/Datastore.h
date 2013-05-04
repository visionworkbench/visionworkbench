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


#ifndef __VW_PLATE_DATASTORE_H__
#define __VW_PLATE_DATASTORE_H__

#include <vw/Plate/FundamentalTypes.h>
#include <vw/Plate/IndexData.pb.h>
#include <vw/Plate/HTTPUtils.h>
#include <vw/Image/PixelTypeInfo.h>
#include <vw/Math/BBox.h>
#include <vw/Core/FundamentalTypes.h>
#include <boost/shared_ptr.hpp>
#include <boost/shared_container_iterator.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/function.hpp>
#include <vector>

namespace vw {
namespace platefile {

typedef boost::shared_ptr<std::vector<vw::uint8> > TileData;
class TileHeader;

// An opaque type that data stores can use to hide their write state
class WriteState {
  public:
    virtual ~WriteState() {}
    // human-readable description
    virtual std::string what() const = 0;
};

struct Tile {
  TileHeader hdr;
  TileData   data;
  Tile() {}
  Tile(const TileHeader& hdr)
    : hdr(hdr) {}
  Tile(const TileHeader& hdr, TileData data)
    : hdr(hdr), data(data) {}
};

class Datastore {
  public:
    typedef std::vector<Tile> TileSearch;
    typedef boost::function<std::ostream&(void)> Logger;

    // without 'header', throws an exception if there is no datastore at that url.
    // with 'header',    if there is no datastore, creates it and uses header
    //                   as the defaults. if there is a datastore, leaves the defaults alone and
    //                   throws if the 'type' parameter does not match the existing format
    static Datastore* open(const Url& url);
    static Datastore* open(const Url& url, const IndexHeader& header);
    virtual ~Datastore() {}

    // TRANSACTION DATA
    virtual Transaction transaction_begin(const std::string& description, TransactionOrNeg transaction_id_override = -1) = 0;
    // Updates the write cursor. Also updates the read cursor if update_read_cursor is true.
    virtual void transaction_end(Transaction transaction_id, bool update_read_cursor) = 0;

    // TILE READ. limit is the limit within a single tile slot
    // head()'s return order will be in transaction order, descending.
    // get()'s  return order will be arbitrary. The return reference is the
    // same as the input reference, it's just there for code flow.
    /// Note: the region is EXCLUSIVE: i.e. BBox2i(0,0,1,1) does not include the point (1,1)
    virtual TileSearch& head(TileSearch& tiles, uint32 level, uint32 row, uint32 col, TransactionRange range, uint32 limit = 0) = 0;
    virtual TileSearch& head(TileSearch& tiles, uint32 level,   const BBox2u& region, TransactionRange range, uint32 limit = 0) = 0;
    virtual TileSearch&  get(TileSearch& tiles, uint32 level, uint32 row, uint32 col, TransactionRange range, uint32 limit = 0);
    virtual TileSearch&  get(TileSearch& tiles, uint32 level,   const BBox2u& region, TransactionRange range, uint32 limit = 0);

    // Return a url that should work to retrieve tile data [might be unimplemented]
    //virtual Url map_to_url(uint32 level, uint32 row, uint32 col, Transaction id, std::string filetype) = 0;
    //virtual Url map_to_url(const TileHeader& t);

    // The return of this is not guaranteed to be the same size as the input;
    // if errors occur on some of the fetches, they will be dropped. Order is
    // not guaranteed, either.
    virtual TileSearch& populate(TileSearch& hdrs) = 0;

    // TILE WRITE
    virtual WriteState* write_request(const Transaction& id) VW_WARN_UNUSED = 0;
    virtual void write_update(WriteState& state, uint32 level, uint32 row, uint32 col, const std::string& filetype, const uint8* data, uint64 size) = 0;
    virtual void write_complete(WriteState& id) = 0;
    virtual void flush() = 0;


    // DEFAULT METADATA

    // These functions are variant, and may cause network IO.
    virtual IndexHeader index_header() const = 0;
    virtual uint32 num_levels() const = 0;
    //virtual Transaction transaction_read_cursor() const = 0;
    //virtual Transaction transaction_write_cursor() const = 0;

    // These functions are invariant once the plate is constructed, and may be cached.
    virtual uint32 id() const = 0;
    virtual uint32 tile_size() const = 0;
    virtual std::string tile_filetype() const = 0;
    virtual PixelFormatEnum pixel_format() const = 0;
    virtual ChannelTypeEnum channel_type() const = 0;

    // LOGGING
    // Use the audit log for events that should be recorded for posterity next to the datastore
    virtual Logger audit_log() const = 0;
    // a sane thing for error_log to do is log to audit_log() and VW_OUT(ErrorMessage, "console")
    virtual Logger error_log() const = 0;
};

}}

#endif
