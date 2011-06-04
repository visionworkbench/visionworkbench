// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
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
};

class Datastore {
  public:
    typedef boost::shared_container_iterator<const std::vector<TileHeader> > meta_iterator;
    typedef boost::iterator_range<meta_iterator> meta_range;
    typedef boost::shared_container_iterator<const std::vector<Tile> > tile_iterator;
    typedef boost::iterator_range<tile_iterator> tile_range;

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

    // Read and Write cursors
    //virtual Transaction transaction_read_cursor() const = 0;
    //virtual Transaction transaction_write_cursor() const = 0;

    // TILE READ. limit is the limit within a single tile slot
    // head()'s return order will be in transaction order, descending.
    // get()'s  return order will be arbitrary.
    virtual meta_range     head(uint32 level, uint32 row, uint32 col, TransactionRange range, uint32 limit = 0) = 0;
    virtual meta_range     head(uint32 level,   const BBox2u& region, TransactionRange range, uint32 limit = 0) = 0;
    virtual tile_range      get(uint32 level, uint32 row, uint32 col, TransactionRange range, uint32 limit = 0);
    virtual tile_range      get(uint32 level,   const BBox2u& region, TransactionRange range, uint32 limit = 0);

    // The return of this is not guaranteed to be 'len' elements; if errors
    // occur on some of the fetches, they will be dropped. Order is not
    // guaranteed, either. 'hdrs' is non-const to allow implementations to sort
    // without copying; hdrs will not be resized, only swapped.
    virtual tile_range populate(TileHeader* hdrs, size_t len) = 0;

    // Return a url that should work to retrieve tile data [might be unimplemented]
    //virtual Url map_to_url(uint32 level, uint32 row, uint32 col, Transaction id, std::string filetype) = 0;
    //virtual Url map_to_url(const TileHeader& t);

    // TILE WRITE
    virtual WriteState* write_request(const Transaction& id) VW_WARN_UNUSED = 0;
    virtual void write_update(WriteState& state, uint32 level, uint32 row, uint32 col, const std::string& filetype, const uint8* data, uint64 size) = 0;
    virtual void write_complete(WriteState& id) = 0;
    virtual void flush() = 0;


    // DEFAULT METADATA
    virtual IndexHeader index_header() const = 0;
    //virtual std::string name() const = 0;

    virtual uint32 id() const;
    virtual uint32 num_levels() const;

    virtual uint32 tile_size() const;
    virtual std::string tile_filetype() const;

    virtual PixelFormatEnum pixel_format() const;
    virtual ChannelTypeEnum channel_type() const;

    // LOGGING
    // Use the audit log for events that should be recorded for posterity next to the datastore
    virtual std::ostream& audit_log() = 0;
    // a sane thing for error_log to do is log to audit_log() and vw_out(ErrorMessage, "console")
    virtual std::ostream& error_log() = 0;
};

}}

#endif
