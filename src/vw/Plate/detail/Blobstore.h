// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_PLATE_DETAIL_BLOBSTORE_H__
#define __VW_PLATE_DETAIL_BLOBSTORE_H__

#include <vw/Plate/Datastore.h>
#include <vw/Core/Log.h>

namespace vw { namespace platefile {
  class Blob;
namespace detail {
  class Index;

class Blobstore : public Datastore {
  private:
    boost::shared_ptr<Index> m_index;
    uint32 m_read_blob_id;
    boost::shared_ptr<Blob> m_read_blob;

    typedef MultiOutputStream<char> log_t;
    log_t m_error_log;

    void init();
    boost::shared_ptr<Blob> open_blob(uint32 blob_id, bool readonly);
  public:
    Blobstore(const Url& u);
    Blobstore(const Url& u, const IndexHeader& d);

    virtual Transaction transaction_begin(const std::string& description, TransactionOrNeg transaction_id_override = -1);
    virtual void transaction_end(Transaction transaction_id, bool update_read_cursor);

    virtual meta_range     head(uint32 level, uint64 row, uint64 col, TransactionRange range, uint32 limit = 1);
    virtual meta_range     head(uint32 level,   const BBox2u& region, TransactionRange range, uint32 limit = 1);
    virtual tile_range      get(uint32 level, uint64 row, uint64 col, TransactionRange range, uint32 limit = 1);
    virtual tile_range      get(uint32 level,   const BBox2u& region, TransactionRange range, uint32 limit = 1);
    virtual tile_range populate(const TileHeader* hdrs, size_t len);

    //virtual Url map_to_url(uint32 level, uint32 row, uint32 col, Transaction id, std::string filetype);
    //virtual Url map_to_url(const TileHeader& t);

    virtual WriteState* write_request(const Transaction& id);
    virtual void write_update(WriteState& state, uint32 level, uint32 row, uint32 col, const std::string& filetype, const uint8* data, uint64 size);
    virtual void write_complete(WriteState& id);
    virtual void flush();

    virtual IndexHeader index_header() const;

    // LOGGING
    virtual std::ostream& audit_log();
    virtual std::ostream& error_log();
};


}}} // vw::platefile::detail

#endif
