// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_PLATE_DETAIL_DIRSTORE_H__
#define __VW_PLATE_DETAIL_DIRSTORE_H__

#include <vw/Plate/Datastore.h>
#include <vw/Core/Log.h>

namespace vw { namespace platefile {
namespace detail {
  class Index;

class Dirstore : public Datastore {
  private:
    void init();
    std::string m_plate_path;
    IndexHeader m_hdr;

    TileSearch& top_id_one_loc(TileSearch& tiles, uint32 level, uint32 row, uint32 col) const;
    TileSearch& one_id_one_loc(TileSearch& tiles, uint32 level, uint32 row, uint32 col, const Transaction& id) const;
    TileSearch& many_id_one_loc(TileSearch& tiles, uint32 level, uint32 row, uint32 col, const TransactionRange& r, uint32 limit) const;

    TileSearch& top_id_many_loc(TileSearch& tiles, uint32 level, const BBox2u& region);
    TileSearch& one_id_many_loc(TileSearch& tiles, uint32 level, const BBox2u& region, const Transaction& id) const;
    TileSearch& many_id_many_loc(TileSearch& tiles, uint32 level, const BBox2u& region, const TransactionRange& r, uint32 limit) const;

    std::string path_by_tid(const platefile::Transaction& id, uint32 level, uint32 row, uint32 col, const std::string& filetype) const;
    std::string path_by_loc(const platefile::Transaction& id, uint32 level, uint32 row, uint32 col, const std::string& filetype) const;
    std::string path_by_loc_no_tid(uint32 level, uint32 row, uint32 col) const;
    std::string path_by_bucket(uint32 level, uint64 bucket) const;
    std::string tid_dir(const Transaction& id) const;

    void add_one_id(Datastore::TileSearch& tiles, uint32 level, uint32 row, uint32 col, const Transaction& id) const;
    void add_top_id(Datastore::TileSearch& tiles, uint32 level, uint32 row, uint32 col) const;
    void add_ids_in_range(Datastore::TileSearch& tiles, uint32 level, uint32 row, uint32 col, const TransactionRange& r, uint32 limit) const;

    // args are row, col
    typedef boost::function<void (uint32, uint32)> add_func_t;
    void bucket_iterate(uint32 level, const BBox2u& region, add_func_t add_func) const;

    std::string get_filetype(uint32 level, uint32 row, uint32 col, const Transaction& id) const;
    std::string get_lockfile(const Transaction& id) const;

    void save_index_file() const;

  public:
    Dirstore(const Url& u);
    Dirstore(const Url& u, const IndexHeader& d);

    virtual Transaction transaction_begin(const std::string& description, TransactionOrNeg transaction_id_override = -1);
    virtual void transaction_end(Transaction transaction_id, bool update_read_cursor);

    virtual TileSearch& head(TileSearch& tiles, uint32 level, uint32 row, uint32 col, TransactionRange range, uint32 limit = 0);
    virtual TileSearch& head(TileSearch& tiles, uint32 level,   const BBox2u& region, TransactionRange range, uint32 limit = 0);
    virtual TileSearch& populate(TileSearch& hdrs);

    virtual WriteState* write_request(const Transaction& id);
    virtual void write_update(WriteState& state, uint32 level, uint32 row, uint32 col, const std::string& filetype, const uint8* data, uint64 size);
    virtual void write_complete(WriteState& id);
    virtual void flush();

    virtual IndexHeader index_header() const;
    virtual vw::uint32 num_levels() const;
    virtual vw::uint32 id() const;
    virtual vw::uint32 tile_size() const;
    virtual std::string tile_filetype() const;
    virtual vw::PixelFormatEnum pixel_format() const;
    virtual vw::ChannelTypeEnum channel_type() const;

    // LOGGING
    virtual Logger audit_log() const;
    virtual Logger error_log() const;
};


}}} // vw::platefile::detail

#endif
