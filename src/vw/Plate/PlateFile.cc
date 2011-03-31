// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Plate/PlateFile.h>
#include <vw/Plate/Datastore.h>
#include <vw/Plate/HTTPUtils.h>
#include <vw/Core/Settings.h>
#include <vw/Core/Debugging.h>
#include <boost/iostreams/tee.hpp>

using namespace vw::platefile;
using namespace vw;

namespace {
  IndexHeader make_hdr(std::string type, std::string description, uint32 tile_size, std::string tile_filetype, PixelFormatEnum pixel_format, ChannelTypeEnum channel_type) {
    IndexHeader hdr;
    hdr.set_type(type);
    hdr.set_description(description);
    hdr.set_tile_size(tile_size);
    hdr.set_tile_filetype(tile_filetype);
    hdr.set_pixel_format(pixel_format);
    hdr.set_channel_type(channel_type);
    return hdr;
  }
}

PlateFile::PlateFile(const Url& url)
  : m_data(Datastore::open(url))
{
  m_error_log.add(vw_out(ErrorMessage, "console"));
  m_error_log.add(m_data->audit_log());
  vw_out(DebugMessage, "platefile") << "Re-opened plate file: \"" << url.string() << "\"\n";
}

PlateFile::PlateFile(const Url& url, std::string type, std::string description, uint32 tile_size, std::string tile_filetype,
                     PixelFormatEnum pixel_format, ChannelTypeEnum channel_type)
  : m_data(Datastore::open(url, make_hdr(type, description, tile_size, tile_filetype, pixel_format, channel_type)))
{
  m_error_log.add(vw_out(ErrorMessage, "console"));
  m_error_log.add(m_data->audit_log());
  vw_out(DebugMessage, "platefile") << "Constructed new platefile: " << url.string() << "\n";
}

//std::string PlateFile::name() const { return m_data->name(); }

IndexHeader PlateFile::index_header() const { return m_data->index_header(); }

std::string PlateFile::default_file_type() const { return m_data->tile_filetype(); }

uint32 PlateFile::default_tile_size() const { return m_data->tile_size(); }

PixelFormatEnum PlateFile::pixel_format() const { return m_data->pixel_format(); }

ChannelTypeEnum PlateFile::channel_type() const { return m_data->channel_type(); }

uint32 PlateFile::num_levels() const { return m_data->num_levels(); }

void PlateFile::sync() const { m_data->flush(); }

void PlateFile::log(std::string message) { m_data->audit_log() << message; }

const Transaction& PlateFile::transaction_begin(const std::string& transaction_description, TransactionOrNeg transaction_id_override) {
  m_transaction.reset(new Transaction(m_data->transaction_begin(transaction_description, transaction_id_override)));
  return *m_transaction;
}

void PlateFile::transaction_resume(const Transaction& tid) {
  m_transaction.reset(new Transaction(tid));
}

void PlateFile::transaction_end(bool update_read_cursor) {
  m_data->transaction_end(*m_transaction, update_read_cursor);
  VW_ASSERT(!m_write_state, LogicErr() << "Must end write before ending transaction");
  m_transaction.reset();
}

const Transaction& PlateFile::transaction_id() const {
  VW_ASSERT(m_transaction, LogicErr() << "Must start transaction before asking for the id");
  return *m_transaction;
}

std::pair<TileHeader, TileData>
PlateFile::read(int col, int row, int level, TransactionOrNeg transaction_id, bool exact_transaction_match) const {
  TransactionRange range(exact_transaction_match ? transaction_id : 0, transaction_id);
  Datastore::tile_range hits = m_data->get(level, row, col, range);
  if (hits.size() == 0)
    vw_throw(TileNotFoundErr() << "No tiles found.");

  return std::make_pair(hits.begin()->hdr, hits.begin()->data);
}

namespace {
  void dump_to_file(const std::string& filename, const vw::uint8* data, size_t size) {
    std::ofstream f(filename.c_str(), std::ios::binary);
    VW_ASSERT(f.is_open(), IOErr() << VW_CURRENT_FUNCTION << ": could not open dst file for writing");
    f.write(reinterpret_cast<const char*>(data), size);
    VW_ASSERT(!f.fail(), IOErr() << VW_CURRENT_FUNCTION << ": could not write to dst file");
    f.close();
  }
}

std::pair<std::string, TileHeader>
PlateFile::read_to_file(std::string const& base_name, int col, int row, int level,
                        TransactionOrNeg transaction_id, bool exact_transaction_match) const
{
  std::pair<TileHeader, TileData> ret = this->read(col, row, level, transaction_id, exact_transaction_match);
  std::string filename = base_name + "." + ret.first.filetype();
  dump_to_file(filename, &(ret.second->operator[](0)), ret.second->size());
  return std::make_pair(filename, ret.first);
}

/// Writing, pt. 1: Locks a blob and returns the blob id that can
/// be used to write tiles.
void PlateFile::write_request() {
  VW_ASSERT(m_transaction, LogicErr() << "Must start transaction before trying to write");
  m_write_state.reset(m_data->write_request(*m_transaction));
}

/// Writing, pt. 3: Signal the completion of the write operation.
void PlateFile::write_complete() {
  VW_ASSERT(m_write_state, LogicErr() << "Must start a transaction before completing it");
  m_data->write_complete(*m_write_state);
  m_write_state.reset();
}

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
//IndexRecord PlateFile::read_record(int col, int row, int level,
//                                                                 TransactionOrNeg transaction_id,
//                                                                 bool exact_transaction_match) {
//  return m_index->read_request(col, row, level, transaction_id, exact_transaction_match);
//}

void PlateFile::write_update(const uint8* data, uint64 data_size, int col, int row, int level, const std::string& type_) {

  std::string type = type_;
  if (type.empty())
    type = this->default_file_type();

  if (!m_write_state)
    vw_throw(LogicErr() << "write_update(): you must first request a write");
  if (type == "auto")
    vw_throw(NoImplErr() << "write_update() does not support filetype 'auto'");

  m_data->write_update(*m_write_state, level, row, col, type, data, data_size);
}

std::list<TileHeader>
PlateFile::search_by_region(int level, vw::BBox2i const& region, const TransactionRange& range) const {
  Datastore::meta_range r = m_data->head(level, region, range, 0);
  std::list<TileHeader> tiles;
  tiles.insert(tiles.begin(), r.begin(), r.end());
  return tiles;
}

std::list<TileHeader>
PlateFile::search_by_location(int col, int row, int level, const TransactionRange& range) {
  Datastore::meta_range r = m_data->head(level, row, col, range, 0);
  std::list<TileHeader> tiles;
  tiles.insert(tiles.begin(), r.begin(), r.end());
  return tiles;
}

std::ostream& PlateFile::audit_log() {
  return m_data->audit_log();
}
std::ostream& PlateFile::error_log() {
  return m_error_log;
}
