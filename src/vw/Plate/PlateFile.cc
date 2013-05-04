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
  TileHeader just_header(const Tile& t) {
    return t.hdr;
  }
}

ReadOnlyPlateFile::ReadOnlyPlateFile(const Url& url)
  : m_data(Datastore::open(url))
{
  VW_OUT(DebugMessage, "platefile") << "Re-opened plate file: \"" << url.string() << "\"\n";
}

ReadOnlyPlateFile::ReadOnlyPlateFile(const Url& url, std::string type, std::string description, uint32 tile_size, std::string tile_filetype,
                     PixelFormatEnum pixel_format, ChannelTypeEnum channel_type)
  : m_data(Datastore::open(url, make_hdr(type, description, tile_size, tile_filetype, pixel_format, channel_type)))
{
  VW_OUT(DebugMessage, "platefile") << "Constructed new platefile: " << url.string() << "\n";
}

PlateFile::PlateFile(const Url& url)
  : ReadOnlyPlateFile(url) {}

PlateFile::PlateFile(const Url& url, std::string type, std::string description, uint32 tile_size, std::string tile_filetype,
                     PixelFormatEnum pixel_format, ChannelTypeEnum channel_type)
  : ReadOnlyPlateFile(url, type, description, tile_size, tile_filetype, pixel_format, channel_type) {}

IndexHeader ReadOnlyPlateFile::index_header() const { return m_data->index_header(); }

std::string ReadOnlyPlateFile::default_file_type() const { return m_data->tile_filetype(); }

uint32 ReadOnlyPlateFile::default_tile_size() const { return m_data->tile_size(); }

PixelFormatEnum ReadOnlyPlateFile::pixel_format() const { return m_data->pixel_format(); }

ChannelTypeEnum ReadOnlyPlateFile::channel_type() const { return m_data->channel_type(); }

uint32 ReadOnlyPlateFile::num_levels() const { return m_data->num_levels(); }

void PlateFile::sync() const { m_data->flush(); }

void PlateFile::log(std::string message) { m_data->audit_log()() << message; }

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
ReadOnlyPlateFile::read(int col, int row, int level, TransactionOrNeg transaction_id, bool exact_transaction_match) const {
  TransactionRange range(exact_transaction_match ? transaction_id : 0, transaction_id);
  Datastore::TileSearch hits;
  m_data->get(hits, level, row, col, range, 1);
  if (hits.size() == 0)
    vw_throw(TileNotFoundErr() << "No tiles found.");

  return std::make_pair(hits.begin()->hdr, hits.begin()->data);
}

Datastore::TileSearch&
ReadOnlyPlateFile::batch_read(Datastore::TileSearch& hdrs) const {
  return m_data->populate(hdrs);
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

boost::tuple<std::string, TileHeader, TileData>
ReadOnlyPlateFile::img_file_name(std::string const& base_name, int col, int row, int level,
                                 TransactionOrNeg transaction_id,
                                 bool exact_transaction_match) const
{
  // Get the image file name and other attributes. 

  std::pair<TileHeader, TileData> ret = this->read(col, row, level, transaction_id, exact_transaction_match);
  std::string filename = base_name + "." + ret.first.filetype();
  return boost::tuple<std::string, TileHeader, TileData>(filename, ret.first, ret.second);
}

std::pair<std::string, TileHeader>
ReadOnlyPlateFile::read_to_file(std::string const& base_name, int col, int row, int level,
                        TransactionOrNeg transaction_id, bool exact_transaction_match) const
{
  // Dump the current image to a file. Return the image file name and extension.
  
  boost::tuple<std::string, TileHeader, TileData> ret = this->img_file_name(base_name, col, row, level,
                                                                            transaction_id, exact_transaction_match);
  std::string filename = ret.get<0>();
  TileHeader  th       = ret.get<1>(); 
  TileData    td       = ret.get<2>(); 
  dump_to_file(filename, &(td->operator[](0)), td->size());
  return std::make_pair(filename, th);
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
ReadOnlyPlateFile::search_by_region(int level, vw::BBox2i const& region, const TransactionRange& range) const {
  Datastore::TileSearch r;
  m_data->head(r, level, region, range, 0);
  std::list<TileHeader> tiles(r.size());
  std::transform(r.begin(), r.end(), tiles.begin(), just_header);
  return tiles;
}

std::list<TileHeader>
ReadOnlyPlateFile::search_by_location(int col, int row, int level, const TransactionRange& range) {
  Datastore::TileSearch r;
  m_data->head(r, level, row, col, range, 0);
  std::list<TileHeader> tiles(r.size());
  std::transform(r.begin(), r.end(), tiles.begin(), just_header);
  return tiles;
}

std::ostream& PlateFile::audit_log() {
  return m_data->audit_log()();
}
std::ostream& PlateFile::error_log() {
  return m_data->error_log()();
}
