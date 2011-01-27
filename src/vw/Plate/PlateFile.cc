// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Plate/PlateFile.h>
#include <vw/Plate/HTTPUtils.h>
#include <vw/Core/Settings.h>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>

namespace fs = boost::filesystem;

using namespace vw::platefile;
using namespace vw;

PlateFile::PlateFile(const Url& url) {
  m_index = Index::construct_open(url);
  vw_out(DebugMessage, "platefile") << "Re-opened plate file: \"" << url.string() << "\"\n";
}

PlateFile::PlateFile(const Url& url, std::string type, std::string description,
                     int tile_size, std::string tile_filetype,
                     PixelFormatEnum pixel_format, ChannelTypeEnum channel_type) {

  IndexHeader hdr;
  hdr.set_type(type);
  hdr.set_description(description);
  hdr.set_tile_size(tile_size);
  hdr.set_tile_filetype(tile_filetype);
  hdr.set_pixel_format(pixel_format);
  hdr.set_channel_type(channel_type);

  vw_out(DebugMessage, "platefile") << "Constructing new platefile: " << url.string() << "\n";
  m_index = Index::construct_create(url, hdr);
}

std::pair<TileHeader, TileData>
PlateFile::read(int col, int row, int level, TransactionOrNeg transaction_id, bool exact_transaction_match) const {
  IndexRecord record = m_index->read_request(col, row, level, transaction_id, exact_transaction_match);

  boost::shared_ptr<Blob> read_blob;
  if (m_write_blob && record.blob_id() == m_write_blob_id) {
    read_blob = m_write_blob;
  } else {
    std::ostringstream blob_filename;
    blob_filename << this->name() << "/plate_" << record.blob_id() << ".blob";
    read_blob.reset(new Blob(blob_filename.str(), true));
  }

  std::pair<TileHeader, TileData> tile;
  tile.first  = read_blob->read_header(record.blob_offset());
  tile.second = read_blob->read_data(record.blob_offset());
  return tile;
}

/// Read the tile header. You supply a base name (without the
/// file's image extension).  The image extension will be appended
/// automatically for you based on the filetype in the TileHeader.
std::pair<std::string, TileHeader>
PlateFile::read_to_file(std::string const& base_name, int col, int row, int level,
                        TransactionOrNeg transaction_id, bool exact_transaction_match) const
{
  // 1. Call index read_request(col,row,level).  Returns IndexRecord.
  IndexRecord record = m_index->read_request(col, row, level, transaction_id, exact_transaction_match);

  // 2. Open the blob file and read the header
  boost::shared_ptr<Blob> read_blob;
  if (m_write_blob && record.blob_id() == m_write_blob_id) {
    read_blob = m_write_blob;
  } else {
    std::ostringstream blob_filename;
    blob_filename << this->name() << "/plate_" << record.blob_id() << ".blob";
    read_blob.reset(new Blob(blob_filename.str(), true));
  }

  // 3. Choose a temporary filename and call BlobIO
  // read_as_file(filename, offset, size) [ offset, size from
  // IndexRecord ]
  std::string filename = base_name + "." + record.filetype();
  read_blob->read_to_file(filename, record.blob_offset());

  TileHeader hdr = read_blob->read_header(record.blob_offset());

  // 4. Return the name of the file
  return std::make_pair(filename, hdr);
}

/// Writing, pt. 1: Locks a blob and returns the blob id that can
/// be used to write tiles.
void vw::platefile::PlateFile::write_request() {

  // Request a blob lock from the index
  uint64 last_size;
  m_write_blob_id = m_index->write_request(last_size);

  // Compute blob filename for writing.
  std::ostringstream blob_filename;
  blob_filename << this->name() << "/plate_" << m_write_blob_id << ".blob";

  // Open the blob for writing.
  m_write_blob.reset( new Blob(blob_filename.str()) );

  if (last_size != 0 && last_size != m_write_blob->size()) {
    std::ostringstream ostr;
    ostr << "WARNING: last close size did not match current size when opening "
         << blob_filename.str()
         << "  ( " << last_size << " != " << m_write_blob->size() << " )\n";
    m_index->log(ostr.str());
  }

  // For debugging:
  std::ostringstream ostr;
  ostr << "Opened blob " << m_write_blob_id << " ( size = " << m_write_blob->size() << " )\n";
  m_index->log(ostr.str());
}

/// Writing, pt. 3: Signal the completion of the write operation.
void vw::platefile::PlateFile::write_complete() {

  // Fetch the size from the blob.
  uint64 new_blob_size = m_write_blob->size();

  // Close the blob for writing.
  m_write_blob.reset();

  // For debugging:
  std::ostringstream ostr;
  ostr << "Closed blob " << m_write_blob_id << " ( size = " << new_blob_size << " )\n";
  m_index->log(ostr.str());

  // Release the blob lock.
  return m_index->write_complete(m_write_blob_id, new_blob_size);
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
vw::platefile::IndexRecord vw::platefile::PlateFile::read_record(int col, int row, int level,
                                                                 TransactionOrNeg transaction_id,
                                                                 bool exact_transaction_match) {
  return m_index->read_request(col, row, level, transaction_id, exact_transaction_match);
}

void PlateFile::write_update(const uint8* data, uint64 data_size, int col, int row,
    int level, Transaction transaction_id, const std::string& type_) {

  std::string type = type_;
  if (type.empty())
    type = this->default_file_type();

  if (type == "auto")
    vw_throw(NoImplErr() << "write_update() does not support filetype 'auto'");
  if (!m_write_blob)
    vw_throw(BlobIoErr() << "write_update(): No blob file open. Are you sure you ran write_request()?");

  TileHeader write_header;
  write_header.set_col(col);
  write_header.set_row(row);
  write_header.set_level(level);
  write_header.set_transaction_id(transaction_id);
  write_header.set_filetype(type);

  // 1. Write the data into the blob
  uint64 blob_offset = m_write_blob->write(write_header, data, data_size);

  // 2. Update the index
  IndexRecord write_record;
  write_record.set_blob_id(m_write_blob_id);
  write_record.set_blob_offset(blob_offset);
  write_record.set_filetype(write_header.filetype());

  m_index->write_update(write_header, write_record);
}
