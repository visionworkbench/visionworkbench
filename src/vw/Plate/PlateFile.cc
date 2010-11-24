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

// -------------------------------------------------------------------------
//                            TEMPORARY TILE FILE
// -------------------------------------------------------------------------

std::string vw::platefile::TemporaryTileFile::unique_tempfile_name(std::string file_extension) {
  std::string base_name = vw_settings().tmp_directory() + "/vw_plate_tile_XXXXXXX";

  boost::scoped_array<char> c_str(new char[base_name.size()+1]);
  strncpy(c_str.get(), base_name.c_str(), base_name.size()+1);
  char* dummy = mktemp(c_str.get());
  std::string ret(dummy);
  return ret + "." + file_extension;
}

vw::platefile::TemporaryTileFile::TemporaryTileFile(std::string filename) :
  m_filename(filename) {
  vw_out(DebugMessage, "plate::tempfile") << "Assumed control of temporary file: "
                                          << m_filename << "\n";
}

vw::platefile::TemporaryTileFile::~TemporaryTileFile() {
  int result = unlink(m_filename.c_str());
  if (result)
    vw_out(ErrorMessage, "plate::tempfile")
      << "WARNING: unlink() failed in ~TemporaryTileFile() for filename \""
      << m_filename << "\"\n";
  vw_out(DebugMessage, "plate::tempfile") << "Destroyed temporary file: "
                                          << m_filename << "\n";
}

/// Opens the temporary file and determines its size in bytes.
vw::int64 vw::platefile::TemporaryTileFile::file_size() const {
  std::ifstream istr(m_filename.c_str(), std::ios::binary);

  if (!istr.is_open())
    vw_throw(IOErr() << "TempPlateFile::file_size() -- could not open \""
             << m_filename << "\".");
  istr.seekg(0, std::ios_base::end);
  return istr.tellg();
}



// -------------------------------------------------------------------------
//                               PLATE FILE
// -------------------------------------------------------------------------

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

  TileHeader hdr = read_blob->read_header<TileHeader>(record.blob_offset());

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

void PlateFile::write_update(const boost::shared_array<uint8> data, uint64 data_size,
                             int col, int row, int level, Transaction transaction_id) {

  // Quick sanity check.
  if (this->default_file_type() == "auto") {
    vw_throw(NoImplErr() << "write_update() does not support writing un-typed "
                         << "data arrays for filetype \'auto\'.\n");
  }


  if (!m_write_blob)
    vw_throw(BlobIoErr() << "Error issuing write_update(). No blob file open. "
                         << "Are you sure your ran write_request()?");

  // 0. Create a write_header
  TileHeader write_header;
  write_header.set_col(col);
  write_header.set_row(row);
  write_header.set_level(level);
  write_header.set_transaction_id(transaction_id);
  write_header.set_filetype(this->default_file_type());

  // 1. Write the data into the blob
  uint64 blob_offset = m_write_blob->write(write_header, data, data_size);

  // 2. Update the index
  IndexRecord write_record;
  write_record.set_blob_id(m_write_blob_id);
  write_record.set_blob_offset(blob_offset);
  write_record.set_filetype(write_header.filetype());

  m_index->write_update(write_header, write_record);
  }
