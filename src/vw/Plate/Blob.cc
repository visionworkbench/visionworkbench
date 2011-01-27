// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Plate/Blob.h>
#include <vw/Plate/Exception.h>

// Vision Workbench
#include <vw/Core/Exception.h>
#include <vw/Core/FundamentalTypes.h>
#include <vw/Core/Log.h>
#include <vw/Core/Debugging.h>

#include <fstream>
#include <string>
#include <boost/shared_array.hpp>
#include <boost/scoped_array.hpp>

#define WHEREAMI (vw::vw_out(VerboseDebugMessage, "platefile.blob") << VW_CURRENT_FUNCTION << ": ")

namespace vw {
namespace platefile {

BlobRecord Blob::read_blob_record(uint16 &blob_record_size) const {

  WHEREAMI << "[Filename: " << m_blob_filename
           << " Offset: " << m_fstream->tellg() << "]\n";

  // Read the blob record
  m_fstream->read(reinterpret_cast<char*>(&blob_record_size), sizeof(blob_record_size));
  WHEREAMI << "[blob_record_size: " << blob_record_size << "]\n";

  boost::shared_array<uint8> blob_rec_data(new uint8[blob_record_size]);
  m_fstream->read(reinterpret_cast<char*>(blob_rec_data.get()), blob_record_size);
  WHEREAMI << "read complete.\n";

  BlobRecord blob_record;
  bool worked = blob_record.ParseFromArray(blob_rec_data.get(),  blob_record_size);
  if (!worked)
    vw_throw(BlobIoErr() << "read_blob_record() failed in " << m_blob_filename
                         << " at offset " << m_fstream->tellg() << "\n");
  return blob_record;
}

bool Blob::iterator::equal (iterator const& iter) const {
  return m_blob->m_blob_filename == iter.m_blob->m_blob_filename
      && m_current_base_offset == iter.m_current_base_offset;
}

void Blob::iterator::increment() {
  m_current_base_offset = m_blob->next_base_offset(m_current_base_offset);
}

TileHeader Blob::iterator::dereference() const {
  return m_blob->read_header(m_current_base_offset);
}

Blob::iterator::iterator( Blob *blob, uint64 base_offset )
  : m_blob(blob), m_current_base_offset(base_offset) {}

uint64 Blob::iterator::current_base_offset() const { return m_current_base_offset; }

TileHeader Blob::read_header(uint64 base_offset64) {

  vw_out(VerboseDebugMessage, "platefile::blob")
    << "Entering read_header() -- " <<" base_offset: " <<  base_offset64 << "\n";

  std::streamoff base_offset = boost::numeric_cast<std::streamoff>(base_offset64);

  // Seek to the requested offset and read the header and data offset
  m_fstream->seekg(base_offset, std::ios_base::beg);

  // Read the blob record
  uint16 blob_record_size;
  BlobRecord blob_record = this->read_blob_record(blob_record_size);

  // The overall blob metadata includes the uint16 of the
  // blob_record_size in addition to the size of the blob_record
  // itself.  The offsets stored in the blob_record are relative to
  // the END of the blob_record.  We compute this offset here.
  // This type-size juggling is to make sure we end up with sane behavior
  // on both 32-bit and 64-bit
  uint32 blob_offset_metadata = boost::numeric_cast<uint32>(sizeof(blob_record_size) + blob_record_size);
  size_t size = boost::numeric_cast<size_t>(blob_record.header_size());
  uint64 offset64 = base_offset + blob_offset_metadata + blob_record.header_offset();

  std::streamoff offset = boost::numeric_cast<std::streamoff>(offset64);

  // Allocate an array of the appropriate size to read the data.
  boost::shared_array<uint8> data(new uint8[size]);

  vw_out(VerboseDebugMessage, "platefile::blob")
    << "\tread_header() -- data offset: " << offset << " size: " << size << "\n";

  m_fstream->seekg(offset, std::ios_base::beg);
  m_fstream->read(reinterpret_cast<char*>(data.get()), size);

  // Throw an exception if the read operation failed (after clearing the error bit)
  if (m_fstream->fail()) {
    m_fstream->clear();
    vw_throw(IOErr() << "Blob::read() -- an error occurred while reading "
             << "data from the blob file.\n");
  }

  // Deserialize the header
  TileHeader header;
  bool worked = header.ParseFromArray(data.get(), boost::numeric_cast<int>(size));
  if (!worked)
    vw_throw(IOErr() << "Blob::read() -- an error occurred while deserializing the header "
             << "from the blob file.\n");

  vw_out(VerboseDebugMessage, "platefile::blob")
    << "\tread_header() -- read " << size << " bytes at " << offset << " from " << m_blob_filename << "\n";

  return header;
}

TileData Blob::read_data(vw::uint64 base_offset) {
  uint64 offset, size;
  std::string dontcare;
  read_sendfile(base_offset, dontcare, offset, size);

  TileData data(new std::vector<uint8>(size));

  m_fstream->seekg(offset, std::ios_base::beg);
  m_fstream->read(reinterpret_cast<char*>(&data->operator[](0)), size);

  // Throw an exception if the read operation failed (after clearing the error bit)
  if (m_fstream->fail()) {
    m_fstream->clear();
    vw_throw(IOErr() << VW_CURRENT_FUNCTION << ": failed to read from blob.");
  }

  return data;
}

boost::shared_array<uint8> Blob::read_data(uint64 base_offset, uint64& data_size) {

  uint64 offset;
  std::string dontcare;

  read_sendfile(base_offset, dontcare, offset, data_size);

  // Allocate an array of the appropriate size to read the data.
  boost::shared_array<uint8> data(new uint8[data_size]);

  m_fstream->seekg(offset, std::ios_base::beg);
  m_fstream->read(reinterpret_cast<char*>(data.get()), data_size);

  // Throw an exception if the read operation failed (after clearing the error bit)
  if (m_fstream->fail()) {
    m_fstream->clear();
    vw_throw(IOErr() << "Blob::read() -- an error occurred while reading "
             << "data from the blob file.\n");
  }

  WHEREAMI << "read " << data_size << " bytes at " << offset
           << " from " << m_blob_filename << "\n";
  return data;
}

uint64 Blob::next_base_offset(uint64 current_base_offset) {

  WHEREAMI << "[current_base_offset: " <<  current_base_offset << "]\n";

  // Seek to the requested offset and read the header and data offset
  m_fstream->seekg(current_base_offset, std::ios_base::beg);

  // Read the blob record
  uint16 blob_record_size;
  BlobRecord blob_record = this->read_blob_record(blob_record_size);

  uint64 blob_offset_metadata = sizeof(blob_record_size) + blob_record_size;
  uint64 next_offset = current_base_offset + blob_offset_metadata + blob_record.data_offset() + blob_record.data_size();

  WHEREAMI << "[next_offset: " <<  next_offset << "]\n";

  return next_offset;
}

/// Returns the data size
uint64 Blob::data_size(uint64 base_offset) const {

  WHEREAMI << "[base_offset: " <<  base_offset << "]\n";

  // Seek to the requested offset and read the header and data offset
  m_fstream->seekg(base_offset, std::ios_base::beg);

  // Read the blob record
  uint16 blob_record_size;
  BlobRecord blob_record = this->read_blob_record(blob_record_size);

  WHEREAMI << "[result size: " <<  blob_record.data_size() << "]\n";

  return blob_record.data_size();
}


// Constructor stores the blob filename for reading & writing
Blob::Blob(std::string filename, bool readonly)
  : m_blob_filename(filename), m_write_count(0)
{

  if (readonly) {
    m_fstream.reset(new std::fstream(m_blob_filename.c_str(),
                                     std::ios::in | std::ios::binary));
    if (!m_fstream->is_open())
        vw_throw(BlobIoErr() << "Could not open blob file \"" << m_blob_filename << "\".");

    WHEREAMI << filename << " (READONLY)\n";
  } else {
    m_fstream.reset(new std::fstream(m_blob_filename.c_str(),
                                     std::ios::in | std::ios::out | std::ios::binary));
    // If the file is not open, then that means that we need to create
    // it.  (Note: the C++ standard does not let you create a file
    // when you specify std::ios::in., hence the gymnastics here.)
    if (!m_fstream->is_open()) {
      m_fstream->clear();                                   // Clear error status
      m_fstream->open(filename.c_str(),                     // Create new file
                      std::ios::out|std::ios::binary);
      if (!m_fstream->is_open())                            // Check for errors
        vw_throw(BlobIoErr() << "Could not create blob file \"" << m_blob_filename << "\".");
      m_end_of_file_ptr = 3 * sizeof(uint64);
      this->write_end_of_file_ptr(m_end_of_file_ptr);       // Initialize EOF pointer
      m_fstream->close();                                   // Close output-only file
      m_fstream->open(filename.c_str(),                     // Reopen as read/write
                      std::ios::out|std::ios::in|std::ios::binary);
    }

    // Set up the fstream so that it throws an exception.
    m_fstream->exceptions ( std::fstream::eofbit | std::fstream::failbit | std::fstream::badbit );

    WHEREAMI << filename << " (READ/WRITE)\n";
  }

  if (!m_fstream->is_open())
    vw_throw(BlobIoErr() << "Could not open blob file \"" << m_blob_filename << "\".");

  // Set the cached copy of the end_of_file_ptr.
  m_end_of_file_ptr = read_end_of_file_ptr();
}

/// Destructor: make sure that we have written the end of file ptr.
Blob::~Blob() {
  this->write_end_of_file_ptr(m_end_of_file_ptr);
  WHEREAMI << m_blob_filename << "\n";
}

void Blob::read_sendfile(uint64 base_offset, std::string& filename, uint64& offset, uint64& size) {
  // Seek to the requested offset and read the header and data offset
  m_fstream->seekg(base_offset, std::ios_base::beg);

  // Read the blob record
  uint16 blob_record_size;
  BlobRecord blob_record = this->read_blob_record(blob_record_size);

  // The overall blob metadata includes the uint16 of the
  // blob_record_size in addition to the size of the blob_record
  // itself.  The offsets stored in the blob_record are relative to
  // the END of the blob_record.  We compute this offset here.
  uint64 blob_offset_metadata = sizeof(blob_record_size) + blob_record_size;

  size     = blob_record.data_size();
  offset   = base_offset + blob_offset_metadata + blob_record.data_offset();
  filename = m_blob_filename;
}

void Blob::write_end_of_file_ptr(uint64 ptr) {

  // We write the end of file pointer three times, because that
  // pretty much gurantees that at least two versions of the
  // pointer will agree if the program terminates for some reason
  // during this write opreration.  (Lazy man's checksum....)
  uint64 data[3];
  data[0] = ptr;
  data[1] = ptr;
  data[2] = ptr;

  // The end of file ptr is stored at the beginning of the blob
  // file.
  m_fstream->seekg(0, std::ios_base::beg);
  m_fstream->write(reinterpret_cast<char*>(&data), 3*sizeof(ptr));
}

uint64 Blob::read_end_of_file_ptr() const {
  uint64 data[3];

  // The end of file ptr is stored at the beginning of the blob
  // file.
  m_fstream->seekg(0, std::ios_base::beg);
  m_fstream->read(reinterpret_cast<char*>(data), 3*sizeof(uint64));

  // Make sure the read ptr is valid by comparing the three
  // entries.
  //
  // If all three agree, then return that value.  If
  // only two agree, then return that entry that two agree on.
  //
  // If none agree, return the end of file but print an error,
  // because this blob file might be corrupt.
  if ((data[0] == data[1]) && (data[0] == data[2]))
    return data[0];
  else if (data[0] == data[1])
    return data[0];
  else if (data[1] == data[2])
    return data[1];
  else if (data[0] == data[2])
    return data[0];
  else {
    vw_out(ErrorMessage) << "\nWARNING: end of file ptr in blobfile " << m_blob_filename
                         << " is inconsistent.  This file may be corrupt.  Proceed with caution.\n";
    m_fstream->seekg(0, std::ios_base::end);
    return m_fstream->tellg();
  }
}

uint64 Blob::write(TileHeader const& header, const uint8* data, uint64 data_size) {

  // Store the current offset of the end of the file.  We'll
  // return that at the end of this function.
  std::streamoff base_offset = boost::numeric_cast<std::streamoff>(m_end_of_file_ptr);
  m_fstream->seekp(base_offset, std::ios_base::beg);

  // Create the blob record and write it to the blob file.
  BlobRecord blob_record;
  blob_record.set_header_offset(0);
  blob_record.set_header_size(header.ByteSize());
  blob_record.set_data_offset(header.ByteSize());
  blob_record.set_data_size(data_size);

  // Write the actual blob record size first.  This will help us
  // read and deserialize this protobuffer later on.
  uint16 blob_record_size = boost::numeric_cast<uint16>(blob_record.ByteSize());
  m_fstream->write(reinterpret_cast<char*>(&blob_record_size), sizeof(blob_record_size));
  blob_record.SerializeToOstream(m_fstream.get());

  // Serialize the header.
  header.SerializeToOstream(m_fstream.get());

  // And write the data.
  m_fstream->write(reinterpret_cast<const char*>(data), boost::numeric_cast<size_t>(data_size));

  // Write the data at the end of the file and return the offset
  // of the beginning of this data file.
  vw_out(VerboseDebugMessage, "platefile::blob") << "Blob::write() -- writing "
                                                     << data_size
                                                     << " bytes to "
                                                     << m_blob_filename << "\n";

  // Update the in-memory copy of the end-of-file pointer
  m_end_of_file_ptr = m_fstream->tellg();

  // The write_count is used to keep track of when we last wrote
  // the end_of_file_ptr to disk.  We don't want to write this too
  // often since this will slow down IO, so we only write it every
  // 10 writes (or when the blob is deconstructed...).
  ++m_write_count;
  if (m_write_count % 10 == 0) {
    this->write_end_of_file_ptr(m_end_of_file_ptr);
  }

  // Return the base_offset
  return base_offset;
}

/// Read data out of the blob and save it as its own file on disk.
void Blob::read_to_file(std::string dest_file, uint64 offset) {
  TileData data = this->read_data(offset);

  // Open the dest_file and write to it.
  std::ofstream ostr(dest_file.c_str(), std::ios::binary);

  if (!ostr.is_open())
    vw_throw(IOErr() << VW_CURRENT_FUNCTION << ": could not open dst file for writing");

  ostr.write(reinterpret_cast<char*>(&data->operator[](0)), data->size());
  ostr.close();
}

void Blob::write_from_file(std::string source_file, TileHeader const& header, uint64& base_offset) {
  // Open the source_file and read data from it.
  std::ifstream istr(source_file.c_str(), std::ios::binary);

  if (!istr.is_open())
    vw_throw(IOErr() << "Blob::write_from_file() -- could not open source file for reading.");

  // Seek to the end and allocate the proper number of bytes of
  // memory, and then seek back to the beginning.
  istr.seekg(0, std::ios_base::end);
  std::streamoff loc = istr.tellg();
  VW_ASSERT(loc > 0, IOErr() << "Failed to identify file length");
  istr.seekg(0, std::ios_base::beg);

  size_t data_size = static_cast<size_t>(loc);

  // Read the data into a temporary memory buffer.
  boost::scoped_array<uint8> data(new uint8[data_size]);
  istr.read(reinterpret_cast<char*>(data.get()), data_size);
  istr.close();

  base_offset = this->write(header, data.get(), data_size);
}

}} // namespace platefile
