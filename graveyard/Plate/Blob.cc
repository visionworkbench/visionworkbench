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

#define WHEREAMI if(::vw::vw_log().is_enabled(VerboseDebugMessage, "platefile.blob")) ::vw::vw_out(VerboseDebugMessage, "platefile.blob")

#if 0
FileSize       = uint64
BlobRecordSize = uint16
TileData       = uint8*
Blob = FileSize FileSize FileSize (BlobRecordSize BlobRecord TileHeader TileData)*
#endif

namespace {
  typedef vw::uint16 BlobRecordSizeType;
}

namespace vw {
namespace platefile {
  using detail::BlobRecord;

void ReadBlob::check_fail(const char* c1, const char* c2) const {
  if (m_fstream->fail()) {
    m_fstream->clear();
    vw_throw(BlobIoErr() << "BlobIoErr occured on blob " << m_blob_filename << " while " << c1 << " " << c2);
  }
}

void ReadBlob::read_at(uint64 offset64, char* dst, uint64 size, const char* context) const {
  if (m_fstream->eof())
    m_fstream->clear();

  std::streamoff offset = boost::numeric_cast<std::streamoff>(offset64);
  m_fstream->seekg(offset, std::ios_base::beg);
  check_fail(context, "while seeking");
  m_fstream->read(dst, boost::numeric_cast<std::streamsize>(size));
  check_fail(context, "while reading data");
}

uint32 tile_header_offset(const uint32& base_offset, const BlobRecord& blob_record, const BlobRecordSizeType& blob_record_size) {
  // The overall blob metadata includes the uint16 of the blob_record_size in
  // addition to the size of the blob_record itself.  The offsets stored in the
  // blob_record are relative to the END of the blob_record.  We compute this
  // offset here.
  uint64 blob_offset_metadata = sizeof(BlobRecordSizeType) + blob_record_size;
  uint64 offset64 = base_offset + blob_offset_metadata + blob_record.header_offset();
  return boost::numeric_cast<uint32>(offset64);
}

uint32 tile_data_offset(const uint32& base_offset, const BlobRecord& blob_record, const BlobRecordSizeType& blob_record_size) {
  uint64 blob_offset_metadata = sizeof(BlobRecordSizeType) + blob_record_size;
  uint64 offset64 = base_offset + blob_offset_metadata + blob_record.data_offset();
  return boost::numeric_cast<uint32>(offset64);
}

BlobRecord ReadBlob::read_blob_record(const uint32& base_offset, BlobRecordSizeType& blob_record_size) const {
  read_at(base_offset, (char*)&blob_record_size, sizeof(BlobRecordSizeType), "reading blob record size");

  boost::shared_array<uint8> blob_rec_data(new uint8[blob_record_size]);
  m_fstream->read(reinterpret_cast<char*>(blob_rec_data.get()), blob_record_size);
  check_fail("reading a blob record");

  BlobRecord blob_record;
  bool worked = blob_record.ParseFromArray(blob_rec_data.get(),  boost::numeric_cast<int>(blob_record_size));
  VW_ASSERT(worked, BlobIoErr() << "failed to parse blob record in " << m_blob_filename << " at base_offset " << base_offset);
  return blob_record;
}

TileHeader ReadBlob::read_tile_header(const uint32& base_offset, const BlobRecord& blob_record, const BlobRecordSizeType& blob_record_size) const {
  uint64 offset = tile_header_offset(base_offset, blob_record, blob_record_size);
  uint64 size   = blob_record.header_size();

  boost::scoped_array<uint8> data(new uint8[size]);
  read_at(offset, (char*)data.get(), size, "reading a tile header");

  TileHeader header;
  bool worked = header.ParseFromArray(data.get(), boost::numeric_cast<int>(size));
  VW_ASSERT(worked, BlobIoErr() << "read_tile_record() failed in " << m_blob_filename << " at offset " << m_fstream->tellg());
  return header;
}

TileData ReadBlob::read_tile_data(const uint32& base_offset, const BlobRecord& blob_record, const BlobRecordSizeType& blob_record_size) const {
  uint64 offset = tile_data_offset(base_offset, blob_record, blob_record_size);
  uint64 size   = blob_record.data_size();

  TileData data(new std::vector<uint8>(size));
  read_at(offset, (char*)(&data->operator[](0)), size, "reading tile data");
  return data;
}

bool ReadBlob::iterator::equal(iterator const& iter) const {
  return m_blob->m_blob_filename == iter.m_blob->m_blob_filename
      && m_current_base_offset == iter.m_current_base_offset;
}

void ReadBlob::iterator::increment() {
  m_current_base_offset = m_blob->next_base_offset(m_current_base_offset);
}

BlobTileRecord ReadBlob::iterator::dereference() const {
  return m_blob->read_record(m_current_base_offset);
}

ReadBlob::iterator::iterator( ReadBlob *blob, uint64 base_offset )
  : m_blob(blob), m_current_base_offset(base_offset) {}

uint64 ReadBlob::iterator::current_base_offset() const { return m_current_base_offset; }

TileHeader ReadBlob::read_header(uint64 base_offset64) {
  VW_OUT(VerboseDebugMessage, "platefile::blob") << "Entering read_header() -- " <<" base_offset: " <<  base_offset64 << "\n";
  std::streamoff base_offset = boost::numeric_cast<std::streamoff>(base_offset64);
  // Read the blob record
  BlobRecordSizeType blob_record_size;
  BlobRecord blob_record = this->read_blob_record(base_offset, blob_record_size);
  return this->read_tile_header(base_offset, blob_record, blob_record_size);
}

TileData ReadBlob::read_data(vw::uint64 base_offset) {
  uint64 offset, size;
  std::string dontcare;
  read_sendfile(base_offset, dontcare, offset, size);

  TileData data(new std::vector<uint8>(size));
  read_at(offset, (char*)(&data->operator[](0)), size, "reading tile data");
  return data;
}

BlobTileRecord ReadBlob::read_record(vw::uint64 base_offset) {
  VW_ASSERT(base_offset >= 24, LogicErr() << "No base_offset will ever be < 24. Something's wrong.");

  BlobTileRecord ret;
  uint16 blob_record_size;
  ret.rec  = read_blob_record(base_offset, blob_record_size);
  ret.hdr  = read_tile_header(base_offset, ret.rec, blob_record_size);
  ret.data = read_tile_data(base_offset, ret.rec, blob_record_size);
  return ret;
}

uint64 ReadBlob::next_base_offset(uint64 current_base_offset) {
  BlobRecordSizeType blob_record_size;
  BlobRecord blob_record = this->read_blob_record(current_base_offset, blob_record_size);

  uint64 blob_offset_metadata = sizeof(BlobRecordSizeType) + blob_record_size;
  uint64 next_offset = current_base_offset + blob_offset_metadata + blob_record.data_offset() + blob_record.data_size();

  WHEREAMI << "[next_offset: " <<  next_offset << "]\n";

  return next_offset;
}

/// Returns the data size
uint64 ReadBlob::data_size(uint64 base_offset) const {

  // Read the blob record
  BlobRecordSizeType blob_record_size;
  BlobRecord blob_record = this->read_blob_record(base_offset, blob_record_size);

  return blob_record.data_size();
}

ReadBlob::ReadBlob(const std::string& filename, bool skip_init)
  : m_blob_filename(filename)
{
  if (!skip_init)
    init();
}

ReadBlob::ReadBlob(const std::string& filename)
  : m_blob_filename(filename)
{ init(); }

void ReadBlob::init() {
  m_fstream.reset(new std::fstream(m_blob_filename.c_str(), std::ios::in | std::ios::binary));
  VW_ASSERT(m_fstream->is_open(), BlobIoErr() << "Could not open blob file " << m_blob_filename);
  m_end_of_file_ptr = read_end_of_file_ptr();
  WHEREAMI << m_blob_filename << std::endl;
}

ReadBlob::~ReadBlob() {
  WHEREAMI << m_blob_filename << "\n";
}

void ReadBlob::read_sendfile(uint64 base_offset, std::string& filename, uint64& offset, uint64& size) {
  // Read the blob record
  BlobRecordSizeType blob_record_size;
  BlobRecord blob_record = this->read_blob_record(boost::numeric_cast<std::streamoff>(base_offset), blob_record_size);

  // The overall blob metadata includes the uint16 of the
  // blob_record_size in addition to the size of the blob_record
  // itself.  The offsets stored in the blob_record are relative to
  // the END of the blob_record.  We compute this offset here.
  uint64 blob_offset_metadata = sizeof(BlobRecordSizeType) + blob_record_size;

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

  if (m_fstream->eof())
    m_fstream->clear();
  // The end of file ptr is stored at the beginning of the blob file.
  m_fstream->seekp(0, std::ios_base::beg);
  check_fail("Failed to seek to beginning of blob for EOF ptr write");
  m_fstream->write(reinterpret_cast<char*>(&data), 3*sizeof(uint64));
  check_fail("Failed to write end of file ptr");
}

uint64 ReadBlob::read_end_of_file_ptr() const {
  uint64 data[3];

  // The end of file ptr is stored at the beginning of the blob file.
  if (m_fstream->eof())
    m_fstream->clear();
  m_fstream->seekg(0, std::ios_base::beg);
  m_fstream->read(reinterpret_cast<char*>(data), 3*sizeof(uint64));
  check_fail("Failed to read end of file ptr");

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
  else if (data[0] == data[1] || data[0] == data[2])
    return data[0];
  else if (data[1] == data[2])
    return data[1];
  else {
    VW_OUT(ErrorMessage) << "end of file ptr in blobfile " << m_blob_filename << " is inconsistent. This file may be corrupt. Proceed with caution.\n";
    if (m_fstream->eof())
      m_fstream->seekg(0, std::ios_base::end);
    return m_fstream->tellg();
  }
}

// Constructor stores the blob filename for reading & writing
Blob::Blob(const std::string& filename_)
  : ReadBlob(filename_, true), m_write_count(0)
{
  m_fstream.reset(new std::fstream(m_blob_filename.c_str(), std::ios::in | std::ios::out | std::ios::binary));

  // If the file is not open, then that means that we need to create it.
  // (Note: the C++ standard does not let you create a file when you specify
  // std::ios::in., hence the gymnastics here.)
  if (!m_fstream->is_open()) {
    m_fstream->clear();
    m_fstream->open(m_blob_filename.c_str(), std::ios::out|std::ios::binary);
    VW_ASSERT(m_fstream->is_open(), BlobIoErr() << "Could not create blob file " << m_blob_filename);

    m_end_of_file_ptr = 3 * sizeof(uint64);
    flush();
    m_fstream->close();
    m_fstream->open(m_blob_filename.c_str(), std::ios::out|std::ios::in|std::ios::binary);
    VW_ASSERT(m_fstream->is_open(), BlobIoErr() << "Could not create blob file " << m_blob_filename);
  }
  m_end_of_file_ptr = read_end_of_file_ptr();
  WHEREAMI << m_blob_filename << std::endl;
}

/// Destructor: make sure that we have written the end of file ptr.
Blob::~Blob() {
  flush();
  WHEREAMI << m_blob_filename << "\n";
}

void Blob::flush() {
  this->write_end_of_file_ptr(m_end_of_file_ptr);
  m_fstream->flush();
  WHEREAMI << m_blob_filename << "\n";
}


uint64 Blob::write(TileHeader const& header, const uint8* data, uint64 data_size) {
  VW_ASSERT(m_end_of_file_ptr >= 24, LogicErr() << "What? This shouldn't happen.");

  // Store the current offset of the end of the file.  We'll
  // return that at the end of this function.
  std::streamoff base_offset = boost::numeric_cast<std::streamoff>(m_end_of_file_ptr);
  if (m_fstream->eof())
    m_fstream->clear();
  m_fstream->seekp(base_offset, std::ios_base::beg);
  check_fail("seeking to write a tile header");

  // Create the blob record and write it to the blob file.
  BlobRecord blob_record;
  blob_record.set_header_offset(0);
  blob_record.set_header_size(header.ByteSize());
  blob_record.set_data_offset(header.ByteSize());
  blob_record.set_data_size(data_size);

  // Write the actual blob record size first.  This will help us
  // read and deserialize this protobuffer later on.
  BlobRecordSizeType blob_record_size = boost::numeric_cast<BlobRecordSizeType>(blob_record.ByteSize());
  m_fstream->write(reinterpret_cast<char*>(&blob_record_size), sizeof(BlobRecordSizeType));
  check_fail("writing a blob record size");
  blob_record.SerializeToOstream(m_fstream.get());
  check_fail("writing a blob record");

  // Serialize the header.
  header.SerializeToOstream(m_fstream.get());
  check_fail("writing a tile header");

  // And write the data.
  m_fstream->write(reinterpret_cast<const char*>(data), boost::numeric_cast<size_t>(data_size));
  check_fail("writing tile data");

  // Write the data at the end of the file and return the offset
  // of the beginning of this data file.
  VW_OUT(VerboseDebugMessage, "platefile::blob") << "Blob::write() -- wrote " << data_size << " bytes to " << m_blob_filename << "\n";

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

const std::string& ReadBlob::filename() const {
  return m_blob_filename;
}

}} // namespace platefile
