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

using namespace vw;

#define WHEREAMI (vw::vw_out(VerboseDebugMessage, "platefile.blob") << VW_CURRENT_FUNCTION << ": ")

// -------------------------------------------------------------------
//                                 BLOB
// -------------------------------------------------------------------

/// read_blob_record()
vw::platefile::BlobRecord vw::platefile::Blob::read_blob_record(uint16 &blob_record_size) const {

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

/// read_data()
boost::shared_array<uint8> vw::platefile::Blob::read_data(vw::uint64 base_offset, vw::uint64& data_size) {

  vw::uint64 offset;
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

vw::uint64 vw::platefile::Blob::next_base_offset(uint64 current_base_offset) {

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
vw::uint64 vw::platefile::Blob::data_size(uint64 base_offset) const {

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
vw::platefile::Blob::Blob(std::string filename, bool readonly) :
  m_blob_filename(filename), m_write_count(0) {

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
vw::platefile::Blob::~Blob() {
  this->write_end_of_file_ptr(m_end_of_file_ptr);
  WHEREAMI << m_blob_filename << "\n";
}

void vw::platefile::Blob::read_sendfile(vw::uint64 base_offset, std::string& filename,
                                        vw::uint64& offset, vw::uint64& size) {
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

void vw::platefile::Blob::write_end_of_file_ptr(uint64 ptr) {

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

uint64 vw::platefile::Blob::read_end_of_file_ptr() const {
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

/// Read data out of the blob and save it as its own file on disk.
void vw::platefile::Blob::read_to_file(std::string dest_file, uint64 offset) {
  uint64 size;
  boost::shared_array<uint8> data = this->read_data(offset, size);

  // Open the dest_file and write to it.
  std::ofstream ostr(dest_file.c_str(), std::ios::binary);

  if (!ostr.is_open())
    vw_throw(IOErr() << "Blob::read_as_file() -- could not open destination "
             << "file for writing..");

  ostr.write(reinterpret_cast<char*>(data.get()), size);
  ostr.close();
}

