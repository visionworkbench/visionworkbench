
// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Plate/Blob.h>

// Vision Workbench
#include <vw/Core/Exception.h>
#include <vw/Core/FundamentalTypes.h>
#include <vw/Core/Log.h>

#include <fstream>
#include <string>
#include <boost/shared_array.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>
namespace fs = boost::filesystem;

using namespace vw;
    
// -------------------------------------------------------------------
//                                 BLOB
// -------------------------------------------------------------------

// Constructor stores the blob filename for reading & writing
vw::platefile::Blob::Blob(std::string filename, bool readonly) : m_blob_filename(filename) {

  // This is an ugly hack, but it appears that with some versions of
  // glibc, it is not possible to CREATE a file at the same time as
  // you open it for reading and writing.  Instead, we open it for
  // writing here, creating it if necessary, then close and reopen it.
  // Blech.
  if (!readonly) {
    std::fstream temporary_stream(m_blob_filename.c_str(),
                                  std::ios::out | std::ios::app);  
    if (!temporary_stream.is_open()) {
      std::cout << "Blob file could not be created!\n";
      vw_throw(IOErr() << "Blob::Blob(): could not create blob file \"" << m_blob_filename << "\".");
    }
    temporary_stream.close();
  }

  // Now we open the file for real.
  if (readonly) {
    m_fstream.reset(new std::fstream(m_blob_filename.c_str(), 
                                     std::ios::in | std::ios::binary));
  } else {
    m_fstream.reset(new std::fstream(m_blob_filename.c_str(), 
                                     std::ios::in | std::ios::out | std::ios::binary));
  }

  if (!m_fstream->is_open()) {
    std::cout << "Blob file is not open!\n";
    vw_throw(IOErr() << "Blob::Blob(): could not open blob file \"" << m_blob_filename << "\".");
  }
}

vw::platefile::Blob::~Blob() {}

void vw::platefile::Blob::read_sendfile(vw::uint64 base_offset, std::string& filename, vw::uint64& offset, vw::uint64& size) {
  // Seek to the requested offset and read the header and data offset
  m_fstream->seekg(base_offset, std::ios_base::beg);

  // Read the blob record
  uint16 blob_record_size;
  BlobRecord blob_record = this->read_blob_record(blob_record_size);

  // The overall blob metadata includes the uint16 of the
  // blob_record_size in addition to the size of the blob_record
  // itself.  The offsets stored in the blob_record are relative to
  // the END of the blob_record.  We compute this offset here.
  uint32 blob_offset_metadata = sizeof(blob_record_size) + blob_record_size;

  size     = blob_record.data_size();
  offset   = base_offset + blob_offset_metadata + blob_record.data_offset();
  filename = m_blob_filename;

}

/// read_data()
boost::shared_array<uint8> vw::platefile::Blob::read_data(vw::uint64 base_offset) {

  vw::uint64 offset, size;
  std::string dontcare;

  read_sendfile(base_offset, dontcare, offset, size);

  // Allocate an array of the appropriate size to read the data.
  boost::shared_array<uint8> data(new uint8[size]);

  m_fstream->seekg(offset, std::ios_base::beg);
  m_fstream->read((char*)(data.get()), size);

  // Throw an exception if the read operation failed (after clearing the error bit)
  if (m_fstream->fail()) {
    m_fstream->clear();
    vw_throw(IOErr() << "Blob::read() -- an error occurred while reading " 
             << "data from the blob file.\n");
  }

  vw::vw_out(vw::VerboseDebugMessage, "plate::blob") << "Blob::read() -- read " 
                                                     << size << " bytes at " << offset 
                                                     << " from " << m_blob_filename << "\n";
  return data;
}

/// Read data out of the blob and save it as its own file on disk.
void vw::platefile::Blob::read_to_file(std::string dest_file, int64 offset) {
  boost::shared_array<uint8> data = this->read_data(offset);
  uint32 size = this->data_size(offset);

  // Open the dest_file and write to it.
  std::ofstream ostr(dest_file.c_str(), std::ios::binary);

  if (!ostr.is_open())
    vw_throw(IOErr() << "Blob::read_as_file() -- could not open destination " 
             << "file for writing..");

  ostr.write((char*)(data.get()), size);
  ostr.close();
}

