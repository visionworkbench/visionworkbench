// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
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
vw::platefile::Blob::Blob(std::string filename) : m_blob_filename(filename) {}

vw::platefile::Blob::~Blob() {}

/// read_data()
boost::shared_array<uint8> vw::platefile::Blob::read_data(vw::uint64 base_offset) {

  
  std::ifstream ifstr(m_blob_filename.c_str(), std::ios::in | std::ios::binary);
  if (!ifstr.is_open()) 
    vw_throw(IOErr() << "Blob::read_data(): could not open blob file \"" << m_blob_filename << "\".");

  // Seek to the requested offset and read the header and data offset
  ifstr.seekg(base_offset, std::ios_base::beg);

  // Read the blob record
  uint16 blob_record_size;
  BlobRecord blob_record = this->read_blob_record(ifstr, blob_record_size);

  // The overall blob metadat includes the uint16 of the
  // blob_record_size in addition to the size of the blob_record
  // itself.  The offsets stored in the blob_record are relative to
  // the END of the blob_record.  We compute this offset here.
  uint32 blob_offset_metadata = sizeof(blob_record_size) + blob_record_size;
  int32 size = blob_record.data_size();
  uint64 offset = base_offset + blob_offset_metadata + blob_record.data_offset();
  
  // Allocate an array of the appropriate size to read the data.
  boost::shared_array<uint8> data(new uint8[size]);
  
  ifstr.seekg(offset, std::ios_base::beg);
  ifstr.read((char*)(data.get()), size);
  
  // Throw an exception if the read operation failed (after clearing the error bit)
  if (ifstr.fail()) {
    ifstr.clear();
    vw_throw(IOErr() << "Blob::read() -- an error occurred while reading " 
             << "data from the blob file.\n");
  }
  
  vw::vw_out(vw::VerboseDebugMessage, "plate::blob") << "Blob::read() -- read " 
                                                     << size << " bytes at " << offset 
                                                     << " from " << m_blob_filename << "\n";
  return data;
}

/// Read data out of the blob and save it as its own file on disk.
void vw::platefile::Blob::read_to_file(std::string dest_file, int64 offset, int64 size) {
  boost::shared_array<uint8> data = this->read_data(offset);

  // Open the dest_file and write to it.
  std::ofstream ostr(dest_file.c_str(), std::ios::binary);

  if (!ostr.is_open())
    vw_throw(IOErr() << "Blob::read_as_file() -- could not open destination " 
             << "file for writing..");

  ostr.write((char*)(data.get()), size);
  ostr.close();
}

