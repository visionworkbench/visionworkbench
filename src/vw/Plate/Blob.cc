#include <vw/Plate/Blob.h>

// Vision Workbench
#include <vw/Core/Exception.h>
#include <vw/Core/FundamentalTypes.h>
#include <vw/Core/Log.h>

#include <fstream>
#include <string>
#include <boost/shared_array.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>
namespace fs = boost::filesystem;
    
// -------------------------------------------------------------------
//                                 BLOB
// -------------------------------------------------------------------

// Constructor stores the blob filename for reading & writing
vw::platefile::Blob::Blob(std::string filename) : m_filename(filename) {
    
  m_blobfile_version = 1;
  
  if( !exists( fs::path( filename, fs::native ) ) ) {
        
    // Check to see if this blob exists.  If not, we create it and
    // write a small blob header.
    std::ofstream ostr(filename.c_str());
    ostr << m_blobfile_version << "\nPlatefile data blob.\n"
         << "Created by the NASA Vision Workbench\n\n";

  } else {

    std::ifstream istr(filename.c_str());
    istr >> m_blobfile_version;
    if (m_blobfile_version != 1) 
      vw_throw(IOErr() << "Could not open plate file blob.  Does not appear to be a valid blob file.");      

  }    
}

int vw::platefile::Blob::version() const { return m_blobfile_version; }

/// Read data out of the blob and save it as its own file on disk.
void vw::platefile::Blob::read_to_file(std::string dest_file, int64 offset, int64 size) {
  boost::shared_array<char> data = this->read<char>(offset, size);

  // Open the dest_file and write to it.
  std::ofstream ostr(dest_file.c_str(), std::ios::binary);

  if (!ostr.is_open())
    vw_throw(IOErr() << "Blob::read_as_file() -- could not open destination " 
             << "file for writing..");

  ostr.write(data.get(), size);
  ostr.close();
}

void vw::platefile::Blob::write_from_file(std::string source_file, int64& offset, int32& size) {

  // Open the source_file and read data from it.
  std::ifstream istr(source_file.c_str(), std::ios::binary);
      
  if (!istr.is_open())
    vw_throw(IOErr() << "Blob::write_from_file() -- could not open source file for reading.");
  // Seek to the end and allocate the proper number of bytes of
  // memory, and then seek back to the beginning.
  istr.seekg(0, std::ios_base::end);
  size = istr.tellg();
  istr.seekg(0, std::ios_base::beg);
      
  // Read the data into a temporary memory buffer.
  boost::shared_array<char> data(new char[size]);
  istr.read(data.get(), size);
  istr.close();

  offset = this->write<char>(data, size);
}

// -------------------------------------------------------------------
//                            BLOB_MANAGER
// -------------------------------------------------------------------

void vw::platefile::BlobManager::BlobManager::next_blob_index() {
  // Move to the next blob
  m_blob_index++;
  if (m_blob_index >= m_blob_locks.size())
    m_blob_index = 0;
}

/// Create a new blob manager.  The max_blob_size is specified in
/// units of megabytes.
vw::platefile::BlobManager::BlobManager(int64 max_blob_size, int nblobs) : 
  m_max_blob_size(max_blob_size * 1024 * 1024), m_blob_index(0) { 
  
  m_blob_locks.resize(nblobs);
  
  for (int i=0; i < m_blob_locks.size(); ++i) {
    m_blob_locks[i] = false;
  }

}

/// Return the number of blobs currently in use.
int vw::platefile::BlobManager::num_blobs() {
  Mutex::Lock lock(m_mutex);
  return m_blob_locks.size();
}

vw::int64 vw::platefile::BlobManager::max_blob_size() { 
  Mutex::Lock lock(m_mutex);
  return m_max_blob_size; 
}

/// Request a blob to write to that has sufficient space to write at
/// least 'size' bytes.  Returns the blob index of a locked blob
/// that you have sole access to write to.
///
/// size is specified in bytes.
//
// TODO: This is pretty simple logic, and would not be very
// efficient because it blocks on write if it catches up to a blob
// that is still locked.  We should add real blob selection logic
// here at a later date.
int vw::platefile::BlobManager::request_lock(int64 size) {
  Mutex::Lock lock(m_mutex);

  // First, we check to see if the next blob is free.  If not, we
  // wait for a release event and recheck.
  while(m_blob_locks[m_blob_index] != false)
    m_blob_release_condition.wait(lock);
  
  // Then we lock it, increment the blob index, and return it.
  m_blob_locks[m_blob_index] = true;      
  int idx = m_blob_index;
  next_blob_index();
  return idx;
}

// Release the blob lock and update its write index (essentially
// "committing" the write to the blob when you are finished with it.).
int vw::platefile::BlobManager::release_lock(int blob_id) {
  Mutex::Lock lock(m_mutex);
  m_blob_locks[blob_id] = false;
  m_blob_release_condition.notify_all();
}

