#ifndef __VW_PLATE_BLOBIO__
#define __VW_PLATE_BLOBIO__

#include <fstream>
#include <string>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/shared_array.hpp>
#include <vw/Core/Exception.h>
#include <vw/Core/FundamentalTypes.h>
#include <vw/Core/Log.h>

namespace fs = boost::filesystem;

namespace vw {
namespace platefile {

  // -------------------------------------------------------------------
  //                                 BLOB
  // -------------------------------------------------------------------

  class Blob {
  
    int m_blobfile_version;
    std::string m_filename;

  public:
    
    // Constructor stores the blob filename for reading & writing
    Blob(std::string filename);

    int version() const;

    // Returns: the binary data in the blob.
    template <class DataT>
    boost::shared_array<DataT> read(vw::int64 offset, vw::int64 size) {

      // Create a new array to store the data that gets read in.
      boost::shared_array<DataT> data(new DataT[size]);

      // Open file, and seek to the very end.
      std::ifstream istr(m_filename.c_str(), std::ios::binary);
      if (!istr.is_open())
        vw_throw(vw::IOErr() << "Blob::read() -- could not open blob file for reading.");

      // Seek to the end and make sure that we haven't requested more
      // bytes than are store in the file!
      istr.seekg(0, std::ios_base::end);
      vw::int64 end_offset = istr.tellg();

      if (end_offset - offset < size) 
        vw::vw_throw(vw::IOErr() << "Blob::read() failed -- read requested more bytes "
                 << " than are availabile in the file.");

      // Seek to the requested offset and start reading.
      istr.seekg(offset, std::ios_base::beg);     
      
      vw::vw_out(vw::VerboseDebugMessage, "plate::blob") << "Blob::reading() -- reading " 
                                                 << size << " bytes at " << offset 
                                                 << " from " << m_filename << "\n";
      istr.read((char*)(data.get()), size);
      istr.close();
      return data;
    }

    // Returns: the offset where the data was written to the blob file.
    template <class DataT>
    vw::int64 write(boost::shared_array<DataT> data, vw::int64 size) {
    
      // Open file, and seek to the very end.
      std::ofstream ostr(m_filename.c_str(), std::ios::app | std::ios::binary);
      if (!ostr.is_open())
        vw::vw_throw(vw::IOErr() << "Blob::write() -- could not open blob file for writing.");

      // Store the current offset of the end of the file.  We'll
      // return that at the end of this function.
      ostr.seekp(0, std::ios_base::end);
      vw::int64 offset = ostr.tellp();
      
      // Write the data at the end of the file and return the offset
      // of the beginning of this data file.
      vw::vw_out(vw::VerboseDebugMessage, "plate::blob") << "Blob::write() -- writing " << size
                                                 << " bytes to " << m_filename << "\n";
      ostr.write((char*)(data.get()), size);
      ostr.close();

      return offset;
    }
    

    /// Read data out of the blob and save it as its own file on disk.
    void read_to_file(std::string dest_file, vw::int64 offset, vw::int64 size);

    /// Write the data file to disk, and the concatenate it into the data blob.
    void write_from_file(std::string source_file, vw::int64& offset, vw::int32& size);

  };

  // -------------------------------------------------------------------
  //                            BLOB_MANAGER
  // -------------------------------------------------------------------

  /// The BlobManager keeps track of how much data has been stored in
  /// each blob so far.  This allows it to pick blobs that have enough
  /// space to store a new block of a given size.  The blob manager also
  /// controls the locking/unlocking of blobs, and can load balance
  /// blobs writes by alternating which blob is offered up for writing
  /// data.
  ///
  /// The BlobManager is thread safe.
  class BlobManager {
    vw::int64 m_max_blob_size;
    std::vector<bool> m_blob_locks;
    int m_blob_index;
    vw::Mutex m_mutex;
    vw::Condition m_blob_release_condition;

    void next_blob_index();

  public:

    /// Create a new blob manager.  The max_blob_size is specified in
    /// units of megabytes.
    BlobManager(vw::int64 max_blob_size = 2048, int nblobs = 2);

    /// Return the number of blobs currently in use.
    int num_blobs();

    vw::int64 max_blob_size();

    /// Request a blob to write to that has sufficient space to write at
    /// least 'size' bytes.  Returns the blob index of a locked blob
    /// that you have sole access to write to.
    ///
    /// size is specified in bytes.
    //
    // TODO: This is pretty simple logic so far, and would not be very
    // efficient because it blocks on write if it catches up to a blob
    // that is still locked.  We should add real blob selection logic
    // here at a later date.
    int request_lock(vw::int64 size);

    // Release the blob lock and update its write index (essentially
    // "committing" the write to the blob when you are finished with it.).
    int release_lock(int blob_id);
  };


}} // namespace vw::platefile

#endif // __VW_PLATE_BLOBIO__
