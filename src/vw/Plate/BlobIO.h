#ifndef __VW_PLATE_BLOBIO__
#define __VW_PLATE_BLOBIO__

#include <fstream>
#include <string>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/shared_array.hpp>
#include <vw/Core/Exception.h>
#include <vw/Core/Log.h>

namespace fs = boost::filesystem;


namespace vw {
namespace platefile {

  class Blob {
  
    int m_blobfile_version;
    std::string m_filename;

  public:
    
    // Constructor stores the blob filename for reading & writing
    Blob(std::string filename) : m_filename(filename) {
    
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

    int version() const { return m_blobfile_version; }

    // Returns: the binary data in the blob.
    template <class DataT>
    boost::shared_array<DataT> read(size_t offset, size_t size) {

      // Create a new array to store the data that gets read in.
      boost::shared_array<DataT> data(new DataT[size]);

      // Open file, and seek to the very end.
      std::ifstream istr(m_filename.c_str(), std::ios::binary);
      if (!istr.is_open())
        vw_throw(IOErr() << "Blob::read() -- could not open blob file for reading.");

      // Seek to the end and make sure that we haven't requested more
      // bytes than are store in the file!
      istr.seekg(0, std::ios_base::end);
      size_t end_offset = istr.tellg();

      if (end_offset - offset < size) 
        vw_throw(IOErr() << "Blob::read() failed -- read requested more bytes "
                 << " than are availabile in the file.");

      // Seek to the requested offset and start reading.
      istr.seekg(offset, std::ios_base::beg);     
      
      vw_out(VerboseDebugMessage, "plate::blob") << "Blob::reading() -- reading " 
                                                 << size << " bytes at " << offset 
                                                 << " from " << m_filename << "\n";
      istr.read((char*)(data.get()), size);
      istr.close();
      return data;
    }

    // Returns: the offset where the data was written to the blob file.
    template <class DataT>
    size_t write(boost::shared_array<DataT> data, size_t size) {

      // Open file, and seek to the very end.
      std::ofstream ostr(m_filename.c_str(), std::ios::app | std::ios::binary);
      ostr.seekp(0, std::ios_base::end);

      // Store the current offset of the end of the file.  We'll
      // return that at the end of this function.
      size_t offset = ostr.tellp();

      if (!ostr.is_open())
        vw_throw(IOErr() << "Blob::write() -- could not open blob file for writing.");

      vw_out(VerboseDebugMessage, "plate::blob") << "Blob::write() -- writing " << size
                                                 << " bytes to " << m_filename << "\n";
      ostr.write((char*)(data.get()), size);
      ostr.close();
      
      return offset;
    }

    // void read_as_file(std::string dest_filename, size_t offset, size_t size) {

    // }

    // void write_from_file(std::string source_file) {

    // }

  };

}} // namespace vw::platefile

#endif // __VW_PLATE_BLOBIO__
