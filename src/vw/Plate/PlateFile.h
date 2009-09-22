#ifndef __VW_PLATE_PLATEFILE_H__
#define __VW_PLATE_PLATEFILE_H__

#include <vw/Math/Vector.h>
#include <vw/Image.h>
#include <vw/Mosaic/ImageComposite.h>

#include <vw/Plate/Index.h>
#include <vw/Plate/Blob.h>

#include <vector>
#include <fstream>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>

// Protocol Buffer
#include <vw/Plate/IndexRecord.pb.h>

namespace fs = boost::filesystem;

namespace vw {
namespace platefile {

  // -------------------------------------------------------------------------
  //                            TEMPORARY TILE FILE
  // -------------------------------------------------------------------------

  // A scoped temporary file object that store a tile under /tmp.  The
  // file is destroyed when this object is deleted.
  class TemporaryTileFile {

    std::string m_filename;

    // Define these as private methods to enforce TemporaryTileFile'
    // non-copyable semantics.
    TemporaryTileFile() {}
    TemporaryTileFile(TemporaryTileFile const&) {}
    TemporaryTileFile& operator=(TemporaryTileFile const&) { return *this; }
  
  public:

    static std::string unique_tempfile_name(std::string file_extension) {
      char base_name[100] = "/tmp/vw_plate_tile_XXXXXXX";
      std::string name = mktemp(base_name);
      return name + "." + file_extension;
    }

    /// This constructor assumes control over an existing file on disk,
    /// and deletes it when the TemporaryTileFile object is de-allocated.
    TemporaryTileFile(std::string filename) : m_filename(filename) {
      vw_out(DebugMessage, "plate::tempfile") << "Assumed control of temporary file: " 
                                              << m_filename << "\n";

    }

    /// This constructor assumes control over an existing file on disk,
    /// and deletes it when the TemporaryTileFile object is de-allocated.
    template <class ViewT>
    TemporaryTileFile(ImageViewBase<ViewT> const& view, std::string file_extension) : 
      m_filename(unique_tempfile_name(file_extension)) {
      write_image(m_filename, view);
      vw_out(DebugMessage, "plate::tempfile") << "Created temporary file: " 
                                              << m_filename << "\n";
    }

    ~TemporaryTileFile() {
      int result = unlink(m_filename.c_str());
      if (result)
        vw_out(ErrorMessage, "plate::tempfile") 
          << "WARNING: unlink() failed in ~TemporaryTileFile() for filename \"" 
          << m_filename << "\"\n";
      vw_out(DebugMessage, "plate::tempfile") << "Destroyed temporary file: " 
                                              << m_filename << "\n";
    }

    std::string file_name() const { return m_filename; }

    /// Opens the temporary file and determines its size in bytes.
    int64 file_size() const {
      std::ifstream istr(m_filename.c_str(), std::ios::binary);
      
      if (!istr.is_open())
        vw_throw(IOErr() << "TempPlateFile::file_size() -- could not open \"" 
                 << m_filename << "\".");
      istr.seekg(0, std::ios_base::end);
      return istr.tellg();
    }

    template <class PixelT> 
    ImageView<PixelT> read() const {
      ImageView<PixelT> img;
      read_image(img, m_filename);
      return img;
    }
          
  };


  // -------------------------------------------------------------------------
  //                            PLATE FILE
  // -------------------------------------------------------------------------

  
  class PlateFile {
    std::string m_plate_name;
    int m_default_block_size;;
    std::string m_default_file_type;
    boost::shared_ptr<IndexBase> m_index;
    FifoWorkQueue m_queue;

  public:
  
    PlateFile(std::string filename, int default_block_size = 256, 
              std::string default_file_type = "jpg") : 
      m_plate_name(filename), m_default_block_size(default_block_size),
      m_default_file_type(default_file_type) {

      // Plate files are stored as an index file and one or more data
      // blob files in a directory.  We create that directory here if
      // it doesn't already exist.
      if( !exists( fs::path( filename, fs::native ) ) ) {
        fs::create_directory(filename);
        m_index = boost::shared_ptr<IndexBase>( new Index(default_block_size, default_file_type) );
        vw_out(DebugMessage, "platefile") << "Creating new plate file: \"" 
                                          << filename << "\"\n";

      // However, if it does exist, then we attempt to open the
      // platefile that is stored there.
      } else {
        try {
          m_index = boost::shared_ptr<IndexBase>( new Index(filename + "/plate.index") );
          vw_out(DebugMessage, "platefile") << "Re-opened plate file: \"" 
                                            << filename << "\"\n";
          m_default_block_size = m_index->default_block_size();
          m_default_file_type = m_index->default_block_filetype();
        } catch (IOErr &e) {
          m_index = boost::shared_ptr<IndexBase>( new Index(default_block_size, default_file_type) );
          vw_out(DebugMessage, "platefile") << "WARNING: Failed to re-open the plate file.  "
                                            << "Creating a new plate file from scratch.";
        }
      }

    }

    /// The destructor saves the platefile to disk. 
    ~PlateFile() { this->save(); }

    /// Force the Platefile to save its index to disk.
    void save() { m_index->save(m_plate_name + "/plate.index"); }

    /// Returns the name of the root directory containing the plate file.
    std::string name() const { return m_plate_name; }

    /// Returns the file type used to store tiles in this plate file.
    std::string default_file_type() const { return m_default_file_type; }

    int default_block_size() const { return m_index->default_block_size(); }

    int depth() const { return m_index->max_depth(); }

    // Vector<int,2>  size_in_tiles() const  {}// return m_index.size(); }

    // Vector<long,2> size_in_pixels() const {}// return size_in_tiles() * tile_size(); }

    /// Read an image from the specified block location in the plate file.
    template <class PixelT>
    IndexRecord read(ImageView<PixelT> &view, int col, int row, int depth) {

      // 1. Call index read_request(col,row,depth).  Returns IndexRecord.
      IndexRecord record = m_index->read_request(col, row, depth);
      if (record.valid()) {
        std::ostringstream blob_filename;
        blob_filename << m_plate_name << "/plate_" << record.blob_id() << ".blob";

        // 2. Choose a temporary filename
        std::string tempfile = TemporaryTileFile::unique_tempfile_name(record.block_filetype());

        // 3. Call BlobIO read_as_file(filename, offset, size) [ offset, size from IndexRecord ]
        Blob blob(blob_filename.str());
        blob.read_to_file(tempfile, record.blob_offset(), record.block_size());
        TemporaryTileFile tile(tempfile);

        // 4. Read data from temporary file.
        view = tile.read<PixelT>();
      }
      return record;
    }

    /// Write an image to the specified block location in the plate file.
    template <class ViewT>
    void write(int col, int row, int depth, ImageViewBase<ViewT> const& view) {      

      // 1. Write data to temporary file. 
      TemporaryTileFile tile(view, m_default_file_type);
      std::string tile_filename = tile.file_name();
      int64 tile_size = tile.file_size();

      // 2. Make write_request(size) to index. Returns blob id.
      int blob_id = m_index->write_request(tile_size);
      std::ostringstream blob_filename;
      blob_filename << m_plate_name << "/plate_" << blob_id << ".blob";

      // 3. Create a blob and call write_from_file(filename).  Returns offset, size.
      Blob blob(blob_filename.str());
      int64 blob_offset;
      int32 block_size;
      blob.write_from_file(tile_filename, blob_offset, block_size);

      // 4. Call write_complete(col, row, depth, record)
      IndexRecord write_record;
      write_record.set_blob_id(blob_id);
      write_record.set_blob_offset(blob_offset);
      write_record.set_block_size(block_size);
      write_record.set_block_filetype(m_default_file_type);
      write_record.set_valid(true);
      m_index->write_complete(col, row, depth, write_record);
    }

    // Returns true if the 
    IndexRecord read_record(int col, int row, int depth) {
      return m_index->read_request(col, row, depth);
    }
  };


}} // namespace vw::plate

#endif // __VW_PLATE_PLATEFILE_H__
