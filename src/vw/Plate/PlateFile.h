#ifndef __VW_PLATE_PLATEFILE_H__
#define __VW_PLATE_PLATEFILE_H__

#include <vw/Math/Vector.h>


namespace vw {
namespace plate {

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

    static std::string unique_tempfile_name(std::string file_extension, 
                                            std::string base_name = /tmp/vw_plate_tile_XXXXXXX) {
      std::string filename = mktemp(base_name.c_str());
      return filename + "." + file_extension;
    }

    /// This constructor assumes control over an existing file on disk,
    /// and deletes it when the TemporaryTileFile object is de-allocated.
    TemporaryFile(std::string filename) : m_filename(filename) {
      vw_out(DebugMessage, "plate::tempfile") << "Assumed control of temporary file: " 
                                              << m_filename << "\n";

    }

    /// This constructor assumes control over an existing file on disk,
    /// and deletes it when the TemporaryTileFile object is de-allocated.
    template <class ViewT>
    TemporaryFile(ImageViewBase<ViewT> const& view, std::string file_extension) : 
      m_filename(unique_tempfile_name(file_extension)) {
      write_image(m_filename, view);
      vw_out(DebugMessage, "plate::tempfile") << "Created temporary file: " 
                                              << m_filename << "\n";
    }

    ~TemporaryFile() {
      int result = unlink(m_filename);
      if (result)
        vw_out(ErrorMessage, "plate::tempfile") 
          << "WARNING: unlink() failed in ~TemporaryTileFile() for filename \"" 
          << m_filename << "\"\n";
      vw_out(DebugMessage, "plate::tempfile") << "Destroyed temporary file: " 
                                              << m_filename << "\n";
    }

    std::string name() const { return m_filename; }

    /// Opens the temporary file and determines its size in bytes.
    size_t file_size() const {
      std::ifstream istr(source_file.c_str(), std::ios::binary);
      
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
    std::string m_tile_filetype;
  

  public:
  
    PlateFile(std::string filename, std::string tile_filetype) : 
      m_plate_name(filename), m_tile_filetype(tile_filetype) {}

    int name() const { return m_plate_name; }
    int tile_filetype() const { return m_tile_filetype; }


    void close() {}
  
    int tile_size() const {}// return m_index.tile_size(); }


    int depth() const {}// return m_index.depth(); }

    Vector<int,2>  size_in_tiles() const  {}// return m_index.size(); }

    Vector<long,2> size_in_pixels() const {}// return size_in_tiles() * tile_size(); }


    template <class PixelT>
    void read(ImageView<PixelT> &view, int col, int row, int depth) {

      // 1. Call index read_request(col,row,depth).  Returns IndexRecord.
      IndexRecord record = m_index.read_request(col, row, depth);
      ostringstream blob_filename;
      blob_filename << m_filename << "/plate_" << record.blob_id() << ".blob";

      // 2. Choose a temporary filename
      std::string tempfile = TemporaryTileFile::unique_tempfile_name(m_tile_filetype);

      // 3. Call BlobIO read_as_file(filename, offset, size) [ offset, size from IndexRecord ]
      Blob blob(blob_filename);
      blob.read_to_file(tempfile, record.offset(), record.size());
      TemporaryTileFile tile(tempfile);

      // 4. Read data from temporary file.
      view = tile.read<PixelT>();
    }

    template <class ViewT>
    void write(int col, int row, int depth, ImageViewBase<ViewT> const& view) {      

      // 1. Write data to temporary file. 
      TemporaryTileFile tile(view, m_tile_filetype);
      std::string tile_filename = tile.filename();
      size_t tile_size = tile.size();

      // 2. Make write_request(size) to index. Returns blob id.
      int blob_id = m_index.write_request(tile_size);
      ostringstream blob_filename;
      blob_filename << m_filename << "/plate_" << blob_id << ".blob";

      // 3. Create a blob and call write_from_file(filename).  Returns offset, size.
      Blob blob(blob_filename.str());
      size_t offset, size;
      blob.write_from_file(tile_filename, offset, size);

      // 4. Call write_complete(col, row, depth, blob_id, offset, size, filetype)
      m_index.write_complete(col, row, depth, blob_id, offset, size, m_tile_filetype.c_str());
    }

  };

}} // namespace vw::plate

#endif // __VW_PLATE_PLATEFILE_H__
