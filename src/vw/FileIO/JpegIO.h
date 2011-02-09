#ifndef __VW_FILEIO_JPEGIO_H__
#define __VW_FILEIO_JPEGIO_H__

#include <vw/FileIO/ScanlineIO.h>

extern "C" {
#include <jpeglib.h>
#include <jerror.h>
}


namespace vw {
namespace fileio {
namespace detail {

// These classes exist to share code between the on-disk and in-memory versions
// of the relevant image resources. They are not intended for use by users
// (thus the detail namespace).
//
// Note: This code does not support suspending data sources! As the libjpeg
// docs say, if you don't know what that means, you don't need it.

class JpegIO {
  private:
    jpeg_error_mgr m_jerr;
  protected:
    void init_base(jpeg_error_mgr** mgr);
};

class JpegIODecompress : public JpegIO, public ScanlineReadBackend {
  protected:
    jpeg_decompress_struct m_ctx;
  public:
    JpegIODecompress();
    virtual ~JpegIODecompress();

    void open();
    bool ready() const;
    void read(uint8* data, size_t bufsize);
};

class JpegIOCompress : public JpegIO, public ScanlineWriteBackend {
  protected:
    jpeg_compress_struct m_ctx;
  public:
    // cols/rows/planes ignored in imageformat
    JpegIOCompress(const ImageFormat& fmt);
    virtual ~JpegIOCompress();

    void open();
    bool ready() const;
    void write(const uint8* buffer, size_t bufsize, size_t rows, size_t cols, size_t planes);
};

void jpeg_ptr_src(j_decompress_ptr cinfo, const uint8* buffer, size_t size);
void jpeg_vector_dest(j_compress_ptr cinfo, std::vector<uint8>* v);

}}} // namespace vw::fileio::detail

#endif
