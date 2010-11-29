#ifndef __VW_FILEIO_JPEGIO_H__
#define __VW_FILEIO_JPEGIO_H__

#include <vw/FileIO/ScanlineIO.h>

extern "C" {
#include <png.h>
}


namespace vw {
namespace fileio {
namespace detail {

// These classes exist to share code between the on-disk and in-memory versions
// of the relevant image resources. They are not intended for use by users
// (thus the detail namespace).

class PngIO {
  protected:
    png_structp m_ctx;
    png_infop m_info;
};

class PngIODecompress : public PngIO, public ScanlineReadBackend {
  private:
    bool m_read;
  public:
    PngIODecompress();
    virtual ~PngIODecompress();

    void open();
    bool ready() const;
    void read(uint8* data, size_t bufsize);
};

class PngIOCompress : public PngIO, public ScanlineWriteBackend {
  private:
    bool m_written;
  public:
    // cols/rows/planes ignored in imageformat
    PngIOCompress(const ImageFormat& fmt);
    virtual ~PngIOCompress();

    void open();
    bool ready() const;
    void write(const uint8* buffer, size_t rows, size_t cols, size_t bufsize);
};

}}} // namespace vw::fileio::detail

#endif
