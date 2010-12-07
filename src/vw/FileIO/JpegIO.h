#ifndef __VW_FILEIO_JPEGIO_H__
#define __VW_FILEIO_JPEGIO_H__

#include <vw/Image/ImageResource.h>

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
    ImageFormat m_fmt;
    size_t m_cstride, m_rstride;

    // attach a data source/sink
    virtual void bind() = 0;

    void init_base(jpeg_error_mgr** mgr);

  public:
    virtual ~JpegIO() {}

    const ImageFormat& fmt() const;

    void reopen();
    size_t chan_bytes();
    size_t line_bytes();

    // After this, protected members should be initialized,  as well as the
    // data source (i.e., your implementation should call bind())
    virtual void open() = 0;
    // return the line the current data source/sink is on
    virtual size_t line() const = 0;
    // resolve pending operations
    virtual void flush() = 0;
    // close data source/sink (paired with open)
    virtual void close() = 0;
};

class JpegIODecompress : public JpegIO {
  protected:
    jpeg_decompress_struct m_ctx;
  public:
    JpegIODecompress();
    virtual ~JpegIODecompress();

    void open();
    size_t line() const;
    void flush();
    void close();

    // this will decode line_bytes() worth of data into 'data'. data must be
    // large enough to contain this.
    void readline(uint8* data);
};

class JpegIOCompress : public JpegIO {
  protected:
    jpeg_compress_struct m_ctx;
  public:
    // cols/rows/planes ignored in imageformat
    JpegIOCompress(const ImageFormat& fmt);
    virtual ~JpegIOCompress();

    void open();
    size_t line() const;
    void flush();
    void close();

    void startwrite(uint32 image_width, uint32 image_height);

    // You must call startwrite before writeline.
    // data must contain line_bytes() worth of data
    void writeline(uint8* data);
};

}}} // namespace vw::fileio::detail

#endif
