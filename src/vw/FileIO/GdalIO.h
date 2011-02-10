#ifndef __VW_FILEIO_GDALIO_H__
#define __VW_FILEIO_GDALIO_H__

#include <vw/FileIO/ScanlineIO.h>

extern "C" {
#include <gdal_priv.h>
#include <cpl_multiproc.h>
}

namespace vw {
  class Mutex;

namespace fileio {
namespace detail {

Mutex& gdal() VW_WARN_UNUSED;

// These classes exist to share code between the on-disk and in-memory versions
// of the relevant image resources. They are not intended for use by users
// (thus the detail namespace).

class GdalIODecompress : public ScanlineReadBackend {
  protected:
    boost::shared_ptr<GDALDataset> m_dataset;
  public:
    GdalIODecompress();
    virtual ~GdalIODecompress();

    void open();
    bool ready() const;
    void read(uint8* data, size_t bufsize);
};

class GdalIOCompress : public ScanlineWriteBackend {
  protected:
    std::string m_fn;
    GDALDriver *m_driver;
    boost::shared_ptr<GDALDataset> m_dataset;

  public:
    // cols/rows/planes ignored in imageformat
    GdalIOCompress(const ImageFormat& fmt);
    virtual ~GdalIOCompress();

    void open();
    bool ready() const;
    void write(const uint8* buffer, size_t bufsize, size_t rows, size_t cols, size_t planes);
};

}}} // namespace vw::fileio::detail

#endif
