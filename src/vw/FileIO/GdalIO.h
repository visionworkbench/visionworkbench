#ifndef __VW_FILEIO_JPEGIO_H__
#define __VW_FILEIO_JPEGIO_H__

#include <vw/FileIO/ScanlineIO.h>

extern "C" {
#include <gdal_priv.h>
#include <cpl_multiproc.h>
}

namespace vw {
  class Mutex;
  class Cache;

namespace fileio {
namespace detail {

Mutex& gdal() VW_WARN_UNUSED;
Cache& gdal_cache() VW_WARN_UNUSED;

class GdalDatasetGenerator {
    std::string m_filename;
  public:
    typedef GDALDataset value_type;
    GdalDatasetGenerator( std::string filename ) : m_filename( filename ) {}

    size_t size() const { return 1; }

    boost::shared_ptr<GDALDataset> generate() const;
};


// These classes exist to share code between the on-disk and in-memory versions
// of the relevant image resources. They are not intended for use by users
// (thus the detail namespace).

class GdalIO {
  protected:
    boost::shared_ptr<GDALDataset> m_dataset;
};

class GdalIODecompress : public GdalIO, public ScanlineReadBackend {
  public:
    GdalIODecompress();
    virtual ~GdalIODecompress();

    void open();
    bool ready() const;
    void read(uint8* data, size_t bufsize);
};

class GdalIOCompress : public GdalIO, public ScanlineWriteBackend {
  protected:
    GDALDriver *m_driver;
    std::string m_fn;

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
