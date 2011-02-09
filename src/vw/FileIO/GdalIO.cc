#include <vw/FileIO/GdalIO.h>
#include <vw/Core/Exception.h>
#include <vw/Core/Log.h>
#include <vw/Core/Settings.h>
#include <vw/Core/Cache.h>

static void CPL_STDCALL gdal_error_handler(CPLErr eErrClass, int nError, const char *pszErrorMsg) {
  vw::MessageLevel lvl;

  switch(eErrClass) {
    case CE_Debug:
    case CE_Warning:
      lvl = vw::WarningMessage;
      break;
    default:
      lvl = vw::ErrorMessage;
      break;
  }

  std::string msg;
  if (pszErrorMsg)
    msg = pszErrorMsg;

  boost::replace_all(msg, "\n", " ");

  if (eErrClass == CE_Fatal)
    vw::vw_throw(vw::IOErr()  << "GdalIO: " << msg << " (code = " << nError << ")");
  else
    vw::vw_out(lvl, "fileio") << "GdalIO: " << msg << " (code = " << nError << ")" << std::endl;
}

// GDAL is not thread-safe, so we keep a global GDAL lock (pointed to
// by gdal_mutex_ptr, below) that we hold anytime we call into the
// GDAL library itself.

namespace {
  vw::RunOnce _gdal_init_once = VW_RUNONCE_INIT;
  vw::Mutex* _gdal_mutex;
  vw::Cache* _gdal_cache;

  void init_gdal() {
    CPLPushErrorHandler(gdal_error_handler);
    // If we run out of handles, GDAL will start closing them behind the
    // scenes. We maintain our own cache so we never hit their limit.
    CPLSetConfigOption("GDAL_MAX_DATASET_POOL_SIZE", "200");
    GDALAllRegister();
    _gdal_cache = new vw::Cache(200);
    _gdal_mutex = new vw::Mutex();
  }
  void kill_gdal() {
    delete _gdal_mutex;
    delete _gdal_cache;
    GDALDumpOpenDatasets(stderr);
    GDALDestroyDriverManager();
    CPLDumpSharedList(0);
    CPLCleanupTLS();
    CPLPopErrorHandler();
  }

  GDALColorInterp bcolor(const boost::shared_ptr<GDALDataset> d, size_t idx /* 1-indexed */) {
    if (ssize_t(idx) > d->GetRasterCount())
      return GCI_Max;
    return d->GetRasterBand(idx)->GetColorInterpretation();
  }

  vw::ChannelTypeEnum channel_gdal_to_vw(GDALDataType gdal_type) {
    switch (gdal_type) {
      case GDT_Byte:    return vw::VW_CHANNEL_UINT8;
      case GDT_Int16:   return vw::VW_CHANNEL_INT16;
      case GDT_UInt16:  return vw::VW_CHANNEL_UINT16;
      case GDT_Int32:   return vw::VW_CHANNEL_INT32;
      case GDT_UInt32:  return vw::VW_CHANNEL_UINT32;
      case GDT_Float32: return vw::VW_CHANNEL_FLOAT32;
      case GDT_Float64: return vw::VW_CHANNEL_FLOAT64;
      default: vw::vw_throw( vw::IOErr() << "Unsupported GDAL channel type (" << gdal_type << ")." );
    }
  }

  GDALDataType channel_vw_to_gdal(vw::ChannelTypeEnum vw_type) {
    switch (vw_type) {
      case vw::VW_CHANNEL_UINT8:   return GDT_Byte;
      case vw::VW_CHANNEL_INT16:   return GDT_Int16;
      case vw::VW_CHANNEL_UINT16:  return GDT_UInt16;
      case vw::VW_CHANNEL_INT32:   return GDT_Int32;
      case vw::VW_CHANNEL_UINT32:  return GDT_UInt32;
      case vw::VW_CHANNEL_FLOAT32: return GDT_Float32;
      case vw::VW_CHANNEL_FLOAT64: return GDT_Float64;
      default: vw::vw_throw( vw::IOErr() << "Unsupported vw channel type (" << vw_type << ")." );
    }
  }
}




namespace vw {
namespace fileio {
namespace detail {

vw::Mutex& gdal() {
  _gdal_init_once.run( init_gdal );
  return *_gdal_mutex;
}
vw::Cache& gdal_cache() {
  _gdal_init_once.run( init_gdal );
  return *_gdal_cache;
}

boost::shared_ptr<GDALDataset> GdalDatasetGenerator::generate() const {
  GDALDataset *dataset = (GDALDataset*) GDALOpen( m_filename.c_str(), GA_ReadOnly );
  if ( !dataset )
    vw_throw(ArgumentErr() << "DiskImageResourceGDAL: Could not open \"" << m_filename << "\"");
  return boost::shared_ptr<GDALDataset>(dataset, GDALClose);
}


////////////////////////////////////////////////////////////////////////////////
// Decompress
////////////////////////////////////////////////////////////////////////////////
GdalIODecompress::GdalIODecompress() { }
GdalIODecompress::~GdalIODecompress() { }

bool GdalIODecompress::ready() const { return true; }

void GdalIODecompress::read(uint8* buffer, size_t bufsize) {
  size_t skip = line_bytes();
  VW_ASSERT(bufsize >= m_fmt.rows * skip, LogicErr() << "Buffer is too small");
  if (m_fmt.pixel_format == VW_PIXEL_SCALAR) {
    // Separate bands
    m_dataset->RasterIO(GF_Read, 0, 0, m_fmt.cols, m_fmt.rows, reinterpret_cast<void*>(buffer), m_fmt.cols, m_fmt.rows,
        channel_vw_to_gdal(m_fmt.channel_type), num_channels(m_fmt.pixel_format), NULL, 0, 0, 0);
  } else {
    // Interleaved pixels
    m_dataset->RasterIO(GF_Read, 0, 0, m_fmt.cols, m_fmt.rows, reinterpret_cast<void*>(buffer), m_fmt.cols, m_fmt.rows,
        channel_vw_to_gdal(m_fmt.channel_type), num_channels(m_fmt.pixel_format), NULL,
        m_cstride, m_rstride, 1);
  }
}

void GdalIODecompress::open() {
  this->bind();

  m_fmt.rows = m_dataset->GetRasterYSize();
  m_fmt.cols = m_dataset->GetRasterXSize();

  size_t chans = m_dataset->GetRasterCount();
  VW_ASSERT(chans > 0, IOErr() << "Cannot read GDAL Image: unknown number of channels");

  if ( bcolor(m_dataset, 1) == GCI_GrayIndex) {
    if (chans == 1) {
      m_fmt.pixel_format = VW_PIXEL_GRAY;
      m_fmt.planes = 1;
    } else if (chans == 2 && bcolor(m_dataset, 2) == GCI_AlphaBand) {
      m_fmt.pixel_format = VW_PIXEL_GRAYA;
      m_fmt.planes = 1;
    }
  } else if ( bcolor(m_dataset, 1) == GCI_RedBand && bcolor(m_dataset, 2) == GCI_GreenBand && bcolor(m_dataset, 3) == GCI_BlueBand) {
    if (chans == 3) {
      m_fmt.pixel_format = VW_PIXEL_RGB;
      m_fmt.planes = 1;
    } else if (chans == 4 && bcolor(m_dataset, 4) == GCI_AlphaBand) {
      m_fmt.pixel_format = VW_PIXEL_RGBA;
      m_fmt.planes = 1;
    }
  }

  if (m_fmt.planes == 0) {
    m_fmt.planes = chans;
    m_fmt.pixel_format = VW_PIXEL_SCALAR;
  }

  m_fmt.channel_type = channel_gdal_to_vw(m_dataset->GetRasterBand(1)->GetRasterDataType());
  m_cstride = num_channels(fmt().pixel_format) * channel_size(fmt().channel_type);
  m_rstride = m_cstride * m_fmt.cols;
}

////////////////////////////////////////////////////////////////////////////////
// Compress
////////////////////////////////////////////////////////////////////////////////
GdalIOCompress::GdalIOCompress(const ImageFormat& fmt) {
  m_fmt = fmt;
}

GdalIOCompress::~GdalIOCompress() { }

bool GdalIOCompress::ready() const {
  return true;
}

void GdalIOCompress::write(const uint8* data, size_t bufsize, size_t rows, size_t cols, size_t planes) {
  size_t skip = cols * chan_bytes();
  VW_ASSERT(bufsize >= rows * skip, LogicErr() << "Buffer is too small");

  VW_ASSERT(planes == 1 || fmt().pixel_format==VW_PIXEL_SCALAR,
      NoImplErr() << "Multi-channel multi-plane images are not supported");

  m_dataset.reset(
    reinterpret_cast<GDALDataset*>(GDALCreate( m_driver, m_fn.c_str(), cols, rows, num_channels(m_fmt.pixel_format), channel_vw_to_gdal(fmt().channel_type), NULL)),
    GDALClose);

  if (fmt().pixel_format == VW_PIXEL_SCALAR) {
    // Separate bands
    m_dataset->RasterIO(GF_Write, 0, 0, cols, rows, const_cast<uint8*>(data), cols, rows,
                        channel_vw_to_gdal(fmt().channel_type), num_channels(m_fmt.pixel_format), NULL, 0, 0, 0);
  } else {
    //Interleaved pixels
    m_dataset->RasterIO(GF_Write, 0, 0, cols, rows, const_cast<uint8*>(data), cols, rows,
                        channel_vw_to_gdal(fmt().channel_type), num_channels(m_fmt.pixel_format), NULL,
                        m_cstride, skip, 1);
  }
  m_dataset.reset();
}

void GdalIOCompress::open() {
  m_driver = (GDALDriver*)GDALGetDriverByName("GTiff");
  VW_ASSERT(GDALGetMetadataItem(m_driver, GDAL_DCAP_VIRTUALIO, NULL ), NoImplErr() << "GDAL's current tiff driver does not support virtual IO");

  this->bind();
  m_cstride = num_channels(fmt().pixel_format) * channel_size(fmt().channel_type);
}

}}} // namespace vw::fileio::detail
