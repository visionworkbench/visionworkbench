#include <vw/FileIO/ScanlineIO.h>

namespace vw { namespace fileio { namespace detail {

const ImageFormat& ScanlineBackend::fmt() const {
  return m_fmt;
}

size_t ScanlineBackend::chan_bytes() const { return m_cstride; }
size_t ScanlineReadBackend::line_bytes() const { return m_rstride; }

}}} // vw::fileio::detail
