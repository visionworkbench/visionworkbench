// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/FileIO/ScanlineIO.h>

namespace vw { namespace fileio { namespace detail {

const ImageFormat& ScanlineBackend::fmt() const {
  return m_fmt;
}

size_t ScanlineBackend::chan_bytes() const { return m_cstride; }
size_t ScanlineReadBackend::line_bytes() const { return m_rstride; }

}}} // vw::fileio::detail
