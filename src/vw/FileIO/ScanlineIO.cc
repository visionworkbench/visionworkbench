// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NASA Vision Workbench is licensed under the Apache License,
//  Version 2.0 (the "License"); you may not use this file except in
//  compliance with the License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__


#include <vw/FileIO/ScanlineIO.h>

namespace vw { namespace fileio { namespace detail {

const ImageFormat& ScanlineBackend::fmt() const {
  return m_fmt;
}

size_t ScanlineBackend::chan_bytes() const { return m_cstride; }
size_t ScanlineReadBackend::line_bytes() const { return m_rstride; }

}}} // vw::fileio::detail
