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


#include <vw/FileIO/TemporaryFile.h>
#include <vw/Core/Exception.h>
#include <vw/Core/Log.h>
#include <vw/Core/Settings.h>
#include <vw/Core/System.h>

#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <unistd.h>

#include <boost/smart_ptr/scoped_array.hpp>
#include <boost/filesystem.hpp>

namespace fs = boost::filesystem;

namespace vw {

void TemporaryFile::init(std::string dir = "", bool delete_on_close = true, const std::string& prefix = "tmp", const std::string& suffix = "", std::ios_base::openmode mode = std::ios_base::binary|std::ios_base::out|std::ios_base::in)
{
  int fd;
  if (dir.empty())
    dir = vw_settings().tmp_directory();

  {
    std::string templ_s = dir + "/" + prefix + "XXXXXX" + suffix;
    boost::scoped_array<char> templ(new char[templ_s.size()+1]);
    ::strcpy(templ.get(), templ_s.c_str());
    fd = ::mkstemps(templ.get(), suffix.size());
    if (fd == -1)
      vw_throw(IOErr() << "Failed to create temporary file from template " << templ_s << ": " << ::strerror(errno));
    m_filename = std::string(templ.get());
  }

    m_buf.reset(new std::filebuf());
    m_buf->open(m_filename.c_str(), mode);
    if (m_buf->is_open())
      ::close(fd);

  if (!m_buf->is_open()) {
    ::close(fd);
    vw_throw(IOErr() << "Failed to create temporary file " << m_filename);
  }

  std::iostream::init(m_buf.get());
  m_delete = delete_on_close;
}

TemporaryFile::TemporaryFile() : std::iostream(0)  {
  init();
}

TemporaryFile::TemporaryFile(const std::string& dir) : std::iostream(0) {
  init(dir);
}

TemporaryFile::TemporaryFile(const std::string& dir, bool delete_on_close) : std::iostream(0) {
  init(dir, delete_on_close);
}

TemporaryFile::TemporaryFile(const std::string& dir, bool delete_on_close, const std::string& prefix) : std::iostream(0) {
  init(dir, delete_on_close, prefix);
}

TemporaryFile::TemporaryFile(const std::string& dir, bool delete_on_close, const std::string& prefix, const std::string& suffix) : std::iostream(0) {
  init(dir, delete_on_close, prefix, suffix);
}

TemporaryFile::TemporaryFile(const std::string& dir, bool delete_on_close, const std::string& prefix, const std::string& suffix, std::ios_base::openmode mode) : std::iostream(0) {
  init(dir, delete_on_close, prefix, suffix, mode);
}

TemporaryFile::~TemporaryFile() {
  std::iostream::init(0);
  m_buf.reset();
  if (m_delete)
    if (::remove(m_filename.c_str()) == -1 && errno != ENOENT)
      VW_OUT(WarningMessage) << "Failed to remove temporary file " << m_filename << ": " << ::strerror(errno) << std::endl;
}

const std::string& TemporaryFile::filename() const {return m_filename;}

void TemporaryDir::init(std::string dir = "", bool delete_on_close = true, const std::string& prefix = "tmp") {
  if (dir.empty())
    dir = vw_settings().tmp_directory();

  {
    std::string templ_s = dir + "/" + prefix + "XXXXXX";
    boost::scoped_array<char> templ(new char[templ_s.size()+1]);
    ::strcpy(templ.get(), templ_s.c_str());
    char *ret = ::mkdtemp(templ.get());
    VW_ASSERT(ret, IOErr() << "Failed to create temporary dir from template " << templ_s << ": " << ::strerror(errno));
    m_filename = std::string(templ.get());
  }

  m_delete = delete_on_close;
}

TemporaryDir::TemporaryDir() {
  init();
}

TemporaryDir::TemporaryDir(const std::string& dir) {
  init(dir);
}

TemporaryDir::TemporaryDir(const std::string& dir, bool delete_on_close) {
  init(dir, delete_on_close);
}

TemporaryDir::TemporaryDir(const std::string& dir, bool delete_on_close, const std::string& prefix) {
  init(dir, delete_on_close, prefix);
}

TemporaryDir::~TemporaryDir() {
  if (m_delete)
    fs::remove_all(m_filename);
}

const std::string& TemporaryDir::filename() const {
  return m_filename;
}

} // namespace vw
