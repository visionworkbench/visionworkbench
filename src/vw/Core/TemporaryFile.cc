#include <vw/Core/TemporaryFile.h>
#include <vw/Core/Exception.h>
#include <vw/Core/Log.h>
#include <vw/Core/Settings.h>
#include <boost/scoped_array.hpp>
#include <cerrno>
#include <cstring>
#include <vw/config.h>

#if defined(VW_HAVE_EXT_STDIO_FILEBUF_H) && VW_HAVE_EXT_STDIO_FILEBUF_H == 1
# include <ext/stdio_filebuf.h>
#endif

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

#if defined(VW_HAVE_EXT_STDIO_FILEBUF_H) && VW_HAVE_EXT_STDIO_FILEBUF_H == 1
    m_buf.reset(new __gnu_cxx::stdio_filebuf<char>(fd, mode));
#else
    m_buf.reset(new std::filebuf());
    m_buf->open(m_filename.c_str(), mode);
    if (m_buf->is_open())
      ::close(fd);
#endif

  if (!m_buf->is_open()) {
    ::close(fd);
    vw_throw(IOErr() << "Failed to create temporary file " << m_filename);
  }

  stream_t::init(m_buf.get());
  m_delete = delete_on_close;
}

TemporaryFile::TemporaryFile() {
  init();
}

TemporaryFile::TemporaryFile(const std::string& dir) {
  init(dir);
}

TemporaryFile::TemporaryFile(const std::string& dir, bool delete_on_close) {
  init(dir, delete_on_close);
}

TemporaryFile::TemporaryFile(const std::string& dir, bool delete_on_close, const std::string& prefix) {
  init(dir, delete_on_close, prefix);
}

TemporaryFile::TemporaryFile(const std::string& dir, bool delete_on_close, const std::string& prefix, const std::string& suffix) {
  init(dir, delete_on_close, prefix, suffix);
}

TemporaryFile::TemporaryFile(const std::string& dir, bool delete_on_close, const std::string& prefix, const std::string& suffix, std::ios_base::openmode mode) {
  init(dir, delete_on_close, prefix, suffix, mode);
}

TemporaryFile::~TemporaryFile() {
  stream_t::init(0);
  m_buf.reset();
  if (m_delete)
    if (::remove(m_filename.c_str()) == -1 && errno != ENOENT)
      vw_out(WarningMessage) << "Failed to remove temporary file " << m_filename << ": " << ::strerror(errno) << std::endl;
}

const std::string& TemporaryFile::filename() const {return m_filename;}

} // namespace vw
