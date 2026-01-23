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
#include <vw/vw_config.h>
#include <vw/Core/Exception.h>
#include <vw/Core/Log.h>
#include <vw/Core/Settings.h>
#include <vw/Core/System.h>

#include <cerrno>
#include <cstring>

#include <boost/smart_ptr/scoped_array.hpp>
#include <boost/filesystem.hpp>

namespace fs = boost::filesystem;

#if defined(VW_HAVE_EXT_STDIO_FILEBUF_H)
# include <ext/stdio_filebuf.h>
#endif

#if defined(VW_HAVE_MKSTEMPS)
#  define vw_mkstemps ::mkstemps
#else

#  include <sys/types.h>
#  include <sys/stat.h>
#  include <fcntl.h>
#  ifdef VW_HAVE_UNISTD_H
#    include <unistd.h>
#  endif

#  ifdef WIN32
#    define stat _stat
     typedef struct _stat struct_stat;
#  else
     typedef struct stat struct_stat;
#  endif

namespace {
// NOTE! arc4random replaced by random, which is not as strong.
// This mkstemps impl is:
/*  $OpenBSD: mktemp.c,v 1.19 2005/08/08 08:05:36 espie Exp $ */
/*
 * Copyright (c) 1987, 1993
 *  The Regents of the University of California.  All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. Neither the name of the University nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 */

static int
_gettemp(char *path, int *doopen, int domkdir, int slen)
{
  char *start, *trv, *suffp;
  struct stat sbuf;
  int rval;
  pid_t pid;

  if (doopen && domkdir) {
    errno = EINVAL;
    return(0);
  }

  for (trv = path; *trv; ++trv)
    ;
  trv -= slen;
  suffp = trv;
  --trv;
  if (trv < path) {
    errno = EINVAL;
    return (0);
  }
  pid = getpid();
  while (trv >= path && *trv == 'X' && pid != 0) {
    *trv-- = (pid % 10) + '0';
    pid /= 10;
  }
  while (trv >= path && *trv == 'X') {
    char c;

    pid = (random() & 0xffff) % (26+26);
    if (pid < 26)
      c = pid + 'A';
    else
      c = (pid - 26) + 'a';
    *trv-- = c;
  }
  start = trv + 1;

  /*
   * check the target directory; if you have six X's and it
   * doesn't exist this runs for a *very* long time.
   */
  if (doopen || domkdir) {
    for (;; --trv) {
      if (trv <= path)
        break;
      if (*trv == '/') {
        *trv = '\0';
        rval = stat(path, &sbuf);
        *trv = '/';
        if (rval != 0)
          return(0);
        if (!S_ISDIR(sbuf.st_mode)) {
          errno = ENOTDIR;
          return(0);
        }
        break;
      }
    }
  }

  for (;;) {
    if (doopen) {
      if ((*doopen =
          open(path, O_CREAT|O_EXCL|O_RDWR, 0600)) >= 0)
        return(1);
      if (errno != EEXIST)
        return(0);
    } else if (domkdir) {
      if (mkdir(path, 0700) == 0)
        return(1);
      if (errno != EEXIST)
        return(0);
    } else if (lstat(path, &sbuf))
      return(errno == ENOENT ? 1 : 0);

    /* tricky little algorithm for backward compatibility */
    for (trv = start;;) {
      if (!*trv)
        return (0);
      if (*trv == 'Z') {
        if (trv == suffp)
          return (0);
        *trv++ = 'a';
      } else {
        if (isdigit(*trv))
          *trv = 'a';
        else if (*trv == 'z') /* inc from z to A */
          *trv = 'A';
        else {
          if (trv == suffp)
            return (0);
          ++*trv;
        }
        break;
      }
    }
  }
  /*NOTREACHED*/
}

int vw_mkstemps(char *path, int slen)
{
  int fd;
  return (_gettemp(path, &fd, 0, slen) ? fd : -1);
}
}
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
    fd = vw_mkstemps(templ.get(), suffix.size());
    if (fd == -1)
      vw_throw(IOErr() << "Failed to create temporary file from template " << templ_s << ": " << ::strerror(errno));
    m_filename = std::string(templ.get());
  }

#if defined(VW_HAVE_EXT_STDIO_FILEBUF_H)
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
