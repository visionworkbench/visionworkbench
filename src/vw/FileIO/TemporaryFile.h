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


#ifndef __VW_CORE_TEMPORARYFILE_H__
#define __VW_CORE_TEMPORARYFILE_H__

#include <iostream>
#include <boost/smart_ptr/shared_ptr.hpp>

/// \file TemporaryFile.h Temporary file and directory classes.

namespace vw {

class TemporaryFile : public std::iostream {
  private:

    boost::shared_ptr<std::filebuf> m_buf;
    std::string m_filename;
    bool m_delete;

    void init(std::string dir, bool delete_on_close, const std::string& prefix, const std::string& suffix, std::ios_base::openmode mode);
  public:
    // defaults:
    //    dir = vw_settings().tmp_directory()
    //    delete_on_close = true
    //    prefix = tmp
    //    suffix = <empty>
    //    mode = binary | out | in

    TemporaryFile();
    TemporaryFile(const std::string& dir);
    TemporaryFile(const std::string& dir, bool delete_on_close);
    TemporaryFile(const std::string& dir, bool delete_on_close, const std::string& prefix);
    TemporaryFile(const std::string& dir, bool delete_on_close, const std::string& prefix, const std::string& suffix);
    TemporaryFile(const std::string& dir, bool delete_on_close, const std::string& prefix, const std::string& suffix, std::ios_base::openmode mode);
    ~TemporaryFile();

    const std::string& filename() const;
};

class TemporaryDir {
  private:
    std::string m_filename;
    bool m_delete;

    void init(std::string dir, bool delete_on_close, const std::string& prefix);
  public:
    TemporaryDir();
    TemporaryDir(const std::string& dir);
    TemporaryDir(const std::string& dir, bool delete_on_close);
    TemporaryDir(const std::string& dir, bool delete_on_close, const std::string& prefix);
    ~TemporaryDir();

    const std::string& filename() const;
};

}

#endif
