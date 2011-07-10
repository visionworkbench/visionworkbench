// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_CORE_TEMPORARYFILE_H__
#define __VW_CORE_TEMPORARYFILE_H__

#include <iostream>
#include <boost/shared_ptr.hpp>

namespace vw {

class TemporaryFile : public std::iostream {
  private:
    typedef std::iostream stream_t;

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
