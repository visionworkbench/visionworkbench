// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_PLATE_MOD_PLATE_UTILS_H__
#define __VW_PLATE_MOD_PLATE_UTILS_H__

#include <vw/Core/Exception.h>
#include <vw/Plate/HTTPUtils.h>

#include <boost/lexical_cast.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/function.hpp>
#include <set>
#include <map>
#include <string>

struct apr_pool_t;
struct request_rec;

namespace vw {
namespace platefile {

class Index;

VW_DEFINE_EXCEPTION(PlateException,   Exception);
VW_DEFINE_EXCEPTION(BadRequest,       PlateException);
VW_DEFINE_EXCEPTION(ServerError,      PlateException);
VW_DEFINE_EXCEPTION(UnknownPlatefile, PlateException);

// Returns the const char* as a string, or an empty() string if it's scary in
// some way.
std::string safe_string_convert(const char* in);

class apache_output : public boost::iostreams::sink {
    request_rec *r;
  public:
    apache_output(request_rec *r) : r(r) {}
    std::streamsize write(const char* s, std::streamsize n);
};

typedef boost::iostreams::stream<apache_output> apache_stream;

struct raii {
  typedef boost::function<void (void)> FuncT;
  FuncT m_leave;
  raii(FuncT enter, FuncT leave) : m_leave(leave) {enter();};
  ~raii() {m_leave();}
};

class WTMLImageSet : public std::map<std::string, std::string> {
  typedef std::map<std::string, std::string> map_t;
  typedef std::set<std::string> child_t;
  child_t child_keys;

  public:
    WTMLImageSet(const std::string& host,
                 const std::string& data_prefix,
                 const std::string& static_prefix,
                 const std::string dem_id,
                 boost::shared_ptr<const Index> index,
                 const std::string& description);

    void serializeToOstream(std::ostream& o) const;
};

class ApacheRequest {
    request_rec *r;
  public:
    std::string url;
    QueryMap args;

    ApacheRequest(request_rec* r);
    bool header_only() const;
    request_rec* writer() const;
};


}} // namespace vw::platefile

#endif
