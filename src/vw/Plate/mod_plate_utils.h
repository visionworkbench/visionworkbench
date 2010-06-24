#ifndef __VW_PLATE_MOD_PLATE_UTILS_H__
#define __VW_PLATE_MOD_PLATE_UTILS_H__

#include <vw/Core/Exception.h>
#include <http_protocol.h>
#include <boost/lexical_cast.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/function.hpp>
#include <set>
#include <map>
#include <string>

namespace vw {
namespace platefile {

class Index;

VW_DEFINE_EXCEPTION(PlateException,   Exception);
VW_DEFINE_EXCEPTION(BadRequest,       PlateException);
VW_DEFINE_EXCEPTION(ServerError,      PlateException);
VW_DEFINE_EXCEPTION(UnknownPlatefile, PlateException);

// Returns a copy of the url with % escapes unquoted, and + replaced by spaces.
std::string url_unquote(const std::string& str);

// Returns the const char* as a string, or an empty() string if it's scary in
// some way.
std::string safe_string_convert(const char* in);

class apache_output : public boost::iostreams::sink {
  public:
    apache_output(request_rec *r) : r(r) {}
    std::streamsize write(const char* s, std::streamsize n) { return ap_rwrite(s, n, r); }
  private:
    request_rec *r;
};

typedef boost::iostreams::stream<apache_output> apache_stream;

struct raii {
  typedef boost::function<void (void)> FuncT;
  FuncT m_leave;
  raii(FuncT enter, FuncT leave) : m_leave(leave) {enter();};
  ~raii() {m_leave();}
};

class QueryMap : private std::map<std::string,std::string> {
    // absolute limit on the number of keys parsed
    static const size_t MAX_KEYS = 10;
    static const size_t MAX_QUERY_LENGTH = 100;

    typedef std::map<std::string,std::string> map_t;

  public:
    QueryMap(const char *query);

    template <typename DataT>
    DataT get(const std::string& key, DataT defaultt) const {
      map_t::const_iterator i = this->find(key);
      if (i == this->end())
        return defaultt;

      try {
        return boost::lexical_cast<DataT>(i->second);
      } catch (const boost::bad_lexical_cast&) {
        vw_throw(BadRequest() << "Illegal query string value");
      }
    }
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


}} // namespace vw::platefile

#endif
