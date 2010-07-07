#ifndef __VW_PLATE_HTTPUTILS_H__
#define __VW_PLATE_HTTPUTILS_H__

#include <vw/Core/Exception.h>
#include <vw/Core/FundamentalTypes.h>
#include <boost/lexical_cast.hpp>
#include <string>
#include <map>

namespace vw {
namespace platefile {

// Returns a copy of the url with % escapes unquoted, and + replaced by spaces.
std::string url_unescape(const std::string& str);

// Returns a copy of the string, url-escaped.
std::string url_escape(const std::string& str, const std::string& safe = "");

class QueryMap {
    typedef std::map<std::string,std::string> map_t;
    map_t m_map;

  public:
    QueryMap() {}
    QueryMap(const std::string& query);

    // serialize querymap to query string (using proper escapes). Only adds
    // prefix if map is nonempty.
    std::string serialize(const std::string& prefix = "?", const std::string& sep = "&") const;

    template <typename DataT>
    DataT get(const std::string& key, DataT defaultt) const {
      map_t::const_iterator i = m_map.find(key);
      if (i == m_map.end())
        return defaultt;

      try {
        return boost::lexical_cast<DataT>(i->second);
      } catch (const boost::bad_lexical_cast&) {
        vw_throw(vw::ArgumentErr() << "Illegal query string value");
      }
    }
};

}} // namespace vw::platefile

#endif
