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

class Url {
  // Refuse to parse a url longer than this.
  static const unsigned MAX_URL_LENGTH = 500;

  std::string m_scheme, m_netloc, m_path, m_fragment;
  QueryMap m_query;

  void parse(const std::string& url, bool parse_query_params);

  public:
    // Parse a URL into 5 components:
    // <scheme>://<netloc>/<path>?<query>#<fragment>
    // pieces are url-decoded.
    // Every URL is parsed as though the scheme is HTTP.
    Url(const char* url, bool parse_query_params=true);
    Url(const std::string& url, bool parse_query_params=true);

    const std::string& scheme()   const { return m_scheme;   }
    const std::string& netloc()   const { return m_netloc;   }
    const std::string& path()     const { return m_path;     }
    const QueryMap& query()       const { return m_query;    }
    const std::string& fragment() const { return m_fragment; }

    std::string& scheme()   { return m_scheme;   }
    std::string& netloc()   { return m_netloc;   }
    std::string& path()     { return m_path;     }
    QueryMap& query()       { return m_query;    }
    std::string& fragment() { return m_fragment; }

    std::string hostname() const {
      size_t f = m_netloc.find(":");
      if (f == std::string::npos)
        return m_netloc;
      return m_netloc.substr(0,f);
    }

    uint16 port() const {
      size_t f = m_netloc.find(":");
      if (f == std::string::npos)
        return 0;
      return boost::lexical_cast<uint16>(m_netloc.substr(f+1));
    }

    std::string url() const {
      return scheme() + "://"
           + m_netloc
           + url_escape((m_path.empty() || m_path[0] != '/' ? "/" : "") + m_path, "/")
           + m_query.serialize()
           + url_escape((m_fragment.empty() ? "" : "#") + m_fragment);
    }
};

}} // namespace vw::platefile

#endif
