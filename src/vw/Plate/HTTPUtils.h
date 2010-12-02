// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_PLATE_HTTPUTILS_H__
#define __VW_PLATE_HTTPUTILS_H__

#include <vw/Core/Exception.h>
#include <vw/Core/FundamentalTypes.h>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/join.hpp>
#include <string>
#include <map>
#include <vector>

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

    // set a value, and return a reference to it for incremental use
    std::string& set(const std::string& key, const std::string& value);

    // is the value set?
    bool has(const std::string& key) const;

    // Dump the map
    void clear();

    template <typename DataT>
    DataT get(const std::string& key) const {
      map_t::const_iterator i = m_map.find(key);
      VW_ASSERT(i != m_map.end(), LogicErr() << "Key " << key << " is not in QueryMap!");
      try {
        return boost::lexical_cast<DataT>(i->second);
      } catch (const boost::bad_lexical_cast&) {
        vw_throw(vw::ArgumentErr() << "Illegal query string value");
      }
    }

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
  static const unsigned MAX_URL_LENGTH;

  std::string m_scheme, m_netloc, m_path, m_fragment;
  QueryMap m_query;

  void parse(const std::string& url, bool parse_query_params);
  bool complete() const;

  public:
    typedef std::vector<std::string> split_t;

    // Parse a URL into 5 components:
    // <scheme>://<netloc>/<path>?<query>#<fragment>
    // pieces are url-decoded.
    // Every URL is parsed as though the scheme is HTTP.
    Url();
    Url(const char* url, bool parse_query_params=true);
    Url(const std::string& url, bool parse_query_params=true);
    virtual ~Url() {}

    const std::string& scheme()   const { return m_scheme;   }
    const std::string& netloc()   const { return m_netloc;   }
    const std::string& path()     const { return m_path;     }
    const QueryMap& query()       const { return m_query;    }
    const std::string& fragment() const { return m_fragment; }

    // "/" -> [""]
    // "/pants" -> ["pants"]
    // "/pants/cheese" -> ["pants", "cheese"]
    split_t path_split() const;

    void scheme(const std::string& s);
    void netloc(const std::string& s);
    void fragment(const std::string& s);
    void path(const std::string& s);

    template <typename SeqSeqT>
    void path_join(const SeqSeqT& s) {
      VW_ASSERT(s.size() > 0, ArgumentErr() << "Cannot join an empty sequence");
      if (s.size() == 1 && s[0].empty())
        this->path("/");
      else
        this->path(boost::join(s, "/"));
    }

    QueryMap& query();

    std::string hostname() const;
    uint16 port() const;

    std::string string() const;

    friend std::istream& operator>>(std::istream& i, Url& val);
    friend std::ostream& operator<<(std::ostream& o, const Url& val);
};

class PlatefileUrl : public Url {
  private:
    void sanity() const;
  public:
    explicit PlatefileUrl(const char* url);
    explicit PlatefileUrl(const std::string& url);
    explicit PlatefileUrl(const Url& u);
    PlatefileUrl(const Url& u, const std::string& name);

    std::string name() const;
    Url base() const;

    void base(const Url& u);
    void name(const std::string& n);


    friend std::istream& operator>>(std::istream& i, PlatefileUrl& val);
};

}} // namespace vw::platefile

#endif
