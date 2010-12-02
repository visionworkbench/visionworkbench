// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Plate/HTTPUtils.h>
#include <vw/Core/Exception.h>
#include <vw/Core/FundamentalTypes.h>
#include <boost/algorithm/string.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/foreach.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/convenience.hpp>
#include <vector>

using std::string;
using std::vector;
using namespace vw;

namespace fs = boost::filesystem;

namespace {
  string snip(const std::string& str, size_t begin, size_t end) {
    return string(str.begin()+begin, str.begin()+end);
  }

  // conditionally safe
  // +&/:;=?@

  const char url_safe[] =
    "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    "abcdefghijklmnopqrstuvwxyz"
    "0123456789"
    "$-_.!*'(),";

  const char scheme_safe[] =
    "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    "abcdefghijklmnopqrstuvwxyz"
    "0123456789"
    "+-.";

  string to_hex_formatter(const boost::iterator_range<std::string::const_iterator>& match) {
    string tmp(match.begin(), match.end());
    if (tmp == " ")
      return "+";
    else {
      std::ostringstream o;
      o << '%' << std::hex << static_cast<uint32>(tmp[0]);
      return o.str();
    }
  }
}

namespace vw {
namespace platefile {

const unsigned Url::MAX_URL_LENGTH = 500;

void Url::parse(const std::string& url, bool parse_query_params) {
  // iterator to the first char that's not part of the scheme
  // That should either be the ":", or end()
  size_t c = url.find("://", 0);
  // iterator to the beginning of the query string (the ?) or end()
  size_t q = url.find('?', c == string::npos ? 0 : c+1);
  // iterator to the beginning of the fragment (the #) or end()
  size_t f = url.find('#',
      q != string::npos ? q+1 :
      c != string::npos ? c+1 : 0);

  // front and back iterators
  size_t begin = 0, end = url.size();

  // did we find a scheme?
  if (c != string::npos) {
    // yes. replace the default "file" with the scheme, and move it up
    m_scheme = boost::to_lower_copy(snip(url, begin, c));
    begin = c+1;
  }

  // did we find a fragment?
  if (f != string::npos) {
    // yes. snip off the fragment and move the end forward
    m_fragment = snip(url, f+1, end);
    end = f;
  }

  // did we find a query string?
  if (q != string::npos) {
    // yes. snip it off and move the end forward
    if (parse_query_params)
      m_query = QueryMap(snip(url, q+1, end));
    end = q;
  }

  // if we had a scheme, snip off the "//" portion of scheme://
  if (url.substr(begin, 2) == "//")
    begin += 2;

  if (c != string::npos) {
    // if we found a scheme
    // [begin, end) now contains just the netloc and the path
    size_t s = url.find('/', begin);

    if (s == string::npos || s > end) {
      m_netloc = snip(url, begin, end);
      m_path   = "";
    } else {
      m_netloc = snip(url, begin, s);
      m_path   = snip(url, s, end);
    }
  } else {
    // if we didn't have a scheme, it's all just path now.
    m_scheme = "file";
    m_netloc = "";
    m_path = snip(url, begin, end);
  }

  if (m_path.empty())
    m_path = "/";
}

Url::Url()
  : m_scheme(""), m_netloc(""), m_path("/"), m_fragment("") {}

Url::Url(const char* url, bool parse_query_params)
  : m_scheme("file"), m_netloc(""), m_path("/"), m_fragment("") {
  VW_ASSERT(url, ArgumentErr() << "Cannot parse empty url");

  size_t len = strlen(url);

  if (len < 1)
    vw_throw(ArgumentErr() << "Refusing to parse an empty url");

  if (len > MAX_URL_LENGTH)
    vw_throw(ArgumentErr() << "Refusing to parse a url longer than " << MAX_URL_LENGTH);

  parse(std::string(url), parse_query_params);
}

Url::Url(const std::string& url, bool parse_query_params)
  : m_scheme("file"), m_netloc(""), m_path("/"), m_fragment("") {
    parse(url, parse_query_params);
}
void Url::scheme(const std::string& s)   { m_scheme = s; }
void Url::netloc(const std::string& s)   { m_netloc = s; }
void Url::fragment(const std::string& s) { m_fragment = s; }
QueryMap& Url::query() { return m_query; }

void Url::path(const std::string& s) {
  VW_ASSERT(!s.empty(), ArgumentErr() << "Invalid path (must be non-empty)");
  VW_ASSERT(m_scheme == "file" || s[0] == '/', ArgumentErr() << "Only file urls can be relative");
  m_path = s;
}

Url::split_t Url::path_split() const {
  VW_ASSERT(m_path.size() > 0, LogicErr() << "Empty Path? Should be impossible.");
  split_t items;

  if (m_path == "/") {
    items.push_back("");
    return items;
  }

  boost::iterator_range<std::string::const_iterator> r(m_path.begin(), m_path.end());
  boost::split(items, r, boost::is_any_of("/"), boost::algorithm::token_compress_on);
  return items;
}

bool Url::complete() const {
  return !(m_scheme.empty() || m_path.empty());
}

std::string Url::hostname() const {
  size_t f = m_netloc.find(":");
  if (f == std::string::npos)
    return m_netloc;
  return m_netloc.substr(0,f);
}

uint16 Url::port() const {
  size_t f = m_netloc.find(":");
  if (f == std::string::npos)
    return 0;
  return boost::lexical_cast<uint16>(m_netloc.substr(f+1));
}

std::string Url::string() const {
  VW_ASSERT(complete(), LogicErr() << "Cannot ask for url string before scheme and path are populated");
  return scheme() + "://"
    + netloc()
    + url_escape(m_path, "/")
    + m_query.serialize()
    + (m_fragment.empty() ? "" : "#") + url_escape(m_fragment, "/=;");
}

std::istream& operator>>(std::istream& i, Url& val) {
  std::string s;
  i >> s;
  val = Url(s);
  return i;
}

std::ostream& operator<<(std::ostream& o, const Url& val) {
  return (o << val.string());
}

void PlatefileUrl::sanity() const {
  Url::split_t sp = path_split();
  VW_ASSERT(sp.size() > 0, ArgumentErr() << "Expected a platefile url (bad path)");
  VW_ASSERT(boost::ends_with(sp.back(), ".plate"), ArgumentErr() << "Expected a platefile url (doesn't end in .plate)");
}

PlatefileUrl::PlatefileUrl(const char* url) : Url(url) {
  sanity();
}

PlatefileUrl::PlatefileUrl(const std::string& url) : Url(url) {
  sanity();
}

PlatefileUrl::PlatefileUrl(const Url& url)
  : Url(url) {
  sanity();
}

PlatefileUrl::PlatefileUrl(const Url& u, const std::string& name)
  : Url(u)
{
  Url::split_t path(this->path_split());
  path.push_back(name);
  this->path_join(path);
  sanity();
}

Url PlatefileUrl::base() const {
  Url::split_t path(this->path_split());
  path.pop_back();
  Url u(*this);
  u.path_join(path);
  return u;
}

void PlatefileUrl::base(const Url& u) {
  PlatefileUrl url(u, this->name());
  *this = url;
}

void PlatefileUrl::name(const std::string& n) {
  VW_ASSERT(boost::ends_with(n, ".plate"), ArgumentErr() << "Expected a platefile name (doesn't end in .plate)");
  Url::split_t path(this->path_split());
  path.back() = n;
  path_join(path);
}

std::string PlatefileUrl::name() const {
  Url::split_t path(this->path_split());
  return path.back();
}

std::istream& operator>>(std::istream& i, PlatefileUrl& val) {
  std::string s;
  i >> s;
  val = PlatefileUrl(s);
  return i;
}

QueryMap::QueryMap(const std::string& query) {
  std::vector<string> items;
  boost::split(items, query, boost::is_any_of(";&"));

  BOOST_FOREACH(const string& item, items) {
    if (item.empty())
      continue;
    size_t eq = item.find('=', 1);
    if (eq == string::npos) {
      m_map[url_unescape(item)] = "";
      continue;
    }
    m_map[url_unescape(item.substr(0, eq))] = url_unescape(item.substr(eq+1));
  }
}

std::string& QueryMap::set(const std::string& key, const std::string& value) {
  return (m_map[key] = value);
}

string vw::platefile::QueryMap::serialize(const std::string& prefix, const std::string& sep) const {
  if (m_map.empty())
    return std::string();

  std::vector<std::string> vals;
  vals.reserve(m_map.size());

  BOOST_FOREACH(const map_t::value_type& p, m_map)
    vals.push_back(url_escape(p.first) + "=" + url_escape(p.second));

  return prefix + boost::join(vals, sep);
}

bool QueryMap::has(const std::string& key) const {
  return m_map.find(key) != m_map.end();
}

void QueryMap::clear() {
  m_map.clear();
}

// This is an adapter so boost::lexical_cast can handle hex.
template <typename ElemT>
struct HexTo {
    ElemT value;
    operator ElemT() const {return value;}
    friend std::istream& operator>>(std::istream& in, HexTo& out) {
      std::ios_base::fmtflags f = in.setf(std::ios_base::hex, std::ios_base::basefield);
      in >> out.value;
      in.setf(f,std::ios_base::basefield);
      return in;
    }
};

unsigned char from_hex(const std::string& str) {
  return boost::numeric_cast<unsigned char>(
      static_cast<uint32_t>(
        boost::lexical_cast<HexTo<uint32_t> >(str)));
}

string url_unescape(const string& str) {
  vector<string> items;
  std::string str2 = boost::replace_all_copy(str, "+", " ");
  boost::split(items, str2, boost::is_any_of("%"));

  for (size_t i = 1; i < items.size(); ++i) {
    unsigned char decoded = from_hex(items[i].substr(0, 2));
    if (decoded == 0)
      items[i] = string(1, '%') + items[i];
    else
      items[i] = string(1, decoded) + items[i].substr(2);
  }
  return boost::join(items, "");
}

string url_escape(const string& str, const string& safe) {
  return boost::find_format_all_copy(
      str,
      boost::token_finder(!boost::is_any_of(url_safe+safe)),
      to_hex_formatter);
}

}} // namespace vw::platefile
