#include <vw/Plate/HTTPUtils.h>
#include <vw/Core/Exception.h>
#include <vw/Core/FundamentalTypes.h>
#include <boost/algorithm/string.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/foreach.hpp>
#include <vector>

using std::string;
using std::vector;
using namespace vw;

namespace {
  string snip(const std::string& str, size_t begin, size_t end) {
    return string(str.begin()+begin, str.begin()+end);
  }
  const char url_safe[] =
    "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    "abcdefghijklmnopqrstuvwxyz"
    "0123456789"
    "_.-";

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
  boost::iterator_range<string::const_iterator> i = boost::find_token(url, !boost::is_any_of(scheme_safe));

  size_t c = bool(i) ? (i.begin() - url.begin()) : 0;
  size_t q = url.find('?', c == 0 ? 0 : c+1);
  size_t f = url.find('#', q == string::npos ? 0 : q+1);

  size_t begin = 0, end = url.size();

  if (c != 0) {
    m_scheme = boost::to_lower_copy(snip(url, begin, c));
    begin = c+1;
  }

  if (f != string::npos) {
    m_fragment = snip(url, f+1, end);
    end = f;
  }

  if (q != string::npos) {
    if (parse_query_params)
      m_query = QueryMap(snip(url, q+1, end));
    end = q;
  }

  if (url.substr(begin, 2) == "//")
    begin += 2;

  size_t s = url.find('/', begin);

  if (s == string::npos || s > end) {
    m_netloc = snip(url, begin, end);
    m_path   = "";
  } else {
    m_netloc = snip(url, begin, s);
    m_path   = snip(url, s, end);
  }
}

Url::Url(const char* url, bool parse_query_params)
  : m_scheme("file"), m_netloc(""), m_path("/"), m_fragment("") {
  VW_ASSERT(url, ArgumentErr() << "Cannot parse empty url");

  if (strlen(url) > MAX_URL_LENGTH)
    vw_throw(ArgumentErr() << "Refusing to parse a url longer than " << MAX_URL_LENGTH);

  parse(string(url), parse_query_params);
}

Url::Url(const string& url, bool parse_query_params)
  : m_scheme("file"), m_netloc(""), m_path("/"), m_fragment("") {
    parse(url, parse_query_params);
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

string vw::platefile::QueryMap::serialize(const std::string& prefix, const std::string& sep) const {
  if (m_map.empty())
    return std::string();

  std::vector<std::string> vals;
  vals.reserve(m_map.size());

  BOOST_FOREACH(const map_t::value_type& p, m_map)
    vals.push_back(url_escape(p.first) + "=" + url_escape(p.second));

  return prefix + boost::join(vals, sep);
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
