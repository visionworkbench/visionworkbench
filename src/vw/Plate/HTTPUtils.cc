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
  const char url_safe[] =
    "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    "abcdefghijklmnopqrstuvwxyz"
    "0123456789"
    "_.-";

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
