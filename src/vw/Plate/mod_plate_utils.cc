#include "mod_plate_utils.h"
#include <vw/Plate/Index.h>

#include <httpd.h>
#include <http_protocol.h>
#include <boost/shared_ptr.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <vector>
#include <map>

using std::string;

namespace vw {
namespace platefile {

string url_unescape(const string& str) {
  // This is sort of ugly... make sure shared_ptr deletes with free
  boost::shared_ptr<char> bytes(::strdup(str.c_str()), free);
  int ret = ap_unescape_url(bytes.get());
  if (ret != OK)
    vw_throw(BadRequest() << "Invalid query string");

  string unquoted(bytes.get());
  boost::replace_all(unquoted, "+", " ");

  return unquoted;
}

string url_escape(const string& str, apr_pool_t *pool) {
  // this char* will be freed by apache (allocated on request pool), so don't
  // try to free it here.
  char *uri_s = ap_escape_uri(pool, str.c_str());
  return std::string(uri_s);
}

QueryMap::QueryMap(const char *query, apr_pool_t *pool) : m_pool(pool) {
  if (!query || strlen(query) > MAX_QUERY_LENGTH)
    return;

  std::vector<string> items;
  boost::split(items, query, boost::is_any_of(";&"));

  size_t found = 0;
  BOOST_FOREACH(const string& item, items) {
    if (item.empty())
      continue;
    if (found++ >= MAX_KEYS)
      break;
    size_t eq = item.find('=', 1);
    if (eq == string::npos) {
      m_map[url_unescape(item)] = "";
      continue;
    }
    m_map[url_unescape(item.substr(0, eq))] = url_unescape(item.substr(eq+1));
  }
}

string QueryMap::serialize(const std::string& prefix, const std::string& sep) const {
  if (m_map.empty())
    return std::string();

  std::vector<std::string> vals;
  vals.reserve(m_map.size());

  BOOST_FOREACH(const map_t::value_type& p, m_map)
    vals.push_back(url_escape(p.first, m_pool) + "=" + url_escape(p.second, m_pool));

  return prefix + boost::join(vals, sep);
}

string safe_string_convert(const char* in) {
  static const size_t MAX_STRING_LENGTH = 255;

  // string blows up if it's handed null, and could allocate a lot of
  // memory if you let it.
  if (!in || strlen(in) > MAX_STRING_LENGTH)
    return string();
  return string(in);
}

std::streamsize apache_output::write(const char* s, std::streamsize n) {
  return ap_rwrite(s, n, r);
}

WTMLImageSet::WTMLImageSet(
    const string& host,
    const string& data_prefix,
    const string& /*static_prefix*/,
    const string dem_id,
    boost::shared_ptr<const Index> index,
    const string& description)
{
  const IndexHeader& hdr = index->index_header();

  (*this)["Generic"]            = "False";
  (*this)["DataSetType"]        = "Planet";
  (*this)["BandPass"]           = "Visible";
  (*this)["BaseTileLevel"]      = "0";
  (*this)["BaseDegreesPerTile"] = "360";
  (*this)["BottomsUp"]          = "False";
  (*this)["Projection"]         = "Toast";
  (*this)["QuadTreeMap"]        = "0123";
  (*this)["CenterX"]            = "0";
  (*this)["CenterY"]            = "0";
  (*this)["OffsetX"]            = "0";
  (*this)["OffsetY"]            = "0";
  (*this)["Rotation"]           = "0";
  (*this)["Sparse"]             = "True";
  (*this)["ElevationModel"]     = "True";
  (*this)["MeanRadius"]         = "3396000";
  (*this)["StockSet"]           = "False";

  (*this)["Name"]         = description;
  (*this)["FileType"]     = string(".") + hdr.tile_filetype();

  // This is the max index, not the number of levels.
  (*this)["TileLevels"]   = boost::lexical_cast<string>(index->num_levels()-1);

  string data_url = host + data_prefix + vw::stringify(hdr.platefile_id());

  (*this)["Url"]          = data_url + "/{1}/{2}/{3}." + hdr.tile_filetype();
  (*this)["ThumbnailUrl"] = data_url + "/0/0/0."       + hdr.tile_filetype();

  // XXX: This is wrong for non-mars!
  (*this)["DemUrl"]       = host + data_prefix + dem_id + "/{0}/{1}/{2}.toast_dem_v1";

  // This probably needs to be in the plate somewhere...
  (*this)["Credits"]      = "";

  child_keys.insert("ThumbnailUrl");
  child_keys.insert("Credits");
}

const string& mapget(const std::map<string,string>& map, const std::string& key) {
  std::map<string,string>::const_iterator i = map.find(key);
  if (i == map.end())
    vw::vw_throw(vw::LogicErr() << "WTMLImageSet set up incorrectly (could not find required key " << key << ")");
  return i->second;
}

void WTMLImageSet::serializeToOstream(std::ostream& o) const {
  o << "<ImageSet";
  for (map_t::const_iterator i = this->begin(), end = this->end(); i != end; ++i) {
    if (child_keys.count(i->first) == 0)
      o << " " << i->first << "='" << i->second << "'";
  }
  o << ">" << std::endl;

  BOOST_FOREACH( const string& key, child_keys ) {
    o << "\t<" << key << "><![CDATA[" << mapget(*this, key) << "]]></" << key << ">" << std::endl;
  }
  o << "</ImageSet>" << std::endl;
}

ApacheRequest::ApacheRequest(request_rec* r)
  : r(r), url(safe_string_convert(r->path_info)), args(QueryMap(r->args, r->pool)) {

  // WWT will append &new to dem urls... even when there's no ?.
  if (url.size() > 4 && url.compare(url.size()-4,4,"&new") == 0)
    url.erase(url.size()-4);
}

bool ApacheRequest::header_only() const {
  return r->header_only;
}

request_rec* ApacheRequest::writer() const {
  return r;
}

}} // namespace vw::platefile
