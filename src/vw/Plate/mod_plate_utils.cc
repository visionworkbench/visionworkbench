#include "mod_plate_utils.h"
#include <vw/Plate/Index.h>

#include <httpd.h>
#include <boost/shared_ptr.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <vector>
#include <map>

using std::string;

string vw::platefile::url_unquote(const string& str) {
  // This is sort of ugly... make sure shared_ptr deletes with free
  boost::shared_ptr<char> bytes(::strdup(str.c_str()), free);
  int ret = ap_unescape_url(bytes.get());
  if (ret != OK)
    vw_throw(BadRequest() << "Invalid query string");

  string unquoted(bytes.get());
  boost::replace_all(unquoted, "+", " ");

  return unquoted;
}

vw::platefile::QueryMap::QueryMap(const char *query) {
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
      (*this)[url_unquote(item)] = "";
      continue;
    }
    (*this)[url_unquote(item.substr(0, eq))] = url_unquote(item.substr(eq+1));
  }
}

string vw::platefile::safe_string_convert(const char* in) {
  static const size_t MAX_STRING_LENGTH = 255;

  // string blows up if it's handed null, and could allocate a lot of
  // memory if you let it.
  if (!in || strlen(in) > MAX_STRING_LENGTH)
    return string();
  return string(in);
}

vw::platefile::WTMLImageSet::WTMLImageSet(
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
  //(*this)["DemUrl"]       = host + static_prefix + "megt128/{0}/{1}/{2}";

  (*this)["DemUrl"]       = host + data_prefix + dem_id + "/{0}/{1}/{2}.toast_dem_v1";


  child_keys.insert("ThumbnailUrl");
  child_keys.insert("Credits");
}

const string& mapget(const std::map<string,string>& map, const std::string& key) {
  std::map<string,string>::const_iterator i = map.find(key);
  if (i == map.end())
    vw::vw_throw(vw::LogicErr() << "WTMLImageSet set up incorrectly");
  return i->second;
}

void vw::platefile::WTMLImageSet::serializeToOstream(std::ostream& o) const {
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
