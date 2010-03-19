// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// Vision Workbench
#include <vw/Core/Stopwatch.h>
#include <vw/Plate/RemoteIndex.h>
#include <vw/Plate/PlateFile.h>

// Apache and the Apache Runtime
#include <httpd.h>
#include <http_config.h>
#include <http_protocol.h>
#include <ap_config.h>
#include <http_core.h>
#include <mod_status.h>
#include <http_log.h>

// STL
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

// Boost
#include <boost/iostreams/concepts.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/regex.hpp>
#include <boost/foreach.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/function.hpp>
#include <boost/algorithm/string/find_iterator.hpp>

#include <cstring>

#include <google/protobuf/service.h>

#include "mod_plate_io.h"
#include "common.h"

using namespace vw;
using namespace vw::platefile;

VW_DEFINE_EXCEPTION(PlateException, Exception);
VW_DEFINE_EXCEPTION(BadRequest,  PlateException);
VW_DEFINE_EXCEPTION(ServerError, PlateException);

typedef std::map<std::string, std::string> QueryMap;

static void null_closure() {}

template <typename MapT, typename DataType>
DataType mapget(MapT& m, const typename MapT::key_type& key, DataType def) {
  if (m.count(key) == 0)
    return def;
  else {
    try {
      return boost::lexical_cast<DataType>(m[key]);
    } catch (const boost::bad_lexical_cast&) {
      vw_throw(BadRequest() << "Illegal query string value");
    }
  }
}

class PlateModule {
  public:
    PlateModule(const server_rec *s);
    ~PlateModule();
    void connect_index();
    int operator()(request_rec *r) const;
    int status(request_rec *r, int flags) const;

    struct IndexCacheEntry {
      std::string shortname;
      std::string filename;
      std::string description;
      boost::shared_ptr<Index> index;
      int read_cursor;
    };

    struct BlobCacheEntry {
      boost::shared_ptr<Blob> blob;
      int platefile_id;
      BlobCacheEntry(boost::shared_ptr<Blob> b, int id) :
        blob(b), platefile_id(id) {}
    };

    typedef std::map<int32, IndexCacheEntry> IndexCache;

    const IndexCache& get_index() const { return index_cache; }
    const boost::shared_ptr<Blob> get_blob(int platefile_id, const std::string& plate_filename, uint32 blob_id) const;
    void sync_index_cache() const;

    std::ostream& logger(MessageLevel level, bool child_id = true) const {
      std::ostream& out = vw_out(level, "plate.apache");

      if (child_id)
        out << "Child[" << m_queue_name << "]: ";
      return out;
    }

    std::string get_dem() const {
      return m_conf->dem_id;
    }

  private:
    boost::shared_ptr<AmqpRpcClient> m_client;
    boost::shared_ptr<IndexService>  m_index_service;

    typedef std::map<std::string, BlobCacheEntry> BlobCache;

    mutable BlobCache  blob_cache;
    mutable IndexCache index_cache;
    bool m_connected;
    // We don't manage the data, and I think apache might modify it behind the scenes.
    // As such, mark it volatile.
    volatile const plate_config *m_conf;
    std::string m_queue_name;
};

namespace {
  boost::shared_ptr<PlateModule> mod_plate_ptr;
}

const PlateModule& mod_plate() {
  VW_ASSERT( mod_plate_ptr, LogicErr() << "Please initialize mod_plate_ptr first" );
  return *mod_plate_ptr;
}

PlateModule& mod_plate_mutable() {
  VW_ASSERT( mod_plate_ptr, LogicErr() << "Please initialize mod_plate_ptr first" );
  return *mod_plate_ptr;
}

std::string url_unquote(const std::string& str) {
  // This is sort of ugly... make sure shared_ptr deletes with free
  boost::shared_ptr<char> bytes(::strdup(str.c_str()), free);
  int ret = ap_unescape_url(bytes.get());
  if (ret != OK)
    vw_throw(BadRequest() << "Invalid query string");

  std::string unquoted(bytes.get());
  boost::replace_all(unquoted, "+", " ");

  return unquoted;
}

void query_to_map(QueryMap& keyval, const char* query) {
  if (!query)
    return;

  //logger(VerboseDebugMessage) << "Parsing query string: " << query << std::endl
  //                            << "Result:" << std::endl;

  std::vector<std::string> items;
  boost::split(items, query, boost::is_any_of(";&"));

  BOOST_FOREACH(const std::string& item, items) {
    if (item.empty())
      continue;
    size_t eq = item.find('=', 1);
    if (eq == std::string::npos) {
      keyval[url_unquote(item)] = "";
      continue;
    }
    keyval[url_unquote(item.substr(0, eq))] = url_unquote(item.substr(eq+1));
  }

  //BOOST_FOREACH( const QueryMap::value_type &i, keyval )
  //  logger(VerboseDebugMessage) << "\t" << i.first << "[" << i.second << "]" << std::endl;
}

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

//static int log_headers(void *null, const char *key, const char *value) {
//  logger(VerboseDebugMessage) << "\t" << key << "[" << value << "]" << std::endl;
//  return 1;
//}


// ---------------------------------------------------
//                 Content Handlers
// ---------------------------------------------------

int handle_image(request_rec *r, const std::string& url) {
  static const boost::regex match_regex("/(\\d+)/(\\d+)/(\\d+)/(\\d+)\\.(\\w+)$");

  boost::smatch match;
  if (!boost::regex_search(url, match, match_regex))
    return DECLINED;

  // We didn't decline. Connect!
  mod_plate_mutable().connect_index();

  QueryMap query;
  query_to_map(query, r->args);

  //logger(VerboseDebugMessage) << "Request Headers: " << std::endl;
  //apr_table_do(log_headers, 0, r->headers_in, NULL);

  int id    = boost::lexical_cast<int>(match[1]),
      level = boost::lexical_cast<int>(match[2]),
      col   = boost::lexical_cast<int>(match[3]),
      row   = boost::lexical_cast<int>(match[4]);
  std::string format = boost::lexical_cast<std::string>(match[5]);

  mod_plate().logger(DebugMessage) << "Request Image: id["  << id
                                   << "] level["  << level
                                   << "] col["    << col
                                   << "] row["    << row
                                   << "] format[" << format << "]" << std::endl;

  PlateModule::IndexCache::const_iterator index_i = mod_plate().get_index().find(id);

  if (index_i == mod_plate().get_index().end()) {
    // If we get an unknown platefile, resync just to make sure
    mod_plate().logger(WarningMessage) << "Platefile not in platefile cache. Resyncing." << std::endl;
    mod_plate().sync_index_cache();

    index_i = mod_plate().get_index().find(id);
    if (index_i == mod_plate().get_index().end())
      vw_throw(BadRequest() << "No such platefile [id = " << id << "]");
  }

  const PlateModule::IndexCacheEntry& index = index_i->second;

  // --------------  Access Plate Index -----------------

  IndexRecord idx_record;
  try {
    int transaction_id = mapget(query, "transaction_id", -1);
    bool exact = mapget(query, "exact", false);

    VW_ASSERT(transaction_id >= -1, BadRequest() << "Illegal transaction_id");

    if (transaction_id == -1) {
      transaction_id = index.index->transaction_cursor();
      exact = false;
    }

    mod_plate().logger(VerboseDebugMessage) << "Sending tile read_request with transaction[" << transaction_id << "] and exact[" << exact << "]" << std::endl;;
    idx_record = index.index->read_request(col,row,level,transaction_id,exact);
  } catch(const TileNotFoundErr &) {
      throw;
  } catch (const BadRequest &) {
    throw;
  } catch(const vw::Exception &e) {
    vw_throw(ServerError() << "Could not read plate index: " << e.what());
  }

  // ---------------- Return the image ------------------

  // Okay, we've gotten this far without error. Set content type now, so HTTP
  // HEAD returns the correct file type
  if (idx_record.filetype() == "png")
    r->content_type = "image/png";
  else if (idx_record.filetype() == "jpg")
    r->content_type = "image/jpg";
  else
    r->content_type = "application/octet-stream";

  if (mapget(query, "nocache", 0u) == 1) {
    apr_table_set(r->headers_out, "Cache-Control", "no-cache");
  }
  else {
    if (level <= 7)
      apr_table_set(r->headers_out, "Cache-Control", "max-age=604800");
    else
      apr_table_set(r->headers_out, "Cache-Control", "max-age=1200");
  }

  // This is as far as we can go without making the request heavyweight. Bail
  // out on a header request now.
  if (r->header_only)
    return OK;

  // These are the sendfile(2) parameters
  std::string filename;
  vw::uint64 offset, size;

  try {
    mod_plate().logger(VerboseDebugMessage) << "Fetching blob" << std::endl;
    // Grab a blob from the blob cache by filename
    boost::shared_ptr<Blob> blob = mod_plate().get_blob(id, index.filename, idx_record.blob_id());

    mod_plate().logger(VerboseDebugMessage) << "Fetching data from blob" << std::endl;
    // And calculate the sendfile(2) parameters
    blob->read_sendfile(idx_record.blob_offset(), filename, offset, size);

  } catch (vw::Exception &e) {
    vw_throw(ServerError() << "Could not load blob data: " << e.what());
  }

  apr_file_t *fd = 0;
  // Open the blob as an apache file with raii (so it goes away when we return)
  raii file_opener(
      boost::bind(apr_file_open, &fd, filename.c_str(), APR_READ|APR_FOPEN_SENDFILE_ENABLED, 0, r->pool),
      boost::bind(apr_file_close, boost::ref(fd)));

  // Use sendfile (if available) to send the proper tile data
  size_t sent;
  apr_status_t ap_ret;

  if ((ap_ret = ap_send_fd(fd, r, offset, size, &sent) != APR_SUCCESS)) {
    char buf[256];
    apr_strerror(ap_ret, buf, 256);
    vw_throw(ServerError() << "ap_send_fd failed: " << buf);
  }
  else if (sent != size)
    vw_throw(ServerError() << "ap_send_fd: short write (expected to send " << size << " bytes, but only sent " << sent);

  //logger(VerboseDebugMessage) << "Reply Headers: " << std::endl;
  //apr_table_do(log_headers, 0, r->headers_out, NULL);

  return OK;
}

class WTMLImageSet : public std::map<std::string, std::string> {
  typedef std::map<std::string, std::string> Map;
  typedef std::set<std::string> Child;
  Child child_keys;

  public:
    WTMLImageSet(const std::string& host,
                 const std::string& data_prefix,
                 const std::string& /*static_prefix*/,
                 const std::string dem_id,
                 const PlateModule::IndexCacheEntry& layer) {
      const IndexHeader& hdr = layer.index->index_header();

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

      (*this)["Name"]         = layer.description;
      (*this)["FileType"]     = std::string(".") + hdr.tile_filetype();
      (*this)["TileLevels"]   = boost::lexical_cast<std::string>(layer.index->num_levels());

      std::string data_url = host + data_prefix + vw::stringify(hdr.platefile_id());

      (*this)["Url"]          = data_url + "/{1}/{2}/{3}." + hdr.tile_filetype();
      (*this)["ThumbnailUrl"] = data_url + "/0/0/0."       + hdr.tile_filetype();
      // XXX: This is wrong for non-mars!
      //(*this)["DemUrl"]       = host + static_prefix + "megt128/{0}/{1}/{2}";

      (*this)["DemUrl"]       = host + data_prefix + dem_id + "/{0}/{1}/{2}.toast_dem_v1";


      child_keys.insert("ThumbnailUrl");
      child_keys.insert("Credits");
    }

    void serializeToOstream(std::ostream& o) {
      o << "<ImageSet";
      for (Map::const_iterator i = this->begin(), end = this->end(); i != end; ++i) {
        if (child_keys.count(i->first) == 0)
          o << " " << i->first << "='" << i->second << "'";
      }
      o << ">" << std::endl;

      BOOST_FOREACH( const std::string& key, child_keys ) {
        o << "\t<" << key << ">" << (*this)[key] << "</" << key << ">" << std::endl;
      }
      o << "</ImageSet>" << std::endl;
    }
};

int handle_wtml(request_rec *r, const std::string& url) {
  static const boost::regex match_regex("/(\\w+\\.wtml)$");

  boost::smatch match;
  if (!boost::regex_match(url, match, match_regex))
    return DECLINED;

  mod_plate_mutable().connect_index();

  std::string filename = boost::lexical_cast<std::string>(match[1]);

  r->content_type = "application/xml";

  if (r->header_only)
    return OK;

  mod_plate().sync_index_cache();

  apache_stream out(r);
  mod_plate().logger(DebugMessage) << "Served WTML[" << filename << "]" << std::endl;

  //WTMLImageSet lalt, mdim2, mdim1, mola;

  //lalt["Name"]       = "Lunar Topography (LALT)";
  //lalt["TileLevels"] = "7";

  //r->headers_out
  //APR_DECLARE(void) apr_table_set(apr_table_t *t, const char *key, const char *val);

  out
    << "<?xml version='1.0' encoding='UTF-8'?>"              << std::endl
    << "<Folder Name='Ames Planetary Content' Group='View'>" << std::endl << std::endl;

  typedef std::pair<int32, PlateModule::IndexCacheEntry> id_cache;

  std::ostringstream prefix;
  prefix << ap_http_scheme(r) << "://" << ap_get_server_name(r);

  unsigned port = ap_get_server_port(r);

  if (!ap_is_default_port(port, r))
      prefix << ":" << port;

  BOOST_FOREACH( const id_cache& e, mod_plate().get_index() ) {

    const std::string filetype = e.second.index->index_header().tile_filetype();
    // WWT can only handle jpg and png
    if (filetype != "jpg" && filetype != "png" && filetype != "auto") {
      vw_out(VerboseDebugMessage) << "Rejecting filetype " << filetype << " for WTML." << std::endl;
      continue;
    }

    WTMLImageSet img(prefix.str(), "/wwt/p/", "/static/", mod_plate().get_dem(), e.second);
    if (r->args) {
      img["Url"]          += std::string("?") + r->args;
      img["ThumbnailUrl"] += std::string("?") + r->args;
      img["DemUrl"]       += std::string("?") + r->args;
    }
    img.serializeToOstream(out);
  }
  out << "</Folder>" << std::endl;

  return OK;
}

PlateModule::PlateModule(const server_rec *s)
  : m_connected(false), m_conf(get_plate_config(s)),
    m_queue_name(AmqpRpcClient::UniqueQueueName("mod_plate"))
{

  // Disable the config file
  vw::vw_settings().set_rc_filename("");

  LogRuleSet rules;
  if (m_conf->rules->nelts == 0)
    rules.add_rule(DebugMessage, "plate.apache");
  else {
    ssize_t i;
    rule_entry *entry = reinterpret_cast<rule_entry*>(m_conf->rules->elts);
    for (i = 0; i < m_conf->rules->nelts; ++i) {
      rules.add_rule(entry->level, entry->name);
      entry++;
    }
  }

  // And log to stderr, which will go to the apache error log
  vw_log().set_console_stream(std::cerr, rules, true);

  logger(DebugMessage) << "child startup" << std::endl;
}

void PlateModule::connect_index() {

  if (m_connected)
      return;

  // Create the necessary services
  boost::shared_ptr<AmqpConnection> conn(new AmqpConnection(m_conf->rabbit_ip));

  std::string exchange = std::string(PLATE_EXCHANGE_NAMESPACE) + "." + m_conf->index_exchange;

  m_client.reset( new AmqpRpcClient(conn, exchange, m_queue_name, "index") );

  // Needs to respond in five seconds
  m_client->timeout(1000);
  m_client->tries(5);

  m_index_service.reset ( new IndexService::Stub(m_client.get() ) );
  m_client->bind_service(m_index_service, m_queue_name);

  m_connected = true;
  logger(DebugMessage) << "child connected to rabbitmq[" << m_conf->rabbit_ip
                       << "] exchange[" << exchange << "]" << std::endl;
}

PlateModule::~PlateModule() {}

int PlateModule::operator()(request_rec *r) const {

    if (!r->path_info)
        return DECLINED;

  if (!m_conf->dem_id) {
    logger(ErrorMessage) << "You forgot to set PlateDemID in the apache config. Refusing to proceed." << std::endl;
    return DECLINED;
  }

  std::string url(r->path_info);

  // WWT will append &new to dem urls... even when there's no ?.
  if (url.size() > 4 && url.compare(url.size()-4,4,"&new") == 0)
      url.erase(url.size()-4);

  typedef boost::function<int (request_rec*, const std::string& url)> Handler;
  static const Handler Handlers[] = {handle_image, handle_wtml};

  BOOST_FOREACH(const Handler h, Handlers) {
    int ret = h(r, url);
    if (ret != DECLINED)
      return ret;
  }
  return DECLINED;
}


int PlateModule::status(request_rec *r, int /*flags*/) const {
  apache_stream out(r);
  out << "Moo!" << std::endl;
  return OK;
}

const boost::shared_ptr<Blob> PlateModule::get_blob(int platefile_id, const std::string& plate_filename, uint32 blob_id) const {
  std::ostringstream ostr;
  ostr << plate_filename << "/plate_" << blob_id << ".blob";
  const std::string& filename = ostr.str();

  BlobCache::const_iterator blob = blob_cache.find(filename);

  // Check the platefile id to make sure the blob wasn't deleted and recreated
  // with a different platefile
  if (blob != blob_cache.end() && blob->second.platefile_id == platefile_id)
    return blob->second.blob;

  boost::shared_ptr<Blob> ret( new Blob(filename, true) );
  blob_cache.insert(std::make_pair(filename, BlobCacheEntry(ret, platefile_id)));
  return ret;
}

void PlateModule::sync_index_cache() const {

  VW_ASSERT(m_connected, LogicErr() << "Must connect before trying to sync cache");

  IndexListRequest request;
  IndexListReply id_list;

  index_cache.clear();

  m_index_service->ListRequest(m_client.get(), &request, &id_list, google::protobuf::NewCallback(&null_closure));

  BOOST_FOREACH( const std::string& name, id_list.platefile_names() ) {

    IndexCacheEntry entry;
    int32 id;

    try {
        std::string index_url = std::string("pf://") + m_conf->rabbit_ip + "/" + m_conf->index_exchange + "/" + name;
        logger(VerboseDebugMessage) << "Trying to load index: " << index_url << std::endl;

        entry.index = Index::construct_open(index_url);
        const IndexHeader& hdr = entry.index->index_header();

        entry.shortname   = name;
        entry.filename    = entry.index->platefile_name();
        entry.read_cursor = entry.index->transaction_cursor();
        entry.description = (hdr.has_description() && ! hdr.description().empty()) ? hdr.description() : entry.shortname + "." + vw::stringify(entry.read_cursor);
        id                = hdr.platefile_id();
    } catch (const vw::Exception& e) {
        logger(ErrorMessage) << "Tried to add " << name << " to the index cache, but failed: " << e.what() << std::endl;
        continue;
    }

    index_cache[id] = entry;
    logger(DebugMessage) << "Adding " << entry.shortname << " to index cache [cursor=" << entry.read_cursor << "]" << std::endl;
  }
}

// --------------------- Apache C++ Entry Points ------------------------

extern "C" int mod_plate_handler(request_rec *r) {
  try {
    return mod_plate()(r);
  } catch (const BadRequest& e) {
    // Client sent a request that was formatted badly
    ap_log_rerror(APLOG_MARK, APLOG_NOTICE, 0, r, "Bad Request: %s", e.what());
    return HTTP_BAD_REQUEST;
  } catch (const TileNotFoundErr& e) {
    // Valid format, but not there
    ap_log_rerror(APLOG_MARK, APLOG_NOTICE, 0, r,  "Tile not found: %s", e.what());
    return HTTP_NOT_FOUND;
  } catch (const ServerError& e) {
    // Something screwed up, but we controlled it
    ap_log_rerror(APLOG_MARK, APLOG_ERR, 0, r,  "Server Error [recovered]: %s", e.what());
    return HTTP_INTERNAL_SERVER_ERROR;
  } catch (const Exception& e) {
    // Something screwed up worse...
    ap_log_rerror(APLOG_MARK, APLOG_CRIT, 0, r, "Server Error [uncaught vw::Exception]: %s", e.what());
    return HTTP_INTERNAL_SERVER_ERROR;
  } catch (const std::exception &e) {
    // Something we don't understand broke. Eek.
    ap_log_rerror(APLOG_MARK, APLOG_ALERT, 0, r, "Server Error [uncaught std::exception]: %s", e.what());
    return HTTP_INTERNAL_SERVER_ERROR;
  }
}

extern "C" int mod_plate_status(request_rec *r, int flags) {
  if (flags & AP_STATUS_SHORT)
    return OK;
  return mod_plate().status(r, flags);
}

extern "C" void mod_plate_child_init(apr_pool_t * /*pchild*/, server_rec *s) {
  try {
    // make sure we create the handler object, and attach the server to it
    mod_plate_ptr.reset(new PlateModule(s));
    mod_plate();
  } catch (const Exception& e) {
    ap_log_error(APLOG_MARK, APLOG_ALERT, 0, s, "Could not start mod_plate child! [uncaught vw::Exception]: %s", e.what());
  } catch (const std::exception& e) {
    ap_log_error(APLOG_MARK, APLOG_ALERT, 0, s, "Could not start mod_plate child! [uncaught std::exception]: %s", e.what());
  } catch (...) {
    ap_log_error(APLOG_MARK, APLOG_ALERT, 0, s, "Could not start mod_plate child! [uncaught unknown exception]");
  }
}
