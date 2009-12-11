
// Vision Workbench
#include <vw/Core/Stopwatch.h>
#include <vw/Plate/RemoteIndex.h>
#include <vw/Plate/PlateFile.h>
#include <vw/Plate/common.h>

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

#include <google/protobuf/service.h>

#include "mod_plate_io.h"

using namespace vw;
using namespace vw::platefile;

VW_DEFINE_EXCEPTION(PlateException, Exception);
VW_DEFINE_EXCEPTION(BadRequest,  PlateException);
VW_DEFINE_EXCEPTION(ServerError, PlateException);

static void null_closure() {}

class PlateModule {
  public:
    PlateModule();
    ~PlateModule();
    int operator()(request_rec *r) const;
    int status(request_rec *r, int flags) const;

    struct IndexCacheEntry {
      std::string shortname;
      std::string filename;
      std::string description;
      boost::shared_ptr<Index> index;
      int read_cursor;
    };

    typedef std::map<int32, IndexCacheEntry> IndexCache;

    const IndexCache& get_index() const { return index_cache; }
    const boost::shared_ptr<Blob> get_blob(const std::string& plate_filename, uint32 blob_id) const;
    void sync_index_cache() const;

  private:
    boost::shared_ptr<AmqpRpcClient> m_client;
    boost::shared_ptr<IndexService>  m_index_service;

    typedef std::map<std::string, boost::shared_ptr<Blob> > BlobCache;

    mutable BlobCache  blob_cache;
    mutable IndexCache index_cache;
};

namespace {
  vw::RunOnce mod_plate_once = VW_RUNONCE_INIT;
  boost::shared_ptr<PlateModule> mod_plate_ptr;
  void init_mod_plate() {
    mod_plate_ptr = boost::shared_ptr<PlateModule>(new PlateModule());
  }
  void kill_mod_plate() {
    // This is safe because it's a shared pointer.
    init_mod_plate();
  }
}

/// Static method to access the singleton instance of the plate module object.
PlateModule& mod_plate() {
  mod_plate_once.run( init_mod_plate );
  return *mod_plate_ptr;
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

// XXX: This is so wrong, and makes me feel dirty. Why doesn't apache have a arg parser?
std::string get_nocache_value(const request_rec* r) {
  boost::smatch match;
  if (r->args) {
    boost::regex nocache_regex("^(.*&)?nocache=(\\d)(&.*)?$");
    if (boost::regex_search(std::string(r->args), match, nocache_regex))
      return match[2];
  }
  return "";
}


// ---------------------------------------------------
//                 Content Handlers
// ---------------------------------------------------

int handle_image(request_rec *r, const std::string& url) {
  static const boost::regex match_regex("/(\\d+)/(\\d+)/(\\d+)/(\\d+)\\.(\\w+)$");

  boost::smatch match;
  if (!boost::regex_search(url, match, match_regex))
    return DECLINED;

  int id    = boost::lexical_cast<int>(match[1]),
      level = boost::lexical_cast<int>(match[2]),
      col   = boost::lexical_cast<int>(match[3]),
      row   = boost::lexical_cast<int>(match[4]);
  std::string format = boost::lexical_cast<std::string>(match[5]);

  vw_out(DebugMessage, "plate.apache") << "Request Image: id["  << id
                                                 << "] level["  << level
                                                 << "] col["    << col
                                                 << "] row["    << row
                                                 << "] format[" << format << "]" << std::endl;

  PlateModule::IndexCache::const_iterator index_i = mod_plate().get_index().find(id);
  if (index_i == mod_plate().get_index().end())
    vw_throw(BadRequest() << "No such platefile [id = " << id << "]");

  const PlateModule::IndexCacheEntry& index = index_i->second;

  // --------------  Access Plate Index -----------------

  IndexRecord idx_record;
  try {
      /// XXX: This doubles the number of amqp messages
    idx_record = index.index->read_request(col,row,level,index.index->transaction_cursor());
  } catch(const TileNotFoundErr &) {
      throw;
  } catch(const vw::Exception &e) {
    vw_throw(ServerError() << "Could not read plate index: " << e.what());
  }

  // ---------------- Return the image ------------------

  // Okay, we've gotten this far without error. Set content type now, so HTTP
  // HEAD returns the correct file type
  r->content_type = "image/png";

  std::string nocache = get_nocache_value(r);
  if (nocache.empty()) {
    if (level <= 7)
      apr_table_set(r->headers_out, "Cache-Control", "max-age=604800");
    else
      apr_table_set(r->headers_out, "Cache-Control", "max-age=1200");
  }
  else
    apr_table_set(r->headers_out, "Cache-Control", "no-cache");

  // This is as far as we can go without making the request heavyweight. Bail
  // out on a header request now.
  if (r->header_only)
    return OK;

  // These are the sendfile(2) parameters
  std::string filename;
  vw::uint64 offset, size;

  if (idx_record.status() != INDEX_RECORD_VALID && idx_record.status() != INDEX_RECORD_STALE) {
      vw_out(DebugMessage, "plate.apache") << "Returning null tile [record=" << idx_record.DebugString() << "]" << std::endl;

      filename = "/big/platefiles/null.png";
      offset = 0;
      apr_finfo_t info;
      apr_stat(&info, filename.c_str(), APR_FINFO_SIZE, r->pool);
      size = info.size;
  } else {

    try {
      // Grab a blob from the blob cache by filename
      boost::shared_ptr<Blob> blob = mod_plate().get_blob(index.filename, idx_record.blob_id());
      // And calculate the sendfile(2) parameters
      blob->read_sendfile(idx_record.blob_offset(), filename, offset, size);

    } catch (vw::Exception &e) {
      vw_throw(ServerError() << "Could not load blob data: " << e.what());
    }
  }

  apr_file_t *fd = 0;
  // Open the blob as an apache file with raii (so it goes away when we return)
  raii file_opener(
      boost::bind(apr_file_open, &fd, filename.c_str(), APR_READ|APR_FOPEN_SENDFILE_ENABLED, 0, r->pool),
      boost::bind(apr_file_close, boost::ref(fd)));

  // Use sendfile (if available) to send the proper tile data
  size_t sent;
  if (ap_send_fd(fd, r, offset, size, &sent) != APR_SUCCESS || sent != size)
    vw_throw(ServerError() << "Did not send tile correctly!");

  return OK;

}

class WTMLImageSet : public std::map<std::string, std::string> {
  typedef std::map<std::string, std::string> Map;
  typedef std::set<std::string> Child;
  Child child_keys;

  public:
    WTMLImageSet(const std::string& url_prefix, const PlateModule::IndexCacheEntry& layer) {
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
      (*this)["ElevationModel"]     = "False";
      (*this)["StockSet"]           = "False";
      // XXX: This is wrong for non-mars!
      (*this)["DemUrl"]             = "http://www.worldwidetelescope.ccom/wwtweb/marsdem.aaspx?q={0},{1},{2},T";

      (*this)["Name"]         = layer.description;
      (*this)["FileType"]     = std::string(".") + hdr.tile_filetype();
      (*this)["TileLevels"]   = boost::lexical_cast<std::string>(layer.index->max_depth());

      const std::string url2 = url_prefix + "p/" + vw::stringify(hdr.platefile_id());

      (*this)["Url"]          = url2 + "/{1}/{2}/{3}." + hdr.tile_filetype();
      (*this)["ThumbnailUrl"] = url2 + "/0/0/0."       + hdr.tile_filetype();

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

  std::string filename = boost::lexical_cast<std::string>(match[1]);

  r->content_type = "application/xml";

  if (r->header_only)
    return OK;

  mod_plate().sync_index_cache();

  apache_stream out(r);
  vw_out(DebugMessage, "plate.apache") << "Served WTML[" << filename << "]" << std::endl;

  //WTMLImageSet lalt, mdim2, mdim1, mola;

  //lalt["Name"]       = "Lunar Topography (LALT)";
  //lalt["TileLevels"] = "7";

  //r->headers_out
  //APR_DECLARE(void) apr_table_set(apr_table_t *t, const char *key, const char *val);

  std::string nocache = get_nocache_value(r);

  out
    << "<?xml version='1.0' encoding='UTF-8'?>"              << std::endl
    << "<Folder Name='Ames Planetary Content' Group='View'>" << std::endl << std::endl;

  typedef std::pair<int32, PlateModule::IndexCacheEntry> id_cache;

  std::ostringstream prefix;
  prefix << ap_http_scheme(r) << "://" << r->server->server_hostname;
  if (r->server->port != 80)
      prefix << ":" << r->server->port;
  prefix << "/wwt/";

  BOOST_FOREACH( const id_cache& e, mod_plate().get_index() ) {
    WTMLImageSet img(prefix.str(), e.second);
    if (!nocache.empty()) {
      img["Url"]          += "?nocache=" + nocache;
      img["ThumbnailUrl"] += "?nocache=" + nocache;
    }
    img.serializeToOstream(out);
  }
  out << "</Folder>" << std::endl;

  return OK;
}

PlateModule::PlateModule() {

  // Disable the config file
  vw::vw_settings().set_rc_filename("");

  LogRuleSet rules;
  rules.add_rule(EveryMessage, "plate.apache");

  // And log to stderr, which will go to the apache error log
  vw_log().set_console_stream(std::cerr, rules, false);

  // Create the necessary services
  std::string queue_name = AmqpRpcClient::UniqueQueueName("mod_plate");

  // XXX: The rabbitmq host needs to be an apache configuration variable or something
  boost::shared_ptr<AmqpConnection> conn(new AmqpConnection("198.10.124.5"));

  AmqpRpcEndpoint::set_default_timeout(1000);

  m_client.reset( new AmqpRpcClient(conn, INDEX_EXCHANGE, queue_name, "index") );
  m_index_service.reset ( new IndexService::Stub(m_client.get() ) );
  m_client->bind_service(m_index_service, queue_name);

  sync_index_cache();

  vw_out(DebugMessage, "plate.apache") << "child startup" << std::endl;
}

PlateModule::~PlateModule() {}

int PlateModule::operator()(request_rec *r) const {

    if (!r->path_info)
        return DECLINED;

  const std::string url(r->path_info);

  typedef boost::function<int (request_rec*, const std::string& url)> Handler;
  static const Handler Handlers[] = {handle_image, handle_wtml};

  BOOST_FOREACH(const Handler h, Handlers) {
    int ret = h(r, url);
    if (ret != DECLINED)
      return ret;
  }
  return DECLINED;
}


int PlateModule::status(request_rec *r, int flags) const {
  apache_stream out(r);
  out << "Moo!" << std::endl;
  return OK;
}

const boost::shared_ptr<Blob> PlateModule::get_blob(const std::string& plate_filename, uint32 blob_id) const {
  std::ostringstream ostr;
  ostr << plate_filename << "/plate_" << blob_id << ".blob";
  const std::string& filename = ostr.str();

  BlobCache::const_iterator blob = blob_cache.find(filename);
  if (blob != blob_cache.end())
    return blob->second;

  boost::shared_ptr<Blob> ret( new Blob(filename, true) );
  blob_cache[filename] = ret;
  return ret;
}

void PlateModule::sync_index_cache() const {
  IndexListRequest request;
  IndexListReply id_list;

  index_cache.clear();

  m_index_service->ListRequest(m_client.get(), &request, &id_list, google::protobuf::NewCallback(&null_closure));

  BOOST_FOREACH( const std::string& name, id_list.platefile_names() ) {

    IndexCacheEntry entry;
    int32 id;

    try {
        // XXX: The rabbitmq host needs to be an apache configuration variable or something
        entry.index = Index::construct_open(std::string("pf://198.10.124.5/index/") + name);
        const IndexHeader& hdr = entry.index->index_header();

        entry.shortname   = name;
        entry.filename    = entry.index->platefile_name();
        entry.read_cursor = entry.index->transaction_cursor();
        entry.description = (hdr.has_description() && ! hdr.description().empty()) ? hdr.description() : entry.shortname + "." + vw::stringify(entry.read_cursor);
        id                = hdr.platefile_id();
    } catch (const vw::Exception& e) {
        vw_out(ErrorMessage, "plate.apache") << "Tried to add " << name << " to the index cache, but failed: " << e.what() << std::endl;
        continue;
    }

    index_cache[id] = entry;
    vw_out(DebugMessage, "plate.apache") << "Adding " << entry.shortname << " to index cache [cursor=" << entry.read_cursor << "]" << std::endl;
  }
}

// --------------------- Apache C++ Entry Points ------------------------

void mod_plate_destroy() {
  kill_mod_plate();
}

extern "C" int mod_plate_handler(request_rec *r) {
  try {
    return mod_plate()(r);
  } catch (const BadRequest& e) {
    // Client sent a request that was formatted badly
    ap_log_rerror(APLOG_MARK, APLOG_NOTICE, 0, r, "Bad Request: %s", e.what());
    return HTTP_BAD_REQUEST;
  } catch (const TileNotFoundErr& e) {
    // Valid format, but not there
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

extern "C" void mod_plate_child_init(apr_pool_t *pchild, server_rec *s) {
  try {
    static_cast<void>(mod_plate());
  } catch (const Exception& e) {
    ap_log_error(APLOG_MARK, APLOG_ALERT, 0, s, "Could not start mod_plate child! [uncaught vw::Exception]: %s", e.what());
  } catch (const std::exception& e) {
    ap_log_error(APLOG_MARK, APLOG_ALERT, 0, s, "Could not start mod_plate child! [uncaught std::exception]: %s", e.what());
  } catch (...) {
    ap_log_error(APLOG_MARK, APLOG_ALERT, 0, s, "Could not start mod_plate child! [uncaught unknown exception]");
  }
}
