
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
VW_DEFINE_EXCEPTION(NoSuchTile,  PlateException);
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
    };

    const IndexCacheEntry get_index(int32 platefile_id) const;
    const boost::shared_ptr<Blob> get_blob(const std::string& plate_filename, uint32 blob_id) const;

  private:
    boost::shared_ptr<AmqpRpcClient>       m_client;
    boost::shared_ptr<IndexService>        m_index_service;
    boost::shared_ptr<IndexManagerService> m_index_mgr_service;

    typedef std::map<std::string, boost::shared_ptr<Blob> > BlobCache;
    typedef std::map<int32, IndexCacheEntry> IndexCache;

    mutable BlobCache  blob_cache;
    mutable IndexCache index_cache;
    void sync_index_cache();
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

// ---------------------------------------------------
//                 Content Handlers
// ---------------------------------------------------

int handle_image(request_rec *r, const std::string& url) {
  static const boost::regex match_regex("/(\\d+)/(\\d+)/(\\d+)/(\\d+)\\.(\\w+)$");

  boost::smatch match;
  if (!boost::regex_match(url, match, match_regex))
    return DECLINED;

  int id    = boost::lexical_cast<int>(match[1]),
      level = boost::lexical_cast<int>(match[2]),
      col   = boost::lexical_cast<int>(match[3]),
      row   = boost::lexical_cast<int>(match[4]);
  std::string format = boost::lexical_cast<std::string>(match[5]);

  ap_log_rerror(APLOG_MARK, APLOG_INFO, 0, r, "Served Image: id[%i] level[%i] col[%i] row[%i] format[%s]", id, level, col, row, format.c_str());

  PlateModule::IndexCacheEntry index = mod_plate().get_index(id);

  // --------------  Access Plate Index -----------------

  IndexRecord idx_record;
  try {
    idx_record = index.index->read_request(col,row,level,-1);
  } catch(vw::Exception &e) {
    vw_throw(ServerError() << "Could not read plate index: " << e.what());
  }

  // ---------------- Return the image ------------------

  if (idx_record.status() != INDEX_RECORD_VALID)
    vw_throw(NoSuchTile() << "The index record was not valid, or just not there.");

  // Okay, we've gotten this far without error. Set content type now, so HTTP
  // HEAD returns the correct file type
  r->content_type = "image/png";

  // This is as far as we can go without making the request heavyweight. Bail
  // out on a header request now.
  if (r->header_only)
    return OK;

  // These are the sendfile(2) parameters
  std::string filename;
  vw::uint64 offset, size;

  try {
    // Grab a blob from the blob cache by filename
    boost::shared_ptr<Blob> blob = mod_plate().get_blob(index.filename, idx_record.blob_id());
    // And calculate the sendfile(2) parameters
    blob->read_sendfile(idx_record.blob_offset(), filename, offset, size);

  } catch (vw::Exception &e) {
    vw_throw(ServerError() << "Could not load blob data: " << e.what());
  }

  apr_file_t *fd = 0;
  // Open the blob as an apache file with raii (so it goes away when we return)
  raii file_opener(
      boost::bind(apr_file_open, &fd, filename.c_str(), APR_READ|APR_FOPEN_SENDFILE_ENABLED, 0, r->pool),
      boost::bind(apr_file_close, fd));

  // Use sendfile (if available) to send the proper tile data
  size_t sent;
  if (ap_send_fd(fd, r, offset, size, &sent) != APR_SUCCESS || sent != size)
    vw_throw(ServerError() << "Did not send tile correctly!");

  return OK;

}

int handle_wtml(request_rec *r, const std::string& url) {
  static const boost::regex match_regex("/(\\w+\\.wtml)$");

  boost::smatch match;
  if (!boost::regex_match(url, match, match_regex))
    return DECLINED;

  std::string filename = boost::lexical_cast<std::string>(match[1]);

  r->content_type = "text/plain";

  if (r->header_only)
    return OK;

  apache_stream out(r);
  ap_log_rerror(APLOG_MARK, APLOG_INFO, 0, r, "Served WTML: %s", filename.c_str());

  out
    << "<?xml version='1.0' encoding='UTF-8'?>" << std::endl
    << "<Folder Name='ARC Test Data' Group='View'>" << std::endl
    << "  <ImageSet Generic='False' DataSetType='Earth' BandPass='Visible' Name='Bluemarble' Url='http://localhost:31337/plate_example/{1}/{2}/{3}.png' BaseTileLevel='0' TileLevels='2' BaseDegreesPerTile='360' FileType='.png' BottomsUp='False' Projection='Toast' QuadTreeMap='0123' CenterX='0' CenterY='0' OffsetX='0' OffsetY='0' Rotation='0' Sparse='False' ElevationModel='False' StockSet='False'>" << std::endl
    << "    <ThumbnailUrl>http://localhost:31337/plate_example/0/0/0.png</ThumbnailUrl>" << std::endl
    << "    <Credits>NASA</Credits>" << std::endl
    << "  </ImageSet>" << std::endl
    << "</Folder>" << std::endl;

  return OK;
}

PlateModule::PlateModule() {

  // Disable the config file
  vw::vw_settings().set_rc_filename("");

  // And log to stderr, which will go to the apache error log
  vw_log().set_console_stream(std::cerr, vw::LogRuleSet(), false);

  // Create the necessary services
  std::string queue_name = AmqpRpcClient::UniqueQueueName("mod_plate");

  m_client.reset  ( new AmqpRpcClient() );
  m_index_service.reset ( new IndexService::Stub(
                            new AmqpRpcChannel(INDEX_EXCHANGE, "index", queue_name),
                            google::protobuf::Service::STUB_OWNS_CHANNEL) );

  m_index_mgr_service.reset ( new IndexManagerService::Stub(
                                new AmqpRpcChannel(INDEX_MGR_EXCHANGE, "index_mgr", queue_name),
                                google::protobuf::Service::STUB_OWNS_CHANNEL) );

  vw_out(0) << "Initializing mod_plate module.";
}

PlateModule::~PlateModule() {}

int PlateModule::operator()(request_rec *r) const {


  const std::string url(r->path_info);

  //apache_stream out(r);
  //r->content_type = "text/plain";
  //out << url << std::endl;
  //out << "path[" << url << "] pathinfo[" << path_tokens.size() << "]" << std::endl;
  //out << "token[";
  //std::copy(path_tokens.begin(), path_tokens.end(), std::ostream_iterator<std::string>(out, "]\ntoken["));
  //out << "]" << std::endl;
  //return OK;

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

void PlateModule::sync_index_cache() {
  IndexListReply id_list;

  {
    IndexListRequest request;
    m_index_mgr_service->ListRequest(m_client.get(), &request, &id_list, google::protobuf::NewCallback(&null_closure));
  }

  BOOST_FOREACH( const std::string& name, id_list.platefile_names() ) {
    boost::shared_ptr<Index> index = Index::construct_open(std::string("pf://index/") + name);

    const IndexHeader& hdr = index->index_header();

    IndexCacheEntry entry;
    entry.shortname   = fs::path(name).leaf();
    entry.filename    = name;
    entry.description = hdr.has_description() ? hdr.description() : "No Description";

    index_cache[hdr.platefile_id()] = entry;

    vw_out(InfoMessage, "plate") << "Adding " << entry.filename << " to index cache" << std::endl;
  }
}

const PlateModule::IndexCacheEntry PlateModule::get_index(int32 platefile_id) const {

  IndexCache::const_iterator index = index_cache.find(platefile_id);
  if (index != index_cache.end())
    return index->second;

  vw_throw(BadRequest() << "No such platefile [id = " << platefile_id << "]");
}

// --------------------- Apache C++ Entry Points ------------------------

extern "C" void mod_plate_init() {
  // call the singleton to make sure it exists
  static_cast<void>(mod_plate());
}

extern "C" void mod_plate_destroy() {
  kill_mod_plate();
}

extern "C" int mod_plate_handler(request_rec *r) {
  try {
    return mod_plate()(r);
  } catch (const BadRequest& e) {
    // Client sent a request that was formatted badly
    ap_log_rerror(APLOG_MARK, APLOG_NOTICE, 0, r, "Bad Request: %s", e.what());
    return HTTP_BAD_REQUEST;
  } catch (const NoSuchTile& e) {
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
