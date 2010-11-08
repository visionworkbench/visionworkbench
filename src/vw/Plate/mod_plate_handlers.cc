// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include "mod_plate_handlers.h"
#include "mod_plate_core.h"
#include "mod_plate_utils.h"
#include "mod_plate.h"

#include <vw/Plate/Blob.h>
#include <vw/Plate/Index.h>
#include <vw/Plate/Exception.h>

#include <httpd.h>
#include <http_core.h>
#include <http_protocol.h>
#include <http_log.h>
#include <mod_status.h>
#include <ap_mpm.h>

#include <boost/regex.hpp>
#include <boost/foreach.hpp>

using namespace vw;
using namespace vw::platefile;

using std::string;

int vw::platefile::handle_image(const ApacheRequest& r) {
  static const boost::regex match_regex("/(\\w+)/(\\d+)/(\\d+)/(\\d+)\\.(\\w+)$");

  boost::smatch match;
  if (!boost::regex_search(r.url, match, match_regex))
    return DECLINED;

  // We didn't decline. Connect!
  mod_plate_mutable().connect_index();

  const string& sid = match[1];

  int level = boost::lexical_cast<int>(match[2]),
      col   = boost::lexical_cast<int>(match[3]),
      row   = boost::lexical_cast<int>(match[4]);
  string format = boost::lexical_cast<string>(match[5]);

  mod_plate().logger(DebugMessage) << "Request Image: id["  << sid
                                   << "] level["  << level
                                   << "] col["    << col
                                   << "] row["    << row
                                   << "] format[" << format << "]" << std::endl;

  const PlateModule::IndexCacheEntry& index = mod_plate().get_index(sid);

  int id = index.index->index_header().platefile_id();

  // --------------  Access Plate Index -----------------

  IndexRecord idx_record;
  try {
    int transaction_id = r.args.get("transaction_id", int(-1));
    bool exact = r.args.get("exact", false);

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

  mod_plate().logger(VerboseDebugMessage) << "Figuring out mime content type from filetype " << idx_record.filetype() << std::endl;
  // Okay, we've gotten this far without error. Set content type now, so HTTP
  // HEAD returns the correct file type
  if (idx_record.filetype() == "png")
    r.writer()->content_type = "image/png";
  else if (idx_record.filetype() == "jpg")
    r.writer()->content_type = "image/jpeg";
  else if (idx_record.filetype() == "tif")
    r.writer()->content_type = "image/tiff";
  else
    r.writer()->content_type = "application/octet-stream";

  if (r.args.get("nocache", 0u) == 1) {
    apr_table_set(r.writer()->headers_out, "Cache-Control", "no-cache");
  }
  else {
    if (level <= 7)
      apr_table_set(r.writer()->headers_out, "Cache-Control", "max-age=604800");
    else
      apr_table_set(r.writer()->headers_out, "Cache-Control", "max-age=1200");
  }

  // This is as far as we can go without making the request heavyweight. Bail
  // out on a header request now.
  if (r.header_only())
    return OK;

  // These are the sendfile(2) parameters
  string filename;
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
      boost::bind(apr_file_open, &fd, filename.c_str(), APR_READ|APR_FOPEN_SENDFILE_ENABLED, 0, r.writer()->pool),
      boost::bind(apr_file_close, boost::ref(fd)));

  // Use sendfile (if available) to send the proper tile data
  size_t sent;
  apr_status_t ap_ret;

  if ((ap_ret = ap_send_fd(fd, r.writer(), offset, size, &sent) != APR_SUCCESS)) {
    char buf[256];
    apr_strerror(ap_ret, buf, 256);
    vw_throw(ServerError() << "ap_send_fd failed: " << buf);
  }
  else if (sent != size)
    vw_throw(ServerError() << "ap_send_fd: short write (expected to send " << size << " bytes, but only sent " << sent);

  return OK;
}

int vw::platefile::handle_wtml(const ApacheRequest& r) {
  static const boost::regex match_regex("/(\\w+\\.wtml)$");

  boost::smatch match;
  if (!boost::regex_match(r.url, match, match_regex))
    return DECLINED;

  mod_plate_mutable().connect_index();

  string filename = boost::lexical_cast<string>(match[1]);

  r.writer()->content_type = "application/xml";

  if (r.header_only())
    return OK;

  if (mod_plate().allow_resync())
    mod_plate().sync_index_cache();

  apache_stream out(r.writer());
  mod_plate().logger(DebugMessage) << "Served WTML[" << filename << "]" << std::endl;

  out
    << "<?xml version='1.0' encoding='UTF-8'?>"              << std::endl
    << "<Folder Name='Ames Planetary Content' Group='View'>" << std::endl << std::endl;

  typedef std::pair<int32, PlateModule::IndexCacheEntry> id_cache;

  string servername = mod_plate().get_servername();

  bool show_all_layers =  r.args.get("all_layers", false);

  BOOST_FOREACH( const id_cache& e, mod_plate().get_index_cache() ) {

    const string filetype = e.second.index->index_header().tile_filetype();
    // WWT can only handle jpg, png, and tif
    if (!show_all_layers) {
      if (filetype != "jpg" && filetype != "png" && filetype != "tif" && filetype != "auto") {
        vw_out(VerboseDebugMessage) << "Rejecting filetype " << filetype << " for WTML." << std::endl;
        continue;
      }
    }

    string dem_id = r.args.get("override_dem", string());
    if (dem_id.empty())
      dem_id = mod_plate().get_dem();

    WTMLImageSet img(servername, "/wwt/p/", "/static/", dem_id, e.second.index, e.second.description);

    // & needs to be escaped in xml attributes. Rather than deal with that,
    // use ; as the query arg sep (yes, ; is in the spec)
    std::string args = r.args.serialize("?", ";");

    if (!args.empty()) {
      img["Url"]          += args;
      img["ThumbnailUrl"] += args;
      img["DemUrl"]       += args;
    }
    img.serializeToOstream(out);
  }
  out << "</Folder>" << std::endl;

  return OK;
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
  } catch (const UnknownPlatefile& e) {
    // Valid format, but not there
    ap_log_rerror(APLOG_MARK, APLOG_NOTICE, 0, r,  "Platefile not found: %s", e.what());
    return HTTP_NOT_FOUND;
  } catch (const ServerError& e) {
    // Something screwed up, but we controlled it
    ap_log_rerror(APLOG_MARK, APLOG_ERR, 0, r,  "Server Error [recovered]: %s", e.what());
    return HTTP_INTERNAL_SERVER_ERROR;
  } catch (const Exception& e) {
    // Something screwed up worse...
    ap_log_rerror(APLOG_MARK, APLOG_CRIT, 0, r, "Server Error [vw::Exception]: %s", e.what());
    return HTTP_INTERNAL_SERVER_ERROR;
  } catch (const std::exception &e) {
    // Something we don't understand broke. Eek.
    ap_log_rerror(APLOG_MARK, APLOG_ALERT, 0, r, "Server Error [std::exception]: %s", e.what());
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
    mod_plate_init(get_plate_config(s));
    mod_plate();
  } catch (const Exception& e) {
    ap_log_error(APLOG_MARK, APLOG_ALERT, 0, s, "Could not start mod_plate child! [uncaught vw::Exception]: %s", e.what());
  } catch (const std::exception& e) {
    ap_log_error(APLOG_MARK, APLOG_ALERT, 0, s, "Could not start mod_plate child! [uncaught std::exception]: %s", e.what());
  } catch (...) {
    ap_log_error(APLOG_MARK, APLOG_ALERT, 0, s, "Could not start mod_plate child! [uncaught unknown exception]");
  }
}

#define CONF_REQUIRE(elt, dir)\
  if (!conf->elt) {\
    fprintf(stderr, "missing required config directive: %s\n", dir);\
    return HTTP_INTERNAL_SERVER_ERROR;\
  }

extern "C" int mod_plate_post_config(apr_pool_t* /*pconf*/, apr_pool_t* /*plog*/, apr_pool_t* /*ptemp*/, server_rec *s) {
  const plate_config* conf = get_plate_config(s);
  CONF_REQUIRE(index_url,      "PlateIndexUrl");
  CONF_REQUIRE(dem_id,         "PlateDemID");
  CONF_REQUIRE(servername,     "PlateServerName");

  int threaded;

  if (ap_mpm_query(AP_MPMQ_IS_THREADED, &threaded) != APR_SUCCESS) {
    fprintf(stderr, "mod_plate could not query mpm: we have a problem. Refusing to register!\n");
    return HTTP_INTERNAL_SERVER_ERROR;
  }

  if (threaded) {
    fprintf(stderr, "Refusing to start mod_plate inside a threaded MPM. We don't know how our threads will interact with apache's.\n");
    return HTTP_INTERNAL_SERVER_ERROR;
  }

  return OK;
}
