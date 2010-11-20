// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include "mod_plate_core.h"
#include "mod_plate_utils.h"
#include "mod_plate_handlers.h"

#include <httpd.h>
#include <apr_tables.h>
#include <vw/Core/Settings.h>
#include <vw/Plate/Index.h>
#include <vw/Plate/Blob.h>
#include <vw/Plate/Rpc.h>
#include <vw/Plate/IndexService.pb.h>

#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <string>

using std::string;

namespace vw {
namespace platefile {

PlateModule::~PlateModule() { }

const PlateModule::IndexCacheEntry& PlateModule::get_index(const string& id_str) const {

  int id;
  PlateModule::IndexCache::const_iterator index_i;

  std::ostringstream msg_name;

  msg_name << "lookup: " << id_str;
  try {
    id = boost::lexical_cast<int>(id_str);
  } catch (const boost::bad_lexical_cast&) {
    id = 0;
  }
  msg_name << " parsed as " << id;

  // try it as an id (unless it's clearly wrong)
  if (id != 0) {
    index_i = index_cache.find(id);
    if (index_i != index_cache.end()) {
      logger(VerboseDebugMessage) << msg_name.str() << std::endl;
      return index_i->second;
    }
  }

  // Try it as an alias
  {
    int id2 = reinterpret_cast<intptr_t>(apr_table_get(m_conf->alias, id_str.c_str()));
    if (id2) {
      msg_name << " resolved as alias to " << id2;
      index_i = index_cache.find(id2);
      if (index_i != index_cache.end()) {
        logger(VerboseDebugMessage) << msg_name.str() << std::endl;
        return index_i->second;
      }
    }
  }

  if (!allow_resync())
    vw_throw(UnknownPlatefile() << "No such platefile (no resync) for " << msg_name.str());

  // If we get an unknown platefile (and we're allowed to do so) resync just to make sure
  logger(WarningMessage) << "Platefile [" << msg_name.str() << "] not in platefile cache. Resyncing." << std::endl;
  sync_index_cache();

  index_i = index_cache.find(id);
  if (index_i == index_cache.end())
    vw_throw(UnknownPlatefile() << "No such platefile (after resync) for " << msg_name.str());

  return index_i->second;
}


std::ostream& PlateModule::logger(MessageLevel level, bool child_id) const {
  std::ostream& out = vw_out(level, "plate.apache");

  if (child_id)
    out << "Child[" << getpid() << "]: ";
  return out;
}

string PlateModule::get_dem() const {
  return m_conf->dem_id;
}

bool PlateModule::allow_resync() const {
  return m_conf->unknown_resync;
}

string PlateModule::get_servername() const {
  return m_conf->servername;
}

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

void mod_plate_init(const plate_config *conf) {
  mod_plate_ptr.reset(new PlateModule(conf));
}

PlateModule::PlateModule(const plate_config* conf)
  : m_connected(false), m_conf(conf), m_base_url(conf->index_url)
{

  // Disable the config file
  vw::vw_settings().set_rc_filename("");

  int i;

  LogRuleSet rules;
  if (m_conf->rules->nelts == 0)
    rules.add_rule(DebugMessage, "plate.apache");
  else {
    for (i = 0; i < m_conf->rules->nelts; ++i) {
      rule_entry *entry = reinterpret_cast<rule_entry*>(&m_conf->rules->elts[i]);
      rules.add_rule(entry->level, entry->name);
    }
  }

  // And log to stderr, which will go to the apache error log
  vw_log().set_console_stream(std::cerr, rules, true);
}

void PlateModule::connect_index() {

  if (m_connected)
    return;

  m_client.reset( new IndexClient(m_base_url) );

  // Total possible wait is timeout*tries
  m_client->set_timeout(m_conf->index_timeout);
  m_client->set_retries(m_conf->index_tries-1);

  m_connected = true;
  logger(DebugMessage) << "child connected to url[" << m_base_url.string() << "]" << std::endl;

  mod_plate().sync_index_cache();
}

int PlateModule::operator()(const ApacheRequest& r) const {

  if (r.url.empty())
    return DECLINED;

  static const Handler Handlers[] = {handle_image, handle_wtml};

  BOOST_FOREACH(const Handler h, Handlers) {
    int ret = h(r);
    if (ret != DECLINED)
      return ret;
  }
  return DECLINED;
}

int PlateModule::status(const ApacheRequest& r, int /*flags*/) const {
  apache_stream out(r.writer());

  out << "IndexCache:<br>" << std::endl;
  BOOST_FOREACH(const IndexCache::value_type& c, get_index_cache())
    out << c.second.shortname << ": " << c.first << "<br>";

  out << "BlobCacheSize: " << get_blob_cache().size() << "<br>";

  return OK;
}

const boost::shared_ptr<Blob> PlateModule::get_blob(int platefile_id, const string& plate_filename, uint32 blob_id) const {
  std::ostringstream ostr;
  ostr << plate_filename << "/plate_" << blob_id << ".blob";
  const string& filename = ostr.str();

  if (m_conf->use_blob_cache) {

    BlobCache::const_iterator blob = blob_cache.find(filename);

    // Check the platefile id to make sure the blob wasn't deleted and recreated
    // with a different platefile
    if (blob != blob_cache.end() && blob->second.platefile_id == platefile_id)
      return blob->second.blob;
  }

  boost::shared_ptr<Blob> ret( new Blob(filename, true) );

  if (m_conf->use_blob_cache)
    blob_cache.insert(std::make_pair(filename, BlobCacheEntry(ret, platefile_id)));

  return ret;
}

void PlateModule::sync_index_cache() const {

  VW_ASSERT(m_connected, LogicErr() << "Must connect before trying to sync cache");

  IndexListRequest request;
  IndexListReply id_list;

  index_cache.clear();

  m_client->ListRequest(m_client.get(), &request, &id_list, null_callback());

  Url base_url = get_base_url();
  base_url.query().set("cache_size", "1");
  Url::split_t path = base_url.path_split();
  path.push_back("");

  BOOST_FOREACH( const string& name, id_list.platefile_names() ) {

    IndexCacheEntry entry;
    int32 id;

    try {
      Url index_url(base_url);
      path.back() = name;
      index_url.path_join(path);
      logger(VerboseDebugMessage) << "Trying to load index: " << index_url.string() << std::endl;

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


}} // namespace vw::platefile
