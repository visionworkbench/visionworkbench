#ifndef __VW_PLATE_MOD_PLATE_CORE_H__
#define __VW_PLATE_MOD_PLATE_CORE_H__

#include "mod_plate_utils.h"
#include "mod_plate.h"

#include <vw/Core/FundamentalTypes.h>
#include <vw/Core/Log.h>
#include <vw/Core/Exception.h>
#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>
#include <string>
#include <map>


namespace vw {
namespace platefile {

class Index;
class Blob;
class AmqpRpcClient;
class IndexService;

typedef boost::function<int (const ApacheRequest&)> Handler;

class PlateModule {
  public:
    PlateModule(const plate_config *conf);
    ~PlateModule();
    void connect_index();
    int operator()(const ApacheRequest& r) const;
    int status(const ApacheRequest& r, int flags) const;

    struct IndexCacheEntry {
      std::string shortname;
      std::string filename;
      std::string description;
      boost::shared_ptr<Index> index;
      int read_cursor;
    };

    typedef std::map<int32, IndexCacheEntry> IndexCache;

    const IndexCache& get_index_cache() const { return index_cache; }

    const IndexCacheEntry& get_index(const std::string& id_str) const;

    const boost::shared_ptr<Blob> get_blob(int platefile_id, const std::string& plate_filename, uint32 blob_id) const;
    void sync_index_cache() const;

    std::ostream& logger(MessageLevel level, bool child_id = true) const;

    std::string get_dem() const;
    bool allow_resync() const;
    std::string get_servername() const;

  private:
    struct BlobCacheEntry {
      boost::shared_ptr<Blob> blob;
      int platefile_id;
      BlobCacheEntry(boost::shared_ptr<Blob> b, int id) :
        blob(b), platefile_id(id) {}
    };

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

const PlateModule& mod_plate();
PlateModule& mod_plate_mutable();
void mod_plate_init(const plate_config *conf);

}} // namespace vw::platefile

#endif
