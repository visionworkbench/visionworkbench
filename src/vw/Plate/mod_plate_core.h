// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_PLATE_MOD_PLATE_CORE_H__
#define __VW_PLATE_MOD_PLATE_CORE_H__

#include "mod_plate_utils.h"
#include "mod_plate.h"

#include <vw/Core/FundamentalTypes.h>
#include <vw/Core/Log.h>
#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>
#include <string>
#include <map>


namespace vw {
namespace platefile {

class Index;
class Blob;
class IndexService;

template <typename ServiceT>
class RpcClient;

typedef RpcClient<IndexService> IndexClient;

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

    struct BlobCacheEntry {
      boost::shared_ptr<Blob> blob;
      int platefile_id;
      BlobCacheEntry(boost::shared_ptr<Blob> b, int id) :
        blob(b), platefile_id(id) {}
    };

    typedef std::map<int32, IndexCacheEntry> IndexCache;
    typedef std::map<std::string, BlobCacheEntry> BlobCache;

    const IndexCache& get_index_cache() const { return index_cache; }
    const BlobCache&  get_blob_cache() const {return blob_cache; }

    const IndexCacheEntry& get_index(const std::string& id_str) const;

    const boost::shared_ptr<Blob> get_blob(int platefile_id, const std::string& plate_filename, uint32 blob_id) const;
    void sync_index_cache() const;

    std::ostream& logger(MessageLevel level, bool child_id = true) const;

    std::string get_dem() const;
    bool allow_resync() const;
    std::string get_servername() const;
    const Url& get_base_url() const;

  private:
    boost::shared_ptr<IndexClient> m_client;

    mutable BlobCache  blob_cache;
    mutable IndexCache index_cache;
    bool m_connected;
    // We don't manage the data, and I think apache might modify it behind the scenes.
    // As such, mark it volatile.
    volatile const plate_config *m_conf;
    Url m_base_url;
};

const PlateModule& mod_plate();
PlateModule& mod_plate_mutable();
void mod_plate_init(const plate_config *conf);

}} // namespace vw::platefile

#endif
