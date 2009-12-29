// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__
#ifndef __VW_PLATEFILE_PAGED_INDEX_H__
#define __VW_PLATEFILE_PAGED_INDEX_H__

#include <vw/Core/Cache.h>
#include <vw/Core/FundamentalTypes.h>
#include <vw/Plate/IndexPage.h>
#include <vw/Plate/Index.h>
#include <vw/Plate/BlobManager.h>
#include <vector>
#include <string>

namespace vw {
namespace platefile {

  // --------------------------------------------------------------------
  //                             INDEX LEVEL
  // --------------------------------------------------------------------
  class IndexLevel {
    int m_level;
    int m_page_width, m_page_height;
    int m_horizontal_pages, m_vertical_pages;
    std::vector<boost::shared_ptr<IndexPageGenerator> > m_cache_generators;
    std::vector<Cache::Handle<IndexPageGenerator> > m_cache_handles;
    vw::Cache m_cache;

  public:
    IndexLevel(std::string base_path, int level, 
               int page_width, int page_height, int cache_size);

    ~IndexLevel();

    /// Fetch the value of an index node at this level.
    IndexRecord get(int32 col, int32 row, int32 transaction_id, bool exact_match = false) const;

    /// Set the value of an index node at this level.
    void set(IndexRecord const& rec, int32 col, int32 row, int32 transaction_id);
  };

  // --------------------------------------------------------------------
  //                             PAGED INDEX
  // --------------------------------------------------------------------

  class PagedIndex : public Index {
    std::string m_base_path;
    int m_page_width, m_page_height;
    int m_default_cache_size;
    boost::shared_ptr<BlobManager> m_blob_manager;
    std::vector<boost::shared_ptr<IndexLevel> > m_levels;
    Mutex m_mutex;
    
  public:

    PagedIndex(std::string base_path, int page_width, int page_height, int default_cache_size);
    virtual ~PagedIndex() {}

    // ----------------------- READ/WRITE REQUESTS  ----------------------

    virtual IndexRecord read_request(int col, int row, int depth, 
                                     int transaction_id, bool exact_transaction_match = false);

    virtual int write_request(int size);

    virtual void write_complete(TileHeader const& header, IndexRecord const& record);

  
    // ----------------------- PROPERTIES  ----------------------

    virtual IndexHeader index_header() const;
  
    virtual int32 version() const;
    virtual int32 max_depth() const;
  
    virtual std::string platefile_name() const;
  
    virtual int32 tile_size() const;
    virtual std::string tile_filetype() const;
  
    virtual PixelFormatEnum pixel_format() const;
    virtual ChannelTypeEnum channel_type() const;

    // --------------------- TRANSACTIONS ------------------------
  
    /// Clients are expected to make a transaction request whenever
    /// they start a self-contained chunk of mosaicking work.  .
    virtual int32 transaction_request(std::string transaction_description,
                                      std::vector<TileHeader> const& tile_headers);
  
    /// Called right before the beginning of the mipmapping pass
    virtual void root_complete(int32 transaction_id,
                               std::vector<TileHeader> const& tile_headers);
  
    /// Once a chunk of work is complete, clients can "commit" their
    /// work to the mosaic by issuding a transaction_complete method.
    virtual void transaction_complete(int32 transaction_id);
  
    // If a transaction fails, we may need to clean up the mosaic.  
    virtual void transaction_failed(int32 transaction_id);
  
    virtual int32 transaction_cursor();

    // --------------------- UTILITIES ------------------------
  
    /// Iterate over all nodes in a tree, calling func for each
    /// location.  Note: this will only be implemented for local
    /// indexes.  This function will throw an error if called on a
    /// remote index.
    virtual void map(boost::shared_ptr<TreeMapFunc> func);

  };

}} // namespace vw::platefile

#endif // __VW_PLATEFILE_PAGED_INDEX_H__
