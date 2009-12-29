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
#include <vw/Plate/LocalIndex.h>
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

  class PagedIndex : public LocalIndex {
    int m_page_width, m_page_height;
    int m_default_cache_size;
    std::vector<boost::shared_ptr<IndexLevel> > m_levels;

    virtual void load_index(std::string plate_filename,
                            std::vector<std::string> const& blob_files);

    void commit_record(IndexRecord const& record, int col, int row, int level, int transaction_id);
    
  public:

    /// Create a new, empty index.
    PagedIndex(std::string plate_filename, IndexHeader new_index_info, 
               int page_width = 256, int page_height = 256, 
               int default_cache_size = 10000);

    /// Open an existing index from a file on disk.
    PagedIndex(std::string plate_filename,
               int page_width = 256, int page_height = 256, 
               int default_cache_size = 10000);

    virtual ~PagedIndex() {}

    // ----------------------- READ/WRITE REQUESTS  ----------------------

    /// Attempt to access a tile in the index.  Throws an
    /// TileNotFoundErr if the tile cannot be found.
    /// 
    /// By default, this call to read will return a tile with the MOST
    /// RECENT transaction_id <= to the transaction_id you specify
    /// here in the function arguments (if a tile exists).  However,
    /// setting exact_transaction_match = true will force the
    /// PlateFile to search for a tile that has the EXACT SAME
    /// transaction_id as the one that you specify.
    ///
    /// A transaction ID of -1 indicates that we should return the
    /// most recent tile, regardless of its transaction id.
    virtual IndexRecord read_request(int col, int row, int level, 
                                     int transaction_id, bool exact_transaction_match = false);

    // Writing, pt. 1: Locks a blob and returns the blob id that can
    // be used to write a tile.
    virtual int write_request(int size);

    // Writing, pt. 2: Supply information to update the index and
    // unlock the blob id.
    virtual void write_complete(TileHeader const& header, IndexRecord const& record);
  
    // ----------------------- PROPERTIES  ----------------------

    virtual int32 num_levels() const;

    void map(boost::shared_ptr<TreeMapFunc> func);
  };

}} // namespace vw::platefile

#endif // __VW_PLATEFILE_PAGED_INDEX_H__
