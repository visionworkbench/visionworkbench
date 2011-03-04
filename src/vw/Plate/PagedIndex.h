// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_PLATEFILE_PAGED_INDEX_H__
#define __VW_PLATEFILE_PAGED_INDEX_H__

#include <vw/Plate/FundamentalTypes.h>
#include <vw/Core/Cache.h>
#include <vw/Plate/Index.h>
#include <vw/Plate/IndexPage.h>

#include <vector>
#include <list>

namespace vw {
namespace platefile {
  class TileHeader;

  // --------------------------------------------------------------------
  //                             INDEX LEVEL
  // --------------------------------------------------------------------
  class IndexLevel {

    typedef Cache::Handle<boost::shared_ptr<PageGeneratorBase> > handle_t;

    boost::shared_ptr<PageGeneratorFactory> m_page_gen_factory;
    int m_level;
    int m_page_width, m_page_height;
    int m_horizontal_pages, m_vertical_pages;
    mutable std::vector<handle_t> m_cache_handles;
    mutable vw::Cache m_cache;
    mutable Mutex m_cache_mutex;

    boost::shared_ptr<IndexPage> fetch_page(int col, int row) const;

  public:
    typedef IndexPage::multi_value_type multi_value_type;

    IndexLevel(boost::shared_ptr<PageGeneratorFactory> page_gen_factory,
               int level, int page_width, int page_height, int cache_size);

    ~IndexLevel();

    /// Sync any unsaved data in the index to disk.
    void sync();

    /// Grab an IndexPage.  Useful if you want to serialize it by hand to disk.
    boost::shared_ptr<IndexPage> get_page(int col, int row) const;

    /// Fetch the value of an index node at this level.
    IndexRecord get(int32 col, int32 row, TransactionOrNeg transaction_id, bool exact_match = false) const;

    /// Set the value of an index node at this level.
    void set(TileHeader const& hdr, IndexRecord const& rec);

    /// Returns a list of valid tiles at this level.
    std::list<TileHeader> search_by_region(BBox2i const& region,
                                           TransactionOrNeg start_transaction_id,
                                           TransactionOrNeg end_transaction_id,
                                           uint32 min_num_matches,
                                           bool fetch_one_additional_entry) const;

    /// Returns a list of valid tiles at this level and specified location
    std::list<TileHeader> search_by_location(int col, int row,
                                             TransactionOrNeg start_transaction_id, TransactionOrNeg end_transaction_id,
                                             bool fetch_one_additional_entry = false) const;
  };

  // --------------------------------------------------------------------
  //                             PAGED INDEX
  // --------------------------------------------------------------------

  class PagedIndex : public Index {

  protected:

    boost::shared_ptr<PageGeneratorFactory> m_page_gen_factory;
    mutable std::vector<boost::shared_ptr<IndexLevel> > m_levels;
    int m_page_width, m_page_height;
    int m_default_cache_size;

    void set_page_generator_factory(boost::shared_ptr<PageGeneratorFactory> page_gen_factory) {
      m_page_gen_factory = page_gen_factory;
    }

  public:
    typedef IndexLevel::multi_value_type multi_value_type;

    /// Open an existing index from a file on disk.
    PagedIndex(boost::shared_ptr<PageGeneratorFactory> page_generator,
               int page_width = 256, int page_height = 256,
               int default_cache_size = 100);

    /// Open an existing index from a file on disk.
    PagedIndex(int page_width = 256, int page_height = 256,
               int default_cache_size = 100);

    virtual ~PagedIndex() {}

    /// Sync any unsaved data in the index to disk.
    virtual void sync();

    virtual void set_default_cache_size(int size) {
      m_default_cache_size = size;
    }

    // ----------------------- READ/WRITE REQUESTS  ----------------------

    /// Grab an IndexPage.  Useful if you want to serialize it by hand to disk.
    virtual boost::shared_ptr<IndexPage> page_request(int col, int row, int level) const;

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
                                     TransactionOrNeg transaction_id, bool exact_transaction_match = false);

    // Writing, pt. 1: Locks a blob and returns the blob id that can
    // be used to write a tile.
    virtual int write_request(uint64 &size) = 0;

    // Writing, pt. 2: Supply information to update the index and
    // unlock the blob id.
    virtual void write_update(TileHeader const& header, IndexRecord const& record);

    /// Writing, pt. 3: Signal the completion
    virtual void write_complete(int blob_id, uint64 blob_offset) = 0;

    // ----------------------- PROPERTIES  ----------------------

    /// Returns a list of valid tiles that match this level, region, and
    /// range of transaction_id's.  Returns a list of TileHeaders with
    /// col/row/level and transaction_id of the most recent tile at each
    /// valid location.  Note: there may be other tiles in the transaction
    /// range at this col/row/level, but valid_tiles() only returns the
    /// first one.
    virtual std::list<TileHeader> search_by_region(int level, BBox2i const& region,
                                                   TransactionOrNeg start_transaction_id,
                                                   TransactionOrNeg end_transaction_id,
                                                   uint32 min_num_matches,
                                                   bool fetch_one_additional_entry) const;

    /// Returns a list of tile headers for a given tile location in
    /// the mosaic, subject to the specified transaction_id range.
    virtual std::list<TileHeader> search_by_location(int col, int row, int level,
                                                     TransactionOrNeg start_transaction_id,
                                                     TransactionOrNeg end_transaction_id,
                                                     bool fetch_one_additional_entry) const;

  };

}} // namespace vw::platefile

#endif // __VW_PLATEFILE_PAGED_INDEX_H__
