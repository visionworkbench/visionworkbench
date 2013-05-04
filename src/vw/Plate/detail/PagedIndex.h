// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NASA Vision Workbench is licensed under the Apache License,
//  Version 2.0 (the "License"); you may not use this file except in
//  compliance with the License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__


#ifndef __VW_PLATEFILE_PAGED_INDEX_H__
#define __VW_PLATEFILE_PAGED_INDEX_H__

#include <vw/Plate/FundamentalTypes.h>
#include <vw/Core/Cache.h>
#include <vw/Plate/detail/Index.h>
#include <vw/Plate/detail/IndexPage.h>

#include <vector>
#include <list>

namespace vw {
namespace platefile {
  class TileHeader;

namespace detail {


  // --------------------------------------------------------------------
  //                             INDEX LEVEL
  // --------------------------------------------------------------------
  class IndexLevel {

    typedef Cache::Handle<boost::shared_ptr<PageGeneratorBase> > handle_t;

    boost::shared_ptr<PageGeneratorFactory> m_page_gen_factory;
    uint32 m_level;
    uint32 m_page_width, m_page_height;
    uint32 m_horizontal_pages, m_vertical_pages;
    mutable std::vector<handle_t> m_cache_handles;
    mutable vw::Cache m_cache;
    mutable Mutex m_cache_mutex;

    boost::shared_ptr<IndexPage> load_page(uint32 col, uint32 row) const;

  public:
    typedef IndexPage::multi_value_type multi_value_type;

    IndexLevel(boost::shared_ptr<PageGeneratorFactory> page_gen_factory,
               uint32 level, uint32 page_width, uint32 page_height, uint32 cache_size);

    ~IndexLevel();

    /// Sync any unsaved data in the index to disk.
    void sync();

    // Returns the page id for a page within the level
    uint32 page_id(uint32 col, uint32 row) const;

    /// Grab an IndexPage.  Useful if you want to serialize it by hand to disk.
    boost::shared_ptr<IndexPage> get_page(uint32 col, uint32 row) const;

    /// Fetch the value of an index node at this level.
    IndexRecord get(int32 col, int32 row, TransactionOrNeg transaction_id, bool exact_match = false) const;

    /// Set the value of an index node at this level.
    void set(TileHeader const& hdr, IndexRecord const& rec);

    /// Returns a list of valid tiles at this level.
    std::list<TileHeader> search_by_region(BBox2i const& region,
                                           TransactionOrNeg start_transaction_id,
                                           TransactionOrNeg end_transaction_id) const;

    /// Returns a list of valid tiles at this level and specified location
    std::list<TileHeader> search_by_location(uint32 col, uint32 row,
                                             TransactionOrNeg start_transaction_id, TransactionOrNeg end_transaction_id) const;
  };

  // --------------------------------------------------------------------
  //                             PAGED INDEX
  // --------------------------------------------------------------------

  class PagedIndex : public Index {

  protected:

    boost::shared_ptr<PageGeneratorFactory> m_page_gen_factory;
    mutable std::vector<boost::shared_ptr<IndexLevel> > m_levels;
    uint32 m_page_width, m_page_height;
    uint32 m_default_cache_size;

    void set_page_generator_factory(boost::shared_ptr<PageGeneratorFactory> page_gen_factory) {
      m_page_gen_factory = page_gen_factory;
    }

  public:
    typedef IndexLevel::multi_value_type multi_value_type;

    /// Open an existing index from a file on disk.
    PagedIndex(boost::shared_ptr<PageGeneratorFactory> page_generator,
               uint32 page_width = 256, uint32 page_height = 256,
               uint32 default_cache_size = 100);

    /// Open an existing index from a file on disk.
    PagedIndex(uint32 page_width = 256, uint32 page_height = 256,
               uint32 default_cache_size = 100);

    virtual ~PagedIndex() {}

    /// Sync any unsaved data in the index to disk.
    virtual void sync();

    virtual void set_default_cache_size(uint32 size) {
      m_default_cache_size = size;
    }

    // ----------------------- READ/WRITE REQUESTS  ----------------------

    /// Grab an IndexPage.  Useful if you want to serialize it by hand to disk.
    virtual boost::shared_ptr<IndexPage> page_request(uint32 col, uint32 row, uint32 level) const;

    // This returns a page id such that two tiles with the same page id will
    // come from the same page (no other ordering is implied)
    virtual uint64 page_id(uint32 col, uint32 row, uint32 level) const;

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
    virtual IndexRecord read_request(uint32 col, uint32 row, uint32 level,
                                     TransactionOrNeg transaction_id, bool exact_transaction_match = false);

    // Writing, pt. 2: Supply information to update the index and
    // unlock the blob id.
    virtual void write_update(TileHeader const& header, IndexRecord const& record);

    // ----------------------- PROPERTIES  ----------------------

    /// Returns a list of valid tiles that match this level, region, and
    /// range of transaction_id's.  Returns a list of TileHeaders with
    /// col/row/level and transaction_id of the most recent tile at each
    /// valid location.  Note: there may be other tiles in the transaction
    /// range at this col/row/level, but valid_tiles() only returns the
    /// first one.
    virtual std::list<TileHeader> search_by_region(uint32 level, BBox2i const& region,
                                                   TransactionOrNeg start_transaction_id,
                                                   TransactionOrNeg end_transaction_id) const;

    /// Returns a list of tile headers for a given tile location in
    /// the mosaic, subject to the specified transaction_id range.
    virtual std::list<TileHeader> search_by_location(uint32 col, uint32 row, uint32 level,
                                                     TransactionOrNeg start_transaction_id,
                                                     TransactionOrNeg end_transaction_id) const;

  };

}}}

#endif // __VW_PLATEFILE_PAGED_INDEX_H__
