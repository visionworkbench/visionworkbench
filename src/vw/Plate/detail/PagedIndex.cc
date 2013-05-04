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


#include <vw/Plate/detail/PagedIndex.h>
#include <vw/Plate/Exception.h>
#include <vw/Core/Debugging.h>

#include <boost/foreach.hpp>

using namespace vw;
using namespace vw::platefile;
using namespace vw::platefile::detail;

#define WHEREAMI (vw::vw_out(VerboseDebugMessage, "platefile.pagedindex") << VW_CURRENT_FUNCTION << ": ")


// --------------------------------------------------------------------
//                             INDEX LEVEL
// --------------------------------------------------------------------

IndexLevel::IndexLevel(boost::shared_ptr<PageGeneratorFactory> page_gen_factory,
                       uint32 level, uint32 page_width, uint32 page_height, uint32 cache_size)
  : m_page_gen_factory(page_gen_factory), m_level(level),
    m_page_width(page_width), m_page_height(page_height), m_cache(cache_size) {

  uint32 tiles_per_side = 1 << level;
  m_horizontal_pages = static_cast<uint32>(ceil(float(tiles_per_side) / float(page_width)));
  m_vertical_pages   = static_cast<uint32>(ceil(float(tiles_per_side) / float(page_height)));

  // Create space for cache handles.  The actual generators are not
  // created until they are needed (because they take enough memory
  // that it's not efficient to allocate the generators ahead of
  // time).  Actual allocation is done automatically by calling
  // load_cache_handle().
  uint32 pages = m_horizontal_pages * m_vertical_pages;
  m_cache_handles.resize(pages);
}

uint32 IndexLevel::page_id(uint32 col, uint32 row) const {
  const size_t MAX_IDX = 1 << m_level;
  VW_ASSERT( col < MAX_IDX && row < MAX_IDX,
             TileNotFoundErr() << "IndexLevel::page_id(" << col << "," << row << ") (level=" << m_level << "): Invalid index");

  uint32 page_col = col / m_page_width;
  uint32 page_row = row / m_page_height;

  return page_row*m_horizontal_pages + page_col;
}

/// Grab an IndexPage.  Useful if you want to serialize it by hand to disk.
boost::shared_ptr<IndexPage> IndexLevel::get_page(uint32 col, uint32 row) const {
  // Access the page.  This will load it into memory if necessary.
  return load_page(col, row);
}

namespace {
  template <typename T, typename U>
  T floorto(T x, U y) {
    BOOST_STATIC_ASSERT(boost::is_integral<T>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<U>::value);
    return x - (x % y);
  }
}

boost::shared_ptr<IndexPage> IndexLevel::load_page(uint32 col, uint32 row) const {
  size_t idx = boost::numeric_cast<size_t>(this->page_id(col,row));

  Mutex::Lock lock(m_cache_mutex);
  // We may need to actually create the page's cache handle if it
  // hasn't been created already.
  if ( !m_cache_handles[idx].attached() ) {
    // the args to create here are the base col and row, so (col / m_page_width) * m_page_width
    boost::shared_ptr<PageGeneratorBase> generator =
      m_page_gen_factory->create(m_level, floorto(col, m_page_width), floorto(row, m_page_height), m_page_width, m_page_height);
      m_cache_handles[idx] = m_cache.insert( generator );
  }

  // WARNING! CACHE MIGHT DELETE YOUR POINTER HERE!
  boost::shared_ptr<IndexPage> result = m_cache_handles[idx];
  m_cache_handles[idx].release();
  return result;
}


IndexLevel::~IndexLevel() {
  Mutex::Lock lock(m_cache_mutex);

  // Make sure the handles drop out of cache
  BOOST_FOREACH( handle_t& h, m_cache_handles )
    h.reset();
}

void IndexLevel::sync() {
  Mutex::Lock lock(m_cache_mutex);

  vw_out(VerboseDebugMessage, "platefile.cache")
    << "Page cache for " << m_page_gen_factory->who() << "@" << m_level << " reports "
    << "hits["   << m_cache.hits()
    << "] misses[" << m_cache.misses()
    << "] evictions[" << m_cache.evictions() << "] since last sync." << std::endl;

  m_cache.clear_stats();

  // Write the index pages to disk by calling their sync() methods.
  BOOST_FOREACH( handle_t& h, m_cache_handles ) {
    if (h.attached() && h.valid())
      h->sync();
  }
}

/// Fetch the value of an index node at this level.
IndexRecord IndexLevel::get(int32 col, int32 row, TransactionOrNeg transaction_id, bool exact_match) const {
  boost::shared_ptr<IndexPage> page = load_page(col, row);
  return page->get(col, row, transaction_id, exact_match);
}

/// Set the value of an index node at this level.
void IndexLevel::set(TileHeader const& header, IndexRecord const& rec) {
  boost::shared_ptr<IndexPage> page = load_page(header.col(), header.row());
  page->set(header, rec);
}

namespace {
  uint32 round_to(uint32 val, uint32 stride) {
    return (val / stride) * stride;
  }
}

/// Returns a list of valid tiles at this level.
std::list<TileHeader>
IndexLevel::search_by_region(BBox2i const& region,
                             TransactionOrNeg start_transaction_id,
                             TransactionOrNeg end_transaction_id) const {

  // Start by computing the search range in pages based on the requested region.
  uint32 min_level_col = round_to(region.min().x(), m_page_width);
  uint32 min_level_row = round_to(region.min().y(), m_page_height);

  uint32 max_level_col = round_to(region.max().x() + m_page_width  - 1,  m_page_width);
  uint32 max_level_row = round_to(region.max().y() + m_page_height - 1, m_page_height);

  WHEREAMI << "[" << min_level_col << " " << min_level_row << "]" << " to [" << max_level_col << " " << max_level_row << "]\n";

  // Iterate over the pages that overlap with the region of interest.
  std::list<TileHeader> result;
  for (uint32 level_row = min_level_row; level_row < max_level_row; level_row += m_page_height) {
    for (uint32 level_col = min_level_col; level_col < max_level_col; level_col += m_page_width) {
      boost::shared_ptr<IndexPage> page = load_page(level_col, level_row);

      // Accumulate valid tiles that overlap with region from this IndexPage.
      std::list<TileHeader> sub_result = page->search_by_region(region, start_transaction_id, end_transaction_id);
      result.splice(result.end(), sub_result);
    }
  }
  return result;
}

/// Fetch the value of an index node at this level.
std::list<TileHeader>
IndexLevel::search_by_location(uint32 col, uint32 row,
                               TransactionOrNeg start_transaction_id,
                               TransactionOrNeg end_transaction_id) const {

  boost::shared_ptr<IndexPage> page = load_page(col, row);
  return page->search_by_location(col, row, start_transaction_id, end_transaction_id);
}


// --------------------------------------------------------------------
//                             PAGED INDEX
// --------------------------------------------------------------------

/// Open an existing index from a file on disk.
PagedIndex::PagedIndex(boost::shared_ptr<PageGeneratorFactory> page_gen_factory,
                       uint32 page_width, uint32 page_height, uint32 default_cache_size)
  : m_page_gen_factory(page_gen_factory), m_page_width(page_width), m_page_height(page_height),
    m_default_cache_size(default_cache_size) {}

/// Open an existing index from a file on disk.
PagedIndex::PagedIndex(uint32 page_width, uint32 page_height, uint32 default_cache_size)
  : m_page_width(page_width), m_page_height(page_height),
    m_default_cache_size(default_cache_size) {}

void PagedIndex::sync() {
  for (unsigned i = 0; i < m_levels.size(); ++i) {
    m_levels[i]->sync();
  }
}

// ----------------------- READ/WRITE REQUESTS  ----------------------

uint64 PagedIndex::page_id(uint32 col, uint32 row, uint32 level) const {
  VW_ASSERT(level < m_levels.size(), ArgumentErr() << "Cannot request page_id for a level that doesn't exist");
  // Use the upper 32 bits for level, and the bottom for level-page-id
  return m_levels[level]->page_id(col, row);
}

/// Grab an IndexPage.  Useful if you want to serialize it by hand to disk.
boost::shared_ptr<IndexPage> PagedIndex::page_request(uint32 col, uint32 row, uint32 level) const {
  if (level >= m_levels.size())
    vw_throw(TileNotFoundErr() << "Requested page at level " << level << " was > than the max index (" << m_levels.size()-1 << ").");
  return m_levels[level]->get_page(col, row);
}

IndexRecord PagedIndex::read_request(uint32 col, uint32 row, uint32 level,
                                     TransactionOrNeg transaction_id, bool exact_transaction_match) {
  if (level >= m_levels.size())
    vw_throw(TileNotFoundErr() << "Requested tile at level " << level << " was greater than the max level (" << m_levels.size() << ").");

  IndexRecord rec = m_levels[level]->get(col, row,  transaction_id, exact_transaction_match);
  if (rec.filetype() == "default_to_index")
      rec.set_filetype(this->tile_filetype());
  return rec;
}

void PagedIndex::write_update(TileHeader const& header, IndexRecord const& record) {
  // First, we check to make sure we have a sufficient number of
  // levels to save the requested data.  If not, we grow the levels
  // vector to the correct size.

  for (uint32 level = m_levels.size(); level <= header.level(); ++level) {
    boost::shared_ptr<IndexLevel> new_level(
        new IndexLevel(m_page_gen_factory, level, m_page_width, m_page_height, m_default_cache_size) );
    m_levels.push_back(new_level);
  }

  m_levels[header.level()]->set(header, record);
}



// ----------------------- PROPERTIES  ----------------------

/// Returns a list of valid tiles that match this level, region, and
/// range of transaction_id's.  Returns a list of TileHeaders with
/// col/row/level and transaction_id of the most recent tile at each
/// valid location.  Note: there may be other tiles in the transaction
/// range at this col/row/level, but search_by_region() only returns the
/// first one.
std::list<TileHeader>
PagedIndex::search_by_region(uint32 level, BBox2i const& region,
                             TransactionOrNeg start_transaction_id, TransactionOrNeg end_transaction_id) const {

  // If the level does not exist, we return an empty list.
  if (level >= m_levels.size())
    return std::list<TileHeader>();

  // Otherwise, we delegate to the search_by_region() method for that level.
  return m_levels[level]->search_by_region(region, start_transaction_id, end_transaction_id);
}

std::list<TileHeader>
PagedIndex::search_by_location(uint32 col, uint32 row, uint32 level,
                               TransactionOrNeg start_transaction_id, TransactionOrNeg end_transaction_id) const {

  // If the level does not exist, we return an empty list.
  if (level >= m_levels.size())
    return std::list<TileHeader>();

  // Otherwise, we delegate to the search_by_location() method for that level.
  return m_levels[level]->search_by_location(col, row, start_transaction_id, end_transaction_id);
}
