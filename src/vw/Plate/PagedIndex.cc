// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Plate/PagedIndex.h>
#include <vw/Plate/Exception.h>
#include <vw/Core/Debugging.h>

#include <boost/foreach.hpp>

using namespace vw;
using namespace vw::platefile;

#define WHEREAMI (vw::vw_out(VerboseDebugMessage, "platefile.pagedindex") << VW_CURRENT_FUNCTION << ": ")


// --------------------------------------------------------------------
//                             INDEX LEVEL
// --------------------------------------------------------------------

IndexLevel::IndexLevel(boost::shared_ptr<PageGeneratorFactory> page_gen_factory,
                       int level, int page_width, int page_height, int cache_size)
  : m_page_gen_factory(page_gen_factory), m_level(level),
    m_page_width(page_width), m_page_height(page_height), m_cache(cache_size) {

  int tiles_per_side = 1 << level;
  m_horizontal_pages = static_cast<int>(ceil(float(tiles_per_side) / float(page_width)));
  m_vertical_pages   = static_cast<int>(ceil(float(tiles_per_side) / float(page_height)));

  // Create space for cache handles.  The actual generators are not
  // created until they are needed (because they take enough memory
  // that it's not efficient to allocate the generators ahead of
  // time).  Actual allocation is done automatically by calling
  // load_cache_handle().
  int pages = m_horizontal_pages * m_vertical_pages;
  m_cache_handles.resize(pages);
}

boost::shared_ptr<IndexPage> IndexLevel::fetch_page(int i, int j) const {
  Mutex::Lock lock(m_cache_mutex);

  size_t idx = j*m_horizontal_pages + i;

  // We may need to actually create the page's cache handle if it
  // hasn't been created already.
  if ( !m_cache_handles[idx].attached() ) {

    WHEREAMI << "Creating cache generator for page " << i << " " << j << " @ " << m_level << "\n";

    boost::shared_ptr<PageGeneratorBase> generator =
      m_page_gen_factory->create(m_level, i * m_page_width, j * m_page_height,
                                 m_page_width, m_page_height);
      m_cache_handles[idx] = m_cache.insert( generator );
  }

  return m_cache_handles[idx];
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

/// Grab an IndexPage.  Useful if you want to serialize it by hand to disk.
boost::shared_ptr<IndexPage> IndexLevel::get_page(int col, int row) const {

  VW_ASSERT( col >= 0 && row >= 0 && col < pow(2,m_level) && row < pow(2,m_level),
             TileNotFoundErr() << "IndexLevel::get_page() failed.  Invalid index [ "
             << col << " " << row << " @ level " << m_level << "]" );

  int32 level_col = col / m_page_width;
  int32 level_row = row / m_page_height;

  WHEREAMI << "(" << level_col << " " << level_row << " @ " << m_level << ")\n";

  // Access the page.  This will load it into memory if necessary.
  return fetch_page(level_col, level_row);
}


/// Fetch the value of an index node at this level.
IndexRecord IndexLevel::get(int32 col, int32 row, TransactionOrNeg transaction_id, bool exact_match) const {

  VW_ASSERT( col >= 0 && row >= 0 && col < pow(2,m_level) && row < pow(2,m_level),
             TileNotFoundErr() << "IndexLevel::get() failed.  Invalid index [ "
             << col << " " << row << " @ level " << m_level << "]" );

  int32 level_col = col / m_page_width;
  int32 level_row = row / m_page_height;

  WHEREAMI << "(" << level_col << " " << level_row << " @ " << m_level << ")\n";

  // Access the page.  This will load it into memory if necessary.
  boost::shared_ptr<IndexPage> page = fetch_page(level_col, level_row);
  return page->get(col, row, transaction_id, exact_match);
}

/// Set the value of an index node at this level.
void IndexLevel::set(TileHeader const& header, IndexRecord const& rec) {

  VW_ASSERT( header.col() >= 0 && header.row() >= 0 &&
             header.col() < pow(2,m_level) && header.row() < pow(2,m_level),
             TileNotFoundErr() << "IndexLevel::set() failed.  Invalid index [ "
             << header.col() << " " << header.row() << " @ level " << m_level << "]" );

  int32 level_col = header.col() / m_page_width;
  int32 level_row = header.row() / m_page_height;

  WHEREAMI << "(" << level_col << " " << level_row << " @ " << m_level << ")\n";

  // Access the page.  This will load it into memory if necessary.
  boost::shared_ptr<IndexPage> page = fetch_page(level_col, level_row);
  page->set(header, rec);
}

/// Returns a list of valid tiles at this level.
std::list<TileHeader>
IndexLevel::search_by_region(BBox2i const& region,
                             TransactionOrNeg start_transaction_id,
                             TransactionOrNeg end_transaction_id,
                             uint32 min_num_matches,
                             bool fetch_one_additional_entry) const {

  // Start by computing the search range in pages based on the requested region.
  int32 min_level_col = region.min().x() / m_page_width;
  int32 min_level_row = region.min().y() / m_page_height;

  int32 max_level_col = static_cast<int32>(ceilf(float(region.max().x()) / float(m_page_width)));
  int32 max_level_row = static_cast<int32>(ceilf(float(region.max().y()) / float(m_page_height)));

  WHEREAMI << "[" << min_level_col << " " << min_level_row << "]"
           << " to [" << max_level_col
           << " " << max_level_row << "]\n";

  // Iterate over the pages that overlap with the region of interest.
  std::list<TileHeader> result;
  for (int32 level_row = min_level_row; level_row < max_level_row; ++level_row) {
    for (int32 level_col = min_level_col; level_col < max_level_col; ++level_col) {
      boost::shared_ptr<IndexPage> page = fetch_page(level_col, level_row);

      // Accumulate valid tiles that overlap with region from this IndexPage.
      std::list<TileHeader> sub_result = page->search_by_region(region,
                                                                start_transaction_id,
                                                                end_transaction_id,
                                                                min_num_matches,
                                                                fetch_one_additional_entry);
      result.splice(result.end(), sub_result);
    }
  }
  return result;
}

/// Fetch the value of an index node at this level.
std::list<TileHeader>
IndexLevel::search_by_location(int32 col, int32 row,
                               TransactionOrNeg start_transaction_id,
                               TransactionOrNeg end_transaction_id,
                               bool fetch_one_additional_entry) const {

  if (col < 0 || row < 0 || col >= 1<<m_level || row >= 1<<m_level)
    return std::list<TileHeader>();

  int32 level_col = col / m_page_width;
  int32 level_row = row / m_page_height;

  // Access the page.  This will load it into memory if necessary.
  boost::shared_ptr<IndexPage> page = fetch_page(level_col, level_row);
  return page->search_by_location(col, row, start_transaction_id,
                                  end_transaction_id, fetch_one_additional_entry);
}


// --------------------------------------------------------------------
//                             PAGED INDEX
// --------------------------------------------------------------------

/// Create a new, empty index.
PagedIndex::PagedIndex(boost::shared_ptr<PageGeneratorFactory> page_gen_factory,
                       IndexHeader /*new_index_info*/,
                       int page_width, int page_height, int default_cache_size)
  : m_page_gen_factory(page_gen_factory), m_page_width(page_width), m_page_height(page_height),
    m_default_cache_size(default_cache_size) {}

/// Open an existing index from a file on disk.
PagedIndex::PagedIndex(boost::shared_ptr<PageGeneratorFactory> page_gen_factory,
                       int page_width, int page_height, int default_cache_size)
  : m_page_gen_factory(page_gen_factory), m_page_width(page_width), m_page_height(page_height),
    m_default_cache_size(default_cache_size) {}

/// Open an existing index from a file on disk.
PagedIndex::PagedIndex(int page_width, int page_height, int default_cache_size)
  : m_page_width(page_width), m_page_height(page_height),
    m_default_cache_size(default_cache_size) {}

void PagedIndex::sync() {
  for (unsigned i = 0; i < m_levels.size(); ++i) {
    m_levels[i]->sync();
  }
}

// ----------------------- READ/WRITE REQUESTS  ----------------------

/// Grab an IndexPage.  Useful if you want to serialize it by hand to disk.
boost::shared_ptr<IndexPage> PagedIndex::page_request(int col, int row, int level) const {
  if (level < 0 || level >= int(m_levels.size()))
    vw_throw(TileNotFoundErr() << "Requested page at " << level
             << " was greater than the max level (" << m_levels.size() << ").");
  return m_levels[level]->get_page(col, row);
}

IndexRecord PagedIndex::read_request(int col, int row, int level,
                                     TransactionOrNeg transaction_id, bool exact_transaction_match) {
  if (level < 0 || level >= int(m_levels.size()))
    vw_throw(TileNotFoundErr() << "Requested tile at level " << level
             << " was greater than the max level (" << m_levels.size() << ").");

  IndexRecord rec = m_levels[level]->get(col, row,  transaction_id, exact_transaction_match);
  if (rec.filetype() == "default_to_index")
      rec.set_filetype(this->tile_filetype());
  return rec;
}

void PagedIndex::write_update(TileHeader const& header, IndexRecord const& record) {
  // First, we check to make sure we have a sufficient number of
  // levels to save the requested data.  If not, we grow the levels
  // vector to the correct size.

  for (int level = boost::numeric_cast<int>(m_levels.size()); level <= header.level(); ++level) {
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
PagedIndex::search_by_region(int level, BBox2i const& region,
                             TransactionOrNeg start_transaction_id, TransactionOrNeg end_transaction_id,
                             uint32 min_num_matches,
                             bool fetch_one_additional_entry) const {

  // If the level does not exist, we return an empty list.
  if (level < 0 || level >= int(m_levels.size()))
    return std::list<TileHeader>();

  // Otherwise, we delegate to the search_by_region() method for that level.
  return m_levels[level]->search_by_region(region,
                                           start_transaction_id, end_transaction_id,
                                           min_num_matches, fetch_one_additional_entry);
}

std::list<TileHeader>
PagedIndex::search_by_location(int col, int row, int level,
                               TransactionOrNeg start_transaction_id, TransactionOrNeg end_transaction_id,
                               bool fetch_one_additional_entry = false) const {

  // If the level does not exist, we return an empty list.
  if (level < 0 || level >= int(m_levels.size()))
    return std::list<TileHeader>();

  // Otherwise, we delegate to the search_by_location() method for that level.
  return m_levels[level]->search_by_location(col, row,
                                             start_transaction_id, end_transaction_id,
                                             fetch_one_additional_entry);

}
