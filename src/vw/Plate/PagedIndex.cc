// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Core/ProgressCallback.h>
#include <vw/Plate/PagedIndex.h>
#include <vw/Plate/Exception.h>
#include <vw/Plate/Blob.h>

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>
namespace fs = boost::filesystem;

using namespace vw;


// --------------------------------------------------------------------
//                             INDEX LEVEL
// --------------------------------------------------------------------

vw::platefile::IndexLevel::IndexLevel(boost::shared_ptr<PageGeneratorFactory> page_gen_factory,
                                      int level, int page_width, int page_height, int cache_size)
    : m_page_gen_factory(page_gen_factory), m_level(level),
      m_page_width(page_width), m_page_height(page_height), m_cache(cache_size) {

  int tiles_per_side = pow(2,level);
  m_horizontal_pages = ceil(float(tiles_per_side) / page_width);
  m_vertical_pages = ceil(float(tiles_per_side) / page_height);

  // Create space for cache handles.  The actual generators are not
  // created until they are needed (because they take enough memory
  // that it's not efficient to allocate the generators ahead of
  // time).  Actual allocation is done automatically by calling
  // load_cache_handle().
  int pages = m_horizontal_pages * m_vertical_pages;
  m_cache_handles.resize(pages);
  m_cache_generators.resize(pages);

}

boost::shared_ptr<vw::platefile::IndexPage> vw::platefile::IndexLevel::fetch_page(int i, int j) const {
  Mutex::Lock lock(m_cache_mutex);

  // We may need to actually create the page's cache handle if it
  // hasn't been created already.
  if ( !(m_cache_generators[j*m_horizontal_pages + i])) {

    vw_out(DebugMessage, "platefile::PagedIndex") << "Creating cache generator for page " 
                                                  << i << " " << j << " @ " << m_level << "\n";
  
    boost::shared_ptr<IndexPageGenerator> generator =
      m_page_gen_factory->create(m_level, i * m_page_width, j * m_page_height, 
                                 m_page_width, m_page_height);
      m_cache_generators[j*m_horizontal_pages + i] = generator;
      m_cache_handles[j*m_horizontal_pages + i] = m_cache.insert( *generator );

  }

  return m_cache_handles[j*m_horizontal_pages + i];
}


vw::platefile::IndexLevel::~IndexLevel() {
  Mutex::Lock lock(m_cache_mutex);

  // We need to free the cache handles first before other things
  // (especially the generators) get unallocated.
  for (unsigned i = 0; i < m_cache_handles.size(); ++i) {
    if (m_cache_generators[i]) 
      m_cache_handles[i].reset();
  }

}

void vw::platefile::IndexLevel::sync() {
  
  // Write the index page to disk by calling it's sync() method.
  //
  // Note: the size of m_cache_handles and m_cache_generators does not
  // change once it's set in the constructor of IndexLevel, so it's
  // safe here to place m_cache_mutex _inside_ the loop.  This allows
  // synchronization operations to be interleaved with fetch()
  // operations, thus preventing the index_server from blocking when
  // it periodically syncs the cache to disk.
  //
  for (unsigned i = 0; i < m_cache_handles.size(); ++i) {
    Mutex::Lock lock(m_cache_mutex);
    if (m_cache_generators[i]) {
      boost::shared_ptr<IndexPage> page = m_cache_handles[i];
      page->sync();
    }
  }

}

/// Grab an IndexPage.  Useful if you want to serialize it by hand to disk.
boost::shared_ptr<vw::platefile::IndexPage> vw::platefile::IndexLevel::get_page(int col, int row) const {
  VW_ASSERT( col >= 0 && row >= 0 && col < pow(2,m_level) && row < pow(2,m_level), 
             TileNotFoundErr() << "IndexLevel::get_page() failed.  Invalid index [ " 
             << col << " " << row << " @ level " << m_level << "]" );
  
  int32 level_col = col / m_page_width;
  int32 level_row = row / m_page_height;

  vw_out(DebugMessage, "platefile::PagedIndex") << "IndexPage::get_page() called " 
                                                << level_col << " " << level_row 
                                                << " @ " << m_level << "\n";

  // Access the page.  This will load it into memory if necessary.
  return fetch_page(level_col, level_row);
}


/// Fetch the value of an index node at this level.
vw::platefile::IndexRecord vw::platefile::IndexLevel::get(int32 col, 
                                                          int32 row, 
                                                          int32 transaction_id,
                                                          bool exact_match) const {

  VW_ASSERT( col >= 0 && row >= 0 && col < pow(2,m_level) && row < pow(2,m_level), 
             TileNotFoundErr() << "IndexLevel::get() failed.  Invalid index [ " 
             << col << " " << row << " @ level " << m_level << "]" );
  
  int32 level_col = col / m_page_width;
  int32 level_row = row / m_page_height;

  vw_out(VerboseDebugMessage, "platefile::PagedIndex") << "IndexPage::get() called " 
                                                       << level_col << " " << level_row 
                                                       << " @ " << m_level << "\n";
  
  // Access the page.  This will load it into memory if necessary.
  boost::shared_ptr<IndexPage> page = fetch_page(level_col, level_row);
  return page->get(col, row, transaction_id, exact_match);
}

/// Set the value of an index node at this level.
void vw::platefile::IndexLevel::set(vw::platefile::TileHeader const& header, 
                                    vw::platefile::IndexRecord const& rec) {

  VW_ASSERT( header.col() >= 0 && header.row() >= 0 && 
             header.col() < pow(2,m_level) && header.row() < pow(2,m_level), 
             TileNotFoundErr() << "IndexLevel::set() failed.  Invalid index [ " 
             << header.col() << " " << header.row() << " @ level " << m_level << "]" );
  
  int32 level_col = header.col() / m_page_width;
  int32 level_row = header.row() / m_page_height;

  vw_out(VerboseDebugMessage, "platefile::PagedIndex") << "IndexPage::set() called " 
                                                       << level_col << " " << level_row 
                                                       << " @ " << m_level << "\n";
  
  // Access the page.  This will load it into memory if necessary.
  boost::shared_ptr<IndexPage> page = fetch_page(level_col, level_row);
  page->set(header, rec);
}

/// Returns a list of valid tiles at this level.
std::list<vw::platefile::TileHeader> 
vw::platefile::IndexLevel::search_by_region(BBox2i const& region,
                                            int start_transaction_id, 
                                            int end_transaction_id, 
                                            int min_num_matches,
                                            bool fetch_one_additional_entry) const {
  
  // Start by computing the search range in pages based on the requested region. 
  int32 min_level_col = region.min().x() / m_page_width;
  int32 min_level_row = region.min().y() / m_page_height;
  
  int32 max_level_col = ceilf(float(region.max().x()) / m_page_width);
  int32 max_level_row = ceilf(float(region.max().y()) / m_page_height);

  vw_out(VerboseDebugMessage, "platefile::PagedIndex") << "IndexPage::search_by_region() called " 
                                                       << "[" << min_level_col << " " 
                                                       << min_level_row << "]"
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
std::list<vw::platefile::TileHeader> 
vw::platefile::IndexLevel::search_by_location(int32 col, int32 row, 
                                              int32 start_transaction_id,
                                              int32 end_transaction_id, 
                                              bool fetch_one_additional_entry) const {

  VW_ASSERT( col >= 0 && row >= 0 && col < pow(2,m_level) && row < pow(2,m_level), 
             TileNotFoundErr() << "IndexLevel::read_headers() failed.  Invalid index [ " 
             << col << " " << row << " @ level " << m_level << "]" );
  
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
vw::platefile::PagedIndex::PagedIndex(boost::shared_ptr<PageGeneratorFactory> page_gen_factory,
                                      IndexHeader new_index_info,
                                      int page_width, int page_height, int default_cache_size) :
  m_page_gen_factory(page_gen_factory),
  m_page_width(page_width), m_page_height(page_height), 
  m_default_cache_size(default_cache_size) {}

/// Open an existing index from a file on disk.
vw::platefile::PagedIndex::PagedIndex(boost::shared_ptr<PageGeneratorFactory> page_gen_factory,
                                      int page_width, int page_height, int default_cache_size) :
  m_page_gen_factory(page_gen_factory),
  m_page_width(page_width), m_page_height(page_height), 
  m_default_cache_size(default_cache_size) {}

void vw::platefile::PagedIndex::sync() {

  // We need to free the cache handles first before other things
  // (especially the generators) get unallocated.
  for (unsigned i = 0; i < m_levels.size(); ++i) {
    m_levels[i]->sync();
  }

}

// ----------------------- READ/WRITE REQUESTS  ----------------------

/// Grab an IndexPage.  Useful if you want to serialize it by hand to disk.
boost::shared_ptr<vw::platefile::IndexPage> vw::platefile::PagedIndex::page_request(int col, int row, int level) const {
  if (level < 0 || level >= int(m_levels.size()))
    vw_throw(TileNotFoundErr() << "Requested page at " << level 
             << " was greater than the max level (" << m_levels.size() << ").");
  return m_levels[level]->get_page(col, row);
}

vw::platefile::IndexRecord vw::platefile::PagedIndex::read_request(int col, int row, int level, 
                                 int transaction_id, bool exact_transaction_match) {
  if (level < 0 || level >= int(m_levels.size()))
    vw_throw(TileNotFoundErr() << "Requested tile at level " << level 
             << " was greater than the max level (" << m_levels.size() << ").");
  
  IndexRecord rec = m_levels[level]->get(col, row,  transaction_id, exact_transaction_match);
  return rec;
}

void vw::platefile::PagedIndex::write_update(TileHeader const& header, IndexRecord const& record) {
  // First, we check to make sure we have a sufficient number of
  // levels to save the requested data.  If not, we grow the levels
  // vector to the correct size.
  while (int(m_levels.size()) <= header.level()) {
    boost::shared_ptr<IndexLevel> new_level( new IndexLevel(m_page_gen_factory,
                                                            m_levels.size(), 
                                                            m_page_width, m_page_height, 
                                                            m_default_cache_size) );
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
std::list<vw::platefile::TileHeader>
vw::platefile::PagedIndex::search_by_region(int level, BBox2i const& region,
                                            int start_transaction_id,
                                            int end_transaction_id,
                                            int min_num_matches, 
                                            bool fetch_one_additional_entry) const {
  
  // If the level does not exist, we return an empty list.
  if (level < 0 || level >= int(m_levels.size()))
    return std::list<TileHeader>();

  // Otherwise, we delegate to the search_by_region() method for that level.
  return m_levels[level]->search_by_region(region, 
                                           start_transaction_id, end_transaction_id,
                                           min_num_matches, fetch_one_additional_entry);
}

std::list<vw::platefile::TileHeader> 
vw::platefile::PagedIndex::search_by_location(int col, int row, int level, 
                                              int start_transaction_id, int end_transaction_id,
                                              bool fetch_one_additional_entry = false) const {

  // If the level does not exist, we return an empty list.
  if (level < 0 || level >= int(m_levels.size()))
    return std::list<TileHeader>();

  // Otherwise, we delegate to the search_by_location() method for that level.
  return m_levels[level]->search_by_location(col, row,
                                             start_transaction_id, end_transaction_id,
                                             fetch_one_additional_entry);

}

