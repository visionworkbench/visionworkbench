// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
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
                                      int level, int page_width, int page_height, int cache_size) : 
  m_level(level), m_page_width(page_width), m_page_height(page_height), m_cache(cache_size) {
  int tiles_per_side = pow(2,level);
  m_horizontal_pages = ceil(float(tiles_per_side) / page_width);
  m_vertical_pages = ceil(float(tiles_per_side) / page_height);
  int pages = m_horizontal_pages * m_vertical_pages;

  // Create the cache handles
  m_cache_handles.resize(pages);
  m_cache_generators.resize(pages);
  for (int j = 0; j < m_vertical_pages; ++j) {
    for (int i = 0; i < m_horizontal_pages; ++i) {
      boost::shared_ptr<IndexPageGenerator> generator = page_gen_factory->create(level, 
                                                                                 i * m_page_width,
                                                                                 j * m_page_height,
                                                                                 page_width, 
                                                                                 page_height) ;
      m_cache_generators[j*m_horizontal_pages + i] = generator;
      m_cache_handles[j*m_horizontal_pages + i] = m_cache.insert( *generator );
    }
  }
}

vw::platefile::IndexLevel::~IndexLevel() {

  // We need to free the cache handles first before other things
  // (especially the generators) get unallocated.
  for (int i = 0; i < m_cache_handles.size(); ++i) {
    m_cache_handles[i].reset();
  }

}

void vw::platefile::IndexLevel::sync() {

  // We need to free the cache handles first before the data gets
  // saved.  
  for (int i = 0; i < m_cache_handles.size(); ++i) {
    m_cache_handles[i].reset();
  }

}

/// Grab an IndexPage.  Useful if you want to serialize it by hand to disk.
boost::shared_ptr<vw::platefile::IndexPage> vw::platefile::IndexLevel::get_page(int col, int row) const {
  VW_ASSERT( col >= 0 && row >= 0 && col < pow(2,m_level) && row < pow(2,m_level), 
             TileNotFoundErr() << "IndexLevel::get_page() failed.  Invalid index [ " 
             << col << " " << row << " @ level " << m_level << "]" );
  
  int32 level_col = col / m_page_width;
  int32 level_row = row / m_page_height;
  
  // Access the page.  This will load it into memory if necessary.
  boost::shared_ptr<IndexPage> page = m_cache_handles[level_row*m_horizontal_pages + level_col];
  return page;
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
  
  // Access the page.  This will load it into memory if necessary.
  boost::shared_ptr<IndexPage> page = m_cache_handles[level_row*m_horizontal_pages + level_col];
  return page->get(col, row, transaction_id, exact_match);
}

/// Fetch the value of an index node at this level.
vw::platefile::IndexLevel::multi_value_type 
vw::platefile::IndexLevel::multi_get(int32 col, 
                                     int32 row, 
                                     int32 start_transaction_id,
                                     int32 end_transaction_id) const {

  VW_ASSERT( col >= 0 && row >= 0 && col < pow(2,m_level) && row < pow(2,m_level), 
             TileNotFoundErr() << "IndexLevel::multi_get() failed.  Invalid index [ " 
             << col << " " << row << " @ level " << m_level << "]" );
  
  int32 level_col = col / m_page_width;
  int32 level_row = row / m_page_height;
  
  // Access the page.  This will load it into memory if necessary.
  boost::shared_ptr<IndexPage> page = m_cache_handles[level_row*m_horizontal_pages + level_col];
  return page->multi_get(col, row, start_transaction_id, end_transaction_id);
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
  
  // Access the page.  This will load it into memory if necessary.
  boost::shared_ptr<IndexPage> page = m_cache_handles[level_row*m_horizontal_pages + level_col];
  page->set(header, rec);
}

/// Returns a list of valid tiles at this level.
std::list<vw::platefile::TileHeader> 
vw::platefile::IndexLevel::valid_tiles(BBox2i const& region,
                                       int start_transaction_id, 
                                       int end_transaction_id, 
                                       int min_num_matches) const {
  
  // Start by computing the search range in pages based on the requested region. 
  int32 min_level_col = region.min().x() / m_page_width;
  int32 min_level_row = region.min().y() / m_page_height;
  
  int32 max_level_col = ceilf(float(region.max().x()) / m_page_width);
  int32 max_level_row = ceilf(float(region.max().y()) / m_page_height);

  // Iterate over the pages that overlap with the region of interest.
  std::list<TileHeader> result;
  for (int32 level_row = min_level_row; level_row < max_level_row; ++level_row) {
    for (int32 level_col = min_level_col; level_col < max_level_col; ++level_col) {
      boost::shared_ptr<IndexPage> page = m_cache_handles[level_row*m_horizontal_pages + level_col];

      // Accumulate valid tiles that overlap with region from this IndexPage.
      std::list<TileHeader> sub_result = page->valid_tiles(region,
                                                           start_transaction_id, 
                                                           end_transaction_id, 
                                                           min_num_matches);
      result.splice(result.end(), sub_result);
    }
  }

  return result;
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
  for (int i = 0; i < m_levels.size(); ++i) {
    m_levels[i]->sync();
  }

}

// ----------------------- READ/WRITE REQUESTS  ----------------------

/// Grab an IndexPage.  Useful if you want to serialize it by hand to disk.
boost::shared_ptr<vw::platefile::IndexPage> vw::platefile::PagedIndex::page_request(int col, int row, int level) const {
  if (level < 0 || level >= m_levels.size())
    vw_throw(TileNotFoundErr() << "Requested page at " << level 
             << " was greater than the max level (" << m_levels.size() << ").");
  return m_levels[level]->get_page(col, row);
}

vw::platefile::IndexRecord vw::platefile::PagedIndex::read_request(int col, int row, int level, 
                                 int transaction_id, bool exact_transaction_match) {
  if (level < 0 || level >= m_levels.size())
    vw_throw(TileNotFoundErr() << "Requested tile at level " << level 
             << " was greater than the max level (" << m_levels.size() << ").");
  
  IndexRecord rec = m_levels[level]->get(col, row,  transaction_id, exact_transaction_match);
  return rec;
}

vw::platefile::PagedIndex::multi_value_type
vw::platefile::PagedIndex::multi_read_request(int col, int row, int level, 
                                              int start_transaction_id, 
                                              int end_transaction_id) {
  if (level < 0 || level >= m_levels.size())
    vw_throw(TileNotFoundErr() << "Requested tile at level " << level 
             << " was greater than the max level (" << m_levels.size() << ").");
  
  return m_levels[level]->multi_get(col, row,  start_transaction_id, end_transaction_id);
}


void vw::platefile::PagedIndex::write_update(TileHeader const& header, IndexRecord const& record) {
  // First, we check to make sure we have a sufficient number of
  // levels to save the requested data.  If not, we grow the levels
  // vector to the correct size.
  while (m_levels.size() <= header.level()) {
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
/// range at this col/row/level, but valid_tiles() only returns the
/// first one.
std::list<vw::platefile::TileHeader> 
vw::platefile::PagedIndex::valid_tiles(int level, BBox2i const& region,
                                       int start_transaction_id,
                                       int end_transaction_id,
                                       int min_num_matches) const {
  
  // If the level does not exist, we return an empty list.
  if (level < 0 || level >= m_levels.size())
    return std::list<TileHeader>();

  // Otherwise, we delegate to the valid_tiles() method for that level.
  return m_levels[level]->valid_tiles(region, 
                                      start_transaction_id, end_transaction_id,
                                      min_num_matches);
}

