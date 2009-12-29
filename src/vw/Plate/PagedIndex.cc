// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/Plate/PagedIndex.h>
using namespace vw;

// --------------------------------------------------------------------
//                             INDEX LEVEL
// --------------------------------------------------------------------

vw::platefile::IndexLevel::IndexLevel(std::string base_path, int level, 
                                      int page_width, int page_height, int cache_size) : 
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
      
      // Generate a filename.
      std::ostringstream filename;
      filename << base_path 
               << "/" << level 
               << "/" << (j * page_height) 
               << "/" << (i * page_width);
      boost::shared_ptr<IndexPageGenerator> generator( new IndexPageGenerator(filename.str(), page_width, page_height) );
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

/// Fetch the value of an index node at this level.
vw::platefile::IndexRecord vw::platefile::IndexLevel::get(int32 col, 
                                                          int32 row, 
                                                          int32 transaction_id,
                                                          bool exact_match) const {

  VW_ASSERT( col >= 0 && row >= 0 && col < pow(2,m_level) && row < pow(2,m_level), 
             ArgumentErr() << "IndexLevel::get() failed.  Invalid index [ " 
             << col << " " << row << " @ level " << m_level << "]" );
  
  int32 level_col = col / m_page_width;
  int32 level_row = row / m_page_height;
  int32 page_col = col % m_page_width;
  int32 page_row = row % m_page_height;
  
  // Access the page.  This will load it into memory if necessary.
  boost::shared_ptr<IndexPage> page = m_cache_handles[level_row*m_horizontal_pages + level_col];
  return page->get(page_col, page_row, transaction_id, exact_match);
}

/// Set the value of an index node at this level.
void vw::platefile::IndexLevel::set(vw::platefile::IndexRecord const& rec, 
                                    int32 col, int32 row, int32 transaction_id) {

  VW_ASSERT( col >= 0 && row >= 0 && col < pow(2,m_level) && row < pow(2,m_level), 
             ArgumentErr() << "IndexLevel::set() failed.  Invalid index [ " 
             << col << " " << row << " @ level " << m_level << "]" );
  
  int32 level_col = col / m_page_width;
  int32 level_row = row / m_page_height;
  int32 page_col = col % m_page_width;
  int32 page_row = row % m_page_height;
  
  // Access the page.  This will load it into memory if necessary.
  boost::shared_ptr<IndexPage> page = m_cache_handles[level_row*m_horizontal_pages + level_col];
  page->set(rec, page_col, page_row, transaction_id);

}



// --------------------------------------------------------------------
//                             PAGED INDEX
// --------------------------------------------------------------------

vw::platefile::PagedIndex::PagedIndex(std::string base_path, int page_width, 
                                      int page_height, int default_cache_size) :
  m_base_path(base_path), m_page_width(page_width), 
  m_page_height(page_height), m_default_cache_size(default_cache_size),
  m_blob_manager(boost::shared_ptr<BlobManager>( new BlobManager() )) {}

// ----------------------- READ/WRITE REQUESTS  ----------------------

vw::platefile::IndexRecord vw::platefile::PagedIndex::read_request(int col, int row, int depth, 
                                 int transaction_id, bool exact_transaction_match) {
  Mutex::Lock lock(m_mutex);
  if (depth < 0 || depth >= m_levels.size())
    vw_throw(TileNotFoundErr() << "Requested tile at level " << depth 
             << " was greater than the max level (" << m_levels.size() << ").");
  
  IndexRecord rec = m_levels[depth]->get(col, row,  transaction_id, exact_transaction_match);
  return rec;
}

int vw::platefile::PagedIndex::write_request(int size) {
  return m_blob_manager->request_lock(size);
}

void vw::platefile::PagedIndex::write_complete(TileHeader const& header, IndexRecord const& record) {
  m_blob_manager->release_lock(record.blob_id(), record.blob_offset());

  Mutex::Lock lock(m_mutex);
  
  // First, we check to make sure we have a sufficient number of
  // levels to save the requested data.  If not, we grow the levels
  // vector to the correct size.
  while (m_levels.size() <= header.depth()) {
    boost::shared_ptr<IndexLevel> new_level( new IndexLevel(m_base_path, m_levels.size(), 
                                                            m_page_width, m_page_height, 
                                                            m_default_cache_size) );
    m_levels.push_back(new_level);
  }

  m_levels[header.depth()]->set(record, header.col(), header.row(), header.transaction_id());
}
  
// ----------------------- PROPERTIES  ----------------------

vw::platefile::IndexHeader vw::platefile::PagedIndex::index_header() const {}
int32 vw::platefile::PagedIndex::version() const {}
int32 vw::platefile::PagedIndex::max_depth() const {}
  
std::string vw::platefile::PagedIndex::platefile_name() const {}
  
int32 vw::platefile::PagedIndex::tile_size() const {}
std::string vw::platefile::PagedIndex::tile_filetype() const {}
  
PixelFormatEnum vw::platefile::PagedIndex::pixel_format() const {}
ChannelTypeEnum vw::platefile::PagedIndex::channel_type() const {}

// --------------------- TRANSACTIONS ------------------------

/// Clients are expected to make a transaction request whenever
/// they start a self-contained chunk of mosaicking work.  .
int32 vw::platefile::PagedIndex::transaction_request(std::string transaction_description,
                                  std::vector<TileHeader> const& tile_headers) {}
  
/// Called right before the beginning of the mipmapping pass
void vw::platefile::PagedIndex::root_complete(int32 transaction_id,
                           std::vector<TileHeader> const& tile_headers) {}
  
/// Once a chunk of work is complete, clients can "commit" their
/// work to the mosaic by issuding a transaction_complete method.
void vw::platefile::PagedIndex::transaction_complete(int32 transaction_id) {}

// If a transaction fails, we may need to clean up the mosaic.  
void vw::platefile::PagedIndex::transaction_failed(int32 transaction_id) {}

int32 vw::platefile::PagedIndex::transaction_cursor() {}

// --------------------- UTILITIES ------------------------

/// Iterate over all nodes in a tree, calling func for each
/// location.  Note: this will only be implemented for local
/// indexes.  This function will throw an error if called on a
/// remote index.
void vw::platefile::PagedIndex::map(boost::shared_ptr<TreeMapFunc> func) { 
  vw_throw(NoImplErr() << "Index::map() not implemented for this index type.");
}
