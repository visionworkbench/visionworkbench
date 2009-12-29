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
  m_level(level), m_cache(cache_size) {
  int tiles_per_side = pow(2,level);
  int horizontal_pages = ceil(float(tiles_per_side) / page_width);
  int vertical_pages = ceil(float(tiles_per_side) / page_height);
  int pages = horizontal_pages * vertical_pages;
  
  // Create the cache handles
  m_cache_handles.resize(pages);
  m_cache_generators.resize(pages);
  for (int j = 0; j < vertical_pages; ++j) {
    for (int i = 0; i < horizontal_pages; ++i) {
      
      // Generate a filename.
      std::ostringstream filename;
      filename << base_path 
               << "/" << level 
               << "/" << (j * page_height) 
               << "/" << (i * page_width);
      boost::shared_ptr<IndexPageGenerator> generator( new IndexPageGenerator(filename.str(), page_width, page_height) );
      m_cache_generators[j*horizontal_pages + i] = generator;
      m_cache_handles[j*horizontal_pages + i] = m_cache.insert( *generator );
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



// --------------------------------------------------------------------
//                             PAGED INDEX
// --------------------------------------------------------------------


// ----------------------- READ/WRITE REQUESTS  ----------------------

vw::platefile::IndexRecord vw::platefile::PagedIndex::read_request(int col, int row, int depth, 
                                 int transaction_id, bool exact_transaction_match) {}

int vw::platefile::PagedIndex::write_request(int size) {}

void vw::platefile::PagedIndex::write_complete(TileHeader const& header, IndexRecord const& record) {}

  
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
