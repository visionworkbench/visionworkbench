// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/Plate/PagedIndex.h>
#include <vw/Core/ProgressCallback.h>

#include <boost/regex.hpp>

using namespace vw;

// --------------------------------------------------------------------
//                             INDEX LEVEL
// --------------------------------------------------------------------

vw::platefile::IndexLevel::IndexLevel(std::string plate_filename, int level, 
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
      filename << plate_filename 
               << "/index"
               << "/" << level 
               << "/" << (j * page_height) 
               << "/" << (i * page_width);
      boost::shared_ptr<IndexPageGenerator> generator( new IndexPageGenerator(filename.str(), 
                                                                              level, 
                                                                              i * m_page_width,
                                                                              j * m_page_height,
                                                                              page_width, 
                                                                              page_height) );
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
             TileNotFoundErr() << "IndexLevel::get() failed.  Invalid index [ " 
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
             TileNotFoundErr() << "IndexLevel::set() failed.  Invalid index [ " 
             << col << " " << row << " @ level " << m_level << "]" );
  
  int32 level_col = col / m_page_width;
  int32 level_row = row / m_page_height;
  int32 page_col = col % m_page_width;
  int32 page_row = row % m_page_height;
  
  // Access the page.  This will load it into memory if necessary.
  boost::shared_ptr<IndexPage> page = m_cache_handles[level_row*m_horizontal_pages + level_col];
  page->set(rec, page_col, page_row, transaction_id);

}

/// Returns a list of valid tiles at this level.
std::list<vw::platefile::TileHeader> vw::platefile::IndexLevel::valid_tiles(int transaction_id, 
                                                                            BBox2i const& region) const {
  
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
      std::list<TileHeader> sub_result = page->valid_tiles(transaction_id, region, false);
      result.splice(result.end(), sub_result);
    }
  }

  return result;
}


// --------------------------------------------------------------------
//                             PAGED INDEX
// --------------------------------------------------------------------

/// Create a new, empty index.
vw::platefile::PagedIndex::PagedIndex(std::string plate_filename, IndexHeader new_index_info,
                                      int page_width, int page_height, int default_cache_size) :
  LocalIndex(plate_filename, new_index_info),       // superclass constructor
  m_page_width(page_width), m_page_height(page_height), 
  m_default_cache_size(default_cache_size) {
  
  std::string base_index_path = plate_filename + "/index";
  if (!fs::exists(base_index_path))
    fs::create_directory(base_index_path);
}

/// Open an existing index from a file on disk.
vw::platefile::PagedIndex::PagedIndex(std::string plate_filename, 
                                      int page_width, int page_height, int default_cache_size) :
  LocalIndex(plate_filename),
  m_page_width(page_width), m_page_height(page_height), 
  m_default_cache_size(default_cache_size) {

  for (int level = 0; level < m_header.num_levels(); ++level) { 
    boost::shared_ptr<IndexLevel> new_level( new IndexLevel(m_plate_filename, level, 
                                                            m_page_width, m_page_height, 
                                                            m_default_cache_size) );
    m_levels.push_back(new_level);
  }
}

// Load index entries by iterating through TileHeaders saved in the
// blob file.  This function essentially rebuilds an index in memory
// using entries that had been previously saved to disk.
void vw::platefile::PagedIndex::rebuild_index(std::string plate_filename) {

  std::cout << "\tRebuilding index: " << plate_filename <<"\n";

  std::vector<std::string> blob_files = this->blob_filenames();
  for (unsigned int i = 0; i < blob_files.size(); ++i) {
    // this->log() << "Loading index entries from blob file: "
    //             << m_plate_filename << "/" << blob_files[i] << "\n";
    
    TerminalProgressCallback tpc(InfoMessage, "\t    " + blob_files[i] + " : ");
    tpc.report_progress(0);
    
    // Extract the current blob id as an integer.
    boost::regex re;
    re.assign("(plate_)(\\d+)(\\.blob)", boost::regex_constants::icase);
    boost::cmatch matches;
    boost::regex_match(blob_files[i].c_str(), matches, re);
    if (matches.size() != 4)
      vw_throw(IOErr() << "LocalIndex::load_index() -- could not parse blob number from blob filename.");
    std::string blob_id_str(matches[2].first, matches[2].second);
    int current_blob_id = atoi(blob_id_str.c_str());
      
    Blob blob(m_plate_filename + "/" + blob_files[i], true);
    Blob::iterator iter = blob.begin();
    while (iter != blob.end()) {
      TileHeader hdr = *iter;
      IndexRecord rec;
      rec.set_blob_id(current_blob_id);
      rec.set_blob_offset(iter.current_base_offset());
      rec.set_status(INDEX_RECORD_VALID);
      this->commit_record(rec, hdr.col(), hdr.row(), hdr.level(), hdr.transaction_id());
      tpc.report_progress(float(iter.current_base_offset()) / blob.size());
      ++iter;
    }
    tpc.report_finished();
  }
}


void vw::platefile::PagedIndex::commit_record(IndexRecord const& record, 
                                              int col, int row, 
                                              int level, int transaction_id) {

  // First, we check to make sure we have a sufficient number of
  // levels to save the requested data.  If not, we grow the levels
  // vector to the correct size.
  int starting_size = m_levels.size();
  while (m_levels.size() <= level) {
    boost::shared_ptr<IndexLevel> new_level( new IndexLevel(m_plate_filename,
                                                            m_levels.size(), 
                                                            m_page_width, m_page_height, 
                                                            m_default_cache_size) );
    m_levels.push_back(new_level);
  }
  
  if (m_levels.size() != starting_size) {
    m_header.set_num_levels(m_levels.size());
    this->save_index_file();
  }

  m_levels[level]->set(record, col, row, transaction_id);
}

// ----------------------- READ/WRITE REQUESTS  ----------------------

vw::platefile::IndexRecord vw::platefile::PagedIndex::read_request(int col, int row, int level, 
                                 int transaction_id, bool exact_transaction_match) {
  if (level < 0 || level >= m_levels.size())
    vw_throw(TileNotFoundErr() << "Requested tile at level " << level 
             << " was greater than the max level (" << m_levels.size() << ").");
  
  IndexRecord rec = m_levels[level]->get(col, row,  transaction_id, exact_transaction_match);
  return rec;
}

int vw::platefile::PagedIndex::write_request(int size) {
  return m_blob_manager->request_lock(size);
}

void vw::platefile::PagedIndex::write_complete(TileHeader const& header, 
                                               IndexRecord const& record) {
  m_blob_manager->release_lock(record.blob_id(), record.blob_offset());
  this->commit_record(record, header.col(), header.row(), 
                      header.level(), header.transaction_id());
}
  
// ----------------------- PROPERTIES  ----------------------

/// Returns a list of tile headers for any valid tiles that exist
/// at a the specified level and transaction_id.  The
/// transaction_id is treated the same as it would be for
/// read_request() above.  The region specifies a tile range of
/// interest.
std::list<vw::platefile::TileHeader> vw::platefile::PagedIndex::valid_tiles(int level, 
                                                                            int transaction_id, 
                                                                            BBox2i const& region) const {
  VW_ASSERT(level >= 0 && level < m_levels.size(),
            ArgumentErr() << "PagedIndex::valid_tiles() failed.  "
            << "Requested tiles from a level that does not exist.");
  return m_levels[level]->valid_tiles(transaction_id, region);
}

int32 vw::platefile::PagedIndex::num_levels() const {
  return m_levels.size();
}

// --------------------- UTILITIES ------------------------

/// Iterate over all nodes in a tree, calling func for each
/// location.  Note: this will only be implemented for local
/// indexes.  This function will throw an error if called on a
/// remote index.
void vw::platefile::PagedIndex::map(boost::shared_ptr<TreeMapFunc> func) { 
  vw_throw(NoImplErr() << "Index::map() not implemented for this index type.");
}
