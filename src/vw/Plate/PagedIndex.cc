// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/Core/ProgressCallback.h>
#include <vw/Plate/PagedIndex.h>
#include <vw/Plate/Exception.h>
#include <vw/Plate/Blob.h>

#include <boost/regex.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>
namespace fs = boost::filesystem;

using namespace vw;


// ----------------------------------------------------------------------
//                            INDEX PAGE
// ----------------------------------------------------------------------

vw::platefile::IndexPage::IndexPage(std::string filename, 
                                    int level, int base_col, int base_row, 
                                    int page_width, int page_height) : 
  m_filename(filename), m_level(level), m_base_col(base_col), m_base_row(base_row),
  m_page_width(page_width), m_page_height(page_height), m_needs_saving(false) {

  vw_out(VerboseDebugMessage, "platefile::index") << "Opening index page [ " 
                                                  << m_base_col << " " << m_base_row 
                                                  << " @ " << m_level << " ]\n";

  if (fs::exists(filename)) {
    this->deserialize();
  } else {
    m_sparse_table.resize(page_width*page_height);
  }
}

vw::platefile::IndexPage::~IndexPage() {
  vw_out(VerboseDebugMessage, "platefile::index") << "Closing index page [ " 
                                                  << m_base_col << " " << m_base_row 
                                                  << " @ " << m_level << " ]\n";
  this->sync();
}

// ----------------------- ACCESSORS  ----------------------

void vw::platefile::IndexPage::set(IndexRecord const& record, 
                                   int col, int row, int transaction_id) {

  // Basic bounds checking
  VW_ASSERT(col >= 0 && col < m_page_width && row >= 0 && row < m_page_height, 
            ArgumentErr() << "IndexPage::set() failed.  Invalid index [" 
            << col << " " << row << "]");

  // Mark this page is 'dirty' so that it gets saved to disk when
  // destroyed.
  m_needs_saving = true;

  std::pair<int32, IndexRecord> p(transaction_id, record);
  int elmnt = row*m_page_width + col;
  if (m_sparse_table.test(elmnt)) {

    // Add to existing entry.
    multi_value_type *entries = m_sparse_table[row*m_page_width + col].operator&();

    // We need to keep this list sorted in decreasing order of
    // transaction ID, do a simple insertion sort here.
    multi_value_type::iterator it = entries->begin();
    while (it != entries->end() && (*it).first >= transaction_id ) {

      // Handle the case where we replace an entry
      if ( (*it).first == transaction_id ) {
        (*it).second = record;
        return;
      }

      // Otherwise, we search forward in the list.
      ++it;
    }
    
    // If we reach this point, we are either at the end of the list,
    // or we have found the correct position.  Either way, we call
    // insert().
    entries->insert(it, p);

  } else {
    
    // Create a new entry
    multi_value_type l;
    l.push_front(p);
    m_sparse_table[elmnt] = l;
    
  }
}

/// Return the IndexRecord for a the given transaction_id at
/// this location.  By default this routine returns the record with
/// the greatest transaction id that is less than or equal to the
/// requested transaction_id.  
///
/// Exceptions to the above behavior:
///
///   - A transaction_id == -1 will return the most recent
///   transaction id.
/// 
///   - Setting exact_match to true forces an exact transaction_id
///   match.
///
vw::platefile::IndexRecord vw::platefile::IndexPage::get(int col, int row, 
                                                             int transaction_id, 
                                                             bool exact_match) const {
  
  // Basic bounds checking
  VW_ASSERT(col >= 0 && col < m_page_width && row >= 0 && row < m_page_height, 
            TileNotFoundErr() << "IndexPage::get() failed.  Invalid index [" 
            << col << " " << row << "]");

  // Interate over entries.
  multi_value_type const& entries = m_sparse_table[row*m_page_width + col];
  multi_value_type::const_iterator it = entries.begin();
  
  // A transaction ID of -1 indicates that we should return the most
  // recent tile (which is the first entry in the list, since it is
  // sorted from most recent to least recent), regardless of its
  // transaction id.
  if (transaction_id == -1 && it != entries.end())
    return (*it).second;
  
  // Otherwise, we search through the entries in the list, looking for
  // the requested t_id.  Note: this search is O(n), so it can be slow
  // if there are a lot of entries and the entry you are looking for
  // is near the end.  However, most pages will contain very few
  // entries, and for those with many entries (i.e. tiles near the
  // root of the mosaic), you will rarely search for old tiles.
  while (it != entries.end()) {
    if (exact_match) {
      if ((*it).first == transaction_id)
        return (*it).second;
    } else {
      if ((*it).first <= transaction_id)
        return (*it).second;
    }
    ++it;
  }
  
  // If we reach this point, then there are no entries before
  // the given transaction_id, so we return an empty (and invalid) record.
  vw_throw(TileNotFoundErr() << "Tiles exist at this location, " 
           << "but none before transaction_id = "  << transaction_id << "\n");
  return IndexRecord(); // never reached
}

vw::platefile::IndexPage::multi_value_type 
vw::platefile::IndexPage::multi_get(int col, int row, 
                                    int start_transaction_id, 
                                    int end_transaction_id) const {

  // Basic bounds checking
  VW_ASSERT(col >= 0 && col < m_page_width && row >= 0 && row < m_page_height, 
            TileNotFoundErr() << "IndexPage::multi_get() failed.  Invalid index [" 
            << col << " " << row << "]");

  // Check first to make sure that there are actually tiles at this location.
  if (!m_sparse_table.test(row*m_page_width + col)) 
    vw_throw(TileNotFoundErr() << "No tiles were found at this location.\n");
    
  // If there are, then we apply the transaction_id filters to select the requested ones.
  multi_value_type results;
  multi_value_type const& entries = m_sparse_table[row*m_page_width + col];
  multi_value_type::const_iterator it = entries.begin();
  while (it != entries.end() && it->first >= start_transaction_id) {
    //    std::cout << "Comparing: " << (it->first) << " and " << start_transaction_id << " <-> " << end_transaction_id << "\n";
    if (it->first >= start_transaction_id && it->first <= end_transaction_id) {
      results.push_back(*it);
    }
    ++it;
  }
  
  if (results.empty()) 
    vw_throw(TileNotFoundErr() << entries.size() << " tiles exist at this location, "
             << "but none match the transaction_id range you specified.\n");

  return results;
}

void vw::platefile::IndexPage::append_if_in_region( std::list<vw::platefile::TileHeader> &results, 
                                                    int transaction_id, IndexRecord const& rec, 
                                                    int col, int row,
                                                    BBox2i const& region) const {

  //  std::cout << "Appending record for level " << m_level << "  at " << col << " " << row << "\n";
        

  // Check to see if the tile is in the specified region.
  Vector2i loc( m_base_col + col, m_base_row + row);
  if ( region.contains( loc ) ) {
    TileHeader hdr;
    hdr.set_col( m_base_col + col );
    hdr.set_row( m_base_row + row );
    hdr.set_level(m_level);
    hdr.set_transaction_id(transaction_id);
    results.push_back(hdr);
  }
}

/// Returns a list of valid tiles in this IndexPage.  Returns a list
/// of TileHeaders with col/row/level and transaction_id of the most
/// recent tile at each valid location.  Note: there may be other
/// tiles in the transaction range at this col/row/level, but
/// valid_tiles() only returns the first one.
std::list<vw::platefile::TileHeader> 
vw::platefile::IndexPage::valid_tiles(vw::BBox2i const& region,
                                      int start_transaction_id,
                                      int end_transaction_id) const {
  std::list<TileHeader> results;

  for (int row = 0; row < m_page_height; ++row) {
    for (int col = 0; col < m_page_width; ++col) {
      if (m_sparse_table.test(row*m_page_width + col)) {

        // Iterate over entries.
        multi_value_type const& entries = m_sparse_table[row*m_page_width + col];
        multi_value_type::const_iterator it = entries.begin();

        // Search through the entries in the list, looking entries
        // that match the requested transaction_id range.  Note: this
        // search is O(n), so it can be slow if there are a lot of
        // entries and the entry you are looking for is near the end.
        // However, most pages will contain very few entries, and for
        // those with many entries (i.e. tiles near the root of the
        // mosaic), you will most likely be searching for recently
        // added tiles, which are sorted to the beginning.
        while (it != entries.end() && it->first >= start_transaction_id) {
          if ( it->first >= start_transaction_id && it->first <= end_transaction_id ) {
            append_if_in_region( results, (*it).first, (*it).second, col, row, region );
            break;
          }
          ++it;
        }

      }
    }
  }

  return results;
}

// ----------------------- DISK I/O  ----------------------

void vw::platefile::IndexPage::sync() {
  if (m_needs_saving) {
    this->serialize();
    m_needs_saving = false;
  }
}

void vw::platefile::IndexPage::serialize() {
  
  // Create the necessary directories if they do not yet exist.
  try {
    fs::path page_path(m_filename);
    fs::create_directories(page_path.parent_path());
  } catch ( fs::basic_filesystem_error<fs::path> &e ) { 
    vw_throw(IOErr() << "Could not create IndexPage.  " << e.what());
  }

  FILE *f = fopen(m_filename.c_str(), "w");
  if (!f)
    vw_throw(IOErr() << "IndexPage::serialize() failed.  Could not open " 
             << m_filename << " for writing.\n");

  // Part 1: Write out the page size
  fwrite(&m_page_width, sizeof(m_page_width), 1, f);
  fwrite(&m_page_height, sizeof(m_page_height), 1, f);

  // Part 2: Write the sparsetable metadata
  m_sparse_table.write_metadata(f);

  // Part 3: Write sparse entries
  for (google::sparsetable<multi_value_type>::nonempty_iterator it = m_sparse_table.nonempty_begin();
       it != m_sparse_table.nonempty_end(); ++it) {

    // Iterate over transaction_id list.
    int32 transaction_list_size = (*it).size();
    fwrite(&transaction_list_size, sizeof(transaction_list_size), 1, f);
    
    multi_value_type::iterator transaction_iter = (*it).begin();
    while (transaction_iter != (*it).end()) {
      
      // Save the transaction id
      int32 t_id = (*transaction_iter).first;
      fwrite(&t_id, sizeof(t_id), 1, f);
      
      // Save the size of each protobuf, and then serialize it to disk.
      uint16 protobuf_size = (*transaction_iter).second.ByteSize();
      boost::shared_array<uint8> protobuf_bytes( new uint8[protobuf_size] );
      (*transaction_iter).second.SerializeToArray(protobuf_bytes.get(), protobuf_size);
      fwrite(&protobuf_size, sizeof(protobuf_size), 1, f);
      fwrite(protobuf_bytes.get(), 1, protobuf_size, f);
    
      ++transaction_iter;
    }
  }

  fclose(f);
  m_needs_saving = false;
}

void vw::platefile::IndexPage::deserialize() {

  FILE *f = fopen(m_filename.c_str(), "r");
  if (!f)
    vw_throw(IOErr() << "IndexPage::deserialize() failed.  Could not open " 
             << m_filename << " for reading.\n");

  // Part 1: Read the page size
  fread(&m_page_width, sizeof(m_page_width), 1, f);
  fread(&m_page_height, sizeof(m_page_height), 1, f);

  // Part 2: Read the sparsetable metadata
  m_sparse_table.read_metadata(f);

  // Part 3: Read sparse entries
  for (google::sparsetable<multi_value_type>::nonempty_iterator it = m_sparse_table.nonempty_begin();
       it != m_sparse_table.nonempty_end(); ++it) {

    // Iterate over transaction_id list.
    int32 transaction_list_size;
    fread(&transaction_list_size, sizeof(transaction_list_size), 1, f);
    
    new (&(*it)) multi_value_type();
    for (int tid = 0; tid < transaction_list_size; ++tid) {

      // Read the transaction id
      int32 t_id;
      fread(&t_id, sizeof(t_id), 1, f);

      // Read the size (in bytes) of this protobuffer and then read
      // the protobuffer and deserialize it.
      uint16 protobuf_size;
      fread(&protobuf_size, sizeof(protobuf_size), 1, f);
      boost::shared_array<uint8> protobuf_bytes( new uint8[protobuf_size] );
      fread(protobuf_bytes.get(), 1, protobuf_size, f);
      IndexRecord rec;
      if (!rec.ParseFromArray(protobuf_bytes.get(), protobuf_size))
        vw_throw(IOErr() << "An error occurred while parsing an IndexEntry in " 
                 << m_filename << ".");
      
      (*it).push_back(std::pair<int32, IndexRecord>(t_id, rec));
    }
  }
  fclose(f);    
  m_needs_saving = false;
}
  

// ----------------------------------------------------------------------
//                         INDEX PAGE GENERATOR
// ----------------------------------------------------------------------

vw::platefile::IndexPageGenerator::IndexPageGenerator( std::string filename, 
                                                       int level, int base_col, int base_row, 
                                                       int page_width, int page_height) : 
  m_filename( filename ), m_level(level), m_base_col(base_col), m_base_row(base_row),
  m_page_width(page_width), m_page_height(page_height) {}  

size_t vw::platefile::IndexPageGenerator::size() const {
  return 1;
}

boost::shared_ptr<vw::platefile::IndexPage> vw::platefile::IndexPageGenerator::generate() const {
  //  vw_out(0) << "Generating cache line for " << m_filename << "\n";
  return boost::shared_ptr<IndexPage>(new IndexPage(m_filename, m_level, m_base_col, m_base_row,
                                                    m_page_width, m_page_height) );
}

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
  int32 page_col = col % m_page_width;
  int32 page_row = row % m_page_height;
  
  // Access the page.  This will load it into memory if necessary.
  boost::shared_ptr<IndexPage> page = m_cache_handles[level_row*m_horizontal_pages + level_col];
  return page->multi_get(page_col, page_row, start_transaction_id, end_transaction_id);
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
std::list<vw::platefile::TileHeader> 
vw::platefile::IndexLevel::valid_tiles(BBox2i const& region,
                                       int start_transaction_id, 
                                       int end_transaction_id) const {
  
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
                                                           end_transaction_id);
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

  //  std::cout << "\tRebuilding index: " << plate_filename <<"\n";

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

vw::platefile::PagedIndex::multi_value_type
vw::platefile::PagedIndex::multi_read_request(int col, int row, int level, 
                                              int start_transaction_id, 
                                              int end_transaction_id) {
  if (level < 0 || level >= m_levels.size())
    vw_throw(TileNotFoundErr() << "Requested tile at level " << level 
             << " was greater than the max level (" << m_levels.size() << ").");
  
  return m_levels[level]->multi_get(col, row,  start_transaction_id, end_transaction_id);
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

/// Returns a list of valid tiles that match this level, region, and
/// range of transaction_id's.  Returns a list of TileHeaders with
/// col/row/level and transaction_id of the most recent tile at each
/// valid location.  Note: there may be other tiles in the transaction
/// range at this col/row/level, but valid_tiles() only returns the
/// first one.
std::list<vw::platefile::TileHeader> 
vw::platefile::PagedIndex::valid_tiles(int level, BBox2i const& region,
                                       int start_transaction_id,
                                       int end_transaction_id) const {
  
  // If the level does not exist, we return an empty list.
  if (level < 0 && level >= m_levels.size())
    return std::list<TileHeader>();

  // Otherwise, we delegate to the valid_tiles() method for that level.
  return m_levels[level]->valid_tiles(region, start_transaction_id, end_transaction_id);
}

int32 vw::platefile::PagedIndex::num_levels() const {
  return m_levels.size();
}


