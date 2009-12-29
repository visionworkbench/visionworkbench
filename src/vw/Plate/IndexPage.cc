// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/Plate/IndexPage.h>
#include <vw/Plate/Exception.h>
#include <vw/Core/FundamentalTypes.h>
#include <vw/Core/Log.h>

#include <boost/shared_array.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>
namespace fs = boost::filesystem;

// ----------------------------------------------------------------------
//                            INDEX PAGE
// ----------------------------------------------------------------------

vw::platefile::IndexPage::IndexPage(std::string filename, 
                                    int level, int base_col, int base_row, 
                                    int page_width, int page_height) : 
  m_filename(filename), m_level(level), m_base_col(base_col), m_base_row(base_row),
  m_page_width(page_width), m_page_height(page_height), m_needs_saving(false) {

  if (fs::exists(filename)) {
    this->deserialize();
  } else {
    m_sparse_table.resize(page_width*page_height);
  }
}

vw::platefile::IndexPage::~IndexPage() {
  this->sync();
}

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
    element_type *entries = m_sparse_table[row*m_page_width + col].operator&();

    // We need to keep this list sorted in decreasing order of
    // transaction ID, do a simple insertion sort here.
    element_type::iterator it = entries->begin();
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
    element_type l;
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
  element_type const& entries = m_sparse_table[row*m_page_width + col];
  element_type::const_iterator it = entries.begin();
  
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

void vw::platefile::IndexPage::append_record( std::list<vw::platefile::TileHeader> &results, 
                                              int transaction_id, IndexRecord const& rec, 
                                              int col, int row,
                                              BBox2i const& region) const {

  std::cout << "Appending record for level " << m_level << "  at " << col << " " << row << "\n";
        

  // Check to see if the region contains a col/row.
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

/// Returns a list of valid tiles in this IndexPage. 
std::list<vw::platefile::TileHeader> vw::platefile::IndexPage::valid_tiles(int transaction_id, 
                                                                           vw::BBox2i const& region,
                                                                           bool exact_match) const {

  std::list<TileHeader> results;

  for (int row = 0; row < m_page_height; ++row) {
    for (int col = 0; col < m_page_width; ++col) {
      if (m_sparse_table.test(row*m_page_width + col)) {

        // Interate over entries.
        element_type const& entries = m_sparse_table[row*m_page_width + col];
        element_type::const_iterator it = entries.begin();

        
        // A transaction ID of -1 indicates that we should return the most
        // recent tile (which is the first entry in the list, since it is
        // sorted from most recent to least recent), regardless of its
        // transaction id.
        if (transaction_id == -1 && it != entries.end())
          append_record( results, (*it).first, (*it).second, col, row, region );

        // Otherwise, we search through the entries in the list, looking for
        // the requested t_id.  Note: this search is O(n), so it can be slow
        // if there are a lot of entries and the entry you are looking for
        // is near the end.  However, most pages will contain very few
        // entries, and for those with many entries (i.e. tiles near the
        // root of the mosaic), you will rarely search for old tiles.
        while (it != entries.end()) {
          if (exact_match) {
            if ((*it).first == transaction_id)
              append_record( results, (*it).first, (*it).second, col, row, region );
          } else {
            if ((*it).first <= transaction_id)
              append_record( results, (*it).first, (*it).second, col, row, region );
          }
          ++it;
        }

      }
    }
  }

  return results;
}

// ----------------------------------------------------------------------
//                            Disk I/O
// ----------------------------------------------------------------------

void vw::platefile::IndexPage::sync() {
  if (m_needs_saving)
    this->serialize();
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
  for (google::sparsetable<element_type>::nonempty_iterator it = m_sparse_table.nonempty_begin();
       it != m_sparse_table.nonempty_end(); ++it) {

    // Iterate over transaction_id list.
    int32 transaction_list_size = (*it).size();
    fwrite(&transaction_list_size, sizeof(transaction_list_size), 1, f);
    
    element_type::iterator transaction_iter = (*it).begin();
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
  for (google::sparsetable<element_type>::nonempty_iterator it = m_sparse_table.nonempty_begin();
       it != m_sparse_table.nonempty_end(); ++it) {

    // Iterate over transaction_id list.
    int32 transaction_list_size;
    fread(&transaction_list_size, sizeof(transaction_list_size), 1, f);
    
    new (&(*it)) element_type();
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
