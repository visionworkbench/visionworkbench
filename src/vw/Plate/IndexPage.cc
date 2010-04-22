// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Plate/IndexPage.h>
#include <vw/Core/Debugging.h>

#include <boost/shared_array.hpp>

#define WHEREAMI (vw::vw_out(VerboseDebugMessage, "platefile.index") << VW_CURRENT_FUNCTION << ": ")

// ----------------------------------------------------------------------
//                            INDEX PAGE
// ----------------------------------------------------------------------

vw::platefile::IndexPage::IndexPage(int level, int base_col, int base_row, 
                                    int page_width, int page_height) : 
  m_level(level), m_base_col(base_col), m_base_row(base_row),
  m_page_width(page_width), m_page_height(page_height) {

  // Set the size of the sparse table.
  m_sparse_table.resize(page_width*page_height);

  WHEREAMI << "[" << m_base_col << " " << m_base_row << " @ " << m_level << "]\n";
}

vw::platefile::IndexPage::~IndexPage() {
  WHEREAMI << "[" << m_base_col << " " << m_base_row << " @ " << m_level << "]\n";
}

void vw::platefile::IndexPage::serialize(std::ostream& ostr) {
  WHEREAMI << "[" << m_base_col << " " << m_base_row << " @ " << m_level << "]\n";

  // Part 1: Write out the page size
  ostr.write((char*)(&m_page_width), sizeof(m_page_width));
  ostr.write((char*)(&m_page_height), sizeof(m_page_height));

  // Part 2: Write the sparsetable metadata
  m_sparse_table.write_metadata(&ostr);

  // Part 3: Write sparse entries
  for (google::sparsetable<multi_value_type>::nonempty_iterator it = m_sparse_table.nonempty_begin();
       it != m_sparse_table.nonempty_end(); ++it) {

    // Iterate over transaction_id list.
    int32 transaction_list_size = (*it).size();
    ostr.write((char*)(&transaction_list_size), sizeof(transaction_list_size));
    
    multi_value_type::iterator transaction_iter = (*it).begin();
    while (transaction_iter != (*it).end()) {
      
      // Save the transaction id
      int32 t_id = (*transaction_iter).first;
      ostr.write((char*)(&t_id), sizeof(t_id));
      
      // Save the size of each protobuf, and then serialize it to disk.
      uint16 protobuf_size = (*transaction_iter).second.ByteSize();
      boost::shared_array<uint8> protobuf_bytes( new uint8[protobuf_size] );
      (*transaction_iter).second.SerializeToArray(protobuf_bytes.get(), protobuf_size);
      ostr.write((char*)(&protobuf_size), sizeof(protobuf_size));
      ostr.write((char*)(protobuf_bytes.get()), protobuf_size);
    
      ++transaction_iter;
    }
  }
}

void vw::platefile::IndexPage::deserialize(std::istream& istr) {

  WHEREAMI << "[" << m_base_col << " " << m_base_row << " @ " << m_level << "]\n";

  // Part 1: Read the page size
  istr.read((char*)(&m_page_width), sizeof(m_page_width));
  istr.read((char*)(&m_page_height), sizeof(m_page_height));

  // Part 2: Read the sparsetable metadata
  m_sparse_table.read_metadata(&istr);

  // Part 3: Read sparse entries
  for (google::sparsetable<multi_value_type>::nonempty_iterator it = m_sparse_table.nonempty_begin();
       it != m_sparse_table.nonempty_end(); ++it) {

    // Iterate over transaction_id list.
    int32 transaction_list_size;
    istr.read((char*)(&transaction_list_size), sizeof(transaction_list_size));
    
    new (&(*it)) multi_value_type();
    for (int tid = 0; tid < transaction_list_size; ++tid) {

      // Read the transaction id
      int32 t_id;
      istr.read((char*)(&t_id), sizeof(t_id));

      // Read the size (in bytes) of this protobuffer and then read
      // the protobuffer and deserialize it.
      uint16 protobuf_size;
      istr.read((char*)(&protobuf_size), sizeof(protobuf_size));
      boost::shared_array<uint8> protobuf_bytes( new uint8[protobuf_size] );
      istr.read((char*)(protobuf_bytes.get()), protobuf_size);
      IndexRecord rec;
      if (!rec.ParseFromArray(protobuf_bytes.get(), protobuf_size))
        vw_throw(IOErr() << "An error occurred while parsing an IndexEntry."); 
      
      (*it).push_back(std::pair<int32, IndexRecord>(t_id, rec));
    }
  }
}

// ----------------------- ACCESSORS  ----------------------

void vw::platefile::IndexPage::set(TileHeader const& header, IndexRecord const& record) {

  VW_ASSERT( header.col() >= 0 && header.row() >= 0,
             ArgumentErr() << "IndexPage::set() failed.  Column and row indices must be positive.");

  int32 page_col = header.col() % m_page_width;
  int32 page_row = header.row() % m_page_height;

  std::pair<int32, IndexRecord> p(header.transaction_id(), record);
  int elmnt = page_row*m_page_width + page_col;
  if (m_sparse_table.test(elmnt)) {

    // Add to existing entry.
    multi_value_type *entries = m_sparse_table[page_row*m_page_width + page_col].operator&();

    // We need to keep this list sorted in decreasing order of
    // transaction ID, do a simple insertion sort here.
    multi_value_type::iterator it = entries->begin();
    while (it != entries->end() && (*it).first >= int64(header.transaction_id()) ) {

      // Handle the case where we replace an entry
      if ( (*it).first == int64(header.transaction_id()) ) {
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

  VW_ASSERT( col >= 0 && row >= 0,
             ArgumentErr() << "IndexPage::get() failed.  Column and row indices must be positive.");

  // Compute page_col and row
  int32 page_col = col % m_page_width;
  int32 page_row = row % m_page_height;

  // Interate over entries.
  multi_value_type const& entries = m_sparse_table[page_row*m_page_width + page_col];
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


void vw::platefile::IndexPage::append_if_in_region( std::list<vw::platefile::TileHeader> &results, 
                                                    multi_value_type const& candidates,
                                                    int col, int row, BBox2i const& region, 
                                                    int min_num_matches) const {

  // Check to see if the tile is in the specified region.
  Vector2i loc( m_base_col + col, m_base_row + row);
  if ( region.contains( loc ) && int(candidates.size()) >= min_num_matches ) {
    TileHeader hdr;
    hdr.set_col( m_base_col + col );
    hdr.set_row( m_base_row + row );
    hdr.set_level(m_level);

    // Return the transaction ID of the first result in the list.
    hdr.set_transaction_id(candidates.begin()->first);
    results.push_back(hdr);
  }
}

/// Returns a list of valid tiles in this IndexPage.  Returns a list
/// of TileHeaders with col/row/level and transaction_id of the most
/// recent tile at each valid location.  Note: there may be other
/// tiles in the transaction range at this col/row/level, but
/// search_by_region() only returns the first one.
std::list<vw::platefile::TileHeader> 
vw::platefile::IndexPage::search_by_region(vw::BBox2i const& region,
                                           int start_transaction_id,
                                           int end_transaction_id,
                                           int min_num_matches,
                                           bool fetch_one_additional_entry) const {
  std::list<TileHeader> results;

  for (int row = 0; row < m_page_height; ++row) {
    for (int col = 0; col < m_page_width; ++col) {
      if (m_sparse_table.test(row*m_page_width + col)) {

        // Iterate over entries.
        multi_value_type const& entries = m_sparse_table[row*m_page_width + col];
        multi_value_type::const_iterator it = entries.begin();
        multi_value_type candidates;

        if (start_transaction_id == -1 && end_transaction_id == -1) {
          
          // If the user has specified a transaction range of [-1, -1],
          // then we only return the last valid tile.
          if (entries.size() > 0)
            candidates.push_back( *(entries.begin()) );

        } else {

          // std::ostringstream ostr;
          // ostr << "Accumulating entries for tile ["
          //      << (m_base_col + col) << " " << (m_base_row + row) <<  " @ " 
          //      << m_level << "] with " << entries.size() << " entries -- ";

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
              //              ostr << it->first << " ";
              candidates.push_back(*it);
            }
            ++it;
          }

          // For snapshotting, we need to fetch one additional entry
          // outside of the specified range.  This next tile
          // represents the "top" tile in the mosaic for entries that
          // may not have been part of the last snapshot.  
          if (fetch_one_additional_entry && it != entries.end()) {
            //ostr << it->first << " ";
            candidates.push_back(*it);
          }
          
          // For debugging
          // if (candidates.size() != 0)// && m_level == 10)
          //   std::cout << ostr.str() << "   ( " << candidates.size() << " )\n";
        }
        
        // Do the region check.
        append_if_in_region( results, candidates, col, row, region, min_num_matches );
      }
    }
  }

  return results;
}

std::list<vw::platefile::TileHeader> 
vw::platefile::IndexPage::search_by_location(int col, int row, 
                                             int start_transaction_id, 
                                             int end_transaction_id, 
                                             bool fetch_one_additional_entry) const {

  int32 page_col = col % m_page_width;
  int32 page_row = row % m_page_height;

  // Basic bounds checking
  VW_ASSERT(page_col >= 0 && page_col < m_page_width && page_row >= 0 && page_row < m_page_height, 
            ArgumentErr() << "IndexPage::read_headers() failed.  Invalid index [" 
            << page_col << " " << page_row << "]");

  // Check first to make sure that there are actually tiles at this location.
  if (!m_sparse_table.test(page_row*m_page_width + page_col))
    return std::list<TileHeader>();

  // If there are, then we apply the transaction_id filters to select the requested ones.
  std::list<vw::platefile::TileHeader> results;
  multi_value_type const& entries = m_sparse_table[page_row*m_page_width + page_col];
  multi_value_type::const_iterator it = entries.begin();
  while (it != entries.end() && it->first >= start_transaction_id) {
    if (it->first >= start_transaction_id && it->first <= end_transaction_id) {
      TileHeader hdr;
      hdr.set_col( m_base_col + page_col );
      hdr.set_row( m_base_row + page_row );
      hdr.set_level(m_level);
      hdr.set_transaction_id(it->first);
      results.push_back(hdr);
    }
    ++it;
  }
  
  // For snapshotting, we need to fetch one additional entry
  // outside of the specified range.  This next tile
  // represents the "top" tile in the mosaic for entries that
  // may not have been part of the last snapshot.  
  if (fetch_one_additional_entry && it != entries.end()) {
    TileHeader hdr;
    hdr.set_col( m_base_col + page_col );
    hdr.set_row( m_base_row + page_row );
    hdr.set_level(m_level);
    hdr.set_transaction_id(it->first);
    results.push_back(hdr);
  }
  
  return results;
}
