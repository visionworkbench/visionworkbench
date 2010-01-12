// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/Plate/IndexPage.h>

#include <boost/shared_array.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>
namespace fs = boost::filesystem;


// ----------------------------------------------------------------------
//                            INDEX PAGE
// ----------------------------------------------------------------------

vw::platefile::IndexPage::IndexPage(int level, int base_col, int base_row, 
                                    int page_width, int page_height) : 
  m_level(level), m_base_col(base_col), m_base_row(base_row),
  m_page_width(page_width), m_page_height(page_height) {

  vw_out(VerboseDebugMessage, "platefile::index") << "Opening index page [ " 
                                                  << m_base_col << " " << m_base_row 
                                                  << " @ " << m_level << " ]\n";
}

vw::platefile::IndexPage::~IndexPage() {
  vw_out(VerboseDebugMessage, "platefile::index") << "Closing index page [ " 
                                                  << m_base_col << " " << m_base_row 
                                                  << " @ " << m_level << " ]\n";
}

// ----------------------- ACCESSORS  ----------------------

void vw::platefile::IndexPage::set(IndexRecord const& record, 
                                   int col, int row, int transaction_id) {

  // Basic bounds checking
  VW_ASSERT(col >= 0 && col < m_page_width && row >= 0 && row < m_page_height, 
            ArgumentErr() << "IndexPage::set() failed.  Invalid index [" 
            << col << " " << row << "]");

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
                                                    multi_value_type const& candidates,
                                                    int col, int row, BBox2i const& region, 
                                                    int min_num_matches) const {

  // Check to see if the tile is in the specified region.
  Vector2i loc( m_base_col + col, m_base_row + row);
  if ( region.contains( loc ) && candidates.size() >= min_num_matches ) {
    TileHeader hdr;
    hdr.set_col( m_base_col + col );
    hdr.set_row( m_base_row + row );
    hdr.set_level(m_level);
    hdr.set_transaction_id(candidates.begin()->first);
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
                                      int end_transaction_id,
                                      int min_num_matches) const {
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
        multi_value_type candidates;
        while (it != entries.end() && it->first >= start_transaction_id) {
          if ( it->first >= start_transaction_id && it->first <= end_transaction_id ) {
            candidates.push_back(*it);
          }
          ++it;
        }
        
        // Do the region check.
        append_if_in_region( results, candidates, col, row, region, min_num_matches );
      }
    }
  }

  return results;
}

