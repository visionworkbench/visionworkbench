// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NASA Vision Workbench is licensed under the Apache License,
//  Version 2.0 (the "License"); you may not use this file except in
//  compliance with the License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__


#include <vw/Plate/detail/IndexPage.h>
#include <vw/Plate/Exception.h>
#include <vw/Math/Vector.h>
#include <vw/Math/BBox.h>
#include <vw/Core/Debugging.h>

#include <boost/shared_array.hpp>
#include <boost/foreach.hpp>

#define WHEREAMI (vw::vw_out(VerboseDebugMessage, "platefile.index") << VW_CURRENT_FUNCTION << ": ")
using namespace vw;
using namespace vw::platefile;
using namespace vw::platefile::detail;

// ----------------------------------------------------------------------
//                            INDEX PAGE
// ----------------------------------------------------------------------

IndexPage::IndexPage(uint32 level, uint32 base_col, uint32 base_row,
                     uint32 page_width, uint32 page_height) :
  m_level(level), m_base_col(base_col), m_base_row(base_row),
  m_page_width(page_width), m_page_height(page_height) {

  // Set the size of the sparse table.
  m_sparse_table.resize(page_width*page_height);

  WHEREAMI << "[" << m_base_col << " " << m_base_row << " @ " << m_level << "]\n";
}

IndexPage::~IndexPage() {
  WHEREAMI << "[" << m_base_col << " " << m_base_row << " @ " << m_level << "]\n";
}

void IndexPage::serialize(std::ostream& ostr) {
  WHEREAMI << "[" << m_base_col << " " << m_base_row << " @ " << m_level << "]\n";

  // Part 1: Write out the page size
  ostr.write(reinterpret_cast<char*>(&m_page_width), sizeof(m_page_width));
  ostr.write(reinterpret_cast<char*>(&m_page_height), sizeof(m_page_height));

  // Part 2: Write the sparsetable metadata
  m_sparse_table.write_metadata(&ostr);

  // Part 3: Write sparse entries
  for (nonempty_iterator it = m_sparse_table.nonempty_begin(); it != m_sparse_table.nonempty_end(); ++it) {

    // Iterate over transaction_id list.
    uint32 transaction_list_size = boost::numeric_cast<uint32>(it->size());
    ostr.write(reinterpret_cast<char*>(&transaction_list_size), sizeof(transaction_list_size));

    multi_value_type::iterator transaction_iter = (*it).begin();
    while (transaction_iter != (*it).end()) {

      // Save the transaction id
      uint32 t_id = (*transaction_iter).first;
      ostr.write(reinterpret_cast<char*>(&t_id), sizeof(t_id));

      // Save the size of each protobuf, and then serialize it to disk.
      uint16 protobuf_size = boost::numeric_cast<uint16>(transaction_iter->second.ByteSize());
      boost::shared_array<uint8> protobuf_bytes( new uint8[protobuf_size] );
      (*transaction_iter).second.SerializeToArray(protobuf_bytes.get(), protobuf_size);
      ostr.write(reinterpret_cast<char*>(&protobuf_size), sizeof(protobuf_size));
      ostr.write(reinterpret_cast<char*>(protobuf_bytes.get()), protobuf_size);

      ++transaction_iter;
    }
  }
}

void IndexPage::deserialize(std::istream& istr) {

  WHEREAMI << "[" << m_base_col << " " << m_base_row << " @ " << m_level << "]\n";

  VW_ASSERT(istr.good(), IOErr() << "while beginning to deserialize.");

  // Part 1: Read the page size
  istr.read(reinterpret_cast<char*>(&m_page_width), sizeof(m_page_width));
  istr.read(reinterpret_cast<char*>(&m_page_height), sizeof(m_page_height));

  VW_ASSERT(istr.good(), IOErr() << "while reading page size.");

  // Part 2: Read the sparsetable metadata
  if (!m_sparse_table.read_metadata(&istr))
    vw_throw(IOErr() << "while reading sparse table metadata.");

  VW_ASSERT(istr.good(), IOErr() << "after reading sparse table metadata.");

  // make sure we initialize everything before we continue
  BOOST_FOREACH(multi_value_type& x, std::make_pair(m_sparse_table.nonempty_begin(), m_sparse_table.nonempty_end()))
    new (&x) multi_value_type();

  size_t count = 0, total = m_sparse_table.num_nonempty();

  // Part 3: Read sparse entries
  BOOST_FOREACH(multi_value_type& x, std::make_pair(m_sparse_table.nonempty_begin(), m_sparse_table.nonempty_end())) {
    if (count++ % 4000 == 0)
      WHEREAMI << "reading tile slot " << count << " of " << total << std::endl;
    // Iterate over transaction_id list.
    uint32 transaction_list_size;
    istr.read(reinterpret_cast<char*>(&transaction_list_size), sizeof(transaction_list_size));

    VW_ASSERT(istr.good(), IOErr() << "while reading transaction list size.");

    for (uint32 tid = 0; tid < transaction_list_size; ++tid) {

      // Read the transaction id
      uint32 t_id;
      istr.read(reinterpret_cast<char*>(&t_id), sizeof(t_id));
      VW_ASSERT(istr.good(), IOErr() << "while reading a transaction id.");

      // Read the size (in bytes) of this protobuffer and then read
      // the protobuffer and deserialize it.
      uint16 protobuf_size;
      istr.read(reinterpret_cast<char*>(&protobuf_size), sizeof(protobuf_size));
      VW_ASSERT(istr.good(), IOErr() << "while reading a message size.");

      boost::shared_array<uint8> protobuf_bytes( new uint8[protobuf_size] );
      istr.read(reinterpret_cast<char*>(protobuf_bytes.get()), protobuf_size);
      VW_ASSERT(istr.good(), IOErr() << "while reading a message.");

      IndexRecord rec;
      if (!rec.ParseFromArray(protobuf_bytes.get(), protobuf_size))
        vw_throw(IOErr() << "while parsing a message.");

      x.push_back(value_type(t_id, rec));
    }
  }

  if (istr.peek() != EOF)
    vw_out(WarningMessage, "platefile.index") << "Unparsed data remaining in index page.\n";
}

// ----------------------- ACCESSORS  ----------------------

void IndexPage::set(TileHeader const& header, IndexRecord const& record) {

  uint32 page_col = header.col() % m_page_width;
  uint32 page_row = header.row() % m_page_height;

  value_type p(header.transaction_id(), record);
  uint32 elmnt = page_row*m_page_width + page_col;
  if (m_sparse_table.test(elmnt)) {

    // Add to existing entry.
    multi_value_type *entries = m_sparse_table[elmnt].operator&();

    // We need to keep this list sorted in decreasing order of
    // transaction ID, do a simple insertion sort here.
    multi_value_type::iterator it = entries->begin();
    while (it != entries->end() && (*it).first >= header.transaction_id() ) {

      // Handle the case where we replace an entry
      if ( (*it).first == header.transaction_id() ) {
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
IndexRecord IndexPage::get(uint32 col, uint32 row, TransactionOrNeg transaction_id_neg, bool exact_match) const {

  // Compute page_col and row
  uint32 page_col = col % m_page_width;
  uint32 page_row = row % m_page_height;

  // Interate over entries.
  multi_value_type const& entries = m_sparse_table[page_row*m_page_width + page_col];
  multi_value_type::const_iterator it = entries.begin();

  if ( it == entries.end() )
    vw_throw(TileNotFoundErr() << "No Tiles exist at this location.");

  // A transaction ID of -1 indicates that we should return the most
  // recent tile (which is the first entry in the list, since it is
  // sorted from most recent to least recent), regardless of its
  // transaction id.
  if (transaction_id_neg.newest())
    return (*it).second;

  Transaction transaction_id = transaction_id_neg.promote();

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
           << "but none before transaction_id = "  << transaction_id);
  return IndexRecord(); // never reached
}

IndexPage::multi_value_type
IndexPage::multi_get(uint32 col, uint32 row,
                     TransactionOrNeg begin_transaction_id, TransactionOrNeg end_transaction_id) const {
  // Compute page_col and row
  uint32 page_col = col % m_page_width;
  uint32 page_row = row % m_page_height;

  // Interate over entries.
  multi_value_type const& entries = m_sparse_table[page_row*m_page_width + page_col];
  multi_value_type::const_iterator it = entries.begin();

  if ( it == entries.end() )
    vw_throw(TileNotFoundErr() << "No Tiles exist at this location.");
  if ( begin_transaction_id <= 1 && end_transaction_id == detail::MAX_TRANSACTION )
    return entries;

  multi_value_type result;
  if ( begin_transaction_id.newest() && end_transaction_id.newest() ) {
    // Pull only the top most transaction
    result.push_back( *it );
    return result;
  } else {
    // Pull only items that our in our transaction range
    BOOST_FOREACH(const value_type& elt, entries ) {
      if ( elt.first < begin_transaction_id )
        break;
      if ( elt.first <= end_transaction_id)
        result.push_back(elt);
    }
  }

  return result;
}

void IndexPage::append_if_in_region( std::list<TileHeader> &results,
                                     multi_value_type const& candidates,
                                     uint32 col, uint32 row, BBox2i const& region) const {

  BOOST_FOREACH(const value_type& t, candidates) {
    Vector2i loc( m_base_col + col, m_base_row + row);
    if ( region.contains( loc ) )
      results.push_back(hdr_from_index(col, row, t));
  }
  // Check to see if the tile is in the specified region.
}

/// Returns a list of valid tiles in this IndexPage.  Returns a list
/// of TileHeaders with col/row/level and transaction_id of the most
/// recent tile at each valid location.  Note: there may be other
/// tiles in the transaction range at this col/row/level, but
/// search_by_region() only returns the first one.
std::list<TileHeader>
IndexPage::search_by_region(BBox2i const& region,
                            TransactionOrNeg start_transaction_id,
                            TransactionOrNeg end_transaction_id) const {

  // empty range means something broke upstream
  VW_ASSERT(start_transaction_id <= end_transaction_id,
            ArgumentErr() << VW_CURRENT_FUNCTION << ": received a null set range ["
                          << start_transaction_id << "," << end_transaction_id << "]");

  if (region.min().x() < 0 || region.min().y() < 0)
    vw_out(WarningMessage) << VW_CURRENT_FUNCTION << ": " << "asked for a region < 0: " << region << std::endl;
  const int32 MAX_IDX = 1 << m_level;
  if (region.max().x() > MAX_IDX || region.max().y() > MAX_IDX)
    vw_out(WarningMessage) << VW_CURRENT_FUNCTION << ": " << "asked for a region outside valid area for level " << m_level << ": " << region << std::endl;

  // Search through the entries in the index page which are in the current region.
  uint32 beg_row = (uint32)std::max((int32)0,             (int32)region.min().y() - (int32)m_base_row);
  uint32 end_row = (uint32)std::min((int32)m_page_height, (int32)region.max().y() - (int32)m_base_row);
  uint32 beg_col = (uint32)std::max((int32)0,             (int32)region.min().x() - (int32)m_base_col);
  uint32 end_col = (uint32)std::min((int32)m_page_width,  (int32)region.max().x() - (int32)m_base_col);

  std::list<TileHeader> results;
  for (uint32 row = beg_row; row < end_row; ++row) {
    for (uint32 col = beg_col; col < end_col; ++col) {
      if (m_sparse_table.test(row*m_page_width + col)) {

        // Iterate over entries.
        multi_value_type const& entries = m_sparse_table[row*m_page_width + col];
        multi_value_type candidates;

        if (start_transaction_id.newest() && end_transaction_id.newest()) {
          // If the user has specified a transaction range of [-1, -1],
          // then we only return the last valid tile.
          if (entries.size() > 0)
            candidates.push_back( *(entries.begin()) );
        } else {
          // Search through the entries in the list, looking entries
          // that match the requested transaction_id range.  Note: this
          // search is O(n), so it can be slow if there are a lot of
          // entries and the entry you are looking for is near the end.
          // However, most pages will contain very few entries, and for
          // those with many entries (i.e. tiles near the root of the
          // mosaic), you will most likely be searching for recently
          // added tiles, which are sorted to the beginning.
          BOOST_FOREACH(const value_type& elt, entries) {
            if (elt.first < start_transaction_id)
              break;
            if (elt.first <= end_transaction_id)
              candidates.push_back(elt);
          }
        }

        append_if_in_region( results, candidates, col, row, region );
      }
    }
  }

  return results;
}

std::list<TileHeader>
IndexPage::search_by_location(uint32 col, uint32 row,
                              TransactionOrNeg start_transaction_id,
                              TransactionOrNeg end_transaction_id) const {

  // empty range means something broke upstream
  VW_ASSERT(start_transaction_id <= end_transaction_id,
            ArgumentErr() << VW_CURRENT_FUNCTION << ": received a null set range ["
                          << start_transaction_id << "," << end_transaction_id << "]");

  uint32 page_col = col % m_page_width;
  uint32 page_row = row % m_page_height;

  std::list<TileHeader> results;

  // Check first to make sure that there are actually tiles at this location.
  if (!m_sparse_table.test(page_row*m_page_width + page_col))
    return results;

  multi_value_type const& entries = m_sparse_table[page_row*m_page_width + page_col];
  if (entries.size() == 0)
    return results;

  if (start_transaction_id.newest() && end_transaction_id.newest())
    results.push_back(hdr_from_index(page_col, page_row, *entries.begin()));
  else {
    BOOST_FOREACH(const value_type& elt, entries) {
      if (elt.first < start_transaction_id)
        break;
      if (elt.first <= end_transaction_id)
        results.push_back(hdr_from_index(page_col, page_row, elt));
    }
  }

  return results;
}
