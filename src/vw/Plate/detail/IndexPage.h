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


#ifndef __VW_PLATEFILE_INDEX_PAGE_H__
#define __VW_PLATEFILE_INDEX_PAGE_H__

#include <vw/Plate/FundamentalTypes.h>
#include <vw/Plate/IndexData.pb.h>
#include <vw/Plate/IndexDataPrivate.pb.h>
#include <vw/Plate/google/sparsetable>
#include <vw/Math/BBox.h>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/shared_ptr.hpp>
#include <string>
#include <list>

namespace vw {
namespace platefile {
namespace detail {

  // ----------------------------------------------------------------------
  //                            INDEX PAGE
  // ----------------------------------------------------------------------

  class IndexPage {

  public:
    typedef std::pair<uint32, IndexRecord> value_type;
    typedef std::list<value_type> multi_value_type;
    typedef google::sparsetable<multi_value_type>::nonempty_iterator nonempty_iterator;

  protected:
    uint32 m_level, m_base_col, m_base_row;
    uint32 m_page_width, m_page_height;
    google::sparsetable<multi_value_type> m_sparse_table;

    void append_if_in_region( std::list<TileHeader> &results,
                              multi_value_type const& candidates,
                              uint32 col, uint32 row, BBox2i const& region) const;

    TileHeader hdr_from_index(uint32 rel_col, uint32 rel_row, const value_type& elt) const {
      TileHeader hdr;
      hdr.set_col( m_base_col + rel_col );
      hdr.set_row( m_base_row + rel_row );
      hdr.set_level(m_level);
      hdr.set_transaction_id(elt.first);
      if (elt.second.has_filetype())
        hdr.set_filetype(elt.second.filetype());
      return hdr;
    }

  public:

    /// Create or open a page file.
    IndexPage(uint32 level, uint32 base_col, uint32 base_row,
              uint32 page_width, uint32 page_height);

    virtual ~IndexPage();

    /// Save any unsaved changes to disk.
    virtual void sync() = 0;

    // For reading/writing to/from disk or a network byte stream.
    void serialize(std::ostream& ostr);
    void deserialize(std::istream& istr);

    // ----------------------- ITERATORS  ----------------------

    nonempty_iterator begin() { return m_sparse_table.nonempty_begin(); }
    nonempty_iterator end() { return m_sparse_table.nonempty_end(); }

    // ----------------------- ACCESSORS  ----------------------

    /// Set the value of an entry in the IndexPage.
    virtual void set(TileHeader const& header, IndexRecord const& record);

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
    IndexRecord get(uint32 col, uint32 row, TransactionOrNeg transaction_id, bool exact_match = false) const;

    /// Return multiple index entries that match the specified
    /// transaction id range.  This range is inclusive of the first
    /// entry AND the last entry: [ begin_transaction_id, end_transaction_id ]
    ///
    /// Results are return as a std::pair<int32, IndexRecord>.  The
    /// first value in the pair is the transaction id for that
    /// IndexRecord.
    ///
    /// Note: this function is mostly used when creating snapshots.
    multi_value_type multi_get(uint32 col, uint32 row,
                               TransactionOrNeg begin_transaction_id, TransactionOrNeg end_transaction_id) const;

    /// Return the number of valid entries in this page.  (Remember
    /// that this is a sparse store of IndexRecords.)
    uint32 sparse_size() { return boost::numeric_cast<uint32>(m_sparse_table.num_nonempty()); }

    /// Returns a list of valid tiles in this IndexPage.
    ///
    /// Note: this function is mostly used when creating snapshots.
    std::list<TileHeader> search_by_region(vw::BBox2i const& region,
                                           TransactionOrNeg start_transaction_id,
                                           TransactionOrNeg end_transaction_id) const;

    /// Return multiple tile headers that match the specified
    /// transaction id range.  This range is inclusive of the first
    /// entry AND the last entry: [ begin_transaction_id, end_transaction_id ]
    ///
    /// Results are return as a std::pair<int32, IndexRecord>.  The
    /// first value in the pair is the transaction id for that
    /// IndexRecord.
    ///
    /// Note: this function is mostly used when creating snapshots.
    std::list<TileHeader> search_by_location(uint32 col, uint32 row,
                                             TransactionOrNeg start_transaction_id,
                                             TransactionOrNeg end_transaction_id) const;

  };

  // ----------------------------------------------------------------------
  //                       INDEX PAGE GENERATOR
  // ----------------------------------------------------------------------

  /// Subclasess of PageGeneratorBase load and unload index pages.
  class PageGeneratorBase {
  public:
    typedef IndexPage value_type;
    virtual size_t size() const { return 1; }
    virtual ~PageGeneratorBase() {};
    virtual boost::shared_ptr<IndexPage> generate() const = 0;
  };

  // PageGeneratorFactory is the base class for entities that can
  // generate PageGenerators.
  class PageGeneratorFactory {

  public:
    virtual ~PageGeneratorFactory() {}
    virtual boost::shared_ptr<PageGeneratorBase>
      create(uint32 level, uint32 base_col, uint32 base_row,
             uint32 page_width, uint32 page_height) = 0;
    // Who is this factory manufacturing pages for? (human-readable)
    virtual std::string who() const = 0;
  };

}}}

#endif

