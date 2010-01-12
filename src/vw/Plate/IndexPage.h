// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__
#ifndef __VW_PLATEFILE_INDEX_PAGE_H__
#define __VW_PLATEFILE_INDEX_PAGE_H__

#include <vw/Core/FundamentalTypes.h>
#include <vw/Core/Log.h>
#include <vw/Math/BBox.h>
#include <vw/Plate/ProtoBuffers.pb.h>
#include <vw/Plate/Exception.h>

#include <google/sparsetable>
#include <string>
#include <list>

namespace vw {
namespace platefile {

  // ----------------------------------------------------------------------
  //                            INDEX PAGE
  // ----------------------------------------------------------------------

  class IndexPage {

  public:
    typedef std::list<std::pair<int32,IndexRecord> > multi_value_type;
    typedef google::sparsetable<multi_value_type>::nonempty_iterator nonempty_iterator;

  protected:
    int m_level, m_base_col, m_base_row;
    int m_page_width, m_page_height;
    google::sparsetable<multi_value_type> m_sparse_table;

    void append_if_in_region( std::list<vw::platefile::TileHeader> &results, 
                              multi_value_type const& candidates,
                              int col, int row, BBox2i const& region, int min_num_matches) const;

  public:
  
    /// Create or open a page file.
    IndexPage(int level, int base_col, int base_row, 
              int page_width, int page_height);

    virtual ~IndexPage();

    /// Save any unsaved changes to disk.
    virtual void sync() = 0;

    // ----------------------- ITERATORS  ----------------------

    nonempty_iterator begin() { return m_sparse_table.nonempty_begin(); }
    nonempty_iterator end() { return m_sparse_table.nonempty_end(); }

    // ----------------------- ACCESSORS  ----------------------

    /// Set the value of an entry in the IndexPage.
    virtual void set(IndexRecord const& record, int col, int row, int transaction_id);

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
    IndexRecord get(int col, int row, int transaction_id, bool exact_match = false) const;

    /// Return multiple index entries that match the specified
    /// transaction id range.  This range is inclusive of the first
    /// entry AND the last entry: [ begin_transaction_id, end_transaction_id ]
    ///
    /// Results are return as a std::pair<int32, IndexRecord>.  The
    /// first value in the pair is the transaction id for that
    /// IndexRecord.
    ///
    /// Note: this function is mostly used when creating snapshots.
    multi_value_type multi_get(int col, int row, 
                               int begin_transaction_id, int end_transaction_id) const;

    /// Return the number of valid entries in this page.  (Remember
    /// that this is a sparse store of IndexRecords.)
    int sparse_size() { return m_sparse_table.num_nonempty(); }

    /// Returns a list of valid tiles in this IndexPage.  
    ///
    /// Note: this function is mostly used when creating snapshots.
    std::list<TileHeader> valid_tiles(vw::BBox2i const& region, 
                                      int start_transaction_id, 
                                      int end_transaction_id, 
                                      int min_num_matches) const;
  };

  // ----------------------------------------------------------------------
  //                       INDEX PAGE GENERATOR
  // ----------------------------------------------------------------------

  /// Subclasess of PageGeneratorBase load and unload index pages.
  class PageGeneratorBase {
  public:
    virtual ~PageGeneratorBase() {};
    virtual boost::shared_ptr<IndexPage> generate() const = 0;
  };

  /// The sole purpose of the IndexPageGenerator class is to hold a
  /// pointer to an instance of the PageGeneratorBase class.  This is
  /// necessary because the VW caching system doesn't store the
  /// generator as a pointer, which prevents us from using the caching
  /// system to store polymorphic index generator types.  So, although
  /// it ain't pretty, this is really necessary for now.
  class IndexPageGenerator {
    boost::shared_ptr<PageGeneratorBase> m_page_gen;

  public:
    typedef IndexPage value_type;

    IndexPageGenerator(boost::shared_ptr<PageGeneratorBase> page_gen) : 
      m_page_gen(page_gen) {}

    size_t size() const { return 1; }

    /// Generate an IndexPage.  If no file exists with the name
    /// m_filename, then an empty IndexPage is generated.
    boost::shared_ptr<IndexPage> generate() const {
      return m_page_gen->generate();
    }
  };

  // PageGeneratorFactory is the base class for entities that can
  // generate PageGenerators.  
  class PageGeneratorFactory {
    
  public:
    virtual ~PageGeneratorFactory() {}
    virtual boost::shared_ptr<IndexPageGenerator> create(int level, int base_col, int base_row, 
                                                         int page_width, int page_height) = 0;
  };

}} // namespace vw::platefile

#endif // __VW_PLATEFILE_INDEX_PAGE_H__

