// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#ifndef __VW_PLATEFILE_INDEX_PAGE_H__
#define __VW_PLATEFILE_INDEX_PAGE_H__

#include <vw/Core/FundamentalTypes.h>
#include <vw/Plate/Protobuffers.pb.h>
#include <google/sparsetable>
#include <list>

#include <boost/shared_ptr.hpp>

namespace vw {
namespace platefile {

  // ----------------------------------------------------------------------
  //                            INDEX PAGE
  // ----------------------------------------------------------------------

  class IndexPage {
    typedef std::list<std::pair<int32,IndexRecord> > element_type;

    std::string m_filename;
    int m_page_width, m_page_height;
    bool m_needs_saving;
    google::sparsetable<element_type> m_sparse_table;

    void serialize();
    void deserialize();

  public:
  
    /// Create or open a page file.
    IndexPage(std::string filename, int page_width, int page_height);

    ~IndexPage();

    // ----------------------------------------------------------------------
    //                          Accessors
    // ----------------------------------------------------------------------

    void set(IndexRecord const& record, int col, int row, int transaction_id);

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

    /// Return the number of valid entries in this page.  (Remember
    /// that this is a sparse store of IndexRecords.)
    int sparse_size() { return m_sparse_table.num_nonempty(); }

    // ----------------------------------------------------------------------
    //                            Disk I/O
    // ----------------------------------------------------------------------

    /// Save any unsaved changes to disk.
    void sync();

  };


  // ----------------------------------------------------------------------
  //                         INDEX PAGE GENERATOR
  // ----------------------------------------------------------------------

  // IndexPageGenerator loads a index page from disk.
  class IndexPageGenerator {
    std::string m_filename;
    int m_page_width;
    int m_page_height;

  public:
    typedef IndexPage value_type;
    IndexPageGenerator( std::string filename, int page_width, int page_height );
    size_t size() const;

    /// Generate an IndexPage.  If no file exists with the name
    /// m_filename, then an empty IndexPage is generated.
    boost::shared_ptr<IndexPage> generate() const;
  };



}} // namespace vw::platefile

#endif // __VW_PLATEFILE_INDEX_PAGE_H__
