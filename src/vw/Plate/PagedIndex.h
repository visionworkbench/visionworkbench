// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__
#ifndef __VW_PLATEFILE_PAGED_INDEX_H__
#define __VW_PLATEFILE_PAGED_INDEX_H__

#include <vw/Core/FundamentalTypes.h>
#include <vw/Core/Cache.h>

#include <vw/Plate/LocalIndex.h>
#include <vw/Plate/BlobManager.h>
#include <vw/Plate/ProtoBuffers.pb.h>

#include <google/sparsetable>

#include <vector>
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

  private:

    std::string m_filename;
    int m_level, m_base_col, m_base_row;
    int m_page_width, m_page_height;
    bool m_needs_saving;
    google::sparsetable<multi_value_type> m_sparse_table;

    void serialize();
    void deserialize();

    void append_record( std::list<vw::platefile::TileHeader> &results, 
                        int transaction_id, 
                        IndexRecord const& rec, int col, int row,
                        BBox2i const& region) const;

  public:
  
    /// Create or open a page file.
    IndexPage(std::string filename, 
              int level, int base_col, int base_row, 
              int page_width, int page_height);

    ~IndexPage();

    // ----------------------- ITERATORS  ----------------------

    nonempty_iterator begin() { return m_sparse_table.nonempty_begin(); }
    nonempty_iterator end() { return m_sparse_table.nonempty_end(); }

    // ----------------------- ACCESSORS  ----------------------

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

    /// Return multiple index entries that match the specified
    /// transaction id range.  This range is inclusive of the first
    /// entry, but not the last entry: [ begin_transaction_id, end_transaction_id )
    ///
    /// Results are return as a std::pair<int32, IndexRecord>.  The
    /// first value in the pair is the transaction id for that
    /// IndexRecord.
    multi_value_type multi_get(int col, int row, 
                               int begin_transaction_id, int end_transaction_id) const;

    /// Return the number of valid entries in this page.  (Remember
    /// that this is a sparse store of IndexRecords.)
    int sparse_size() { return m_sparse_table.num_nonempty(); }

    /// Returns a list of valid tiles in this IndexPage. 
    std::list<TileHeader> valid_tiles(int transaction_id, 
                                      vw::BBox2i const& region,
                                      bool exact_match) const;

    // ----------------------- DISK I/O  ----------------------

    /// Save any unsaved changes to disk.
    void sync();
  };


  // ----------------------------------------------------------------------
  //                         INDEX PAGE GENERATOR
  // ----------------------------------------------------------------------

  // IndexPageGenerator loads a index page from disk.
  class IndexPageGenerator {
    std::string m_filename;
    int m_level, m_base_col, m_base_row;
    int m_page_width, m_page_height;

  public:
    typedef IndexPage value_type;
    IndexPageGenerator( std::string filename, int level, int base_col, int base_row, 
                        int page_width, int page_height );
    size_t size() const;

    /// Generate an IndexPage.  If no file exists with the name
    /// m_filename, then an empty IndexPage is generated.
    boost::shared_ptr<IndexPage> generate() const;
  };

  // --------------------------------------------------------------------
  //                             INDEX LEVEL
  // --------------------------------------------------------------------
  class IndexLevel {
    int m_level;
    int m_page_width, m_page_height;
    int m_horizontal_pages, m_vertical_pages;
    std::vector<boost::shared_ptr<IndexPageGenerator> > m_cache_generators;
    std::vector<Cache::Handle<IndexPageGenerator> > m_cache_handles;
    vw::Cache m_cache;

  public:
    typedef IndexPage::multi_value_type multi_value_type;

    IndexLevel(std::string base_path, int level, 
               int page_width, int page_height, int cache_size);

    ~IndexLevel();

    /// Fetch the value of an index node at this level.
    IndexRecord get(int32 col, int32 row, int32 transaction_id, bool exact_match = false) const;

    /// Fetch the value of an index node at this level.
    multi_value_type multi_get(int32 col, int32 row, 
                               int32 begin_transaction_id, 
                               int32 end_transaction_id) const; 

    /// Set the value of an index node at this level.
    void set(IndexRecord const& rec, int32 col, int32 row, int32 transaction_id);

    /// Returns a list of valid tiles at this level.
    std::list<TileHeader> valid_tiles(int transaction_id, BBox2i const& region) const;
  };

  // --------------------------------------------------------------------
  //                             PAGED INDEX
  // --------------------------------------------------------------------

  class PagedIndex : public LocalIndex {
    int m_page_width, m_page_height;
    int m_default_cache_size;
    std::vector<boost::shared_ptr<IndexLevel> > m_levels;

    virtual void rebuild_index(std::string plate_filename);

    void commit_record(IndexRecord const& record, int col, int row, int level, int transaction_id);
    
  public:
    typedef IndexLevel::multi_value_type multi_value_type;

    /// Create a new, empty index.
    PagedIndex(std::string plate_filename, IndexHeader new_index_info, 
               int page_width = 256, int page_height = 256, 
               int default_cache_size = 10000);

    /// Open an existing index from a file on disk.
    PagedIndex(std::string plate_filename,
               int page_width = 256, int page_height = 256, 
               int default_cache_size = 10000);

    virtual ~PagedIndex() {}

    // ----------------------- READ/WRITE REQUESTS  ----------------------

    /// Attempt to access a tile in the index.  Throws an
    /// TileNotFoundErr if the tile cannot be found.
    /// 
    /// By default, this call to read will return a tile with the MOST
    /// RECENT transaction_id <= to the transaction_id you specify
    /// here in the function arguments (if a tile exists).  However,
    /// setting exact_transaction_match = true will force the
    /// PlateFile to search for a tile that has the EXACT SAME
    /// transaction_id as the one that you specify.
    ///
    /// A transaction ID of -1 indicates that we should return the
    /// most recent tile, regardless of its transaction id.
    virtual IndexRecord read_request(int col, int row, int level, 
                                     int transaction_id, bool exact_transaction_match = false);

    /// Return multiple index entries that match the specified
    /// transaction id range.  This range is inclusive of the first
    /// entry, but not the last entry: [ begin_transaction_id, end_transaction_id )
    ///
    /// Results are return as a std::pair<int32, IndexRecord>.  The
    /// first value in the pair is the transaction id for that
    /// IndexRecord.
    virtual std::list<std::pair<int32, IndexRecord> > multi_read_request(int col, int row, int level, 
                                                                         int begin_transaction_id, 
                                                                         int end_transaction_id);

    // Writing, pt. 1: Locks a blob and returns the blob id that can
    // be used to write a tile.
    virtual int write_request(int size);

    // Writing, pt. 2: Supply information to update the index and
    // unlock the blob id.
    virtual void write_complete(TileHeader const& header, IndexRecord const& record);
  
    // ----------------------- PROPERTIES  ----------------------

    /// Returns a list of tile headers for any valid tiles that exist
    /// at a the specified level and transaction_id.  The
    /// transaction_id is treated the same as it would be for
    /// read_request() above.  The region specifies a tile range of
    /// interest.
    virtual std::list<TileHeader> valid_tiles(int level, int transaction_id, 
                                              BBox2i const& region) const;

    virtual int32 num_levels() const;
    void map(boost::shared_ptr<TreeMapFunc> func);
  };

}} // namespace vw::platefile

#endif // __VW_PLATEFILE_PAGED_INDEX_H__
