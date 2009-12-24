// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_PLATEFILE_PAGED_INDEX_H__
#define __VW_PLATEFILE_PAGED_INDEX_H__

namespace vw {
namespace platefile {

  // --------------------------------------------------------------------
  //                             INDEX LEVEL
  // --------------------------------------------------------------------
  class IndexLevel {
    int m_level;
    int m_cache_size;

  public:
    IndexLevel(int level, int cache_size) : m_level(level), m_cache_size(cache_size) {
      int num_pages = pow(2,level*2);
    }
  };

  // --------------------------------------------------------------------
  //                             PAGED INDEX
  // --------------------------------------------------------------------

  class PagedIndex : public Index {
  public:
    virtual ~PagedIndex() {}

    // ----------------------- READ/WRITE REQUESTS  ----------------------

    virtual IndexRecord read_request(int col, int row, int depth, 
                                     int transaction_id, bool exact_transaction_match = false) {}

    virtual int write_request(int size) {}

    virtual void write_complete(TileHeader const& header, IndexRecord const& record) {}

  
    // ----------------------- PROPERTIES  ----------------------

    virtual IndexHeader index_header() const = 0;
  
    virtual int32 version() const = 0;
    virtual int32 max_depth() const = 0;
  
    virtual std::string platefile_name() const = 0;
  
    virtual int32 tile_size() const = 0;
    virtual std::string tile_filetype() const = 0;
  
    virtual PixelFormatEnum pixel_format() const = 0;
    virtual ChannelTypeEnum channel_type() const = 0;

    // --------------------- TRANSACTIONS ------------------------
  
    /// Clients are expected to make a transaction request whenever
    /// they start a self-contained chunk of mosaicking work.  .
    virtual int32 transaction_request(std::string transaction_description,
                                      std::vector<TileHeader> const& tile_headers) = 0;
  
    /// Called right before the beginning of the mipmapping pass
    virtual void root_complete(int32 transaction_id,
                               std::vector<TileHeader> const& tile_headers) = 0;
  
    /// Once a chunk of work is complete, clients can "commit" their
    /// work to the mosaic by issuding a transaction_complete method.
    virtual void transaction_complete(int32 transaction_id) = 0;
  
    // If a transaction fails, we may need to clean up the mosaic.  
    virtual void transaction_failed(int32 transaction_id) = 0;
  
    virtual int32 transaction_cursor() = 0;

    // --------------------- UTILITIES ------------------------
  
    /// Iterate over all nodes in a tree, calling func for each
    /// location.  Note: this will only be implemented for local
    /// indexes.  This function will throw an error if called on a
    /// remote index.
    virtual void map(boost::shared_ptr<TreeMapFunc> func) { 
      vw_throw(NoImplErr() << "Index::map() not implemented for this index type.");
    }
  };

}} // namespace vw::platefile

#endif // __VW_PLATEFILE_PAGED_INDEX_H__
