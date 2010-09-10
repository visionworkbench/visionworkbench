// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_PLATE_REMOTE_INDEX_H__
#define __VW_PLATE_REMOTE_INDEX_H__

#include <vw/Plate/PagedIndex.h>
#include <vw/Plate/IndexPage.h>
#include <queue>

namespace vw {
namespace platefile {

  class AmqpRpcClient;
  class IndexService;

  // ----------------------------------------------------------------------
  //                         LOCAL INDEX PAGE
  // ----------------------------------------------------------------------

  class RemoteIndexPage : public IndexPage {
    int m_platefile_id;
    boost::shared_ptr<AmqpRpcClient> m_rpc_controller;
    boost::shared_ptr<IndexService> m_index_service;

    // For packetizing write requests.
    std::queue<IndexWriteUpdate> m_write_queue;
    void flush_write_queue();

  public:

    RemoteIndexPage(int platefile_id,
                    boost::shared_ptr<AmqpRpcClient> rpc_controller,
                    boost::shared_ptr<IndexService> index_service,
                    int level, int base_col, int base_row,
                    int page_width, int page_height);

    virtual ~RemoteIndexPage();

    /// Set the value of an entry in the RemoteIndexPage.
    virtual void set(TileHeader const& header, IndexRecord const& record);

    /// Save any unsaved changes to disk.
    virtual void sync();
  };

  // ----------------------------------------------------------------------
  //                       REMOTE INDEX PAGE GENERATOR
  // ----------------------------------------------------------------------

  class RemotePageGenerator : public PageGeneratorBase {
    int m_platefile_id;
    boost::shared_ptr<AmqpRpcClient> m_rpc_controller;
    boost::shared_ptr<IndexService> m_index_service;
    int m_level, m_base_col, m_base_row;
    int m_page_width, m_page_height;

  public:
    RemotePageGenerator( int platefile_id,
                         boost::shared_ptr<AmqpRpcClient> rpc_controller,
                         boost::shared_ptr<IndexService> index_service,
                         int level, int base_col, int base_row,
                         int page_width, int page_height );
    virtual ~RemotePageGenerator() {}

    /// Generate an IndexPage.  If no IndexPage exists on the Remote
    /// Server by this name, an empty IndexPage is generated.
    virtual boost::shared_ptr<IndexPage> generate() const;
  };

  /// The RemotePageGeneratorFactory creates a generator that can
  /// produce pages from a file on disk.
  class RemotePageGeneratorFactory : public PageGeneratorFactory {
    int m_platefile_id;
    boost::shared_ptr<AmqpRpcClient> m_rpc_controller;
    boost::shared_ptr<IndexService> m_index_service;

  public:
    RemotePageGeneratorFactory() : m_platefile_id(-1) {}
    RemotePageGeneratorFactory(int platefile_id,
                               boost::shared_ptr<AmqpRpcClient> rpc_controller,
                               boost::shared_ptr<IndexService> index_service) :
      m_platefile_id(platefile_id), m_rpc_controller(rpc_controller),
      m_index_service(index_service) {}
    virtual ~RemotePageGeneratorFactory() {}
    virtual boost::shared_ptr<PageGeneratorBase> create(int level, int base_col, int base_row,
                                                        int page_width, int page_height);
    virtual std::string who() const;
  };

  // -------------------------------------------------------------------
  //                            REMOTE INDEX
  // -------------------------------------------------------------------

  class RemoteIndex : public PagedIndex {

    int m_platefile_id;
    IndexHeader m_index_header;
    std::string m_short_plate_filename;
    std::string m_full_plate_filename;

    // Remote connection
    boost::shared_ptr<AmqpRpcClient> m_rpc_controller;
    boost::shared_ptr<IndexService> m_index_service;

  public:
    /// Constructor (for opening an existing index)
    RemoteIndex(std::string const& url);

    /// Constructor (for creating a new index)
    RemoteIndex(std::string const& url, IndexHeader new_index_info);

    /// destructor
    virtual ~RemoteIndex();

    // Writing, pt. 1: Locks a blob and returns the blob id that can
    // be used to write a tile.
    virtual int write_request(uint64 &size);

    /// Writing, pt. 3: Signal the completion
    virtual void write_complete(int blob_id, uint64 blob_offset);

    /// Log a message to the platefile log.
    virtual void log(std::string message);

    // /// Returns a list of valid tiles at this level.
    // virtual std::list<TileHeader> valid_tiles(int level, BBox2i const& region,
    //                                           int start_transaction_id,
    //                                           int end_transaction_id,
    //                                           int min_num_matches) const;

    // ----------------------- PROPERTIES  ----------------------

    // /// Returns a list of tile headers for any valid tiles that exist
    // /// at a the specified level and transaction_id.  The
    // /// transaction_id is treated the same as it would be for
    // /// read_request() above.  The region specifies a tile range of
    // /// interest.
    // virtual std::list<TileHeader> valid_tiles(int level, BBox2i const& region,
    //                                           int begin_transaction_id,
    //                                           int end_transaction_id,
    //                                           int min_num_matches) const;

    virtual IndexHeader index_header() const;

    virtual int32 version() const;
    virtual int32 num_levels() const;

    virtual std::string platefile_name() const;

    virtual int32 tile_size() const;
    virtual std::string tile_filetype() const;

    virtual PixelFormatEnum pixel_format() const;
    virtual ChannelTypeEnum channel_type() const;

    // --------------------- TRANSACTIONS ------------------------

    // Clients are expected to make a transaction request whenever
    // they start a self-contained chunk of mosaicking work.  .
    virtual int32 transaction_request(std::string transaction_description,
                                      int transaction_id_override);

    // Once a chunk of work is complete, clients can "commit" their
    // work to the mosaic by issuding a transaction_complete method.
    virtual void transaction_complete(int32 transaction_id, bool update_read_cursor);

    // If a transaction fails, we may need to clean up the mosaic.
    virtual void transaction_failed(int32 transaction_id);

    virtual int32 transaction_cursor();

  };

}} // namespace vw::plate

#endif // __VW_PLATE_REMOTE_INDEX_H__
