// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_PLATE_REMOTE_INDEX_H__
#define __VW_PLATE_REMOTE_INDEX_H__

#include <vw/Plate/PagedIndex.h>
#include <vw/Plate/IndexPage.h>
#include <vw/Plate/HTTPUtils.h>
#include <queue>

namespace vw {
namespace platefile {

  template <typename ServiceT>
  class RpcClient;

  class IndexService;
  class IndexWriteUpdate;

  typedef RpcClient<IndexService> IndexClient;

  // ----------------------------------------------------------------------
  //                         LOCAL INDEX PAGE
  // ----------------------------------------------------------------------

  class RemoteIndexPage : public IndexPage {
    int m_platefile_id;
    boost::shared_ptr<IndexClient> m_client;

    // For packetizing write requests.
    std::queue<IndexWriteUpdate> m_write_queue;
    void flush_write_queue();

  public:

    RemoteIndexPage(int platefile_id,
                    boost::shared_ptr<IndexClient> client,
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
    boost::shared_ptr<IndexClient> m_client;
    int m_level, m_base_col, m_base_row;
    int m_page_width, m_page_height;

  public:
    RemotePageGenerator( int platefile_id,
                         boost::shared_ptr<IndexClient> client,
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
    boost::shared_ptr<IndexClient> m_client;

  public:
    RemotePageGeneratorFactory(int platefile_id, boost::shared_ptr<IndexClient> client)
      : m_platefile_id(platefile_id), m_client(client) {}
    virtual ~RemotePageGeneratorFactory() {}

    virtual boost::shared_ptr<PageGeneratorBase>
    create(int level, int base_col, int base_row, int page_width, int page_height);

    virtual std::string who() const;
  };

  // -------------------------------------------------------------------
  //                            REMOTE INDEX
  // -------------------------------------------------------------------

  class RemoteIndex : public PagedIndex {
    Url m_url;
    int m_platefile_id;
    std::string m_short_plate_filename;
    std::string m_full_plate_filename;
    IndexHeader m_index_header;

    // Remote connection
    boost::shared_ptr<IndexClient> m_client;

  public:
    /// Constructor (for opening an existing index)
    RemoteIndex(const Url& url);

    /// Constructor (for creating a new index)
    RemoteIndex(const Url& url, IndexHeader new_index_info);

    /// destructor
    virtual ~RemoteIndex();

    // Writing, pt. 1: Locks a blob and returns the blob id that can
    // be used to write a tile.
    virtual int write_request(uint64 &size);

    /// Writing, pt. 3: Signal the completion
    virtual void write_complete(int blob_id, uint64 blob_offset);

    /// Log a message to the platefile log.
    virtual void log(std::string message);

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
    virtual Transaction transaction_request(std::string transaction_description,
                                            TransactionOrNeg transaction_id_override);

    // Once a chunk of work is complete, clients can "commit" their
    // work to the mosaic by issuding a transaction_complete method.
    virtual void transaction_complete(Transaction transaction_id, bool update_read_cursor);

    // If a transaction fails, we may need to clean up the mosaic.
    virtual void transaction_failed(Transaction transaction_id);

    virtual Transaction transaction_cursor();

  };

}} // namespace vw::plate

#endif // __VW_PLATE_REMOTE_INDEX_H__
