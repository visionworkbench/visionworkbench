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


#ifndef __VW_PLATE_REMOTE_INDEX_H__
#define __VW_PLATE_REMOTE_INDEX_H__

#include <vw/Plate/detail/PagedIndex.h>
#include <vw/Plate/detail/IndexPage.h>
#include <vw/Plate/HTTPUtils.h>
#include <queue>

namespace vw {
namespace platefile {

  template <typename ServiceT>
  class RpcClient;

  class IndexService;
  class IndexWriteUpdate;

  typedef RpcClient<IndexService> IndexClient;

namespace detail {

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
                    uint32 level, uint32 base_col, uint32 base_row,
                    uint32 page_width, uint32 page_height);

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
    uint32 m_level, m_base_col, m_base_row;
    uint32 m_page_width, m_page_height;

  public:
    RemotePageGenerator( int platefile_id,
                         boost::shared_ptr<IndexClient> client,
                         uint32 level, uint32 base_col, uint32 base_row,
                         uint32 page_width, uint32 page_height );
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
    create(uint32 level, uint32 base_col, uint32 base_row, uint32 page_width, uint32 page_height);

    virtual std::string who() const;
  };

  // -------------------------------------------------------------------
  //                            REMOTE INDEX
  // -------------------------------------------------------------------

  class RemoteIndex : public PagedIndex {
    vw::platefile::Url m_url;
    int m_platefile_id;
    std::string m_short_plate_filename;
    std::string m_full_plate_filename;
    mutable IndexHeader m_index_header;

    // Remote connection
    boost::shared_ptr<IndexClient> m_client;

    // Log streamer
    struct LogRequestSink;
    friend class RemoteIndex::LogRequestSink;
    boost::shared_ptr<std::ostream> m_logger;

    void update_header() const;

  public:
    /// Constructor (for opening an existing index)
    RemoteIndex(const Url& url);

    /// Constructor (for creating a new index)
    RemoteIndex(const Url& url, IndexHeader new_index_info);

    /// destructor
    virtual ~RemoteIndex();

    // Writing, pt. 1: Locks a blob and returns the blob id that can
    // be used to write a tile.
    virtual uint32 write_request();

    /// Writing, pt. 3: Signal the completion
    virtual void write_complete(uint32 blob_id);

    /// Log a message to the platefile log.
    virtual std::ostream& log();

    virtual IndexHeader index_header() const;

    virtual uint32 platefile_id() const;
    virtual uint32 version() const;
    virtual uint32 num_levels() const;

    virtual std::string platefile_name() const;

    virtual uint32 tile_size() const;
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

}}}

#endif // __VW_PLATE_REMOTE_INDEX_H__
