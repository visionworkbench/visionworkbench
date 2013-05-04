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


#ifndef __VW_PLATEFILE_LOCAL_INDEX_H__
#define __VW_PLATEFILE_LOCAL_INDEX_H__

#include <vw/Plate/FundamentalTypes.h>
#include <vw/Image/PixelTypeInfo.h>
#include <vw/Plate/detail/PagedIndex.h>

namespace vw {
namespace platefile {

  class BlobManager;
  class TileHeader;

namespace detail {

  // ----------------------------------------------------------------------
  //                         LOCAL INDEX PAGE
  // ----------------------------------------------------------------------

  class LocalIndexPage : public IndexPage {
    std::string m_filename;
    bool m_needs_saving;

    // For reading/writing to/from disk.
    void serialize();
    void deserialize();

  public:
    /// Create or open a page file.
    LocalIndexPage(std::string filename,
                   uint32 level, uint32 base_col, uint32 base_row,
                   uint32 page_width, uint32 page_height);

    virtual ~LocalIndexPage();

    /// Set the value of an entry in the IndexPage.
    virtual void set(TileHeader const& header, IndexRecord const& record);

    /// Save any unsaved changes to disk.
    virtual void sync();

  };

  // ----------------------------------------------------------------------
  //                       LOCAL INDEX PAGE GENERATOR
  // ----------------------------------------------------------------------

  // loads a index page from disk.
  class LocalPageGenerator : public PageGeneratorBase {
    std::string m_filename;
    uint32 m_level, m_base_col, m_base_row;
    uint32 m_page_width, m_page_height;

  public:
    LocalPageGenerator( std::string filename, uint32 level, uint32 base_col, uint32 base_row,
                        uint32 page_width, uint32 page_height );
    virtual ~LocalPageGenerator() {}

    /// Generate an IndexPage.  If no file exists with the name
    /// m_filename, then an empty IndexPage is generated.
    virtual boost::shared_ptr<IndexPage> generate() const;
  };

  /// The LocalPageGeneratorFactory creates a generator that can
  /// produce pages from a file on disk.
  class LocalPageGeneratorFactory : public PageGeneratorFactory {
    std::string m_plate_filename;

  public:
    LocalPageGeneratorFactory(std::string plate_filename) :
      m_plate_filename(plate_filename) {}
    virtual ~LocalPageGeneratorFactory() {}

    virtual boost::shared_ptr<PageGeneratorBase> create(uint32 level, uint32 base_col, uint32 base_row,
                                                        uint32 page_width, uint32 page_height);

    virtual std::string who() const;
  };

  // -------------------------------------------------------------------
  //                            LOCAL INDEX
  // -------------------------------------------------------------------

  class LocalIndex : public PagedIndex {
    std::string m_plate_filename;
    IndexHeader m_header;
    boost::shared_ptr<BlobManager> m_blob_manager;
    boost::shared_ptr<vw::LogInstance> m_log;

    void save_index_file() const;
    std::string index_filename() const;
    std::string log_filename() const;
    std::vector<std::string> blob_filenames() const;

    void open_impl();

  public:

    /// Create a new, empty index.
    LocalIndex( std::string plate_filename, IndexHeader new_index_info );

    /// Open an existing index from a file on disk.
    LocalIndex( std::string plate_filename );

    /// Destructor
    virtual ~LocalIndex() {}

    // Rebuild an index from blob file entries.  You should only do
    // this if you lose or corrupt an index.  This may take a long
    // time.
    void rebuild_index();

    /// Use this to send data to the index's logfile like this:
    ///
    ///   index_instance.log() << "some text for the log...\n";
    ///
    std::ostream& log();

    // -----------------------    I/O      ----------------------

    // Writing, pt. 1: Locks a blob and returns the blob id that can
    // be used to write a tile.
    virtual uint32 write_request();

    // Writing, pt. 2: Supply information to update the index and
    // unlock the blob id.
    virtual void write_update(TileHeader const& header, IndexRecord const& record);

    /// Writing, pt. 3: Signal the completion
    virtual void write_complete(uint32 blob_id);

    // ----------------------- PROPERTIES  ----------------------

    virtual uint32 version() const { return m_header.version(); }
    virtual std::string platefile_name() const { return m_plate_filename; }
    virtual IndexHeader index_header() const { return m_header; }
    virtual uint32 tile_size() const { return m_header.tile_size(); }
    virtual std::string tile_filetype() const { return m_header.tile_filetype(); }
    virtual uint32 num_levels() const { return m_header.num_levels(); }
    virtual uint32 platefile_id() const {return m_header.platefile_id(); }

    virtual PixelFormatEnum pixel_format() const {
      return PixelFormatEnum(m_header.pixel_format());
    }
    virtual ChannelTypeEnum channel_type() const {
      return ChannelTypeEnum(m_header.channel_type());
    }

    // Clients are expected to make a transaction request whenever
    // they start a self-contained chunk of mosaicking work.  Use
    // transaction_id_override to force the use of a transaction ID
    // for an upcoming transaction.  Setting transaction_id_override
    // to -1 lets the platefile choose its own transaction_id.
    virtual Transaction transaction_request(std::string transaction_description,
                                            TransactionOrNeg transaction_id_override);

    // Once a chunk of work is complete, clients can "commit" their
    // work to the mosaic by issuding a transaction_complete method.
    virtual void transaction_complete(Transaction transaction_id, bool update_read_cursor);

    // If a transaction fails, we may need to clean up the mosaic.
    virtual void transaction_failed(Transaction transaction_id);

    // Return the current location of the transaction cursor.  This
    // will be the last transaction id that refers to a coherent
    // version of the mosaic.
    virtual Transaction transaction_cursor();
  };

}}}

#endif
