// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_PLATEFILE_LOCAL_INDEX_H__
#define __VW_PLATEFILE_LOCAL_INDEX_H__

#include <vw/Core/FundamentalTypes.h>
#include <vw/Core/Log.h>
#include <vw/Image/PixelTypeInfo.h>

#include <vw/Plate/PagedIndex.h>
#include <vw/Plate/BlobManager.h>
#include <vw/Plate/ProtoBuffers.pb.h>

namespace vw {
namespace platefile {

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

  protected:
    virtual void commit_record(IndexRecord const& record, int col, int row, 
                               int level, int transaction_id);
    

  public:

    /// Create a new, empty index.
    LocalIndex( std::string plate_filename, IndexHeader new_index_info );

    /// Open an existing index from a file on disk.
    LocalIndex( std::string plate_filename );

    /// Destructor
    virtual ~LocalIndex() {}
    
    // Rebuild an index from blob file entries
    virtual void rebuild_index(std::string plate_filename);

    /// Use this to send data to the index's logfile like this:
    ///
    ///   index_instance.log() << "some text for the log...\n";
    ///
    std::ostream& log();

    // -----------------------    I/O      ----------------------

    // Writing, pt. 1: Locks a blob and returns the blob id that can
    // be used to write a tile.
    virtual int write_request(int size);

    /// Writing, pt. 3: Signal the completion 
    virtual void write_complete(int blob_id, uint64 blob_offset);

    // ----------------------- PROPERTIES  ----------------------

    virtual int version() const { return m_header.version(); }
    virtual std::string platefile_name() const { return m_plate_filename; }
    virtual IndexHeader index_header() const { return m_header; }
    virtual int32 tile_size() const { return m_header.tile_size(); }
    virtual std::string tile_filetype() const { return m_header.tile_filetype(); }
    virtual int32 num_levels() const { return m_header.num_levels(); }

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
    virtual int32 transaction_request(std::string transaction_description,
                                      int transaction_id_override);

    // Once a chunk of work is complete, clients can "commit" their
    // work to the mosaic by issuding a transaction_complete method.
    virtual void transaction_complete(int32 transaction_id, bool update_read_cursor);

    // If a transaction fails, we may need to clean up the mosaic.  
    virtual void transaction_failed(int32 transaction_id);

    // Return the current location of the transaction cursor.  This
    // will be the last transaction id that refers to a coherent
    // version of the mosaic.
    virtual int32 transaction_cursor();
  };

}} // namespace vw::plate

#endif // __VW_PLATEFILE_LOCAL_INDEX_H__
