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

#include <vw/Plate/Index.h>
#include <vw/Plate/BlobManager.h>
#include <vw/Plate/ProtoBuffers.pb.h>

namespace vw {
namespace platefile {

  // -------------------------------------------------------------------
  //                            LOCAL INDEX
  // -------------------------------------------------------------------

  class LocalIndex : public Index { 
  
  protected:
    std::string m_plate_filename;
    IndexHeader m_header;
    boost::shared_ptr<BlobManager> m_blob_manager;
    boost::shared_ptr<vw::LogInstance> m_log;

    void save_index_file() const;
    std::string index_filename() const;
    std::string log_filename() const;
    std::vector<std::string> blob_filenames() const;
    virtual void rebuild_index(std::string plate_filename) = 0;

  public:

    /// Create a new, empty index.
    LocalIndex( std::string plate_filename, IndexHeader new_index_info);

    /// Open an existing index from a file on disk.
    LocalIndex(std::string plate_filename);

    /// Destructor
    virtual ~LocalIndex() {}

    /// Use this to send data to the index's logfile like this:
    ///
    ///   index_instance.log() << "some text for the log...\n";
    ///
    std::ostream& log ();

    virtual IndexHeader index_header() const { return m_header; }

    // /// Save an index out to a file on disk.  This serializes the
    // /// tree.
    // virtual void save(std::string const& filename);

    virtual int version() const { return m_header.version(); }
    
    virtual std::string platefile_name() const { return m_plate_filename; }

    virtual int32 tile_size() const { return m_header.tile_size(); }
    virtual std::string tile_filetype() const { return m_header.tile_filetype(); }

    virtual PixelFormatEnum pixel_format() const { 
      return PixelFormatEnum(m_header.pixel_format()); 
    }

    virtual ChannelTypeEnum channel_type() const {
      return ChannelTypeEnum(m_header.channel_type());
    }

    // Clients are expected to make a transaction request whenever
    // they start a self-contained chunk of mosaicking work.  .
    virtual int32 transaction_request(std::string transaction_description,
                                      std::vector<TileHeader> const& tile_headers);

    // Once a chunk of work is complete, clients can "commit" their
    // work to the mosaic by issuding a transaction_complete method.
    virtual void transaction_complete(int32 transaction_id);

    // If a transaction fails, we may need to clean up the mosaic.  
    virtual void transaction_failed(int32 transaction_id);

    // Return the current location of the transaction cursor.  This
    // will be the last transaction id that refers to a coherent
    // version of the mosaic.
    virtual int32 transaction_cursor();
  };

}} // namespace vw::plate

#endif // __VW_PLATEFILE_LOCAL_INDEX_H__
