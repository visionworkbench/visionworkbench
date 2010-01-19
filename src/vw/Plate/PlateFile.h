// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

/// \file PlateFile.h
///
/// Once written, tiles in the plate file are never deleted.  New tiles
/// can supercede old ones, but all of the data remains in the blob files,
/// and the index maintains a sorted list of pointers to previous versions
/// of a tile.  The sorting key for tiles is a field in each index entry
/// called 'transaction_id', which works much like a timestamp to keep
/// track of the order with which tiles were added to the platefile
/// system.
///
/// transaction_id is not simply a timestamp.  Instead, it is a unique
/// integer assigned to each new image that is added to the mosaic.
/// transaction_ids are incremented as each new image is added.
///
/// READING TILES using the transaction_id
/// --------------------------------------
///
/// Read requests can come from the mod_plate apache module (in response
/// to a web request) or from a mosaicking client (in the process of
/// compositing new tiles in the mosaic).  In both cases, it is important
/// for a read request to return tiles from a consistent version of the
/// mosaic.  In other words, we don't want to read tiles from a region
/// where tiles are being actively written -- doing so would result in a
/// patchwork of new and old tiles, and that would look very bad(tm).
///
/// transaction_ids allow you to go back to the latest consistent version
/// (or any previous version) of the mosaic.  Each read request can
/// specify an transaction_id, thereby forcing the server to access the
/// index as it would have appeared ON OR BEFORE the requested
/// transaction_id.  The index server keeps track of the latest "complete"
/// transaction ID (lets call it the transaction_cursor), and will furnish
/// that to clients if they are simply interested in accessing the latest
/// complete version of the mosaic.
///
///
/// WRITING TILES using the transaction_id
/// --------------------------------------
///
/// transaction_id's are assigned on a per-image (i.e. per transaction)
/// basis in the mosaic.  Prior to writing any tiles, the mosaicking
/// client is expected to send a TRANSACTION_REQUEST message to obtain a
/// transaction_id from the index server.  (Note to self: that request
/// should contain some indentifying information that is written to the
/// index log so that we can reconstruct what happened later if the
/// transaction fails.)  The index server doles these out in increasing
/// order so that each client gets a unique transaction id, and that these
/// ids are issued in the order that they were requested.
///
/// Mosaicking clients should also request the transaction_id pointed to
/// by the transaction_cursor (see above) so that it has a consistent
/// version of the mosaic to use for its read requests during compositing.
/// Both the read_transaction_id and write_transaction_id are stored for
/// use in subsequent mosaicking steps.
///
/// When the mosaicking client completes its work for a given transaction,
/// it is expected to send a TRANSACTION_COMPLETE message to inform the
/// index server of its succesful completion.  The index server uses this
/// information to update the transaction_cursor that can be used for read
/// requests.
///
///
/// CAVEATS
/// -------
///
/// - What happens if a mosaicking client fails in mid-mosaic and never
///   issues a TRANSACTION_COMPLETE message?
///
///   This is bad!!  The transaction_cursor will get stuck waiting for the
///   half-finished transaction ID.
///
///   Possbile solutions: keep track of the number of completed
///  transactions, and force the transaction ID to complete if more than
///   a certain number of transaction have elapsed.  (i.e. assume it has
///   failed and note this failure in the index server log.)
///
///   (But then we should make sure that a TRANSACTION_COMPLETE message
///   that comes back from the dead doesn't somehow mess up the
///   transaction_cursor!)
///
/// - What happens to the transaction_cursor if the index server crashes?
///
///   Transactions updates don't happen very often.  We can afford to
///   write the new transaction cursor location to the index file on
///   disk before we actually update it in the live index.
///

#ifndef __VW_PLATE_PLATEFILE_H__
#define __VW_PLATE_PLATEFILE_H__

#include <vw/Math/Vector.h>
#include <vw/Image/ImageView.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/Core/ThreadPool.h>

#include <vw/Plate/Index.h>
#include <vw/Plate/Blob.h>
#include <vw/Plate/Exception.h>

#include <vector>
#include <fstream>
#include <stdlib.h>

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/scoped_array.hpp>
namespace fs = boost::filesystem;

namespace vw {
namespace platefile {

  // -------------------------------------------------------------------------
  //                            TEMPORARY TILE FILE
  // -------------------------------------------------------------------------

  // A scoped temporary file object that store a tile under /tmp.  The
  // file is destroyed when this object is deleted.
  class TemporaryTileFile {

    std::string m_filename;

    // Define these as private methods to enforce TemporaryTileFile'
    // non-copyable semantics.
    TemporaryTileFile() {}
    TemporaryTileFile(TemporaryTileFile const&) {}
    TemporaryTileFile& operator=(TemporaryTileFile const&) { return *this; }
  
  public:
    
    /// Generate a unique filename ( usually in /tmp, though this can
    /// be overridden using vw_settings().tmp_directory() ).
    static std::string unique_tempfile_name(std::string file_extension);

    /// This constructor assumes control over an existing file on disk,
    /// and deletes it when the TemporaryTileFile object is de-allocated.
    TemporaryTileFile(std::string filename);

    /// This constructor assumes control over an existing file on disk,
    /// and deletes it when the TemporaryTileFile object is de-allocated.
    template <class ViewT>
    TemporaryTileFile(ImageViewBase<ViewT> const& view, std::string file_extension) : 
      m_filename(unique_tempfile_name(file_extension)) {
      write_image(m_filename, view);
      vw_out(DebugMessage, "plate::tempfile") << "Created temporary file: " 
                                              << m_filename << "\n";
    }

    ~TemporaryTileFile();

    std::string file_name() const { return m_filename; }

    /// Opens the temporary file and determines its size in bytes.
    int64 file_size() const;

    // Read an image from the temporary tile file.
    template <class PixelT> 
    ImageView<PixelT> read() const {
      ImageView<PixelT> img;
      read_image(img, m_filename);
      return img;
    }
          
  };


  // -------------------------------------------------------------------------
  //                            PLATE FILE
  // -------------------------------------------------------------------------

  
  class PlateFile {
    boost::shared_ptr<Index> m_index;
    FifoWorkQueue m_queue;
    boost::shared_ptr<Blob> m_write_blob;
    int m_write_blob_id;

  public:
    PlateFile(std::string url);

    PlateFile(std::string url, std::string type, std::string description,
              int tile_size, std::string tile_filetype, 
              PixelFormatEnum pixel_format, ChannelTypeEnum channel_type);

    /// The destructor saves the platefile to disk. 
    ~PlateFile() {}

    /// Returns the name of the root directory containing the plate file.
    std::string name() const { return m_index->platefile_name(); }

    /// Returns the name of the root directory containing the plate file.
    IndexHeader index_header() const { return m_index->index_header(); }

    /// Returns the file type used to store tiles in this plate file.
    std::string default_file_type() const { return m_index->tile_filetype(); }

    int default_tile_size() const { return m_index->tile_size(); }

    PixelFormatEnum pixel_format() const { return m_index->pixel_format(); }

    ChannelTypeEnum channel_type() const { return m_index->channel_type(); }

    int num_levels() const { return m_index->num_levels(); }

    /// Read the tile header. You supply a base name (without the
    /// file's image extension).  The image extension will be appended
    /// automatically for you based on the filetype in the TileHeader.
    std::string read_to_file(std::string const& base_name, 
                             int col, int row, int level, int transaction_id);

    /// Read an image from the specified tile location in the plate file.  
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
    template <class ViewT>
    TileHeader read(ViewT &view, int col, int row, int level, 
                    int transaction_id, bool exact_transaction_match = false) {

      TileHeader result;
      
      // 1. Call index read_request(col,row,level).  Returns IndexRecord.
      IndexRecord record = m_index->read_request(col, row, level, 
                                                 transaction_id, exact_transaction_match);

      // 2. Open the blob file and read the header.  If we are reading
      // from the same blob as we already have open for writing, we go
      // ahead and use that already-open file pointer.  Otherwise, we
      // open the new blob for reading.
      boost::shared_ptr<Blob> read_blob;
      // if (m_write_blob && record.blob_id() == m_write_blob_id) {
      //   read_blob = m_write_blob;
      // } else {
        std::ostringstream blob_filename;
        blob_filename << this->name() << "/plate_" << record.blob_id() << ".blob";
        read_blob.reset(new Blob(blob_filename.str(), true));
      // }
      TileHeader header = read_blob->read_header<TileHeader>(record.blob_offset());
      
      // 3. Choose a temporary filename and call BlobIO
      // read_as_file(filename, offset, size) [ offset, size from
      // IndexRecord ]
      std::string tempfile = TemporaryTileFile::unique_tempfile_name(header.filetype());
      read_blob->read_to_file(tempfile, record.blob_offset());
      TemporaryTileFile tile(tempfile);
      
      // 4. Read data from temporary file.
      view = tile.read<typename ViewT::pixel_type>();
      
      // 5. Return the tile header.
      return header;
    }

    /// Read one ore more images at a specified location in the
    /// platefile by specifying a range of transaction ids of
    /// interest.  This range is inclusive of the first entry, but not
    /// the last entry: [ begin_transaction_id, end_transaction_id )
    ///
    /// This is mostly useful when compositing tiles during mipmapping.
    ///
    template <class ViewT>
    std::list<TileHeader> multi_read(std::list<ViewT> &tiles, int col, int row, int level, 
                                     int begin_transaction_id, int end_transaction_id) {

      std::list<TileHeader> results;
      
      // 1. Call index read_request(col,row,level).  Returns IndexRecord.
      std::list<std::pair<int32, IndexRecord> > records = m_index->multi_read_request(col, row, 
                                                                                      level, 
                                                                                      begin_transaction_id, 
                                                                                      end_transaction_id);

      std::list<std::pair<int32, IndexRecord> >::iterator iter = records.begin();
      while (iter != records.end()) {
        IndexRecord &record = iter->second;
          
        // 2. Open the blob file and read the header
        boost::shared_ptr<Blob> read_blob;
        // if (m_write_blob && record.blob_id() == m_write_blob_id) {
        //   read_blob = m_write_blob;
        // } else {
          std::ostringstream blob_filename;
          blob_filename << this->name() << "/plate_" << record.blob_id() << ".blob";
          read_blob.reset(new Blob(blob_filename.str(), true));
          //        }
        TileHeader header = read_blob->read_header<TileHeader>(record.blob_offset());
          
        // 3. Choose a temporary filename and call BlobIO
        // read_as_file(filename, offset, size) [ offset, size from
        // IndexRecord ]
        std::string tempfile = TemporaryTileFile::unique_tempfile_name(header.filetype());
        read_blob->read_to_file(tempfile, record.blob_offset());
        TemporaryTileFile tile_file(tempfile);
        
        // 4. Read data from temporary file.
        ViewT tile = tile_file.read<typename ViewT::pixel_type>();
        tiles.push_back(tile);
        
        // 5. Access the tile header and return it.
        results.push_back( header );
        
        ++iter;
      }
      return results;
    }

    
    /// Writing, pt. 1: Locks a blob and returns the blob id that can
    /// be used to write tiles.
    void write_request();

    /// Writing, pt. 2: Write an image to the specified tile location
    /// in the plate file.  Returns the size (in bytes) written by
    /// this write_update.
    template <class ViewT>
    void write_update(ImageViewBase<ViewT> const& view, 
                      int col, int row, int level, int transaction_id) {      

      if (!m_write_blob)
        vw_throw(BlobIoErr() << "Error issueing write_update().  No blob file open.  "
                 << "Are you sure your ran write_request()?");

      // 0. Create a write_header
      TileHeader write_header;
      write_header.set_col(col);
      write_header.set_row(row);
      write_header.set_level(level);
      write_header.set_transaction_id(transaction_id);
      write_header.set_filetype(this->default_file_type());

      // 1. Write data to temporary file. 
      TemporaryTileFile tile(view, this->default_file_type());
      std::string tile_filename = tile.file_name();

      // 3. Create a blob and call write_from_file(filename).  Returns
      // offset, size.  
      int64 blob_offset;
      m_write_blob->write_from_file(tile_filename, write_header, blob_offset);

      // 4. Call write_update(col, row, level, record) to update the
      // index with the new data.
      IndexRecord write_record;
      write_record.set_blob_id(m_write_blob_id);
      write_record.set_blob_offset(blob_offset);
      
      m_index->write_update(write_header, write_record);
    }

    /// Writing, pt. 3: Signal the completion of the write operation.
    void write_complete();

    /// Read a record out of the platefile.  
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
    IndexRecord read_record(int col, int row, int level, 
                            int transaction_id, bool exact_transaction_match = false);

    // --------------------- TRANSACTIONS ------------------------

    // Clients are expected to make a transaction request whenever
    // they start a self-contained chunk of mosaicking work.  .
    virtual int32 transaction_request(std::string transaction_description,
                                      int transaction_id_override) {
      return m_index->transaction_request(transaction_description, transaction_id_override);
    }

    // Once a chunk of work is complete, clients can "commit" their
    // work to the mosaic by issuding a transaction_complete method.
    virtual void transaction_complete(int32 transaction_id, bool update_read_cursor) {
      m_index->transaction_complete(transaction_id, update_read_cursor);
    }

    // If a transaction fails, we may need to clean up the mosaic.  
    virtual void transaction_failed(int32 transaction_id) {
      m_index->transaction_failed(transaction_id);
    }

    virtual int32 transaction_cursor() {
      return m_index->transaction_cursor();
    }
    
    // ----------------------- UTILITIES --------------------------

    /// Returns a list of valid tiles that match this level, region, and
    /// range of transaction_id's.  Returns a list of TileHeaders with
    /// col/row/level and transaction_id of the most recent tile at each
    /// valid location.  Note: there may be other tiles in the transaction
    /// range at this col/row/level, but valid_tiles() only returns the
    /// first one.
    std::list<TileHeader> valid_tiles(int level, vw::BBox2i const& region,
                                      int start_transaction_id, 
                                      int end_transaction_id, 
                                      int min_num_matches) const {
      return m_index->valid_tiles(level, region, 
                                  start_transaction_id, end_transaction_id,
                                  min_num_matches);
    }


  };



}} // namespace vw::plate

#endif // __VW_PLATE_PLATEFILE_H__
