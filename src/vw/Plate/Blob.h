// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_PLATE_BLOBIO__
#define __VW_PLATE_BLOBIO__

/// \file Blob.h
///
/// Platefile data is stored in data "blobs" that exist as files on the
/// filesystem.  New data is added to the data blobs whenever a new tile
/// is written to the platefile.  To ensure data consistency (as one would
/// have in a journaling filesystem), the new data is first written to a
/// journal file (a sidecar that exists alongside the blob), and then it
/// is copied from the journal file to the blob file.  In this way, the
/// data will always be written to the blob if it is succesfully written
/// to the journal, and if the data is not successfully written to the
/// journal, then it does not end up corrupting the blob. (The actual
/// process is considerably more complicated and optimized, but this is
/// the general idea.)
///
/// Data in the blob is stored in stanzas with the following way:
///
///   [ BLOB HEADER_SIZE ]  [ uint16 ]
///   [ BLOB HEADER ]       [ uint8 - serialized BlobHeader protobuffer ]
///     (contains HEADER_OFFSET, DATA_OFFSET, HEADER_SIZE, DATA_SIZE)
///
///   [ HEADER ]            [ uint8 - serialized IndexRecord protobuffer ]
///
///   [ DATA ]              [ uint8 - N raw bytes of data ]
///
///

#include <vw/Plate/IndexData.pb.h>
#include <vw/Core/Exception.h>
#include <vw/Core/FundamentalTypes.h>
#include <vw/Core/Log.h>
#include <boost/shared_array.hpp>
#include <fstream>
#include <string>

namespace vw {
namespace platefile {

  // -------------------------------------------------------------------
  //                                 BLOB
  // -------------------------------------------------------------------

  class Blob : boost::noncopyable {

    std::string m_blob_filename;
    boost::shared_ptr<std::fstream> m_fstream;
    uint64 m_end_of_file_ptr;
    uint64 m_write_count;

    /// Returns the metadata (i.e. BlobRecord) for a blob entry.
    BlobRecord read_blob_record(uint16 &blob_record_size) const;

    // End-of-file point manipulation.
    void write_end_of_file_ptr(uint64 ptr);
    uint64 read_end_of_file_ptr() const;

  public:

    // ------------------------ Blob Iterator ------------------------------

    /// An STL-compliant iterator for iterating over the index entries
    /// in a blob.  This can be used by Index.h to read in and rebuild
    /// the Index tree.
    class iterator : public boost::iterator_facade<Blob::iterator, TileHeader,
                                                   boost::forward_traversal_tag,
                                                   TileHeader, int64> {

      // This is required for boost::iterator_facade
      friend class boost::iterator_core_access;

      // Private variables
      Blob& m_blob;
      uint64 m_current_base_offset;

      // Iterator methods.  The boost iterator facade takes these and
      // uses them to construct normal iterator methods.
      bool equal (iterator const& iter) const {
        return (m_current_base_offset == iter.m_current_base_offset);
      }

      void increment() {
        m_current_base_offset = m_blob.next_base_offset(m_current_base_offset);
      }

      TileHeader const dereference() const {
        return m_blob.read_header<TileHeader>(m_current_base_offset);
      }

    public:

      // Constructors
      iterator( Blob &blob, uint64 base_offset )
        : m_blob(blob), m_current_base_offset(base_offset) {}

      uint64 current_base_offset() const { return m_current_base_offset; }
      uint64 current_data_size() const { return m_blob.data_size(m_current_base_offset); }

    };

    // -----------------------------------------------------------------------

    /// Constructor
    Blob(std::string filename, bool readonly = false);

    /// The destructor flushes any unwritten journal entries and
    /// closes the blob and journal files.
    ~Blob();

    /// Returns the size of the blob in bytes.  Note: only counts
    /// valid entries.  (Invalid data may exist beyond the end of the
    /// end_of_file_ptr)
    uint64 size() const { return m_end_of_file_ptr; }

    /// Returns an iterator pointing to the first TileHeader in the blob.
    ///
    /// 3*sizeof(uint64) is the very first byte in the file after the
    /// end-of-file pointer.  (See the *_end_of_file_ptr() routines
    /// above for more info...)
    iterator begin() { return iterator(*this, 3*sizeof(uint64) ); }

    /// Returns an iterator pointing one past the last TileHeader in the blob.
    iterator end() { return iterator(*this, m_end_of_file_ptr ); }

    /// Seek to the next base offset given the current base offset.
    uint64 next_base_offset(uint64 current_base_offset);

    /// Returns binary index record (a serialized protobuffer) for an
    /// entry starting at base_offset.
    template <class ProtoBufT>
    ProtoBufT read_header(vw::uint64 base_offset64) {

      vw_out(VerboseDebugMessage, "platefile::blob")
        << "Entering read_header() -- " <<" base_offset: " <<  base_offset64 << "\n";

      std::streamoff base_offset = boost::numeric_cast<std::streamoff>(base_offset64);

      // Seek to the requested offset and read the header and data offset
      m_fstream->seekg(base_offset, std::ios_base::beg);

      // Read the blob record
      uint16 blob_record_size;
      BlobRecord blob_record = this->read_blob_record(blob_record_size);

      // The overall blob metadata includes the uint16 of the
      // blob_record_size in addition to the size of the blob_record
      // itself.  The offsets stored in the blob_record are relative to
      // the END of the blob_record.  We compute this offset here.
      // This type-size juggling is to make sure we end up with sane behavior
      // on both 32-bit and 64-bit
      uint32 blob_offset_metadata = boost::numeric_cast<uint32>(sizeof(blob_record_size) + blob_record_size);
      size_t size = boost::numeric_cast<size_t>(blob_record.header_size());
      uint64 offset64 = base_offset + blob_offset_metadata + blob_record.header_offset();

      std::streamoff offset = boost::numeric_cast<std::streamoff>(offset64);

      // Allocate an array of the appropriate size to read the data.
      boost::shared_array<uint8> data(new uint8[size]);

      vw_out(VerboseDebugMessage, "platefile::blob")
        << "\tread_header() -- data offset: " << offset << " size: " << size << "\n";

      m_fstream->seekg(offset, std::ios_base::beg);
      m_fstream->read(reinterpret_cast<char*>(data.get()), size);

      // Throw an exception if the read operation failed (after clearing the error bit)
      if (m_fstream->fail()) {
        m_fstream->clear();
        vw_throw(IOErr() << "Blob::read() -- an error occurred while reading "
                 << "data from the blob file.\n");
      }

      // Deserialize the header
      ProtoBufT header;
      bool worked = header.ParseFromArray(data.get(), boost::numeric_cast<int>(size));
      if (!worked)
        vw_throw(IOErr() << "Blob::read() -- an error occurred while deserializing the header "
                 << "from the blob file.\n");

      vw_out(vw::VerboseDebugMessage, "platefile::blob")
        << "\tread_header() -- read " << size << " bytes at " << offset << " from " << m_blob_filename << "\n";

      return header;
    }

    /// Returns the binary data for an entry starting at base_offset.
    boost::shared_array<uint8> read_data(vw::uint64 base_offset, vw::uint64& data_size);

    /// Returns the parameters necessary to call sendfile(2)
    void read_sendfile(vw::uint64 base_offset, std::string& filename, vw::uint64& offset, vw::uint64& size);

    /// Returns the data size
    uint64 data_size(uint64 base_offset) const;

    /// Write a tile to the blob file. You must supply the header
    /// (e.g. a serialized TileHeader protobuffer) and the data as
    /// shared_arrays.  Returns the base_offset where the data was
    /// written to the blob file.
    template <class ProtoBufT>
    vw::uint64 write(ProtoBufT const& header, boost::shared_array<uint8> data, uint64 data_size) {

      // Store the current offset of the end of the file.  We'll
      // return that at the end of this function.
      std::streamoff base_offset = boost::numeric_cast<std::streamoff>(m_end_of_file_ptr);
      m_fstream->seekp(base_offset, std::ios_base::beg);

      // Create the blob record and write it to the blob file.
      BlobRecord blob_record;
      blob_record.set_header_offset(0);
      blob_record.set_header_size(header.ByteSize());
      blob_record.set_data_offset(header.ByteSize());
      blob_record.set_data_size(data_size);

      // Write the actual blob record size first.  This will help us
      // read and deserialize this protobuffer later on.
      uint16 blob_record_size = boost::numeric_cast<uint16>(blob_record.ByteSize());
      m_fstream->write(reinterpret_cast<char*>(&blob_record_size), sizeof(blob_record_size));
      blob_record.SerializeToOstream(m_fstream.get());

      // Serialize the header.
      header.SerializeToOstream(m_fstream.get());

      // And write the data.
      m_fstream->write(reinterpret_cast<char*>(data.get()), boost::numeric_cast<size_t>(data_size));

      // Write the data at the end of the file and return the offset
      // of the beginning of this data file.
      vw::vw_out(vw::VerboseDebugMessage, "platefile::blob") << "Blob::write() -- writing "
                                                         << data_size
                                                         << " bytes to "
                                                         << m_blob_filename << "\n";

      // Update the in-memory copy of the end-of-file pointer
      m_end_of_file_ptr = m_fstream->tellg();

      // The write_count is used to keep track of when we last wrote
      // the end_of_file_ptr to disk.  We don't want to write this too
      // often since this will slow down IO, so we only write it every
      // 10 writes (or when the blob is deconstructed...).
      ++m_write_count;
      if (m_write_count % 10 == 0) {
        this->write_end_of_file_ptr(m_end_of_file_ptr);
      }

      // Return the base_offset
      return base_offset;
    }


    /// Read data out of the blob and save it as its own file on disk.
    void read_to_file(std::string dest_file, vw::uint64 offset);

    /// Write the data file to disk, and the concatenate it into the data blob.
    template <class ProtoBufT>
    void write_from_file(std::string source_file, ProtoBufT const& header,
                         uint64& base_offset) {

      // Open the source_file and read data from it.
      std::ifstream istr(source_file.c_str(), std::ios::binary);

      if (!istr.is_open())
        vw_throw(IOErr() << "Blob::write_from_file() -- could not open source file for reading.");

      // Seek to the end and allocate the proper number of bytes of
      // memory, and then seek back to the beginning.
      istr.seekg(0, std::ios_base::end);
      std::streamoff loc = istr.tellg();
      VW_ASSERT(loc > 0, IOErr() << "Failed to identify file length");
      istr.seekg(0, std::ios_base::beg);

      size_t data_size = static_cast<size_t>(loc);

      // Read the data into a temporary memory buffer.
      boost::shared_array<uint8> data(new uint8[data_size]);
      istr.read(reinterpret_cast<char*>(data.get()), data_size);
      istr.close();

      base_offset = this->write(header, data, data_size);
    }

  };

}} // namespace vw::platefile

#endif // __VW_PLATE_BLOBIO__
