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
#include <vw/Plate/IndexDataPrivate.pb.h>
#include <vw/Core/Exception.h>
#include <vw/Core/FundamentalTypes.h>
#include <vw/Core/Log.h>
#include <boost/shared_array.hpp>
#include <fstream>
#include <string>

namespace vw {
namespace platefile {

  typedef boost::shared_ptr<std::vector<vw::uint8> > TileData;

  struct BlobTileRecord {
    detail::BlobRecord rec;
    TileHeader hdr;
    TileData data;
  };

  class ReadBlob : boost::noncopyable {
    protected:
      std::string m_blob_filename;
      uint64 m_end_of_file_ptr;
      boost::shared_ptr<std::fstream> m_fstream;

      void read_at(uint64 offset, char* dst, uint64 size, const char* context) const;
      void check_fail(const char* context, const char* context2 = "") const;

      typedef vw::uint16 BlobRecordSizeType;
      /// Returns the metadata (i.e. BlobRecord) for a blob entry.
      detail::BlobRecord read_blob_record(const uint32& base_offset, BlobRecordSizeType &blob_record_size) const;
      TileHeader         read_tile_header(const uint32& base_offset, const detail::BlobRecord& blob_record, const BlobRecordSizeType& blob_record_size) const;
      TileData           read_tile_data  (const uint32& base_offset, const detail::BlobRecord& blob_record, const BlobRecordSizeType& blob_record_size) const;

      uint64 read_end_of_file_ptr() const;

      void init();

      // protected constructor so WriteBlob can do its own initialization
      ReadBlob(const std::string& filename, bool skip_init);
    public:
      class iterator : public boost::iterator_facade<iterator, const BlobTileRecord, boost::forward_traversal_tag, const BlobTileRecord> {
        private:
          // This is required for boost::iterator_facade
          friend class boost::iterator_core_access;
          // Private variables
          ReadBlob* m_blob;
          uint64 m_current_base_offset;

          bool equal (iterator const& iter) const;
          void increment();
          BlobTileRecord dereference() const;
        public:
          iterator( ReadBlob *blob, uint64 base_offset );
          uint64 current_base_offset() const;
      };

      typedef iterator const_iterator;

      explicit ReadBlob(const std::string& filename);
      ~ReadBlob();

      /// Returns the size of the blob in bytes.  Note: only counts
      /// valid entries.  (Invalid data may exist beyond the end of the
      /// end_of_file_ptr)
      uint64 size() const { return m_end_of_file_ptr; }

      const std::string& filename() const;

      /// Returns an iterator pointing to the first TileHeader in the blob.
      ///
      /// 3*sizeof(uint64) is the very first byte in the file after the
      /// end-of-file pointer.  (See the *_end_of_file_ptr() routines
      /// above for more info...)
      iterator begin() { return iterator(this, 3*sizeof(uint64) ); }

      /// Returns an iterator pointing one past the last TileHeader in the blob.
      iterator end() { return iterator(this, m_end_of_file_ptr ); }

      /// Seek to the next base offset given the current base offset.
      uint64 next_base_offset(uint64 current_base_offset);

      /// Returns binary index record (a serialized protobuffer) for an
      /// entry starting at base_offset.
      TileHeader read_header(vw::uint64 base_offset64);

      /// Returns the binary data for an entry starting at base_offset.
      TileData read_data(vw::uint64 base_offset);

      /// Returns the whole blob record (this is faster than calling read_header then real_tile_data)
      BlobTileRecord read_record(vw::uint64 base_offset);

      /// Returns the parameters necessary to call sendfile(2)
      void read_sendfile(vw::uint64 base_offset, std::string& filename, vw::uint64& offset, vw::uint64& size);

      /// Returns the data size
      uint64 data_size(uint64 base_offset) const;
  };

  class Blob : public ReadBlob {
    private:
      uint64 m_write_count;
      void write_end_of_file_ptr(uint64 ptr);
    public:
      explicit Blob(const std::string& filename);

      /// The destructor flushes any unwritten journal entries and
      /// closes the blob and journal files.
      ~Blob();

      /// Write a tile to the blob file. You must supply the header
      /// (e.g. a serialized TileHeader protobuffer) and the data as
      /// shared_arrays.  Returns the base_offset where the data was
      /// written to the blob file.
      vw::uint64 write(TileHeader const& header, const uint8* data, uint64 data_size);

      // Flush all pending changes
      void flush();
  };

}} // namespace vw::platefile

#endif // __VW_PLATE_BLOBIO__
