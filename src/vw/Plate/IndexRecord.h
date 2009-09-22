#ifndef __VW_PLATEFILE_INDEXRECORD_H__
#define __VW_PLATEFILE_INDEXRECORD_H__

#include <vw/Core/FundamentalTypes.h>
#include <vw/Core/Exception.h>

#define VW_PLATE_INDEXRECORD_FILETYPE_SIZE 5
#define VW_PLATE_INDEX_VERSION 2

namespace vw {
namespace platefile {

  // -------------------------------------------------------------------
  //                          INDEX_RECORD
  // 
  // The IndexRecord stores all of the metadata needed to recover an
  // image from the blob.  It will be marked as valid if it contains
  // good data, or as invalid if it does not.  Records may be marked
  // as invalid if they represent low-res tiles that need to be
  // re-renedered from higher-res tiles that have been updated.
  // -------------------------------------------------------------------

  class IndexRecord {

    int32 m_blob_id;
    int64 m_blob_offset;
    int32 m_block_size;
    uint8 m_block_filetype[VW_PLATE_INDEXRECORD_FILETYPE_SIZE];
    uint8 m_valid;

  public:

    IndexRecord() : m_valid(false) {}

    IndexRecord(int32 blob_id, int64 blob_offset, 
                int32 block_size, std::string block_filetype) :
      m_blob_id(blob_id), m_blob_offset(blob_offset), 
      m_block_size(block_size), m_valid(true) {

      if (block_filetype.size() > 4) 
        vw_throw(ArgumentErr() << "IndexRecord: filetype argument must be 4 characters or fewer.");
      strncpy((char*)m_block_filetype, block_filetype.c_str(), 5);

    }

    IndexRecord(std::istream &istr) {
      this->deserialize(istr);
    }


    int32 blob_id() const { return m_blob_id; }
    void set_blob_id(int32 blob_id) { m_blob_id = blob_id; }
    
    int64 blob_offset() const { return m_blob_offset; }
    void set_blob_offset(int64 blob_offset) { m_blob_offset = blob_offset; }

    int32 block_size() const { return m_block_size; }
    void set_block_size(int32 block_size) { m_block_size = block_size; }

    std::string block_filetype() const { return (char*)m_block_filetype; }
    void set_block_filetype(std::string block_filetype) {
      if (block_filetype.size() > 4) 
        vw_throw(ArgumentErr() << "IndexRecord: filetype argument must be 4 characters or fewer.");
      strncpy((char*)m_block_filetype, block_filetype.c_str(), 5);
    }

    bool valid() const { return m_valid; }
    void validate() { m_valid = true; }
    void invalidate() { m_valid = false; }

    /// Serialize the index record as a series of bytes.
    void serialize(std::ostream &ostr) const {
      ostr.write( (char*)&m_blob_id, sizeof(m_blob_id));
      ostr.write( (char*)&m_blob_offset, sizeof(m_blob_offset));
      ostr.write( (char*)&m_block_size, sizeof(m_block_size));
      uint8 filetype_size = VW_PLATE_INDEXRECORD_FILETYPE_SIZE;
      ostr.write( (char*)&filetype_size , sizeof(filetype_size));
      for (int i = 0; i < filetype_size; ++i) 
        ostr.write( (char*)&(m_block_filetype[i]), sizeof(*m_block_filetype));
      ostr.write( (char*)&m_valid, sizeof(m_valid));
    }

    /// Deserialize the index record from a stream of bytes
    void deserialize(std::istream &istr) {
      istr.read( (char*)&m_blob_id, sizeof(m_blob_id) );
      istr.read( (char*)&m_blob_offset, sizeof(m_blob_offset) );
      istr.read( (char*)&m_block_size, sizeof(m_block_size) );
      uint8 filetype_size;
      istr.read( (char*)&filetype_size, sizeof(filetype_size) );
      for (int i = 0; i < filetype_size; ++i) 
        istr.read( (char*)&(m_block_filetype[i]), sizeof(*m_block_filetype));
      istr.read( (char*)&m_valid, sizeof(m_valid));
    }

  };

}} // namespace vw::platefile

#endif // __VW_PLATE_INDEXRECORD__
