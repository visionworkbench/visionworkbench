#ifndef __VW_FILEIO_MEMORYIMAGERESOURCEGDAL_H__
#define __VW_FILEIO_MEMORYIMAGERESOURCEGDAL_H__

#include <vw/FileIO/MemoryImageResource.h>
#include <boost/noncopyable.hpp>

namespace vw {

  class SrcMemoryImageResourceGDAL : public SrcMemoryImageResource, private boost::noncopyable {
      struct Data;
      mutable boost::shared_ptr<Data> m_data;

    public:
      SrcMemoryImageResourceGDAL(boost::shared_array<const uint8> buffer, size_t len);

      virtual void read( ImageBuffer const& buf, BBox2i const& bbox ) const;

      virtual ImageFormat format() const;

      virtual bool has_block_read() const  {return false;}
      virtual bool has_nodata_read() const;
      virtual double nodata_read() const;
  };

  class DstMemoryImageResourceGDAL : public DstMemoryImageResource {
      struct Data;
      boost::shared_ptr<Data> m_data;

    public:
      DstMemoryImageResourceGDAL(const ImageFormat& fmt);

      virtual void write( ImageBuffer const& buf, BBox2i const& bbox );
      virtual void flush() {}

      virtual bool has_block_write()  const {return false;}
      virtual bool has_nodata_write() const {return true;}

      virtual void set_nodata_write(double value);

      virtual const uint8* data() const;
      virtual size_t size() const;
  };

}

#endif
