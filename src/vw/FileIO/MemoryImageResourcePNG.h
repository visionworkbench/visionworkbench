#ifndef __VW_FILEIO_MEMORYIMAGERESOURCEPNG_H__
#define __VW_FILEIO_MEMORYIMAGERESOURCEPNG_H__

#include <vw/FileIO/MemoryImageResource.h>
#include <boost/noncopyable.hpp>

namespace vw {

  class SrcMemoryImageResourcePNG : public SrcMemoryImageResource, private boost::noncopyable {
      struct Data;
      mutable boost::shared_ptr<Data> m_data;

    public:
      SrcMemoryImageResourcePNG(const uint8* buffer, size_t len);

      virtual void read( ImageBuffer const& buf, BBox2i const& bbox ) const;

      virtual int32 cols() const;
      virtual int32 rows() const;
      virtual int32 planes() const;
      virtual PixelFormatEnum pixel_format() const;
      virtual ChannelTypeEnum channel_type() const;

      virtual bool has_block_read() const  {return false;}
      virtual bool has_nodata_read() const {return false;}
  };

  class DstMemoryImageResourcePNG : public DstMemoryImageResource {
      struct Data;
      boost::shared_ptr<Data> m_data;

    public:
      DstMemoryImageResourcePNG(std::vector<uint8>* buffer, const ImageFormat& fmt);

      virtual void write( ImageBuffer const& buf, BBox2i const& bbox );
      virtual void flush() {}

      virtual bool has_block_write()  const {return false;}
      virtual bool has_nodata_write() const {return false;}
  };

}

#endif
