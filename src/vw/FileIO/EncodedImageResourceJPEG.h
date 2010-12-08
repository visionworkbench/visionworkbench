#ifndef __VW_FILEIO_ENCODEDIMAGERESOURCEJPEG_H__
#define __VW_FILEIO_ENCODEDIMAGERESOURCEJPEG_H__

#include <vw/FileIO/EncodedImageResource.h>
#include <boost/noncopyable.hpp>

namespace vw {

  class SrcEncodedImageResourceJPEG : public SrcEncodedImageResource, private boost::noncopyable {
      struct Data;
      boost::shared_ptr<Data> m_data;
      const uint8* m_encoded;
      size_t m_size;

    public:
      SrcEncodedImageResourceJPEG(const uint8* data, size_t len) { this->set_encoded_data(data, len); }

      virtual void set_encoded_data( const uint8* data, size_t len);
      virtual void read( ImageBuffer const& buf, BBox2i const& bbox ) const;

      virtual int32 cols() const;
      virtual int32 rows() const;
      virtual int32 planes() const;
      virtual PixelFormatEnum pixel_format() const;
      virtual ChannelTypeEnum channel_type() const;

      virtual bool has_block_read() const  {return false;}
      virtual bool has_nodata_read() const {return false;}
  };

  class DstEncodedImageResourceJPEG : public DstEncodedImageResource {
      struct Data;
      boost::shared_ptr<Data> m_data;
      const uint8* m_encoded;
      size_t m_size;

    public:

      DstEncodedImageResourceJPEG(VarArray<uint8>& data) { this->set_decoded_data(data); }

      virtual void set_decoded_data( VarArray<uint8>& data);
      virtual void write( ImageBuffer const& buf, BBox2i const& bbox );

      virtual bool has_block_write()  const {return false;}
      virtual bool has_nodata_write() const {return false;}

      virtual void flush();
  };

}

#endif
