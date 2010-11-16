// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// ImageResource to read/write OpenCV in-memory images. Do not include this
// header directly.

#ifndef __VW_IMAGE_IMAGERESOURCEOPENCV_H__
#define __VW_IMAGE_IMAGERESOURCEOPENCV_H__

#ifndef __INSIDE_VW_IMAGE_IMAGERESOURCEIMPL_H__
# error Do not include this file directly. Instead, include <vw/Image/ImageResourceImpl.h>
#endif

#include <vw/Image/ImageResource.h>

namespace cv {class Mat;}

namespace vw {

class ImageResourceOpenCV : public ImageResource {
  private:
    ImageFormat m_format;
    boost::shared_ptr<cv::Mat> m_matrix;

  protected:
    // is the given bbox roi stored contiguously in memory?
    bool contiguous_roi(const BBox2i& bbox) const;
    // identify the current matrix type
    ImageFormat identify() const;

  public:

    ImageResourceOpenCV(boost::shared_ptr<cv::Mat> matrix);

    virtual ~ImageResourceOpenCV() {};

    virtual int32 cols() const {
      return m_format.cols;
    }

    virtual int32 rows() const {
      return m_format.rows;
    }

    virtual int32 planes() const {
      return m_format.planes;
    }

    virtual PixelFormatEnum pixel_format() const {
      return m_format.pixel_format;
    }

    virtual ChannelTypeEnum channel_type() const {
      return m_format.channel_type;
    }

    // Read the image resource at the given location into the given buffer.
    virtual void read( ImageBuffer const& dst_buf, BBox2i const& bbox ) const;

    // Write the given buffer to the image resource at the given location.
    virtual void write( ImageBuffer const& buf, BBox2i const& bbox );

    virtual void flush() {}

    virtual bool has_block_write() const  {return false;}
    virtual bool has_nodata_write() const {return false;}
    virtual bool has_block_read() const   {return false;}
    virtual bool has_nodata_read() const  {return false;}
};

}

#endif
