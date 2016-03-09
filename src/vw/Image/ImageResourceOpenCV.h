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

/// ImageResource to read/write OpenCV in-memory images.
class ImageResourceOpenCV : public ImageResource {
  private:
    ImageFormat                m_format;
    boost::shared_ptr<cv::Mat> m_matrix;

  protected:
    /// is the given bbox roi stored contiguously in memory?
    bool contiguous_roi(const BBox2i& bbox) const;
    /// identify the current matrix type
    ImageFormat identify() const;

  public:

    ImageResourceOpenCV(boost::shared_ptr<cv::Mat> matrix);

    virtual ~ImageResourceOpenCV() {};

    virtual ImageFormat format() const {
      return m_format;
    }

    // Read the image resource at the given location into the given buffer.
    virtual void read( ImageBuffer const& dst_buf, BBox2i const& bbox ) const;

    // Write the given buffer to the image resource at the given location.
    virtual void write( ImageBuffer const& buf, BBox2i const& bbox );

    virtual void flush() {}

    virtual bool has_block_write () const {return false;}
    virtual bool has_nodata_write() const {return false;}
    virtual bool has_block_read  () const {return false;}
    virtual bool has_nodata_read () const {return false;}
};

}

#endif
