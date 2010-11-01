// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#define __INSIDE_VW_IMAGE_IMAGERESOURCEIMPL_H__ 1
#include <vw/Image/ImageResourceOpenCV.h>
#undef  __INSIDE_VW_IMAGE_IMAGERESOURCEIMPL_H__

#include <vw/Core/Debugging.h>
#include <opencv/cxcore.h>

namespace {
  cv::Rect makeRect(const vw::BBox2i& bbox) {
    return cv::Rect(bbox.min().x(), bbox.min().y(), bbox.width(), bbox.height());
  }
}

namespace vw {

bool ImageResourceOpenCV::contiguous_roi(const BBox2i& bbox) const {
  // if the opencv data is contiguous, do it the easy way. It's contiguous if:
  //  1) The data is from a single row (opencv guarantees this, only rows can have gaps), or
  //  2) Rows are contiguous (isContinuous) and the roi spans all columns
  if (bbox.height() == 1)
    return true;

  return m_matrix->isContinuous()
         && bbox.min().x() == 0
         && uint32(bbox.width()) == m_format.cols;
}

ImageFormat ImageResourceOpenCV::identify() const {
  ImageFormat fmt;
  fmt.cols = m_matrix->cols;
  fmt.rows = m_matrix->rows;
  fmt.planes = 1;

  switch(m_matrix->channels()) {
    case 1:  fmt.pixel_format = VW_PIXEL_SCALAR;            break;
    case 2:  fmt.pixel_format = VW_PIXEL_GRAYA;             break;
    case 3:  fmt.pixel_format = VW_PIXEL_RGB;               break;
    case 4:  fmt.pixel_format = VW_PIXEL_RGBA;              break;
    case 5:  fmt.pixel_format = VW_PIXEL_GENERIC_5_CHANNEL; break;
    case 6:  fmt.pixel_format = VW_PIXEL_GENERIC_6_CHANNEL; break;
    default:
      vw_throw(ArgumentErr() << "ImageResourceOpenCV: Too many channels ("
                             << m_matrix->channels() << ") in input matrix");
  }

  switch(m_matrix->depth()) {
    case CV_8U:  fmt.channel_type = VW_CHANNEL_UINT8;   break;
    case CV_8S:  fmt.channel_type = VW_CHANNEL_INT8;    break;
    case CV_16U: fmt.channel_type = VW_CHANNEL_UINT16;  break;
    case CV_16S: fmt.channel_type = VW_CHANNEL_INT16;   break;
    case CV_32S: fmt.channel_type = VW_CHANNEL_INT32;   break;
    case CV_32F: fmt.channel_type = VW_CHANNEL_FLOAT32; break;
    case CV_64F: fmt.channel_type = VW_CHANNEL_FLOAT64; break;
    default:
      vw_throw(ArgumentErr() << "ImageResourceOpenCV: Unknown channel type "
                             << m_matrix->depth() << " in input matrix");
  }
  return fmt;
}

ImageResourceOpenCV::ImageResourceOpenCV(boost::shared_ptr<cv::Mat> matrix)
  : m_matrix(matrix)
{
  VW_ASSERT(!m_matrix->empty(), ArgumentErr() << VW_CURRENT_FUNCTION << ": Matrix must be allocated already.");
  m_format = identify();
}

void ImageResourceOpenCV::read( ImageBuffer const& dst_buf, BBox2i const& bbox ) const {
  VW_ASSERT(dst_buf.format.cols == uint32(bbox.width()) && dst_buf.format.rows == uint32(bbox.height()),
      LogicErr()    << VW_CURRENT_FUNCTION << ": Destination buffer has wrong dimensions!" );
  VW_ASSERT(bbox.min().x() >= 0 && bbox.min().y() >= 0 && uint32(bbox.max().x()) <= m_format.cols && uint32(bbox.max().y()) <= m_format.rows,
      ArgumentErr() << VW_CURRENT_FUNCTION << ": Bounding box must be inside matrix.");

  ImageFormat src_format = m_format;
  src_format.rows = bbox.height();
  src_format.cols = bbox.width();

  // make sure we don't accidentially delete our matrix
  boost::shared_ptr<cv::Mat> src_matrix = m_matrix;
  Vector2i origin = bbox.min();

  if (!contiguous_roi(bbox)) {
    src_matrix.reset(new cv::Mat());
    // copyTo handles resizing
    m_matrix->operator()(makeRect(bbox)).copyTo(*src_matrix);
    VW_ASSERT(src_matrix->isContinuous(), LogicErr() << "copied matrix is still not continuous!");
    origin = Vector2i(0,0);
  }

  const ImageBuffer src_buf(src_format, src_matrix->ptr(origin.y()));
  convert(dst_buf, src_buf, false);
}

void ImageResourceOpenCV::write( ImageBuffer const& src_buf, BBox2i const& bbox ) {
  VW_ASSERT(src_buf.format.cols == uint32(bbox.width()) && src_buf.format.rows == uint32(bbox.height()),
      LogicErr()    << VW_CURRENT_FUNCTION << ": Source buffer has wrong dimensions!" );
  VW_ASSERT(bbox.min().x() >= 0 && bbox.min().y() >= 0 && uint32(bbox.max().x()) <= m_format.cols && uint32(bbox.max().y()) <= m_format.rows,
      ArgumentErr() << VW_CURRENT_FUNCTION << ": Bounding box must be inside matrix.");

  ImageFormat dst_format = m_format;
  dst_format.rows = bbox.height();
  dst_format.cols = bbox.width();

  if (contiguous_roi(bbox)) {
    ImageBuffer dst_buf(dst_format, m_matrix->ptr(bbox.min().y()));
    convert(dst_buf, src_buf, false);
  } else {
    // Set up a buffer that's actually backed by an opencv matrix
    ImageBuffer tmp_buf(dst_format, 0);
    cv::Mat tmp_mat(bbox.width(), bbox.height(), m_matrix->type());
    VW_ASSERT(tmp_mat.isContinuous(), LogicErr() << "tmp matrix is still not continuous!");
    tmp_buf.data = tmp_mat.ptr();

    // perform conversion as usual
    convert(tmp_buf, src_buf, false);

    // copy the newly-converted data to the necessary roi
    cv::Mat submatrix(*m_matrix, makeRect(bbox));
    tmp_mat.copyTo(submatrix);
  }
}

} // namespace vw
