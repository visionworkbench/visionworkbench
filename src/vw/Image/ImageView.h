// __BEGIN_LICENSE__
// 
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
// 
// Copyright 2006 Carnegie Mellon University. All rights reserved.
// 
// This software is distributed under the NASA Open Source Agreement
// (NOSA), version 1.3.  The NOSA has been approved by the Open Source
// Initiative.  See the file COPYING at the top of the distribution
// directory tree for the complete NOSA document.
// 
// THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
// KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
// LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
// SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
// A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
// THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
// DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
// 
// __END_LICENSE__

/// \file ImageView.h
/// 
/// Defines the core in-memory image view type.
///
#ifndef __VW_IMAGE_IMAGEVIEW_H__
#define __VW_IMAGE_IMAGEVIEW_H__

#include <string.h> // For memset()

#include <boost/smart_ptr.hpp>
#include <boost/type_traits.hpp>

#include <vw/Core/Exception.h>

#include <vw/Image/ImageViewBase.h>
#include <vw/Image/ImageResource.h>
#include <vw/Image/PixelAccessors.h>

namespace vw {

  /// The standard image container for in-memory image data.
  ///
  /// This class represents an image stored in memory or, more
  /// precisely, a view onto such an image.  That is, the ImageView
  /// object itself does not contain the image data itself but rather
  /// a pointer to it.  More than one ImageView object can point to
  /// the same data, and they can even choose to interpret that data
  /// differently.  In particular, copying an ImageView object is a
  /// shallow, lightweight operation.  The underlying image data is
  /// reference counted, so the user does not usually need to be
  /// concerned with memory allocation and deallocation.
  /// 
  /// A more complete name for this class might be something like
  /// MemoryImageView, or StandardImageView, but it is so ubiquitous
  /// that we decided to keep the name short.
  ///
  template <class PixelT>
  class ImageView : public ImageViewBase<ImageView<PixelT> >
  {
    boost::shared_array<PixelT> m_data;
    unsigned m_cols, m_rows, m_planes;
    PixelT *m_origin;
    ptrdiff_t m_cstride, m_rstride, m_pstride;

  public:
    /// The base type of the image.
    typedef ImageViewBase<ImageView<PixelT> > base_type;

    /// The pixel type of the image.
    typedef PixelT pixel_type;

    /// The data type returned when accessing the image.
    typedef PixelT& result_type;

    /// The image's %pixel_accessor type.
    typedef MemoryStridingPixelAccessor<PixelT> pixel_accessor;

    /// Constructs an empty image with zero size.
    ImageView()
      : m_cols(0), m_rows(0), m_planes(0), m_origin(0), m_cstride(0), 
        m_rstride(0), m_pstride(0) {}

    /// Copy-constructs a view pointing to the same data.
    /// Provided explicitly to clarify its precedence over 
    /// the templatized generalized copy constructor.
    ImageView( ImageView const& other )
      : m_data(other.m_data), m_cols(other.m_cols), 
        m_rows(other.m_rows), m_planes(other.m_planes),
        m_origin(other.m_origin), m_cstride(other.m_cstride), 
        m_rstride(other.m_rstride), m_pstride(other.m_pstride) {}

    /// Constructs an empty image with the given dimensions.
    ImageView( unsigned cols, unsigned rows, unsigned planes=1 )
      : m_cols(0), m_rows(0), m_planes(0), m_origin(0), m_cstride(0), 
        m_rstride(0), m_pstride(0) {
      set_size( cols, rows, planes );
    }

    /// Constructs an image view and rasterizes the given view into it.
    template <class ViewT>
    ImageView( ViewT const& view )
      : m_cols(0), m_rows(0), m_planes(0), m_origin(0), m_cstride(0),
        m_rstride(0), m_pstride(0) {
      set_size( view.cols(), view.rows(), view.planes() );
      view.rasterize( *this, BBox2i(0,0,view.cols(),view.rows()) );
    }

    /// Rasterizes the given view into the image, adjusting the size if needed.
    template <class SrcT>
    inline ImageView& operator=( ImageViewBase<SrcT> const& view ) {
      set_size( view.impl().cols(), view.impl().rows(), view.impl().planes() );
      view.impl().rasterize( *this, BBox2i(0,0,view.impl().cols(),view.impl().rows()) );
      return *this;
    }

    /// Rasterizes the given view into the image.
    template <class SrcT>
    inline ImageView const& operator=( ImageViewBase<SrcT> const& view ) const {
      view.impl().rasterize( *this, BBox2i(0,0,view.impl().cols(),view.impl().rows()) );
      return *this;
    }

    /// Returns the number of columns in the image.
    inline unsigned cols() const { return m_cols; }

    /// Returns the number of rows in the image.
    inline unsigned rows() const { return m_rows; }

    /// Returns the number of planes in the image.
    inline unsigned planes() const { return m_planes; }

    /// Returns a pixel_accessor pointing to the top-left corner of the first plane.
    inline pixel_accessor origin() const {
      return pixel_accessor( m_origin, m_cstride, m_rstride, m_pstride );
    }

    /// Returns the pixel at the given position in the given plane.
    inline result_type operator()( int col, int row, int plane=0 ) const {
      return *(m_origin + col*m_cstride + row*m_rstride + plane*m_pstride);
    }
  
    /// Adjusts the size of the image, allocating a new buffer if the size has changed.
    void set_size( unsigned cols, unsigned rows, unsigned planes = 1 ) {
      if( cols==m_cols && rows==m_rows && planes==m_planes ) return;
        
      unsigned size = cols*rows*planes;
      if( size==0 ) {
        m_data.reset();
      }
      else {
        boost::shared_array<PixelT> data( new PixelT[size] );
        m_data = data;
      }

      m_cols = cols;
      m_rows = rows;
      m_planes = planes;
      m_origin = m_data.get();
      m_cstride = 1;
      m_rstride = cols;
      m_pstride = rows*cols;

      // Fundamental types might not be initialized.  Really this is
      // true of all POD types, but there's no good way to detect
      // those.  We can only hope that the user will never use a
      // custom POD pixel type.
      // 
      // Note that this is a copy of the fill algorithm that resides
      // in ImageAlgorithms.h, however including ImageAlgorithms.h
      // directly causes an include file cycle.
      if( boost::is_fundamental<pixel_type>::value ) {
        memset( m_data.get(), 0, m_rows*m_cols*m_planes*sizeof(PixelT) );
      }
    }

    /// Adjusts the size of the image to match the dimensions of another image.
    template <class ImageT>
    void set_size( ImageViewBase<ImageT> &img ) {
      this->set_size(img.impl().cols(), img.impl().rows(), img.impl().planes());
    }

    /// Resets to an empty image with zero size.
    void reset() {
      m_data.reset();
      m_cols = m_rows = m_planes = 0;
      m_origin = 0;
      m_cstride = m_rstride = m_pstride = 0;
    }

    /// Returns a pointer to the origin of the image in memory.
    pixel_type *data() const {
      return m_origin;
    }

    /// A safe bool conversion intermediate type.
    typedef typename boost::shared_array<PixelT>::unspecified_bool_type unspecified_bool_type;
    /// Evaluates to true in a bool context if this ImageView points
    /// to a valid block of memory.  (It is false if e.g. the object is
    /// default-constructed, or if reset() is called, or if set_size()
    /// is called with zero dimensions.)
    operator unspecified_bool_type() const {
      return m_data;
    }

    /// Returns true if no other ImageView object is sharing 
    /// this block of memory.
    bool unique() const {
      return (!m_data) || m_data.unique();
    }

    /// Returns an ImageBuffer describing the image data.
    ImageBuffer buffer() const {
      ImageBuffer buffer;
      buffer.data = data();
      buffer.format = base_type::format();
      buffer.cstride = sizeof(PixelT);
      buffer.rstride = sizeof(PixelT)*cols();
      buffer.pstride = sizeof(PixelT)*cols()*rows();
      buffer.unpremultiplied = false;
      return buffer;
    }

    /// The return type of prerasterize().
    typedef ImageView prerasterize_type;

    /// Prepare an ImageView to be rasterized.  Simply returns the 
    /// original image view.
    inline prerasterize_type prerasterize( BBox2i bbox ) const { return *this; }

    /// Rasterize the image view.  Simply invokes the default 
    /// rasterization function.
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i bbox ) const {
      vw::rasterize( prerasterize(bbox), dest, bbox );
    }
  };

  // Image view traits
  /// Specifies that ImageView objects are resizable.
  template <class PixelT>
  struct IsResizable<ImageView<PixelT> > : public true_type {};

  /// Specifies that ImageView objects are fast to access.
  template <class PixelT>
  struct IsMultiplyAccessible<ImageView<PixelT> > : public true_type {};


  // *******************************************************************
  // Image view reading and writing functions.
  // *******************************************************************

  template <class PixelT>
  inline void read_image( ImageView<PixelT>& dst, ImageResource const& src, BBox2i const& bbox ) {
    unsigned planes = 1;
    if( ! IsCompound<PixelT>::value ) {
      // The image has a fundamental pixel type
      if( src.planes()>1 && src.channels()>1 )
        vw_throw( ArgumentErr() << "Cannot read a multi-plane multi-channel image resource into a single-channel view." );
      planes = (std::max)( src.planes(), src.channels() );
    }
    dst.set_size( bbox.width(), bbox.height(), planes );
    src.read( dst.buffer(), bbox );
  }

  template <class PixelT>
  inline void read_image( ImageView<PixelT>& dst, ImageResource const& src ) {
    read_image( dst, src, BBox2i(0,0,src.cols(),src.rows()) );
  }

  template <class PixelT>
  inline void read_image( ImageView<PixelT> const& dst, ImageResource const& src, BBox2i const& bbox ) {
    src.read( dst.buffer(), bbox );
  }

  template <class PixelT>
  inline void read_image( ImageView<PixelT> const& dst, ImageResource const& src ) {
    read_image( dst, src, BBox2i(0,0,src.cols(),src.rows()) );
  }

  template <class ImageT>
  inline void read_image( ImageViewBase<ImageT> const& dst, ImageResource const& src, BBox2i const& bbox ) {
    ImageView<typename ImageT::pixel_type> intermediate;
    read_image( intermediate, src, bbox );
    dst = intermediate;
  }

  template <class ImageT>
  inline void read_image( ImageViewBase<ImageT> const& dst, ImageResource const& src ) {
    read_image( dst, src, BBox2i(0,0,src.cols(),src.rows()) );
  }

  template <class PixelT>
  inline void write_image( ImageResource &dst, ImageView<PixelT> const& src, BBox2i const& bbox ) {
    dst.write( src.buffer(), bbox );
  }

  template <class PixelT>
  inline void write_image( ImageResource &dst, ImageView<PixelT> const& src ) {
    write_image( dst, src, BBox2i(0,0,dst.cols(),dst.rows()) );
  }

  template <class ImageT>
  inline void write_image( ImageResource &dst, ImageViewBase<ImageT> const& src, BBox2i const& bbox ) {
    ImageView<typename ImageT::pixel_type> intermediate = src;
    write_image( dst, intermediate, bbox );
  }

  template <class ImageT>
  inline void write_image( ImageResource &dst, ImageViewBase<ImageT> const& src ) {
    ImageView<typename ImageT::pixel_type> intermediate = src;
    write_image( dst, intermediate );
  }

} // namespace vw

#endif // __VW_IMAGE_IMAGEVIEW_H__
