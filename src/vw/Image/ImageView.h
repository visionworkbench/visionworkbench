// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file ImageView.h
///
/// Defines the core in-memory image view type.
///
#ifndef __VW_IMAGE_IMAGEVIEW_H__
#define __VW_IMAGE_IMAGEVIEW_H__

#include <cstring> // For memset()

#include <boost/smart_ptr.hpp>
#include <boost/type_traits.hpp>

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
    // Protect user from themselves. ImageView can never contain
    // another ImageView.
    BOOST_STATIC_ASSERT( !IsImageView<PixelT>::value );

    boost::shared_array<PixelT> m_data;
    int32 m_cols, m_rows, m_planes;
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
      : ImageViewBase<ImageView<PixelT> >(other),
        m_data(other.m_data), m_cols(other.m_cols),
        m_rows(other.m_rows), m_planes(other.m_planes),
        m_origin(other.m_origin), m_cstride(other.m_cstride),
        m_rstride(other.m_rstride), m_pstride(other.m_pstride) {}

    /// Constructs an empty image with the given dimensions.
    ImageView( int32 cols, int32 rows, int32 planes=1 )
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
    inline int32 cols() const { return m_cols; }

    /// Returns the number of rows in the image.
    inline int32 rows() const { return m_rows; }

    /// Returns the number of planes in the image.
    inline int32 planes() const { return m_planes; }

    /// Returns a pixel_accessor pointing to the top-left corner of the first plane.
    inline pixel_accessor origin() const {
#if defined(VW_ENABLE_BOUNDS_CHECK) && (VW_ENABLE_BOUNDS_CHECK==1)
      return pixel_accessor( m_origin, m_cstride, m_rstride, m_pstride,
                             cols(), rows(), planes() );
#else
      return pixel_accessor( m_origin, m_cstride, m_rstride, m_pstride );
#endif
    }

    /// Returns the pixel at the given position in the given plane.
    inline result_type operator()( int32 col, int32 row, int32 plane=0 ) const {
#if defined(VW_ENABLE_BOUNDS_CHECK) && (VW_ENABLE_BOUNDS_CHECK==1)
      if (col < 0 || col >= cols() || row < 0 || row >= rows() || plane < 0 || plane >= planes())
        vw_throw(ArgumentErr() << "ImageView::operator() - invalid index [" << col << " " << row << " " << plane << "] for ImageView with dimensions [" << cols() << " " << rows() << " " << planes() << "]");
#endif
      return *(m_origin + col*m_cstride + row*m_rstride + plane*m_pstride);
    }

    /// Adjusts the size of the image, allocating a new buffer if the size has changed.
    void set_size( int32 cols, int32 rows, int32 planes = 1 ) {
      if( cols==m_cols && rows==m_rows && planes==m_planes ) return;

      // FIXME: We can do better than this on 64-bit machines.
      int32 size = cols*rows*planes;
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
    void set_size( const ImageViewBase<ImageT> &img ) {
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
    inline prerasterize_type prerasterize( BBox2i /*bbox*/ ) const { return *this; }

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

} // namespace vw

#endif // __VW_IMAGE_IMAGEVIEW_H__
