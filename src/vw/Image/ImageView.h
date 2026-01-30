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
#include <vw/Image/PixelTypes.h>
#include <vw/Image/PixelMask.h>
#include <vw/Math/Vector.h>

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
  /// - WARNING: Never refer to these objects by reference!  The
  ///            behaviour is undefined.
  /// - The user manual discourages users from having multiple PLANES.
  ///   This is distinct from multiple CHANNELS which is achieved by
  ///   using a pixel type container object with multiple values, such
  ///   as Vector3.  The manual does not address why this is but it is
  ///   safe to say that using planes is not well supported.
  /// - Because of the above, image data is stored internally in 
  ///   INTERLEAVED (BIP) format.
  template <class PixelT>
  class ImageView : public ImageViewBase<ImageView<PixelT>> {
    // Protect user from themselves. ImageView can never contain
    // another ImageView.
    BOOST_STATIC_ASSERT(!IsImageView<PixelT>::value);

    boost::shared_array<PixelT> m_data; ///< Underlying data pointer
    int32 m_cols, m_rows, m_planes;     ///< Data dimensions
    PixelT *m_origin;                   ///< Generally points to m_data.get()
    ssize_t m_rstride, m_pstride;       ///< Row stride and plane stride in PixelT counts.

  public:
    /// The base type of the image.
    typedef ImageViewBase<ImageView<PixelT>> base_type;

    /// The pixel type of the image.
    typedef PixelT pixel_type;

    /// The data type returned when accessing the image.
    typedef PixelT& result_type;

    /// The image's %pixel_accessor type.
    typedef MemoryStridingPixelAccessor<PixelT> pixel_accessor;

    /// Constructs an empty image with zero size.
    ImageView()
      : m_cols(0), m_rows(0), m_planes(0), m_origin(0),
        m_rstride(0), m_pstride(0) {}

    /// Copy-constructs a view pointing to the same data.
    /// Provided explicitly to clarify its precedence over
    /// the templatized generalized copy constructor.
    ImageView(ImageView const& other)
      : ImageViewBase<ImageView<PixelT>>(other),
        m_data(other.m_data), m_cols(other.m_cols),
        m_rows(other.m_rows), m_planes(other.m_planes),
        m_origin(other.m_origin),
        m_rstride(other.m_rstride), m_pstride(other.m_pstride) {}

    /// Constructs an empty image with the given dimensions.
    ImageView(int32 cols, int32 rows, int32 planes=1)
      : m_cols(0), m_rows(0), m_planes(0), m_origin(0),
        m_rstride(0), m_pstride(0) {
      set_size(cols, rows, planes);
    }

    /// Constructs an image view and rasterizes the given view into it.
    template <class ViewT>
    ImageView(ViewT const& view)
      : m_cols(0), m_rows(0), m_planes(0), m_origin(0),
        m_rstride(0), m_pstride(0) {
      set_size(view.cols(), view.rows(), view.planes());
      view.rasterize(*this, BBox2i(0,0,view.cols(),view.rows()));
    }

    /// Note that this is almost a copy of read_image in ImageIO, but actually
    /// including that is a circular dependency.
    explicit ImageView(const SrcImageResource& src) {
      int32 planes = 1;
      if(! IsCompound<PixelT>::value) {
        // The image has a fundamental pixel type
        if(src.planes()>1 && src.channels()>1)
          vw_throw(ArgumentErr() << "Cannot read a multi-plane multi-channel image resource into a single-channel view.");
        planes = (std::max)(src.planes(), src.channels());
      }
      set_size(src.cols(), src.rows(), planes);
      src.read(this->buffer(), BBox2i(0,0,src.cols(),src.rows()));
    }

    /// Rasterizes the given view into the image, adjusting the size if needed.
    template <class SrcT>
    inline ImageView& operator=(ImageViewBase<SrcT> const& view) {
      set_size(view.impl().cols(), view.impl().rows(), view.impl().planes());
      view.impl().rasterize(*this, BBox2i(0,0,view.impl().cols(),view.impl().rows()));
      return *this;
    }

    /// Rasterizes the given view into the image.
    template <class SrcT>
    inline ImageView const& operator=(ImageViewBase<SrcT> const& view) const {
      view.impl().rasterize(*this, BBox2i(0,0,view.impl().cols(),view.impl().rows()));
      return *this;
    }

    inline int32 cols  () const { return m_cols;   } ///< Returns the number of columns in the image.
    inline int32 rows  () const { return m_rows;   } ///< Returns the number of rows in the image.
    inline int32 planes() const { return m_planes; } ///< Returns the number of planes in the image.
    
    /// Returns a pixel_accessor pointing to the top-left corner of the first plane.
    inline pixel_accessor origin() const {
#if defined(VW_ENABLE_BOUNDS_CHECK) && (VW_ENABLE_BOUNDS_CHECK==1)
      return pixel_accessor(m_origin, m_rstride, m_pstride,
                             cols(), rows(), planes());
#else
      return pixel_accessor(m_origin, m_rstride, m_pstride);
#endif
    }

    /// Returns the pixel at the given position in the given plane.
    /// - See the warning about using planes at the top of the class.
    inline result_type operator()(int32 col, int32 row, int32 plane=0) const {
#if defined(VW_ENABLE_BOUNDS_CHECK) && (VW_ENABLE_BOUNDS_CHECK==1)
      if (col < 0 || col >= cols() || row < 0 || row >= rows() || plane < 0 || plane >= planes())
        vw_throw(ArgumentErr() << "ImageView::operator() - invalid index [" << col << " " << row << " " << plane << "] for ImageView with dimensions [" << cols() << " " << rows() << " " << planes() << "]");
#endif
      return *(m_origin + col + row*m_rstride + plane*m_pstride);
    }

    /// Adjusts the size of the image, allocating a new buffer if the size has changed.
    void set_size(int32 cols, int32 rows, int32 planes = 1) {
      // These sizes are pretty large for in-memory images and should only come up
      //  in the case of bugs in the code.
      static const int32  MAX_PIXEL_SIZE   = 80000;
      static const int32  MAX_PLANE_COUNT  = 1024; // Really should never be using these anyways
      static const uint64 MAX_TOTAL_PIXELS = 6400000000;
      
      // Check if we already have the correct size
      if(cols==m_cols && rows==m_rows && planes==m_planes)
          return;

      VW_ASSERT(cols >= 0 && rows >= 0 && planes >= 0, // No negative sizes!
                ArgumentErr() << "Cannot allocate image with negative pixel count (you requested " 
                              << cols << " x "  << rows << " x " << planes << ")");

      // Make sure the image is not too big.
      VW_ASSERT(cols < MAX_PIXEL_SIZE && rows < MAX_PIXEL_SIZE,
          ArgumentErr() << "Refusing to allocate an image larger than " << MAX_PIXEL_SIZE-1 
                        << " pixels on a side (you requested " << cols << " x " << rows << ")");
      VW_ASSERT(planes < MAX_PLANE_COUNT,
          ArgumentErr() << "Refusing to allocate an image with more than " << MAX_PLANE_COUNT-1 
                        << " planes (you requested " << planes << ")");

      uint64 size64 = uint64(cols) * uint64(rows) * uint64(planes);

      VW_ASSERT(size64 < MAX_TOTAL_PIXELS,
          ArgumentErr() << "Refusing to allocate an image with more than " << MAX_TOTAL_PIXELS 
                        << " pixels (you requested " << size64 << ")");
      
      size_t size = size64;

      if(size==0)
        m_data.reset();
      else {
        boost::shared_array<PixelT> data(new (std::nothrow) PixelT[size]);
        if (!data) {
          // print it and throw it for the benefit of OSX, which
          // doesn't print the exception what() on terminate()
          VW_OUT(ErrorMessage) << "Cannot allocate enough memory for a " 
                               << cols << "x" << rows << "x" << planes
                               << " image: too many bytes!" << std::endl;
          vw_throw(ArgumentErr() << "Cannot allocate enough memory for a " 
                   << cols << "x" << rows << "x" << planes << " image: too many bytes!");
        }
        m_data = data;
      }

      m_cols    = cols;
      m_rows    = rows;
      m_planes  = planes;
      m_origin  = m_data.get();
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
      if(boost::is_fundamental<pixel_type>::value) {
        memset(m_data.get(), 0, m_rows*m_cols*m_planes*sizeof(PixelT));
      }
    }

    /// Adjusts the size of the image to match the dimensions of another image.
    template <class ImageT>
    void set_size(const ImageViewBase<ImageT> &img) {
      this->set_size(img.impl().cols(), img.impl().rows(), img.impl().planes());
    }

    /// Resets to an empty image with zero size.
    void reset() {
      m_data.reset();
      m_cols    = m_rows = m_planes = 0;
      m_origin  = 0;
      m_rstride = m_pstride = 0;
    }

    /// Returns a pointer to the origin of the image in memory.
    pixel_type *data() const {
      return m_origin;
    }

    bool is_valid_image() const {
      return !(!m_data);
    }

    /// Returns true if no other ImageView object is sharing this block of memory.
    bool unique() const {
      return (!m_data) || m_data.unique();
    }

    /// Returns an ImageBuffer describing the image data.
    ImageBuffer buffer() const {
      ImageBuffer buffer;
      buffer.data    = data();
      buffer.format  = base_type::format();
      buffer.cstride = sizeof(PixelT);
      buffer.rstride = sizeof(PixelT)*cols();
      buffer.pstride = sizeof(PixelT)*cols()*rows();
      return buffer;
    }

    /// The return type of prerasterize().
    typedef ImageView prerasterize_type;

    /// Prepare an ImageView to be rasterized.  Simply returns the
    /// original image view.
    inline prerasterize_type prerasterize(BBox2i /*bbox*/) const { return *this; }

    /// Rasterize the image view.  Simply invokes the default
    /// rasterization function.
    template <class DestT> inline void rasterize(DestT const& dest, BBox2i bbox) const {
      vw::rasterize(prerasterize(bbox), dest, bbox);
    }
  };

  // Image view traits
  /// Specifies that ImageView objects are resizable.
  template <class PixelT>
  struct IsResizable<ImageView<PixelT>> : public true_type {};

  /// Specifies that ImageView objects are fast to access.
  template <class PixelT>
  struct IsMultiplyAccessible<ImageView<PixelT>> : public true_type {};

  // Explicit template instantiation declarations for common types.
  // These are instantiated once in ImageView.cc to reduce
  // compilation time and object file size.
  extern template class ImageView<float>;
  extern template class ImageView<double>;
  extern template class ImageView<uint8>;
  extern template class ImageView<int16>;
  extern template class ImageView<uint16>;
  extern template class ImageView<int32>;
  extern template class ImageView<uint32>;
  extern template class ImageView<PixelGray<float>>;
  extern template class ImageView<PixelGray<double>>;
  extern template class ImageView<PixelGray<uint8>>;
  extern template class ImageView<PixelRGB<uint8>>;
  extern template class ImageView<PixelRGB<float>>;
  extern template class ImageView<PixelRGB<double>>;
  extern template class ImageView<PixelRGBA<uint8>>;
  extern template class ImageView<PixelMask<uint8>>;
  extern template class ImageView<PixelMask<float>>;
  extern template class ImageView<PixelMask<double>>;
  extern template class ImageView<PixelMask<PixelGray<uint8>>>;
  extern template class ImageView<PixelMask<PixelGray<float>>>;
  extern template class ImageView<PixelMask<Vector2f>>;
  extern template class ImageView<PixelMask<Vector2i>>;
  extern template class ImageView<Vector2>;
  extern template class ImageView<Vector3>;

} // namespace vw

#endif // __VW_IMAGE_IMAGEVIEW_H__
