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

/// \file ImageViewRef.h
///
/// A generic image view reference class.
///
/// Most Vision Workbench image processing functions simply return image view
/// objects that lazily represent the processed data. Under some circumstances
/// it is helpful to be able to hold onto such a processed view without
/// rasterizing it.  Ordinarily this requires knowing the full type of the view.
/// When this is not acceptable, the \ref vw::ImageViewRef class allows you to
/// hold a virtualized reference to an arbitrary image view with a given pixel
/// type.
///
/// - WARNING: Never refer to these objects by reference!  The
//             behaviour is undefined.
#ifndef __VW_IMAGE_IMAGEVIEWREF_H__
#define __VW_IMAGE_IMAGEVIEWREF_H__

#include <boost/type_traits.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>

#include <vw/Image/ImageView.h>
#include <vw/Image/Manipulation.h>
#include <vw/Image/SparseImageCheck.h>

namespace vw {

  /// \cond INTERNAL
  template <class PixelT>
  class ImageViewRefAccessorBase {
  public:
    virtual ~ImageViewRefAccessorBase() {}
    virtual ImageViewRefAccessorBase* copy() const = 0;
    virtual void next_col  () = 0;
    virtual void prev_col  () = 0;
    virtual void next_row  () = 0;
    virtual void prev_row  () = 0;
    virtual void next_plane() = 0;
    virtual void prev_plane() = 0;
    virtual void advance(ssize_t di, ssize_t dj, ssize_t dp=0) = 0;
    virtual PixelT operator*() const = 0;
  };

  template <class IterT>
  class ImageViewRefAccessorImpl: public ImageViewRefAccessorBase<typename IterT::pixel_type> {
  private:
    IterT m_iter;
    // We know what this is, but the base class doesn't. So just cast.
    typedef typename IterT::offset_type iter_offset_type;
  public:
    typedef typename IterT::pixel_type pixel_type;

    ImageViewRefAccessorImpl(IterT const& iter): m_iter(iter) {}
    virtual ~ImageViewRefAccessorImpl() {}

    virtual ImageViewRefAccessorBase<pixel_type>* copy() const { return new ImageViewRefAccessorImpl(m_iter); }

    virtual void next_col  () { m_iter.next_col();   }
    virtual void prev_col  () { m_iter.prev_col();   }
    virtual void next_row  () { m_iter.next_row();   }
    virtual void prev_row  () { m_iter.prev_row();   }
    virtual void next_plane() { m_iter.next_plane(); }
    virtual void prev_plane() { m_iter.prev_plane(); }
    virtual void advance(ssize_t di, ssize_t dj, ssize_t dp=0) {
      m_iter.advance((iter_offset_type)di,(iter_offset_type)dj,dp);
    }
    virtual pixel_type operator*() const { return *m_iter; }
  };
  /// \endcond

  /// A special virtualized accessor adaptor.
  ///
  /// This accessor adaptor is used by the \ref vw::ImageViewRef class.
  template <class PixelT>
  class ImageViewRefAccessor {
  private:
    boost::scoped_ptr<ImageViewRefAccessorBase<PixelT>> m_iter;
  public:
    typedef PixelT  pixel_type;
    typedef PixelT  result_type;
    typedef ssize_t offset_type;

    template <class IterT> ImageViewRefAccessor(IterT const& iter): m_iter(new ImageViewRefAccessorImpl<IterT>(iter)) {}
    ~ImageViewRefAccessor() {}

    ImageViewRefAccessor(ImageViewRefAccessor const& other): m_iter(other.m_iter->copy()) {}
    ImageViewRefAccessor& operator=(ImageViewRefAccessor const& other) { 
      m_iter = other.m_iter->copy(); return *this; 
    }

    inline ImageViewRefAccessor& next_col  () { m_iter->next_col();   return *this; }
    inline ImageViewRefAccessor& prev_col  () { m_iter->prev_col();   return *this; }
    inline ImageViewRefAccessor& next_row  () { m_iter->next_row();   return *this; }
    inline ImageViewRefAccessor& prev_row  () { m_iter->prev_row();   return *this; }
    inline ImageViewRefAccessor& next_plane() { m_iter->next_plane(); return *this; }
    inline ImageViewRefAccessor& prev_plane() { m_iter->prev_plane(); return *this; }
    inline ImageViewRefAccessor& advance(ssize_t di, ssize_t dj, ssize_t dp=0) { 
      m_iter->advance(di,dj,dp=0); return *this; 
    }
    inline pixel_type operator*() const { return *(*m_iter); }
  };


  /// \cond INTERNAL
  // Base class definition
  template <class PixelT>
  class ImageViewRefBase {
  public:
    typedef PixelT pixel_type;
    typedef ImageViewRefAccessor<PixelT> pixel_accessor;

    virtual ~ImageViewRefBase() {}

    virtual int32 cols  () const = 0;
    virtual int32 rows  () const = 0;
    virtual int32 planes() const = 0;
    virtual pixel_type operator()(int32 i,  int32 j,  int32 p) const = 0;
    virtual pixel_type operator()(double i, double j, int32 p) const = 0;
    virtual pixel_accessor origin() const = 0;

    virtual bool sparse_check(BBox2i const& bbox) const = 0;
    virtual void rasterize(ImageView<pixel_type> const& dest, BBox2i const& bbox) const = 0;
  };

  // ImageViewRef class implementation
  template <class ViewT>
  class ImageViewRefImpl: public ImageViewRefBase<typename ViewT::pixel_type> {
  private:
    ViewT m_view;
  public:
    typedef typename ViewT::pixel_type pixel_type;
    typedef ImageViewRefAccessor<typename ViewT::pixel_type> pixel_accessor;

    ImageViewRefImpl(ImageViewBase<ViewT> const& view): m_view(view.impl()) {}
    virtual ~ImageViewRefImpl() {}

    virtual int32          cols  () const { return m_view.cols();   }
    virtual int32          rows  () const { return m_view.rows();   }
    virtual int32          planes() const { return m_view.planes(); }
    virtual pixel_accessor origin() const { return m_view.origin(); }
    
    virtual pixel_type     operator()(int32  i, int32  j, int32 p) const { return m_view(i,j,p); }
    virtual pixel_type     operator()(double i, double j, int32 p) const { return m_view(i,j,p); }

    virtual bool sparse_check(BBox2i const& bbox) const { return vw::sparse_check(m_view, bbox); }
    virtual void rasterize(ImageView<pixel_type> const& dest, BBox2i const& bbox) const { m_view.rasterize(dest, bbox); }

    ViewT const& child() const { return m_view; }
  };
  /// \endcond


  /// A virtualized image view reference object.
  ///
  /// This class behaves as a reference to an arbitrary image view
  /// with the given pixel type.  The purpose of this class is to
  /// hide the full type of a view behind a veil of abstraction,
  /// making things like run-time polymorphic behavior possible.
  /// The inevitable cost of this flexibility is one virtual
  /// function call per method invocation.  In many cases there
  /// are additional costs associated with not being able to
  /// perform template-based optimizations at compile time.
  ///
  /// Like any C++ reference, you bind an ImageViewRef to a view
  /// using a constructor and future operations act on the bound
  /// object instead of the reference object itself.
  ///
  /// The current implementation of ImageViewRef is read-only.
  template <class PixelT>
  class ImageViewRef: public ImageViewBase<ImageViewRef<PixelT>> {
  private:
    boost::shared_ptr<ImageViewRefBase<PixelT>> m_view;
  public:
    typedef PixelT pixel_type;
    typedef PixelT result_type;
    typedef ImageViewRefAccessor<PixelT> pixel_accessor;

    // Default construction conjures up an empty memory image as a placeholder.
    // This allows an ImageViewRef to be created without any arguments, which
    // makes it suitable for use in situations where creation and assignment
    // must happen as separate steps, such as in STL containers.
    ImageViewRef(): m_view(new ImageViewRefImpl<ImageView<PixelT>>(ImageView<PixelT>())) {}

    // Assignment constructor creates an ImageViewRef from another ImageView.
    template <class ViewT> ImageViewRef(ImageViewBase<ViewT> const& view):
     m_view(new ImageViewRefImpl<ViewT>(view)) {}
    ~ImageViewRef() {}

    template <class ViewT> void reset(ImageViewBase<ViewT> const& view) {
      m_view.reset(new ImageViewRefImpl<ViewT>(view)); 
    }

    inline int32 cols  () const { return m_view->cols();   }
    inline int32 rows  () const { return m_view->rows();   }
    inline int32 planes() const { return m_view->planes(); }

    // These difficult enable-ifs are to avoid ambagious operator overload. The
    // rule below is, if the user passes anything float like, we'll cast all of
    // the input to double. This is done with out consideration if the
    // underlining type really is floating point accessible.
    template <class T1, class T2>
    inline typename boost::enable_if<boost::mpl::and_<boost::is_integral<T1>,boost::is_integral<T2>>,pixel_type>::type
    operator()(T1 i, T2 j, int32 p=0) const {
      return m_view->operator()(int32(i),int32(j),p);
    }

    template <class T1, class T2>
    inline typename boost::disable_if<boost::mpl::and_<boost::is_integral<T1>,boost::is_integral<T2>>,pixel_type>::type
    operator()(T1 i, T2 j, int32 p=0) const {
      return m_view->operator()(double(i),double(j),p);
    }

    inline pixel_accessor origin() const {
      return m_view->origin();
    }

    inline bool sparse_check(BBox2i const& bbox) const { 
      return m_view->sparse_check(bbox);
    }

    /// \cond INTERNAL
    typedef CropView<ImageView<PixelT>> prerasterize_type;

    inline prerasterize_type prerasterize(BBox2i const& bbox) const {
      // If we're wrapping a plain ImageView, we can avoid copying the data.
      ImageViewRefImpl<ImageView<PixelT>> *image_ptr = dynamic_cast<ImageViewRefImpl<ImageView<PixelT>>*>(m_view.get());
      if (image_ptr)
       return CropView<ImageView<PixelT>>(image_ptr->child(), 0, 0, cols(), rows());
       
      // Otherwise, we must rasterize ourselves.
      ImageView<PixelT> buf(bbox.width(), bbox.height(), planes());
      m_view->rasterize(buf, bbox);
      return CropView<ImageView<PixelT>>(buf, BBox2i(-bbox.min().x(),-bbox.min().y(),cols(),rows()));
    }

    template <class DestT> inline void rasterize(DestT const& dest, BBox2i const& bbox) const {
      vw::rasterize(prerasterize(bbox), dest, bbox);
    }

    // A special performance-enhancing overload for rasterizing directly into
    // an ImageView with the proper pixel type.  This cannot be templatized
    // or otherwise generalized because it calls m_view's virtual rasterize method.
    inline void rasterize(ImageView<PixelT> const& dest, BBox2i const& bbox) const {
      m_view->rasterize(dest, bbox);
    }
    /// \endcond
  };

  template <class PixelT>
  class SparseImageCheck<ImageViewRef<PixelT>> {
    ImageViewRef<PixelT> const& image;
  public:
    SparseImageCheck(ImageViewRef<PixelT> const& image): image(image) {}
    bool operator() (BBox2i const& bbox) {
      return image.sparse_check(bbox);
    }
  };

  // Explicit template instantiation declarations for common types.
  // These are instantiated once in ImageViewRef.cc to reduce
  // compilation time and object file size.
  extern template class ImageViewRef<float>;
  extern template class ImageViewRef<double>;
  extern template class ImageViewRef<uint8>;
  extern template class ImageViewRef<PixelGray<float>>;
  extern template class ImageViewRef<PixelGray<uint8>>;
  extern template class ImageViewRef<PixelRGB<uint8>>;
  extern template class ImageViewRef<PixelMask<float>>;
  extern template class ImageViewRef<PixelMask<double>>;
  extern template class ImageViewRef<PixelMask<uint8>>;
  extern template class ImageViewRef<PixelMask<Vector2f>>;
  extern template class ImageViewRef<PixelMask<PixelRGB<uint8>>>;
  extern template class ImageViewRef<Vector2>;
  extern template class ImageViewRef<Vector3>;

} // namespace vw

#endif // __VW_IMAGE_IMAGEVIEWREF_H__
