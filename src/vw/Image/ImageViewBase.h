// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file ImageViewBase.h
///
/// Image view base class and default rasterization code.
///
/// This file provides the core image view functionality.  You should not
/// need to care about this file unless you are writing your own view.
/// First there is a templatized base class for image views, \ref vw::ImageViewBase.
/// Then there are the default templates for several image view traits
/// classes.  Finally, there is the default rasterization function,
/// \ref vw::rasterize, which iterates over source and destination views
/// copying pixels from one into the other.
///
#ifndef __VW_IMAGE_IMAGEVIEWBASE_H__
#define __VW_IMAGE_IMAGEVIEWBASE_H__

#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>

#include <vw/Core/ProgressCallback.h>
#include <vw/Image/ImageResource.h>
#include <vw/Image/PixelIterator.h>

namespace vw {

  // *******************************************************************
  // The core image view CRTP base class
  // *******************************************************************

  /// A CRTP image view base class.
  // Consider a function that takes a single argument, which should be
  // an arbitrary image view.  You may not want to define it like this:
  //    template<class ImageViewT> void foo(ImageViewT& v) {...}
  // because that will accept any type, not just an image view type.
  // Instead, to resolve ambiguities, you can define it like this:
  //    template<class T> void foo(ImageViewBase<T>& v) {...}
  // which affords you better compile-time type checking and greater
  // function overload flexibility.
  template <class ImplT>
  struct ImageViewBase {

    /// \cond INTERNAL
    // Methods to access the derived type
    inline ImplT& impl() { return static_cast<ImplT&>(*this); }
    inline ImplT const& impl() const { return static_cast<ImplT const&>(*this); }
    /// \endcond

    /// An STL-compatible iterator type.
    typedef PixelIterator<ImplT> iterator;

    /// Returns an iterator pointing to the first pixel in the image.
    iterator begin() const { return iterator(impl(),0,0,0); }

    /// Returns an iterator pointing one past the last pixel in the image.
    iterator end() const { return iterator(impl(),0,0,impl().planes()); }

    /// Returns the number of channels in the image's pixel type.
    inline int32 channels() const { return CompoundNumChannels<typename ImplT::pixel_type>::value; }

    /// Returns the format ID of the image's pixel type.
    inline PixelFormatEnum pixel_format() const { return PixelFormatID<typename ImplT::pixel_type>::value; }

    /// Returns the channel type ID of the image's pixel type.
    inline ChannelTypeEnum channel_type() const { return ChannelTypeID<typename CompoundChannelType<typename ImplT::pixel_type>::type>::value; }

    /// Returns an ImageFormat object describing the image format.
    ImageFormat format() const {
      ImageFormat format;
      format.cols = impl().cols();
      format.rows = impl().rows();
      format.planes = impl().planes();
      format.pixel_format = pixel_format();
      format.channel_type = channel_type();
      return format;
    }

  /// \cond INTERNAL
    protected:
    // These are defined here to prevent the user from accidentally
    // copy constructing an ImageViewBase.
    ImageViewBase() {}
    ImageViewBase(ImageViewBase const&) {}
    ImageViewBase& operator=(ImageViewBase const&) { return *this; }
    /// \endcond
  };


  // *******************************************************************
  // Default image traits templates
  // *******************************************************************

  /// Indicates whether a view can be resized via <B>set_size()</B>.
  template <class ImplT>
  struct IsResizable : public false_type {};

  /// Indicates whether a view type can be accessed at floating-point positions.
  template <class ImplT>
  struct IsFloatingPointIndexable : public false_type {};

  /// Indicates whether or not a type is an image view type.
  template <class ImageT>
  struct IsImageView : public boost::is_base_of<ImageViewBase<ImageT>,ImageT>::type {};

  /// Indicates whether or not a view can be accessed multiple times
  /// just as efficiently as a locally-cached version.
  template <class ImplT>
  struct IsMultiplyAccessible : public false_type {};


  // *******************************************************************
  // Pixel iteration functions
  // *******************************************************************

  template <class ViewT, class FuncT>
  void for_each_pixel_( const ImageViewBase<ViewT> &view_, FuncT &func, const ProgressCallback &progress ) {
    const ViewT& view = view_.impl();
    typedef typename ViewT::pixel_accessor pixel_accessor;
    pixel_accessor plane_acc = view.origin();
    for( int32 plane = view.planes(); plane; --plane ) {
      pixel_accessor row_acc = plane_acc;
      for( int32 row = 0; row<view.rows(); ++row ) {
        progress.report_fractional_progress(row,view.rows());
        pixel_accessor col_acc = row_acc;
        for( int32 col = view.cols(); col; --col ) {
          func( *col_acc );
          col_acc.next_col();
        }
        row_acc.next_row();
      }
      plane_acc.next_plane();
    }
    progress.report_finished();
  }

  template <class ViewT, class FuncT>
  void for_each_pixel( const ImageViewBase<ViewT> &view, FuncT &func, const ProgressCallback &progress = ProgressCallback::dummy_instance() ) {
    for_each_pixel_<ViewT,FuncT>(view,func,progress);
  }

  template <class ViewT, class FuncT>
  void for_each_pixel( const ImageViewBase<ViewT> &view, const FuncT &func, const ProgressCallback &progress = ProgressCallback::dummy_instance() ) {
    for_each_pixel_<ViewT,const FuncT>(view,func,progress);
  }

  template <class View1T, class View2T, class FuncT>
  void for_each_pixel_( const ImageViewBase<View1T> &view1_, const ImageViewBase<View2T> &view2_, FuncT &func ) {
    const View1T& view1 = view1_.impl();
    const View2T& view2 = view2_.impl();
    VW_ASSERT( view1.cols()==view2.cols() && view1.rows()==view2.rows() && view1.planes()==view2.planes(),
               ArgumentErr() << "for_each_pixel_: Image arguments must have the same dimensions." );
    typedef typename View1T::pixel_accessor pixel_accessor_1;
    typedef typename View2T::pixel_accessor pixel_accessor_2;
    pixel_accessor_1 plane_acc_1 = view1.origin();
    pixel_accessor_2 plane_acc_2 = view2.origin();
    for( int32 plane = view1.planes(); plane; --plane ) {
      pixel_accessor_1 row_acc_1 = plane_acc_1;
      pixel_accessor_2 row_acc_2 = plane_acc_2;
      for( int32 row = view1.rows(); row; --row ) {
        pixel_accessor_1 col_acc_1 = row_acc_1;
        pixel_accessor_2 col_acc_2 = row_acc_2;
        for( int32 col = view1.cols(); col; --col ) {
          func( *col_acc_1, *col_acc_2 );
          col_acc_1.next_col();
          col_acc_2.next_col();
        }
        row_acc_1.next_row();
        row_acc_2.next_row();
      }
      plane_acc_1.next_plane();
      plane_acc_2.next_plane();
    }
  }

  template <class View1T, class View2T, class FuncT>
  void for_each_pixel( const ImageViewBase<View1T> &view1, const ImageViewBase<View2T> &view2, FuncT &func ) {
    for_each_pixel_<View1T,View2T,FuncT>(view1,view2,func);
  }

  template <class View1T, class View2T, class FuncT>
  void for_each_pixel( const ImageViewBase<View1T> &view1, const ImageViewBase<View2T> &view2, const FuncT &func ) {
    for_each_pixel_<View1T,View2T,const FuncT>(view1,view2,func);
  }

  template <class View1T, class View2T, class View3T, class FuncT>
  void for_each_pixel_( const ImageViewBase<View1T> &view1_, const ImageViewBase<View2T> &view2_, const ImageViewBase<View3T> &view3_, FuncT &func ) {
    const View1T& view1 = view1_.impl();
    const View2T& view2 = view2_.impl();
    const View3T& view3 = view3_.impl();
    VW_ASSERT( view1.cols()==view2.cols() && view1.rows()==view2.rows() && view1.planes()==view2.planes() &&
               view1.cols()==view3.cols() && view1.rows()==view3.rows() && view1.planes()==view3.planes(),
               ArgumentErr() << "for_each_pixel_: Image arguments must have the same dimensions." );
    typedef typename View1T::pixel_accessor pixel_accessor_1;
    typedef typename View2T::pixel_accessor pixel_accessor_2;
    typedef typename View3T::pixel_accessor pixel_accessor_3;
    pixel_accessor_1 plane_acc_1 = view1.origin();
    pixel_accessor_2 plane_acc_2 = view2.origin();
    pixel_accessor_3 plane_acc_3 = view3.origin();
    for( int32 plane = view1.planes(); plane; --plane ) {
      pixel_accessor_1 row_acc_1 = plane_acc_1;
      pixel_accessor_2 row_acc_2 = plane_acc_2;
      pixel_accessor_3 row_acc_3 = plane_acc_3;
      for( int32 row = view1.rows(); row; --row ) {
        pixel_accessor_1 col_acc_1 = row_acc_1;
        pixel_accessor_2 col_acc_2 = row_acc_2;
        pixel_accessor_3 col_acc_3 = row_acc_3;
        for( int32 col = view1.cols(); col; --col ) {
          func( *col_acc_1, *col_acc_2, *col_acc_3 );
          col_acc_1.next_col();
          col_acc_2.next_col();
          col_acc_3.next_col();
        }
        row_acc_1.next_row();
        row_acc_2.next_row();
        row_acc_3.next_row();
      }
      plane_acc_1.next_plane();
      plane_acc_2.next_plane();
      plane_acc_3.next_plane();
    }
  }

  template <class View1T, class View2T, class View3T, class FuncT>
  void for_each_pixel( const ImageViewBase<View1T> &view1, const ImageViewBase<View2T> &view2, const ImageViewBase<View3T> &view3, FuncT &func ) {
    for_each_pixel_<View1T,View2T,View3T,FuncT>(view1,view2,view3,func);
  }

  template <class View1T, class View2T, class View3T, class FuncT>
  void for_each_pixel( const ImageViewBase<View1T> &view1, const ImageViewBase<View2T> &view2, const ImageViewBase<View3T> &view3, const FuncT &func ) {
    for_each_pixel_<View1T,View2T,View3T,const FuncT>(view1,view2,view3,func);
  }

  // *******************************************************************
  // The master rasterization function
  // *******************************************************************

  /// This function is called by image views that do not have special
  /// optimized rasterization methods.  The user can also call it
  /// explicitly when pixel-by-pixel rasterization is preferred to
  /// the default optimized rasterization behavior.  This can be
  /// useful in some cases, such as when the views are heavily
  /// subsampled.
  template <class SrcT, class DestT>
  inline void rasterize( SrcT const& src, DestT const& dest, BBox2i bbox ) {
    typedef typename DestT::pixel_type DestPixelT;
    typedef typename SrcT::pixel_accessor SrcAccT;
    typedef typename DestT::pixel_accessor DestAccT;
    VW_ASSERT( int(dest.cols())==bbox.width() && int(dest.rows())==bbox.height() && dest.planes()==src.planes(),
               ArgumentErr() << "rasterize: Source and destination must have same dimensions." );
    SrcAccT splane = src.origin().advance(bbox.min().x(),bbox.min().y());
    DestAccT dplane = dest.origin();
    for( int32 plane=src.planes(); plane; --plane ) {
      SrcAccT srow = splane;
      DestAccT drow = dplane;
      for( int32 row=bbox.height(); row; --row ) {
        SrcAccT scol = srow;
        DestAccT dcol = drow;
        for( int32 col=bbox.width(); col; --col ) {
          *dcol = DestPixelT(*scol);
          scol.next_col();
          dcol.next_col();
        }
        srow.next_row();
        drow.next_row();
      }
      splane.next_plane();
      dplane.next_plane();
    }
  }

  /// A convenience overload to rasterize the entire source.
  template <class SrcT, class DestT>
  inline void rasterize( SrcT const& src, DestT const& dest ) {
    rasterize( src, dest, BBox2i(0,0,src.cols(),src.rows()) );
  }

  /// A specialization for resizable destination views.
  ///
  /// This function resizes the destination view prior to
  /// rasterization.
  /// \see vw::rasterize
  template <class SrcT, class DestT>
  typename boost::enable_if<IsResizable<DestT>, void>::type
  inline rasterize( SrcT const& src, DestT& dest, BBox2i bbox ) {
    dest.set_size( bbox.width(), bbox.height(), src.planes() );
    rasterize( src, const_cast<DestT const&>(dest), bbox );
  }

  /// A convenience overload to rasterize the entire source.
  template <class SrcT, class DestT>
  typename boost::enable_if<IsResizable<DestT>, void>::type
  inline rasterize( SrcT const& src, DestT& dest ) {
    dest.set_size( src.cols(), src.rows(), src.planes() );
    rasterize( src, const_cast<DestT const&>(dest), BBox2i(0,0,src.cols(),src.rows()) );
  }

} // namespace vw

#endif // __VW_IMAGE_IMAGEVIEWBASE_H__
