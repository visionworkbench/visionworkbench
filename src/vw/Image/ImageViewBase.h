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

#include <vw/Core/Exception.h>
#include <vw/Core/FundamentalTypes.h>
#include <vw/Core/CompoundTypes.h>
#include <vw/Math/BBox.h>
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
  void for_each_pixel_( const ImageViewBase<ViewT> &view_, FuncT &func ) {
    const ViewT& view = view_.impl();
    typedef typename ViewT::pixel_accessor pixel_accessor;
    pixel_accessor plane_acc = view.origin();
    for (int32 p = 0; p < view.planes(); p++, plane_acc.next_plane()) { 
      pixel_accessor row_acc = plane_acc;
      for (int32 j = 0; j < view.rows(); j++, row_acc.next_row()) {
        pixel_accessor col_acc = row_acc;
        for (int32 i = 0; i < view.cols(); i++, col_acc.next_col()) {
          func( *col_acc );
        }
      }  
    }
  }

  template <class ViewT, class FuncT>
  void for_each_pixel( const ImageViewBase<ViewT> &view, FuncT &func ) {
    for_each_pixel_<ViewT,FuncT>(view,func);
  }

  template <class ViewT, class FuncT>
  void for_each_pixel( const ImageViewBase<ViewT> &view, const FuncT &func ) {
    for_each_pixel_<ViewT,const FuncT>(view,func);
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
    rasterize( src, const_cast<DestT const&>(dest), BBox2i(0,0,src.cols(),src.rows()) );
  }

} // namespace vw

#endif // __VW_IMAGE_IMAGEVIEWBASE_H__
