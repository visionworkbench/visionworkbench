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

/// \file GenericImageBuffer.h
/// 
/// Describes a run-type-specified image buffer.
///
#ifndef __VW_IMAGE_GENERIC_IMAGE_BUFFER_H__
#define __VW_IMAGE_GENERIC_IMAGE_BUFFER_H__

#include <vw/Image/PixelTypes.h>
#include <vw/Image/ImageView.h>

namespace vw {

  /// Describes the format of an image, i.e. its dimensions, pixel
  /// structure, and channel type.  Chiefly used by the
  /// GenericImageBuffer class.
  struct GenericImageFormat {
    unsigned cols, rows, planes;
    PixelFormatEnum pixel_format;
    ChannelTypeEnum channel_type;
    
    GenericImageFormat()
      : cols(0), rows(0), planes(0),
        pixel_format(VW_PIXEL_UNKNOWN),
        channel_type(VW_CHANNEL_UNKNOWN)
    {}

    template <class ImageT>
    GenericImageFormat( ImageViewBase<ImageT> const& image )
      : cols( image.impl().cols() ),
        rows( image.impl().rows() ),
        planes( image.impl().planes() ),
        pixel_format( PixelFormatID<typename ImageT::pixel_type>::value ),
        channel_type( ChannelTypeID<typename CompoundChannelType<typename ImageT::pixel_type>::type>::value )
    {}
  };


  /// Represents a generic image buffer in memory, with dimensions and
  /// pixel format specified at run time.  This class does not
  /// allocate any memory, but rather provides a common format for
  /// describing an existing in-memory buffer of pixels.  The primary
  /// purpose of this class is to provide some common ground for
  /// converting between image formats using the convert() function.
  /// To allocate a fresh buffer for an image, see ImageView.
  struct GenericImageBuffer {
    void* data;
    GenericImageFormat format;
    ptrdiff_t cstride, rstride, pstride;
    bool unpremultiplied;
    
    /// Return a pointer to the pixel at (u,v,p)
    void* operator()( unsigned i, unsigned j, unsigned p = 0 ) const {
      // Cast 
      return ((uint8*)data) + (i*cstride + j*rstride + p*pstride);
    }
    
    /// Default constructor; constructs an undefined buffer
    GenericImageBuffer()
      : data(0), format(),
        cstride(0), rstride(0), pstride(0),
        unpremultiplied(false)
    {}

    /// Constructs a GenericImageBuffer pointing to an ImageView's data
    template <class PixelT>
    GenericImageBuffer( ImageView<PixelT> const& image )
      : data( image.data() ),
        format( image ),
        cstride( sizeof(PixelT) ),
        rstride( sizeof(PixelT)*format.cols ),
        pstride( sizeof(PixelT)*format.cols*format.rows ),
        unpremultiplied( false )
    {}

  };


  /// Copies image pixel data from the source buffer to the destination 
  /// buffer, converting the pixel format and channel type as required.
  void convert( GenericImageBuffer const& dst, GenericImageBuffer const& src );

} // namespace vw

#endif // __VW_IMAGE_GENERIC_IMAGE_BUFFER_H__
