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

/// \file FileIO/DiskImageResource.h
/// 
/// An abstract base class referring to an image on disk.
/// 
#ifndef __VW_FILEIO_DISK_IMAGE_RESOURCE_H__
#define __VW_FILEIO_DISK_IMAGE_RESOURCE_H__

#include <string>
#include <boost/type_traits.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>

#include <vw/Core/Debugging.h>
#include <vw/Math/BBox.h>
#include <vw/Math/Vector.h>
#include <vw/Image/ImageResource.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/PixelTypes.h>

namespace vw {

  // *******************************************************************
  // The DiskImageResource abstract base class
  // *******************************************************************

  /// Base class from which specific file handlers derive.
  class DiskImageResource : public ImageResource {
  public:

    virtual ~DiskImageResource() {};

    /// Returns the number of columns in an image on disk.
    unsigned cols() const { return m_format.cols; }

    /// Returns the number of rows in an image on disk.
    unsigned rows() const { return m_format.rows; }

    /// Returns the number of planes in an image on disk.
    unsigned planes() const { return m_format.planes; }

    /// Returns the pixel format of an image on disk.
    PixelFormatEnum pixel_format() const { return m_format.pixel_format; }

    /// Returns the channel type of an image on disk.
    ChannelTypeEnum channel_type() const { return m_format.channel_type; }

    /// Return the filename of the disk image file.
    std::string filename() const { return m_filename; }

    /// Create a new DiskImageResource of the appropriate type
    /// pointing to an existing file on disk.
    ///
    /// Don't forget to delete the DiskImageResource object when
    /// you're finished with it!
    static DiskImageResource* open( std::string const& filename );

    /// Create a new DiskImageResource of the appropriate type 
    /// pointing to a newly-created empty file on disk.
    ///
    /// The underlying driver chooses a file format that it supports
    /// that matches the one you requested as closely as possible.  
    /// If you care exactly what format it chose, you can ask the
    /// resource after it's created.  If you need finer-grained
    /// control you must manually create a resource of the appropraite
    /// type.  Don't forget to delete this DiskImageResource object
    /// when you're finished with it!
    static DiskImageResource* create( std::string const& filename, ImageFormat const& format );

    typedef DiskImageResource* (*construct_open_func)( std::string const& filename );
    
    typedef DiskImageResource* (*construct_create_func)( std::string const& filename,
                                                         ImageFormat const& format );

    static void register_file_type( std::string const& extension,
                                    construct_open_func open_func,
                                    construct_create_func create_func );

  protected:
    DiskImageResource( std::string const& filename ) : m_filename(filename) {}
    ImageFormat m_format;
    std::string m_filename;
  };


  // *******************************************************************
  // Free functions using the DiskImageResource interface
  // *******************************************************************

  /// Read an image on disk into a vw::ImageView<T> object.
  template <class PixelT>
  void read_image( ImageView<PixelT>& in_image, const std::string &filename ) {

    vw_out(InfoMessage+1) << "\tLoading image: " << filename << "\t";

    // Open the file for reading
    DiskImageResource *r = DiskImageResource::open( filename );

    vw_out(InfoMessage+1) << r->cols() << "x" << r->rows() << "x" << r->planes() << "  " << r->channels() << " channel(s)\n";

    // Read it in and wrap up
    read_image( in_image, *r );
    delete r;
  }


  /// Write an vw::ImageView<T> to disk.  If you supply a filename
  /// with an asterisk ('*'), each plane of the image will be saved as
  /// a seperate file on disk and the asterisk will be replaced with
  /// the plane number.
  template <class ImageT>
  void write_image( const std::string &filename, ImageViewBase<ImageT> const& out_image ) {

    VW_ASSERT( out_image.impl().cols() != 0 && out_image.impl().rows() != 0 && out_image.impl().planes() != 0,
               ArgumentErr() << "write_image: cannot write empty image to disk" );

    // Rasterize the image if needed
    ImageView<typename ImageT::pixel_type> image( out_image.impl() );
    ImageBuffer buf = image.buffer();

    unsigned files = 1;
    // If there's an asterisk, save one file per plane
    if( boost::find_last(filename,"*") ) {
      files = buf.format.planes;
      buf.format.planes = 1;
    }
    
    for( unsigned i=0; i<files; ++i ) {
      std::string name = filename;
      if( files > 1 ) boost::replace_last( name, "*",  str( boost::format("%1%") % i ) );
      vw_out(InfoMessage+1) << "\tSaving image: " << name << "\t";
      DiskImageResource *r = DiskImageResource::create( name, buf.format );
      vw_out(InfoMessage+1) << r->cols() << "x" << r->rows() << "x" << r->planes() << "  " << r->channels() << " channel(s)\n";
      r->write( buf, BBox2i(0,0,r->cols(),r->rows()) );
      delete r;
      buf.data = (uint8*)buf.data + buf.pstride;
    }
  }

  /// Write an vw::ImageView<T> to disk.  If you supply a filename
  /// with an asterisk ('*'), each plane of the image will be saved as
  /// a seperate file on disk and the asterisk will be replaced with
  /// the plane number.
  template <class ElemT>
  void write_image( const std::string &filename, std::vector<ElemT> const& out_image_vector ) {

    // If there's an asterisk, save one file per plane
    if( ! boost::find_last(filename,"*") ) {
      vw_throw( vw::ArgumentErr() << "write_image: filename must contain * when writing a vector of image views\n" );
    }

    for (unsigned i=0; i<out_image_vector.size(); i++){
      std::string name = filename;
      boost::replace_last( name, "*",  str( boost::format("%1%") % i ) );
      write_image( name, out_image_vector[i] );
    }
  }

} // namespace vw

#endif // __VW_FILEIO_DISK_IMAGE_RESOURCE_H__
