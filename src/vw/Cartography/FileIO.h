// __BEGIN_LICENSE__
//
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
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

/// \file Cartography/FileIO.h
/// 
/// Function for reading and writing georeferenced images to/from disk.
/// The file format is inferred from the suffix of the filename.
///
#ifndef __VW_CARTOGRAPHY_FILE_IO_H__
#define __VW_CARTOGRAPHY_FILE_IO_H__

#include <boost/type_traits.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>

#include <vw/config.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/GenericImageBuffer.h>
#include <vw/Cartography/DiskImageResourceGDAL.h>

#include <iostream>

namespace vw {
namespace cartography {

  /// Read an image on disk into a vw::ImageView<T> object.  If you
  /// supply a filename with an asterisk ('*'), each plane of the
  /// image will be read from a seperate file on disk with the asterisk
  /// replaced by the plane number.
  template <class PixelT>
  void read_georeferenced_image( ImageView<PixelT>& in_image, GeoReference& georef, const std::string &filename ) {

    // XXX If there's an asterisk in the filename, handle it here.
    if( boost::find_last(filename,"*") )
      throw NoImplErr() << "File-per-plane reading is not yet supported.";

    std::cout << "\tLoading georeferenced image: " << filename << "\t";

    // Open the file for reading
    DiskImageResource *r = DiskImageResourceGDAL::construct_open( filename );

    std::cout << r->cols() << "x" << r->rows() << "x" << r->planes() << "  " << r->channels() << " channel(s)\n";

    // Read it in and wrap up
    static_cast<DiskImageResourceGDAL*>(r)->read_georeference(georef);
    r->read( in_image );
    delete r;
  }


  /// Write an vw::ImageView<T> to disk.  If you supply a filename
  /// with an asterisk ('*'), each plane of the image will be saved as
  /// a seperate file on disk and the asterisk will be replaced with
  /// the plane number.
  template <class ImageT>
  void write_georeferenced_image(const std::string &filename, ImageViewBase<ImageT> const& out_image, GeoReference const& georef) {

    // Rasterize the image if needed
    ImageView<typename ImageT::pixel_type> image( out_image.impl() );
    GenericImageBuffer buf(image);

    int files = 1;
    // If there's an asterisk, save one file per plane
    if( boost::find_last(filename,"*") ) {
      files = buf.format.planes;
      buf.format.planes = 1;
    }
    
    for( int i=0; i<files; ++i ) {
      std::string name = filename;
      if( files > 1 ) boost::replace_last( name, "*",  str( boost::format("%1%") % i ) );
      std::cout << "\tSaving image: " << name << "\t";
      DiskImageResource *r = DiskImageResourceGDAL::construct_create( name, buf.format );
      std::cout << r->cols() << "x" << r->rows() << "x" << r->planes() << "  " << r->channels() << " channel(s)\n";
      static_cast<DiskImageResourceGDAL*>(r)->write_georeference(georef);
      r->write( buf );
      delete r;
      buf.data += buf.pstride;
    }
  }

}} // namespace vw::cartography

#endif // __VW_CARTOGRAPHY_FILEIO_H__  
