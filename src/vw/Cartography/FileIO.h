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
#include <vw/FileIO/DiskImageResource.h>
#include <vw/FileIO/FileMetadata.h>

#include <iostream>

namespace vw {
namespace cartography {

  /// Read an image on disk into a vw::ImageView<T> object.
  template <class PixelT>
  void read_georeferenced_image( ImageView<PixelT>& in_image, GeoReference& georef, const std::string &filename ) {

    FileMetadataCollection fmeta;
    fmeta.associate_file_metadata( &georef );
    read_image( in_image, fmeta, filename );
  }


  /// Read in only the georeference object.
  inline void read_georeference( GeoReference& georef, const std::string &filename ) {
    FileMetadataCollection fmeta;
    fmeta.associate_file_metadata( &georef );
    read_metadata( fmeta, filename );
  }


  /// Write an vw::ImageView<T> to disk.  If you supply a filename
  /// with an asterisk ('*'), each plane of the image will be saved as
  /// a seperate file on disk and the asterisk will be replaced with
  /// the plane number.
  template <class ImageT>
  void write_georeferenced_image(const std::string &filename, ImageViewBase<ImageT> const& out_image, GeoReference const& georef, const ProgressCallback &progress_callback = ProgressCallback::dummy_instance() ) {

    FileMetadataCollection fmeta;
    fmeta.associate_file_metadata_const( &georef );
    write_image( filename, out_image, fmeta, progress_callback);
  }

}} // namespace vw::cartography

#endif // __VW_CARTOGRAPHY_FILEIO_H__  
