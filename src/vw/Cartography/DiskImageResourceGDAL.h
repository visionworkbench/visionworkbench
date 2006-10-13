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

/// \file DiskImageResourceGDAL.h
/// 
/// Provides support for georeferenced files via the GDAL library.
///
#ifndef __VW_DISK_IMAGE_RESOUCE_GDAL_H__
#define __VW_DISK_IMAGE_RESOUCE_GDAL_H__

#include <vw/config.h>
#include <string>

// VW Headers
#include <vw/Image/PixelTypes.h>
#include <vw/Image/GenericImageBuffer.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/Math/Matrix.h>
#include <vw/Cartography/GeoReference.h>

namespace vw {
namespace cartography {

  class DiskImageResourceGDAL : public DiskImageResource {
  public:

    DiskImageResourceGDAL( std::string const& filename )
    {
      m_dataset = NULL;
      open( filename );
    }

    DiskImageResourceGDAL( std::string const& filename, 
                           GenericImageFormat const& format )
    {
      m_dataset = NULL;
      create( filename, format );
    }
    
    virtual ~DiskImageResourceGDAL() {
      flush();
    }
    
    virtual void read( GenericImageBuffer const& dest ) const;
    virtual void write( GenericImageBuffer const& dest );
    virtual void flush();

    void read_georeference( GeoReference& georef );
    void write_georeference( GeoReference const& georef );

    void open( std::string const& filename );    
    void create( std::string const& filename,
                 GenericImageFormat const& format );
    
    static DiskImageResource* construct_open( std::string const& filename );

    static DiskImageResource* construct_create( std::string const& filename,
                                                GenericImageFormat const& format );

  private:
    std::string m_filename;
    void* m_dataset;
    Matrix<double,3,3> m_geo_transform;
  };

}} // namespace vw::cartography

#endif // __VW_DISK_IMAGE_RESOUCE_GDAL_H__
