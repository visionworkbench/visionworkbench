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
#ifndef __VW_FILEIO_DISKIMAGERESOUCEGDAL_H__
#define __VW_FILEIO_DISKIMAGERESOUCEGDAL_H__

#include <vw/config.h>
#include <string>

// VW Headers
#include <vw/Image/PixelTypes.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/Math/Matrix.h>

namespace vw {

  class DiskImageResourceGDAL : public DiskImageResource {
  public:

    DiskImageResourceGDAL( std::string const& filename )
      : DiskImageResource( filename )
    {
      m_dataset = NULL;
      m_convert_jp2 = false;
      open( filename );
    }

    DiskImageResourceGDAL( std::string const& filename, 
                           ImageFormat const& format )
      : DiskImageResource( filename )
    {
      m_dataset = NULL;
      m_convert_jp2 = false;
      create( filename, format );
    }
    
    virtual ~DiskImageResourceGDAL() {
      flush();
    }

    /// Returns the type of disk image resource.
    static std::string type_static() { return "GDAL"; }
    
    /// Returns the type of disk image resource.
    virtual std::string type() { return type_static(); }
    
    virtual void read( ImageBuffer const& dest, BBox2i const& bbox ) const;
    virtual void write( ImageBuffer const& dest, BBox2i const& bbox );
    virtual Vector2i native_block_size() const;
    virtual void flush();

    void* dataset() { return m_dataset; }

    void open( std::string const& filename );    
    void create( std::string const& filename,
                 ImageFormat const& format );
    
    static DiskImageResource* construct_open( std::string const& filename );

    static DiskImageResource* construct_create( std::string const& filename,
                                                ImageFormat const& format );

  private:
    std::string m_filename;
    void* m_dataset;
    Matrix<double,3,3> m_geo_transform;
    bool m_convert_jp2;
  };

} // namespace vw

#endif // __VW_FILEIO_DISKIMAGERESOUCEGDAL_H__
