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

/// \file FileIO/DiskImageResourceTIFF.h
/// 
/// Provides support for TIFF image files.
///
#ifndef __VW_FILEIO_DISKIMAGERESOUCETIFF_H__
#define __VW_FILEIO_DISKIMAGERESOUCETIFF_H__

#include <string>

#include <vw/Image/PixelTypes.h>
#include <vw/FileIO/DiskImageResource.h>

namespace vw {

  class DiskImageResourceInfoTIFF;

  class DiskImageResourceTIFF : public DiskImageResource {
  public:

    DiskImageResourceTIFF( std::string const& filename );

    DiskImageResourceTIFF( std::string const& filename, 
                           ImageFormat const& format,
                           bool use_compression = false );
    
    virtual ~DiskImageResourceTIFF() {}
    
    /// Returns the type of disk image resource.
    static std::string type_static() { return "TIFF"; }

    /// Returns the type of disk image resource.
    virtual std::string type() { return type_static(); }
    
    virtual Vector2i native_block_size() const;

    virtual void read( ImageBuffer const& buf, BBox2i const& bbox ) const;

    virtual void write( ImageBuffer const& dest, BBox2i const& bbox );

    void open( std::string const& filename );

    void create( std::string const& filename,
                 ImageFormat const& format );

    static DiskImageResource* construct_open( std::string const& filename );

    static DiskImageResource* construct_create( std::string const& filename,
                                                ImageFormat const& format );

    void use_lzw_compression(bool state) { m_use_compression = state; }

  protected:
    void check_retval(const int retval, const int error_val) const;

  private:
    boost::shared_ptr<DiskImageResourceInfoTIFF> m_info;
    bool m_use_compression;
  };

} // namespace vw

#endif // __VW_FILEIO_DISK_IMAGE_RESOUCE_TIFF_H__
