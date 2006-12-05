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
#ifndef __VW_FILEIO_DISK_IMAGE_RESOUCE_TIFF_H__
#define __VW_FILEIO_DISK_IMAGE_RESOUCE_TIFF_H__

#include <string>

#include <vw/Image/PixelTypes.h>
#include <vw/Image/GenericImageBuffer.h>
#include <vw/FileIO/DiskImageResource.h>

namespace vw {

  struct DiskImageResourceInfoTIFF;

  class DiskImageResourceTIFF : public DiskImageResource {
  public:

    DiskImageResourceTIFF( std::string const& filename );

    DiskImageResourceTIFF( std::string const& filename, 
                           GenericImageFormat const& format );
    
    virtual ~DiskImageResourceTIFF() {}
    
    virtual Vector2i native_read_block_size() const;
    virtual void read_generic( GenericImageBuffer const& buf, BBox2i bbox ) const;
    virtual void read_generic( GenericImageBuffer const& dest ) const {
      read_generic( dest, BBox2i(0,0,cols(),rows()) );
    }

    virtual void write_generic( GenericImageBuffer const& dest );
    virtual void flush() {}

    void open( std::string const& filename );

    void create( std::string const& filename,
                 GenericImageFormat const& format );

    static DiskImageResource* construct_open( std::string const& filename );

    static DiskImageResource* construct_create( std::string const& filename,
                                                GenericImageFormat const& format );

  private:
    boost::shared_ptr<DiskImageResourceInfoTIFF> m_info;
  };

} // namespace vw

#endif // __VW_FILEIO_DISK_IMAGE_RESOUCE_TIFF_H__
