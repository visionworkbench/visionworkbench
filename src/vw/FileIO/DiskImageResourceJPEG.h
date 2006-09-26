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

/// \file DiskImageResourceJPEG.h
/// 
/// Provides support for file formats via libJPEG.
///
#ifndef __VW_FILEIO_DISK_IMAGE_RESOUCE_JPEG_H__
#define __VW_FILEIO_DISK_IMAGE_RESOUCE_JPEG_H__

#include <vw/config.h>

#if defined(VW_HAVE_PKG_JPEG) && VW_HAVE_PKG_JPEG==1

#include <string>

#include <vw/Image/PixelTypes.h>
#include <vw/Image/GenericImageBuffer.h>
#include <vw/FileIO/DiskImageResource.h>

namespace vw {

  class DiskImageResourceJPEG : public DiskImageResource {
  public:

    DiskImageResourceJPEG( std::string const& filename )
    {
      m_file_ptr = NULL;
      m_jpg_compress_header = NULL;
      m_jpg_decompress_header = NULL;
      open( filename );
    }
    
    DiskImageResourceJPEG( std::string const& filename, 
                           GenericImageFormat const& format )
    {
      m_quality = 0.85;
      m_file_ptr = NULL;
      m_jpg_compress_header = NULL;
      m_jpg_decompress_header = NULL;
      create( filename, format );
    }
    
    virtual ~DiskImageResourceJPEG();
    
    virtual void read( GenericImageBuffer const& dest ) const;
    virtual void write( GenericImageBuffer const& dest );
    virtual void flush();

    /// Set the compression quality of the jpeg image.  The quality is
    /// a value between 0.0 and 1.0.  The lower the quality, the more
    /// lossy the compression.
    void set_quality(float quality) { m_quality = quality; }

    void open( std::string const& filename );

    void create( std::string const& filename,
                 GenericImageFormat const& format );

    static DiskImageResource* construct_open( std::string const& filename );

    static DiskImageResource* construct_create( std::string const& filename,
                                                GenericImageFormat const& format );

  private:
    
    std::string m_filename;
    float m_quality;
    void* m_jpg_decompress_header;
    void* m_jpg_compress_header;
    void* m_file_ptr;
  };

} // namespace vw

#endif // HAVE_PGK_JPEG

#endif // __VW_FILEIO_DISK_IMAGE_RESOUCE_JPEG_H__
