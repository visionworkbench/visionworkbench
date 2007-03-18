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

/// \file FileIO/DiskImageResourceOpenEXR.h
/// 
/// Provides support for the OpenEXR file format.
///
#ifndef __VW_FILEIO_DISKIMAGERESOUCEOPENEXR_H__
#define __VW_FILEIO_DISKIMAGERESOUCEOPENEXR_H__

#include <string>

#include <vw/Image/PixelTypes.h>
#include <vw/FileIO/DiskImageResource.h>

// Forward-declare the Imf::InputFile class

/// \cond INTERNAL
namespace Imf { class InputFile; }
/// \endcond

namespace vw {

  class DiskImageResourceOpenEXR : public DiskImageResource {
  public:

    DiskImageResourceOpenEXR( std::string const& filename )
      : DiskImageResource( filename )
    {
      m_file_ptr = 0;
      open( filename );
    }

    DiskImageResourceOpenEXR( std::string const& filename, 
                              ImageFormat const& format )
      : DiskImageResource( filename )
    {
      m_file_ptr = 0;
      create( filename, format );
    }
    
    virtual ~DiskImageResourceOpenEXR();
    
    virtual Vector2i native_block_size() const;

    virtual void read( ImageBuffer const& dest, BBox2i const& bbox ) const;

    virtual void write( ImageBuffer const& dest, BBox2i const& bbox );

    void open( std::string const& filename );

    void create( std::string const& filename,
                 ImageFormat const& format );

    static DiskImageResource* construct_open( std::string const& filename );

    static DiskImageResource* construct_create( std::string const& filename,
                                                ImageFormat const& format );

  private:
    std::string m_filename;
    Vector2i m_block_size;
    Imf::InputFile* m_file_ptr;
  };

} // namespace vw

#endif // __VW_FILEIO_DISKIMAGERESOUCEOPENEXR_H__
