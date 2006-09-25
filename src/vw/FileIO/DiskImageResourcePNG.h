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

/// \file DiskImageResourcePNG.h
/// 
/// Provides support for the PNG file format.
///
#ifndef __VW_FILEIO_DISK_IMAGE_RESOUCE_PNG_H__
#define __VW_FILEIO_DISK_IMAGE_RESOUCE_PNG_H__

#include <vw/config.h>

#if defined(VW_HAVE_PKG_PNG) && VW_HAVE_PKG_PNG==1

#include <string>
#include <fstream>
#include <boost/shared_ptr.hpp>

#include <vw/Image/PixelTypes.h>
#include <vw/Image/GenericImageBuffer.h>
#include <vw/FileIO/DiskImageResource.h>

namespace vw {

  struct DiskImageResourceInfoPNG;

  class DiskImageResourcePNG : public DiskImageResource {
  public:

    DiskImageResourcePNG( std::string const& filename );

    DiskImageResourcePNG( std::string const& filename, 
                          GenericImageFormat const& format );
    
    virtual ~DiskImageResourcePNG();
    
    virtual void read( GenericImageBuffer const& dest ) const;
    virtual void write( GenericImageBuffer const& dest );
    virtual void flush() {}

    void open( std::string const& filename );

    void create( std::string const& filename,
                 GenericImageFormat const& format );

    static DiskImageResource* construct_open( std::string const& filename );

    static DiskImageResource* construct_create( std::string const& filename,
                                                GenericImageFormat const& format );

  private:
    boost::shared_ptr<DiskImageResourceInfoPNG> m_info;
  };

} // namespace vw

#endif // HAVE_PKG_PNG

#endif // __VW_FILEIO_DISK_IMAGE_RESOUCE_PNG_H__
