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

/// \file DiskImageView.h
///
/// A read-only disk image view.  This is now just a thin 
/// wrapper around the more general ImageResourceView.
///
#ifndef __VW_FILEIO_DISKIMAGEVIEW_H__
#define __VW_FILEIO_DISKIMAGEVIEW_H__

#include <string>
#include <map>

#include <vw/Core/Cache.h>
#include <vw/Core/Debugging.h>
#include <vw/Image/ImageResourceView.h>
#include <vw/FileIO/DiskImageResource.h>

namespace vw {

  /// A view of an image on disk.
  template <class PixelT>
  class DiskImageView : public ImageResourceView<PixelT>
  {
    typedef ImageResourceView<PixelT> base_type;

  public:
    /// Constructs a DiskImageView of the given file on disk.
    DiskImageView( std::string const& filename, bool cache=true )
      : base_type( DiskImageResource::open( filename ), cache ) {}

    /// Constructs a DiskImageView of the given file on disk 
    /// using the specified cache area.
    DiskImageView( std::string const& filename, Cache& cache )
      : base_type( DiskImageResource::open( filename ), cache ) {}

    /// Constructs a DiskImageView of the given resource.  Takes
    /// ownership of the resouqce object (i.e. deletes it when it's
    /// done using it).
    DiskImageView( DiskImageResource *resource, bool cache=true )
      : base_type( resource, cache ) {}

    /// Constructs a DiskImageView of the given resource using the
    /// specified cache area.  Takes ownership of the resource object
    /// (i.e. deletes it when it's done using it).
    DiskImageView( DiskImageResource *resource, Cache& cache )
      : base_type( resource, cache ) {}

    virtual ~DiskImageView() {}
    
    std::string filename() const { return dynamic_cast<DiskImageResource const*>(base_type::resource())->filename(); }

  };

} // namespace vw

#endif // __VW_FILEIO_DISKIMAGEVIEW_H__
