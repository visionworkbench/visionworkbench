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
#ifndef __VW_CARTOGRAPHY_DISKIMAGERESOURCEGEOREFERENCEHELPERPDS_H__
#define __VW_CARTOGRAPHY_DISKIMAGERESOURCEGEOREFERENCEHELPERPDS_H__

#include <vw/FileIO/FileMetadata.h>
#include <vw/FileIO/DiskImageResource.h>

// Boost
#include <boost/algorithm/string.hpp>

namespace vw {
namespace cartography {

  // This implements some very preliminary support for reading
  // georeferencing info from a PDS image.  So far, this is only
  // guranteed to work on the MOLA MEGDR data products.
  struct DiskImageResourceGeoReferenceHelperPDS {
    static void read_georeference( FileMetadata* georef_, DiskImageResource* r_ );
    static void write_georeference( DiskImageResource* r_, FileMetadata const* georef_ );
  };

}} // namespace vw::cartography

#endif // __VW_CARTOGRAPHY_DISKIMAGERESOURCEGEOREFERENCEHELPERPDS_H__
