// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_CARTOGRAPHY_GEOREFERENCEHELPERPDS_H__
#define __VW_CARTOGRAPHY_GEOREFERENCEHELPERPDS_H__

#include <vw/FileIO/DiskImageResourcePDS.h>
#include <vw/Cartography/GeoReference.h>

// Boost
#include <boost/algorithm/string.hpp>

namespace vw {
namespace cartography {

  void read_pds_georeference( GeoReference& georef, DiskImageResourcePDS const& resource );
  // We do not support writing PDS images at this time.
  // void write_pds_georeference( DiskImageResourcePDS& resource, GeoReference const& georef );

}} // namespace vw::cartography

#endif // __VW_CARTOGRAPHY_GEOREFERENCEHELPERPDS_H__
