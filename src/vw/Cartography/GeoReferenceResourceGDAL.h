// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_CARTOGRAPHY_GEOREFERENCEHELPERGDAL_H__
#define __VW_CARTOGRAPHY_GEOREFERENCEHELPERGDAL_H__

#include <vw/FileIO/DiskImageResourceGDAL.h>
#include <vw/Cartography/GeoReference.h>

// Boost
#include <boost/algorithm/string.hpp>

namespace vw {
namespace cartography {

  bool read_gdal_georeference( GeoReference& georef, DiskImageResourceGDAL const& resource );
  void write_gdal_georeference( DiskImageResourceGDAL& resource, GeoReference const& georef );

}} // namespace vw::cartography

#endif // __VW_CARTOGRAPHY_GEOREFERENCEHELPERGDAL_H__
