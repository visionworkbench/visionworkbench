// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NASA Vision Workbench is licensed under the Apache License,
//  Version 2.0 (the "License"); you may not use this file except in
//  compliance with the License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__


#ifndef __VW_CARTOGRAPHY_GEOREFERENCEHELPERGDAL_H__
#define __VW_CARTOGRAPHY_GEOREFERENCEHELPERGDAL_H__

#include <vw/FileIO/DiskImageResourceGDAL.h>
#include <vw/Cartography/GeoReference.h>

// Boost
#include <boost/algorithm/string.hpp>

/// \file GeoReferenceResourceGDAL.h
/// Functions to read/write a \ref vw::cartography::GeoReference from GDAL.

namespace vw {
namespace cartography {

  bool read_gdal_georeference( GeoReference& georef, DiskImageResourceGDAL const& resource );
  void write_gdal_georeference( DiskImageResourceGDAL& resource, GeoReference const& georef );

  // Read an arbitrary name = value pair from the geoheader.
  bool read_gdal_string( DiskImageResourceGDAL const& resource, std::string const& str_name,
                         std::string & str_val );

  // Read all name = value pairs from the geoheader.
  bool read_gdal_strings( DiskImageResourceGDAL const& resource, 
                          std::map<std::string, std::string> & value_pairs);

  // Write an arbitrary name = value pair in the geoheader.
  void write_gdal_string( DiskImageResourceGDAL& resource, std::string const& str_name,
                          std::string const& str_val );

}} // namespace vw::cartography

#endif // __VW_CARTOGRAPHY_GEOREFERENCEHELPERGDAL_H__
