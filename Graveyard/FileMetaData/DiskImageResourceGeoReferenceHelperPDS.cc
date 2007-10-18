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
#include <vw/Math/Matrix.h>
#include <vw/Core/Exception.h>

#include <vw/FileIO/DiskImageResourcePDS.h>
#include <vw/Cartography/DiskImageResourceGeoReferenceHelperPDS.h>
#include <vw/Cartography/Datum.h>
#include <vw/Cartography/GeoReference.h>

// For boost::split
#include <boost/algorithm/string.hpp>

namespace vw {
namespace cartography {

  /*static*/ void DiskImageResourceGeoReferenceHelperPDS::read_georeference( FileMetadata* georef_, DiskImageResource* r_ ) {
    DiskImageResourcePDS* r = (DiskImageResourcePDS*)r_;
    GeoReference* georef = (GeoReference*)georef_;
    
    std::string key, value;
    bool failed = false;
    
    Datum datum;
    datum.name() = "PDS Datum";
    datum.spheroid_name() = "PDS Spheroid";
    // Collect datum information
    key = "A_AXIS_RADIUS";
    if ( r->query(key,value) ) { datum.set_semi_major_axis(atof(value.c_str()) * 1000); } else {failed = true; }
    key = "B_AXIS_RADIUS";
    if ( r->query(key,value) ) { datum.set_semi_minor_axis(atof(value.c_str()) * 1000); } else {failed = true; }
    georef->set_datum(datum);

    // Collect projection information
    key = "MAP_PROJECTION_TYPE";
    if ( r->query(key,value) ) { 
      if (value != "\"SIMPLE CYLINDRICAL\"") {
        vw_out(ErrorMessage) << "Unsupported map projection type in PDS header."; 
        failed = true;
      }
      datum.set_semi_minor_axis(atof(value.c_str()) * 1000); 
    } else {failed = true; }


    // Set affine transform
    Matrix<double,3,3> transform;
    transform.set_identity();

    key = "WESTERNMOST_LONGITUDE";
    if ( r->query(key,value) ) { transform(0,2) = atof(value.c_str()); } else {failed = true; }
    key = "MAXIMUM_LATITUDE";
    if ( r->query(key,value) ) { transform(1,2) = atof(value.c_str()); } else {failed = true; }
    key = "MAP_RESOLUTION";
    if ( r->query(key,value) ) { 
      transform(0,0) = 1/atof(value.c_str()); 
      transform(1,1) = -1/atof(value.c_str()); 
    } else {failed = true; }

    georef->set_transform(transform);

    if (failed) 
      vw_throw(IOErr() << "DiskImageResourcePDS: Error reading georeferencing information");
  }
  
  /*static*/ void DiskImageResourceGeoReferenceHelperPDS::write_georeference( DiskImageResource* r_, FileMetadata const* georef_ ) {
    vw_throw(IOErr() << "DiskImageResourcePDS: does not support writing of georeferncing information");
  }

}} // namespace vw::cartography
