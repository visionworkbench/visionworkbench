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
#include <vw/Cartography/GeoReference.h>

// GDAL
#include "ogr_spatialref.h"

namespace vw {
namespace cartography {

  // Utility function for creating a gdal spatial reference object
  // from a vw GeoReference object.
  OGRSpatialReference gdal_spatial_ref_from_georef(GeoReference const* georef) {
    OGRSpatialReference gdal_spatial_ref;
    char* wkt_copy = new char[2048];
    strncpy(wkt_copy, georef->wkt_str().c_str(), 2048);
    gdal_spatial_ref.importFromWkt(&wkt_copy);

    // This causes the c runtime to complain about freeing an object
    // that wasn't malloc'd, but it seems like this is a potential
    // (small) memory leak.
    //
    //    delete [] wkt_copy;

    return gdal_spatial_ref;
  }

  /// Construct a default georeference.  This georeference will use
  /// the identity matrix as the initial transformation matrix, and
  /// select the default datum (???) and projection (geographic).
  GeoReference::GeoReference() {
    m_transform.set_identity();
    m_is_projected = false;
    OGRSpatialReference oSRS;
    this->set_spatial_ref(&oSRS);
  }
  
  /// Takes a void pointer to an OGRSpatialReference. The affine transform defaults to the identity matrix.
  GeoReference::GeoReference(void* spatial_ref_ptr) {
    m_transform.set_identity();
    this->set_spatial_ref(spatial_ref_ptr);
  }

  /// Takes a void pointer to an OGRSpatialReference and an affine transformation matrix.
  GeoReference::GeoReference(void* spatial_ref_ptr, Matrix<double,3,3> const& transform) {
    m_transform = transform;
    this->set_spatial_ref(spatial_ref_ptr);
  }
  
  /// Takes a string in proj.4 format. The affine transform defaults to the identity matrix.
  GeoReference::GeoReference(std::string const proj4_str) {
    m_transform.set_identity();
    OGRSpatialReference oSRS;
    oSRS.importFromProj4(proj4_str.c_str());
    this->set_spatial_ref(&oSRS);
  }

  /// Takes a string in proj.4 format and an affine transformation matrix.
  GeoReference::GeoReference(std::string const proj4_str, Matrix<double,3,3> const& transform) {
    m_transform = transform;
    OGRSpatialReference oSRS;
    oSRS.importFromProj4(proj4_str.c_str());
    this->set_spatial_ref(&oSRS);
  }    
  
  /// Takes a geodetic datum.  The affine transform defaults to the identity matrix.
  GeoReference::GeoReference(GeoDatum const& datum) {
    m_transform.set_identity();
    OGRSpatialReference oSRS;
    oSRS.SetGeogCS("Vision Workbench GeoReference", 
                   datum.name.c_str(),
                   datum.spheroid_name.c_str(),
                   double(datum.semi_major_axis), double(datum.semi_minor_axis),
                   datum.meridian_name.c_str(), double(datum.meridian_offset),
                   "degree", atof(SRS_UA_DEGREE_CONV) );  /// XXX Fixme!!
    this->set_spatial_ref(&oSRS);
  }
  
  /// Takes a geodetic datum and an affine transformation matrix
  GeoReference::GeoReference(GeoDatum const& datum, Matrix<double,3,3> const& transform) {
    m_transform = transform;
    OGRSpatialReference oSRS;
    oSRS.SetGeogCS("Vision Workbench GeoReference", 
                   datum.name.c_str(),
                   datum.spheroid_name.c_str(),
                   double(datum.semi_major_axis), double(datum.semi_minor_axis),
                   datum.meridian_name.c_str(), double(datum.meridian_offset),
                   "degree", atof(SRS_UA_DEGREE_CONV) );  /// XXX FIXME!!
    this->set_spatial_ref(&oSRS);
  }

  /// Re-initialize this spatial reference using a string in
  /// Well-Known Text (WKT) format.
  void GeoReference::set_wkt_str(std::string const& wkt_str) {
    OGRSpatialReference gdal_spatial_ref;
    char* wkt_copy = new char[2048];
    strncpy(wkt_copy, wkt_str.c_str(), 2048);
    gdal_spatial_ref.importFromWkt(&wkt_copy);
    this->set_spatial_ref(&gdal_spatial_ref);
  }

  /// Re-initialize this spatial reference using a string in Proj.4 format.
  void GeoReference::set_proj4_str(std::string const& proj4_str) {
    OGRSpatialReference gdal_spatial_ref;
    gdal_spatial_ref.importFromProj4(proj4_str.c_str());
    this->set_spatial_ref(&gdal_spatial_ref);
  }
  
  void GeoReference::set_spatial_ref(void* spatial_ref_ptr) { 

    OGRSpatialReference* gdal_spatial_ref_ptr = (OGRSpatialReference*) spatial_ref_ptr;
    
    // Grab the parameters for the georeference itself.
    char *proj4_str = NULL , *wkt_str = NULL;
    gdal_spatial_ref_ptr->exportToProj4( &proj4_str );
    gdal_spatial_ref_ptr->exportToWkt( &wkt_str );
    m_proj4_str = proj4_str;
    m_wkt_str = wkt_str;
    delete proj4_str;
    delete wkt_str;

    const char* georef_name = gdal_spatial_ref_ptr->GetAttrValue("GEOGCS");
    if (georef_name) { 
      m_name = georef_name; 
    } else { 
      m_name = "Unknown Geospatial Reference Frame"; 
    }
    m_is_projected = gdal_spatial_ref_ptr->IsProjected();
  }

  /// Returns a datum object for the current spatial reference
  GeoDatum GeoReference::datum() const {
    OGRSpatialReference gdal_spatial_ref = gdal_spatial_ref_from_georef(this);

    GeoDatum datum;
    // Set up the parameters in the geodetic datum.
    const char* datum_name = gdal_spatial_ref.GetAttrValue("DATUM");
    if (datum_name) { datum.name = datum_name; }
    const char* spheroid_name = gdal_spatial_ref.GetAttrValue("SPHEROID");
    if (spheroid_name) { datum.spheroid_name = spheroid_name; }
    const char* meridian_name = gdal_spatial_ref.GetAttrValue("PRIMEM");
    if (meridian_name) { datum.meridian_name = meridian_name; }
    OGRErr e1, e2;
    double semi_major = gdal_spatial_ref.GetSemiMajor(&e1);
    double semi_minor = gdal_spatial_ref.GetSemiMinor(&e2);
    if (e1 != OGRERR_FAILURE && e2 != OGRERR_FAILURE) { 
      datum.semi_major_axis = semi_major;
      datum.semi_minor_axis = semi_minor;
    }
    datum.meridian_offset = gdal_spatial_ref.GetPrimeMeridian();
    return datum;
  }

  /// Returns true or false depending on whether the spatial reference
  /// is projected.
  bool GeoReference::is_projected() const {
    return m_is_projected;
  }

  /// Returns a GeoProjection object.
  std::string GeoReference::projection_name() const {
    OGRSpatialReference gdal_spatial_ref = gdal_spatial_ref_from_georef(this);

    // Set up the parameters for this mapping projection
    const char* projection_name = gdal_spatial_ref.GetAttrValue("PROJECTION");
    if (projection_name) { return projection_name; }
    else { return "NONE"; }
  }

  void GeoReference::set_well_known_geogcs(std::string name) {
    OGRSpatialReference gdal_spatial_ref = gdal_spatial_ref_from_georef(this);
    gdal_spatial_ref.SetWellKnownGeogCS(name.c_str());
    set_spatial_ref(&gdal_spatial_ref);
  }

  void GeoReference::set_sinusoidal(double center_longitude, double false_easting, double false_northing) {
    OGRSpatialReference gdal_spatial_ref = gdal_spatial_ref_from_georef(this);
    gdal_spatial_ref.SetSinusoidal(center_longitude, false_easting, false_northing);
    set_spatial_ref(&gdal_spatial_ref);
  }
  
  void GeoReference::set_mercator(double center_latitude, double center_longitude, double scale, double false_easting, double false_northing) {
    OGRSpatialReference gdal_spatial_ref = gdal_spatial_ref_from_georef(this);
    gdal_spatial_ref.SetMercator(center_latitude, center_longitude, scale, false_easting, false_northing);
    set_spatial_ref(&gdal_spatial_ref);
  }

  void GeoReference::set_orthographic(double center_latitude, double center_longitude, double false_easting, double false_northing) {
    OGRSpatialReference gdal_spatial_ref = gdal_spatial_ref_from_georef(this);
    gdal_spatial_ref.SetOrthographic(center_latitude, center_longitude, false_easting, false_northing);
    set_spatial_ref(&gdal_spatial_ref);  
  }

  void GeoReference::set_stereographic(double center_latitude, double center_longitude, double scale, double false_easting, double false_northing) {
    OGRSpatialReference gdal_spatial_ref = gdal_spatial_ref_from_georef(this);
    gdal_spatial_ref.SetStereographic(center_latitude, center_longitude, scale, false_easting, false_northing);
    set_spatial_ref(&gdal_spatial_ref);
  }
  
  void GeoReference::set_UTM(int zone, int north) {
    OGRSpatialReference gdal_spatial_ref = gdal_spatial_ref_from_georef(this);
    gdal_spatial_ref.SetUTM(zone, north);
    set_spatial_ref(&gdal_spatial_ref);    
  }

}} // vw::cartography
