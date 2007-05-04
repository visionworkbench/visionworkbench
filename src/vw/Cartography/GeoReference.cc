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
#include <vw/FileIO/PropertyMultiMap.h>

//NOTE: conditional is for when GeoReference becomes independent of GDAL
#if defined(VW_HAVE_PKG_GDAL) && VW_HAVE_PKG_GDAL
#include <vw/Cartography/DiskImageResourceGeoReferenceHelperGDAL.h>
#include <vw/FileIO/DiskImageResourceGDAL.h>
#endif

#if defined(VW_HAVE_PKG_JPEG2K) && VW_HAVE_PKG_JPEG2K
#include <vw/Cartography/DiskImageResourceGeoReferenceHelperJP2.h>
#include <vw/FileIO/DiskImageResourceJP2.h>
#endif

// xmlParser
#include <xmlParser.h>

// GDAL
#include "ogr_spatialref.h"

namespace {
  vw::FileMetadata::ReadMetadataMapType *read_map = 0;
  vw::FileMetadata::WriteMetadataMapType *write_map = 0;
}

void vw::cartography::GeoReference::register_disk_image_resource(std::string const& disk_image_resource_type,
                                                                 vw::FileMetadata::read_metadata_func read_func,
                                                                 vw::FileMetadata::write_metadata_func write_func) {
  if(!read_map) read_map = new vw::FileMetadata::ReadMetadataMapType();
  if(!write_map) write_map = new vw::FileMetadata::WriteMetadataMapType();
  read_map->insert(std::make_pair(disk_image_resource_type, read_func));
  write_map->insert(std::make_pair(disk_image_resource_type, write_func));
  vw::DiskImageResource::register_metadata_type(disk_image_resource_type, vw::cartography::GeoReference::metadata_type_static());
}

static void register_default_disk_image_resources() {
  static bool already = false;
  if( already ) return;
  already = true;

#if defined(VW_HAVE_PKG_GDAL) && VW_HAVE_PKG_GDAL==1
  vw::cartography::GeoReference::register_disk_image_resource(vw::DiskImageResourceGDAL::type_static(), &vw::cartography::DiskImageResourceGeoReferenceHelperGDAL::read_georeference, &vw::cartography::DiskImageResourceGeoReferenceHelperGDAL::write_georeference);
#endif

#if defined(VW_HAVE_PKG_JPEG2K) && VW_HAVE_PKG_JPEG2K==1
  vw::cartography::GeoReference::register_disk_image_resource(vw::DiskImageResourceJP2::type_static(), &vw::cartography::DiskImageResourceGeoReferenceHelperJP2::read_georeference, &vw::cartography::DiskImageResourceGeoReferenceHelperJP2::write_georeference);
#endif
}

namespace {

  // Make GML from a GeoReference
  //NOTE: caller must free returned string
  //FIXME: needs to take a GeoReference instead of just making boilerplate GML
  //FIXME: remember that affine transform is [lat; lon] = dehom(A * [u; v; 1])
  char* make_gml(vw::cartography::GeoReference const& georef) {
    char* s;
    XMLNode t, n1, n2, n3, n4, n5, n6, n7, n8, n9, n10;
    
    t = XMLNode::createXMLTopNode("xml", TRUE);
    t.addAttribute("version", "1.0");
    t.addAttribute("encoding", "UTF-8");
    
    //NOTE: this was initially generated from a sample gml file using gml_print_code("/home/ttemplet/gmljp2/gmlsamples/minimalgml.xml")
    
    n1 = t.addChild("gml:FeatureCollection");
    n1.addAttribute("xmlns", "http://www.opengis.net/gml");
    n1.addAttribute("xmlns:gml", "http://www.opengis.net/gml");
    n1.addAttribute("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance");
    n1.addAttribute("xsi:schemaLocation", "http://www.opengis.net/gml gmlJP2Profile.xsd"); //FIXME: see p. 41 of GMLJP2 standard, which appears to say that gmlJP2Profile.xsd must be included in the jp2 file and cannot be linked like in the commented-out line below
    //n1.addAttribute("xsi:schemaLocation", "http://www.opengis.net/gml http://schemas.opengis.net/gml/3.1.1/profiles/gmlJP2Profile/1.0.0/gmlJP2Profile.xsd");
    
    n2 = n1.addChild("gml:boundedBy");
    n3 = n2.addChild("gml:Envelope");
    n4 = n3.addChild("gml:lowerCorner");
    n4.addText("270379.500 3942462.000");
    n4 = n3.addChild("gml:upperCorner");
    n4.addText("518842.500 3942462.000");
    
    n2 = n1.addChild("gml:featureMember");
    n3 = n2.addChild("gml:FeatureCollection");
    
    n4 = n3.addChild("gml:boundedBy");
    n5 = n4.addChild("gml:Envelope");
    n6 = n5.addChild("gml:lowerCorner");
    n6.addText("270379.500 3942462.000");
    n6 = n5.addChild("gml:upperCorner");
    n6.addText("518842.500 3942462.000");
    
    n4 = n3.addChild("gml:featureMember");
    n5 = n4.addChild("gml:RectifiedGridCoverage");
    n5.addAttribute("dimension", "2");
    n5.addAttribute("gml:id", "RGC0001");
    n6 = n5.addChild("gml:description");
    n6.addText("This GMLJP2 Minimal Root Instance contains a GML Rectified Grid. The rectified grid is embedded in a RectifiedGridCoverage with generic range parameters (to be ignored).");
    n6 = n5.addChild("gml:rectifiedGridDomain");
    n7 = n6.addChild("gml:RectifiedGrid");
    n7.addAttribute("dimension", "2");
    
    n8 = n7.addChild("gml:limits");
    n9 = n8.addChild("gml:GridEnvelope");
    n10 = n9.addChild("gml:low");
    n10.addText("0 0");
    n10 = n9.addChild("gml:high");
    n10.addText("8718 7812");
    
    n8 = n7.addChild("gml:axisName");
    n8.addText("x");
    n8 = n7.addChild("gml:axisName");
    n8.addText("y");
    
    n8 = n7.addChild("gml:origin");
    n9 = n8.addChild("gml:Point");
    n9.addAttribute("gml:id", "Pt001");
    n9.addAttribute("srsName", "urn:ogc:def:crs:EPSG:6.6:32612");
    n10 = n9.addChild("gml:description");
    n10.addText("\"Upper-left\" image origin");
    n10 = n9.addChild("gml:coordinates");
    n10.addText("270379.500000, 3942462.000000");
    
    n8 = n7.addChild("gml:offsetVector");
    n8.addAttribute("srsName", "urn:ogc:def:crs:EPSG:6.6:32612");
    n8.addText("28.5 0");
    n8 = n7.addChild("gml:offsetVector");
    n8.addAttribute("srsName", "urn:ogc:def:crs:EPSG:6.6:32612");
    n8.addText("0 28.5");
    
    n6 = n5.addChild("gml:rangeSet");
    n7 = n6.addChild("gml:File");
    n8 = n7.addChild("gml:rangeParameters");
    n9 = n8.addChild("gml:QuantityList");
    n9.addAttribute("uom", "urn:ogc:def:crs:EPSG:6.6:32612");
    n9.addText("inapplicable");
    n8 = n7.addChild("gml:fileName");
    n8.addText("Not Applicable");
    n8 = n7.addChild("gml:fileStructure");
    n8.addText("Record Interleaved");
    
    n6 = n5.addChild("gml:coverageFunction");
    n7 = n6.addChild("gml:GridFunction");
    n8 = n7.addChild("gml:sequenceRule");
    n8.addAttribute("order", "+x+y");
    n8.addText("Linear");
    n8 = n7.addChild("gml:startPoint");
    n8.addText("0 0");
    
    s = t.createXMLString(false);
    return s;
  }
  
  #if 0
  // Recursive implementation of gml_print_code()
  void gml_print_code_recursive(XMLNode const& n, int d) {
    int num_attributes, num_text, num_children;
    int i;
    
    std::cout << "n" << d << " = n" << (d - 1) << ".addChild(\"" << (strchr(n.getName(), ':') ? "" : "gml:") << n.getName() << "\");" << std::endl;
    
    num_attributes = n.nAttribute();
    for(i = 0; i < num_attributes; i++)
      std::cout << "n" << d << ".addAttribute(\"" << n.getAttributeName(i) << "\", \"" << n.getAttributeValue(i) << "\");" << std::endl;
      
    num_text = n.nText();
    for(i = 0; i < num_text; i++)
      std::cout << "n" << d << ".addText(\"" << n.getText(i) << "\");" << std::endl;
  
    num_children = n.nChildNode();
    for(i = 0; i < num_children; i++)
      gml_print_code_recursive(n.getChildNode(i), d + 1);
  }
  
  // Print C++ code to create the XML file fn
  void gml_print_code(const char* fn) {
    XMLNode t;
     
    t = XMLNode::openFileHelper(fn, "FeatureCollection");
    gml_print_code_recursive(t, 1);
  }
  #endif
  
}

namespace vw {
namespace cartography {

  // Utility function for creating a gdal spatial reference object
  // from a vw GeoReference object.
  OGRSpatialReference gdal_spatial_ref_from_georef(GeoReference const* georef) {
    OGRSpatialReference gdal_spatial_ref;
    char wkt_copy[2048];
    strncpy(wkt_copy, georef->wkt_str().c_str(), 2048);
    char* wkt_ptr = &(wkt_copy[0]);
    gdal_spatial_ref.importFromWkt(&wkt_ptr);
    return gdal_spatial_ref;
  }

  /// Construct a default georeference.  This georeference will use
  /// the identity matrix as the initial transformation matrix, and
  /// select the default datum (WGS84) and projection (geographic).
  GeoReference::GeoReference() {
    m_transform.set_identity();
    m_is_projected = false;
    OGRSpatialReference oSRS;
    oSRS.SetWellKnownGeogCS("WGS84");
    this->set_spatial_ref(&oSRS);
    register_default_disk_image_resources();
  }
  
  /// Takes a void pointer to an OGRSpatialReference. The affine transform defaults to the identity matrix.
  GeoReference::GeoReference(void* spatial_ref_ptr) {
    m_transform.set_identity();
    this->set_spatial_ref(spatial_ref_ptr);
    register_default_disk_image_resources();
  }

  /// Takes a void pointer to an OGRSpatialReference and an affine transformation matrix.
  GeoReference::GeoReference(void* spatial_ref_ptr, Matrix<double,3,3> const& transform) {
    m_transform = transform;
    this->set_spatial_ref(spatial_ref_ptr);
    register_default_disk_image_resources();
  }
  
  /// Takes a string in proj.4 format. The affine transform defaults to the identity matrix.
  GeoReference::GeoReference(std::string const proj4_str) {
    m_transform.set_identity();
    OGRSpatialReference oSRS;
    oSRS.importFromProj4(proj4_str.c_str());
    this->set_spatial_ref(&oSRS);
    register_default_disk_image_resources();
  }

  /// Takes a string in proj.4 format and an affine transformation matrix.
  GeoReference::GeoReference(std::string const proj4_str, Matrix<double,3,3> const& transform) {
    m_transform = transform;
    OGRSpatialReference oSRS;
    oSRS.importFromProj4(proj4_str.c_str());
    this->set_spatial_ref(&oSRS);
    register_default_disk_image_resources();
  }    
  
  /// Takes a geodetic datum.  The affine transform defaults to the identity matrix.
  GeoReference::GeoReference(GeoDatum const& datum) {
    m_transform.set_identity();
    OGRSpatialReference oSRS;
    oSRS.SetGeogCS("Vision Workbench GeoReference", 
                   datum.name().c_str(),
                   datum.spheroid_name().c_str(),
                   double(datum.semi_major_axis()), double(datum.semi_minor_axis()),
                   datum.meridian_name().c_str(), double(datum.meridian_offset()),
                   "degree", atof(SRS_UA_DEGREE_CONV) );
    this->set_spatial_ref(&oSRS);
    register_default_disk_image_resources();
  }
  
  /// Takes a geodetic datum and an affine transformation matrix
  GeoReference::GeoReference(GeoDatum const& datum, Matrix<double,3,3> const& transform) {
    m_transform = transform;
    OGRSpatialReference oSRS;
    oSRS.SetGeogCS("Vision Workbench GeoReference", 
                   datum.name().c_str(),
                   datum.spheroid_name().c_str(),
                   double(datum.semi_major_axis()), double(datum.semi_minor_axis()),
                   datum.meridian_name().c_str(), double(datum.meridian_offset()),
                   "degree", atof(SRS_UA_DEGREE_CONV) );
    this->set_spatial_ref(&oSRS);
    register_default_disk_image_resources();
  }

  /// Re-initialize this spatial reference using a string in
  /// Well-Known Text (WKT) format.
  void GeoReference::set_wkt_str(std::string const& wkt_str) {
    OGRSpatialReference gdal_spatial_ref;
    char wkt_copy[2048];
    strncpy(wkt_copy, wkt_str.c_str(), 2048);
    char* wkt_ptr = &(wkt_copy[0]);
    gdal_spatial_ref.importFromWkt(&wkt_ptr);
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
    char *proj4_str = NULL , *wkt_str = NULL , *gml_str = NULL;
    gdal_spatial_ref_ptr->exportToProj4( &proj4_str );
    gdal_spatial_ref_ptr->exportToWkt( &wkt_str );
    gml_str = make_gml( *this );
    m_proj4_str = proj4_str;
    m_wkt_str = wkt_str;
    m_gml_str = gml_str;
    delete proj4_str;
    delete wkt_str;
    free(gml_str);

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
    if (datum_name) { datum.name() = datum_name; }
    const char* spheroid_name = gdal_spatial_ref.GetAttrValue("SPHEROID");
    if (spheroid_name) { datum.spheroid_name() = spheroid_name; }
    const char* meridian_name = gdal_spatial_ref.GetAttrValue("PRIMEM");
    if (meridian_name) { datum.meridian_name() = meridian_name; }
    OGRErr e1, e2;
    double semi_major = gdal_spatial_ref.GetSemiMajor(&e1);
    double semi_minor = gdal_spatial_ref.GetSemiMinor(&e2);
    if (e1 != OGRERR_FAILURE && e2 != OGRERR_FAILURE) { 
      datum.semi_major_axis() = semi_major;
      datum.semi_minor_axis() = semi_minor;
    }
    datum.meridian_offset() = gdal_spatial_ref.GetPrimeMeridian();
    return datum;
  }

  /// Return the box that bounds the area represented by the
  /// geotransform for an image of the given dimensions.
  BBox2 GeoReference::bounding_box(int width, int height) const {
    
    Vector3 upper_left = m_transform * Vector3(0,0,1);
    Vector3 lower_right = m_transform * Vector3(width,height,1);
    Vector3 upper_right = m_transform * Vector3(width,0,1);
    Vector3 lower_left = m_transform * Vector3(0,height,1);

    // Renormalize in case these transforms are not homgeneous.  This
    // would be very unusual but it's better to be rigorous, I think.
    upper_left /= upper_left(2);
    lower_right /= lower_right(2);
    upper_right /= upper_right(2);
    lower_left /= lower_left(2);

    BBox2 final_result;
    final_result.grow(subvector(upper_left,0,2));
    final_result.grow(subvector(lower_left,0,2));
    final_result.grow(subvector(upper_right,0,2));
    final_result.grow(subvector(lower_right,0,2));

    return final_result;
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

  void GeoReference::set_transverse_mercator(double center_latitude, double center_longitude, double scale, double false_easting, double false_northing) {
    OGRSpatialReference gdal_spatial_ref = gdal_spatial_ref_from_georef(this);
    gdal_spatial_ref.SetTM(center_latitude, center_longitude, scale, false_easting, false_northing);
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

  void GeoReference::set_polar_stereographic(double center_latitude, double center_longitude, double scale, double false_easting, double false_northing) {
    OGRSpatialReference gdal_spatial_ref = gdal_spatial_ref_from_georef(this);
    gdal_spatial_ref.SetPS(center_latitude, center_longitude, scale, false_easting, false_northing);
    set_spatial_ref(&gdal_spatial_ref);
  }

  void GeoReference::set_lambert_azimuthal(double center_latitude, double center_longitude, double false_easting, double false_northing) {
    OGRSpatialReference gdal_spatial_ref = gdal_spatial_ref_from_georef(this);
    gdal_spatial_ref.SetLAEA(center_latitude, center_longitude, false_easting, false_northing);
    set_spatial_ref(&gdal_spatial_ref);
  }
  
  void GeoReference::set_UTM(int zone, int north) {
    OGRSpatialReference gdal_spatial_ref = gdal_spatial_ref_from_georef(this);
    gdal_spatial_ref.SetUTM(zone, north);
    set_spatial_ref(&gdal_spatial_ref);    
  }
  
  /*virtual*/ void GeoReference::read_file_metadata(DiskImageResource* r) {
    if(read_map) {
      std::pair<std::string,read_metadata_func> item;
      bool found_func = read_map->find(item, r->type());
      if(found_func) {
        item.second(this, r);
        return;
      }
    }
    vw_throw(NoImplErr() << "Unsuppported DiskImageResource: " << r->type());
  }
  
  /*virtual*/ void GeoReference::write_file_metadata(DiskImageResource* r) const {
    if(write_map) {
      std::pair<std::string,write_metadata_func> item;
      bool found_func = write_map->find(item, r->type());
      if(found_func) {
        item.second(r, this);
        return;
      }
    }
    vw_throw(NoImplErr() << "Unsuppported DiskImageResource: " << r->type());
  }

}} // vw::cartography
