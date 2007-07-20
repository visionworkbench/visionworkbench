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

#include <vw/Cartography/DiskImageResourceGeoReferenceHelperPDS.h>
#include <vw/FileIO/DiskImageResourcePDS.h>

// xmlParser
#include <xmlParser.h>

// libProj.4
#include <projects.h>

// Boost
#include <boost/algorithm/string.hpp>



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
  vw::cartography::GeoReference::register_disk_image_resource(vw::DiskImageResourcePDS::type_static(), &vw::cartography::DiskImageResourceGeoReferenceHelperPDS::read_georeference, &vw::cartography::DiskImageResourceGeoReferenceHelperPDS::write_georeference);
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

  // Here is some machinery to keep track of an initialized proj.4
  // projection context using a smart pointer.  Using a smart pointer
  // here simplies the rest of the GeoReference class considerably, and
  // reduces the possibility of a memory related bug.
  class ProjContext {
    PJ* m_proj_ptr;

    char** split_proj4_string(std::string const& proj4_str, int &num_strings) {
      std::vector<std::string> arg_strings;
      std::string trimmed_proj4_str = boost::trim_copy(proj4_str);
      boost::split( arg_strings, trimmed_proj4_str, boost::is_any_of(" ") ); 
    
      char** strings = new char*[arg_strings.size()]; 
      for ( unsigned i = 0; i < arg_strings.size(); ++i ) {
        strings[i] = new char[2048];
        strncpy(strings[i], arg_strings[i].c_str(), 2048);
      }
      num_strings = arg_strings.size();
      return strings;
    }

    // These are private so that we don't accidentally call them.
    // Copying a ProjContext is bad.
    ProjContext(ProjContext const& ctx) {}
    ProjContext& operator=(ProjContext const& ctx) {}

  public:
    ProjContext(std::string const& proj4_str) {

      // proj.4 is expecting the parameters to be split up into seperate
      // c-style strings.
      int num;
      char** proj_strings = split_proj4_string(proj4_str, num);
      m_proj_ptr = pj_init(num, proj_strings);
                                                
      if (!m_proj_ptr) 
        vw::vw_throw(vw::LogicErr() << "GeoReference: an error occured while initializing proj.4 with string: " << proj4_str);
      
      for (int i = 0; i < num; i++) 
        delete [] proj_strings[i];
      delete [] proj_strings;
    }

    ~ProjContext() { pj_free(m_proj_ptr); }
    inline PJ* proj_ptr() { return m_proj_ptr; }
  };



  std::string GeoReference::proj4_str() const { 
    std::string proj4_str = m_proj_projection_str + " " + m_datum.proj4_str() + " +no_defs";
    return proj4_str;
  }
  
  void GeoReference::init_proj() {
    m_proj_context = boost::shared_ptr<ProjContext>(new ProjContext(proj4_str()));
  }

  /// Construct a default georeference.  This georeference will use
  /// the identity matrix as the initial transformation matrix, and
  /// select the default datum (WGS84) and projection (geographic).
  GeoReference::GeoReference() {
    set_transform(vw::math::identity_matrix<3>());
    register_default_disk_image_resources();
    set_geographic();
    init_proj();
  }

  /// Takes a geodetic datum.  The affine transform defaults to the identity matrix.
  GeoReference::GeoReference(Datum const& datum) : GeoReferenceBase(datum){
    set_transform(vw::math::identity_matrix<3>());
    register_default_disk_image_resources();
    set_geographic();
    init_proj();
  }
  
  /// Takes a geodetic datum and an affine transformation matrix
  GeoReference::GeoReference(Datum const& datum, Matrix<double,3,3> const& transform) : GeoReferenceBase(datum) {
    set_transform(transform);
    register_default_disk_image_resources();
    set_geographic();
    init_proj();
  }

  void GeoReference::set_transform(Matrix3x3 transform) {
    m_transform = transform;
    m_shifted_transform = m_transform;
    m_shifted_transform(0,2) += 0.5*m_transform(0,0);
    m_shifted_transform(1,2) += 0.5*m_transform(1,1);
    m_inv_transform = vw::math::inverse(m_transform);
    m_inv_shifted_transform = vw::math::inverse(m_shifted_transform);
  }

  // We override the base classes method here so that we have the
  // opportunity to call init_proj()
  void GeoReference::set_datum(Datum const& datum) {
    m_datum = datum;
    init_proj();
  }
  
  // Adjust the affine transform to the VW convention ( [0,0] is at
  // the center of upper left pixel) if file is georeferenced
  // according to the convention that [0,0] is the upper left hand
  // corner of the upper left pixel.
  inline Matrix<double,3,3> GeoReference::vw_native_transform() const {
    if (pixel_interpretation() == GeoReference::PixelAsArea)
      return m_shifted_transform;
    else
      return m_transform;
  }

  inline Matrix<double,3,3> GeoReference::vw_native_inverse_transform() const {
    if (pixel_interpretation() == GeoReference::PixelAsArea) 
      return m_inv_shifted_transform;
    else
      return m_inv_transform;
  }

  void GeoReference::set_well_known_geogcs(std::string name) {
    m_datum.set_well_known_datum(name);
    init_proj();
  }

  void GeoReference::set_geographic() {
    m_is_projected = false;
    m_proj_projection_str = "+proj=longlat";
    init_proj();
  }

  void GeoReference::set_sinusoidal(double center_longitude, double false_easting, double false_northing) {
    std::ostringstream strm;
    strm << "+proj=sinu +lon_0=" << center_longitude << " +x_0=" << false_easting << " +y_0=" << false_northing << " +units=m";
    m_proj_projection_str = strm.str();
    m_is_projected = true;
    init_proj();
  }
  
  void GeoReference::set_mercator(double center_latitude, double center_longitude, double latitude_of_true_scale, double false_easting, double false_northing) {
    std::ostringstream strm;
    strm << "+proj=merc +lon_0=" << center_longitude << " +lat_0=" << center_latitude << " +lat_ts=" << latitude_of_true_scale << " +x_0=" << false_easting << " +y_0=" << false_northing << " +units=m";
    m_proj_projection_str = strm.str();
    m_is_projected = true;
    init_proj();
  }

  void GeoReference::set_transverse_mercator(double center_latitude, double center_longitude, double scale, double false_easting, double false_northing) {
    std::ostringstream strm;
    strm << "+proj=tmerc +lon_0=" << center_longitude << " +lat_0=" << center_latitude << " +k=" << scale << " +x_0=" << false_easting << " +y_0=" << false_northing << " +units=m";
    m_proj_projection_str = strm.str();
    m_is_projected = true;
    init_proj();
  }

  void GeoReference::set_orthographic(double center_latitude, double center_longitude, double false_easting, double false_northing) {
    std::ostringstream strm;
    strm << "+proj=ortho +lon_0=" << center_longitude << " +lat_0=" << center_latitude << " +x_0=" << false_easting << " +y_0=" << false_northing << " +units=m";
    m_proj_projection_str = strm.str();
    m_is_projected = true;
    init_proj();
  }

  void GeoReference::set_oblique_stereographic(double center_latitude, double center_longitude, double scale, double false_easting, double false_northing) {
    std::ostringstream strm;
    strm << "+proj=sterea +lon_0=" << center_longitude << " +lat_0=" << center_latitude << " +k=" << scale << " +x_0=" << false_easting << " +y_0=" << false_northing << " +units=m";
    m_proj_projection_str = strm.str();
    m_is_projected = true;
    init_proj();
  }

  void GeoReference::set_stereographic(double center_latitude, double center_longitude, double scale, double false_easting, double false_northing) {
    std::ostringstream strm;
    strm << "+proj=stere +lon_0=" << center_longitude << " +lat_0=" << center_latitude << " +k=" << scale << " +x_0=" << false_easting << " +y_0=" << false_northing << " +units=m";
    m_proj_projection_str = strm.str();
    m_is_projected = true;
    init_proj();
  }

  void GeoReference::set_lambert_azimuthal(double center_latitude, double center_longitude, double false_easting, double false_northing) {
    std::ostringstream strm;
    strm << "+proj=laea +lon_0=" << center_longitude << " +lat_0=" << center_latitude << " +x_0=" << false_easting << " +y_0=" << false_northing << " +units=m";
    m_proj_projection_str = strm.str();
    m_is_projected = true;
    init_proj();
  }
  
  void GeoReference::set_UTM(int zone, int north) {
    std::ostringstream strm;
    strm << "+proj=utm +zone=" << zone;
    if (!north) strm << " +south";
    strm << " +units=m";
    m_proj_projection_str = strm.str();
    m_is_projected = true;
    init_proj();
  }

  void GeoReference::set_proj4_projection_str(std::string const& s) { 
    m_proj_projection_str = s; 
    if (s.find("+proj=longlat") == 0)
      m_is_projected = false;
    else
      m_is_projected = true;
    init_proj();
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


  /// For a given pixel coordinate, compute the position of that
  /// pixel in this georeferenced space.
  Vector2 GeoReference::pixel_to_point(Vector2 pix) const {
    Vector2 loc;
    Matrix<double,3,3> M = this->vw_native_transform();
    loc[0] = pix[0] * M(0,0) + pix[1] * M(0,1) + M(0,2) / (M(2,0) + M(2,1) + M(2,2));
    loc[1] = pix[0] * M(1,0) + pix[1] * M(1,1) + M(1,2) / (M(2,0) + M(2,1) + M(2,2));
    return loc;
  }
  
  /// For a given location 'loc' in projected space, compute the
  /// corresponding pixel coordinates in the image.
  Vector2 GeoReference::point_to_pixel(Vector2 loc) const {
    Vector2 pix;
    Matrix<double,3,3> M = this->vw_native_inverse_transform();
    pix[0] = loc[0] * M(0,0) + loc[1] * M(0,1) + M(0,2) / (M(2,0) + M(2,1) + M(2,2));
    pix[1] = loc[0] * M(1,0) + loc[1] * M(1,1) + M(1,2) / (M(2,0) + M(2,1) + M(2,2));
    return pix;
  }

  
  /// For a point in the projected space, compute the position of
  /// that point in unprojected (Geographic) coordinates (lat,lon).
  Vector2 GeoReference::point_to_lonlat(Vector2 loc) const {
    XY projected;  
    LP unprojected;

    projected.u = loc[0];
    projected.v = loc[1];

    // Proj.4 expects the (lon,lat) pair to be in radians, so we
    // must make a conversion if the CS in geographic (lat/lon).
    if ( !m_is_projected ) {
      projected.u *= DEG_TO_RAD;
      projected.v *= DEG_TO_RAD;
    }
    unprojected = pj_inv(projected, m_proj_context->proj_ptr());

    // Convert from radians to degrees.
    return Vector2(unprojected.u * RAD_TO_DEG, unprojected.v * RAD_TO_DEG);
  }
  
  /// Given a position in geographic coordinates (lat,lon), compute
  /// the location in the projected coordinate system.
  Vector2 GeoReference::lonlat_to_point(Vector2 lon_lat) const {
    XY projected;  
    LP unprojected;

    // Proj.4 expects the (lon,lat) pair to be in radians, so we
    // must make a conversion if the CS in geographic (lat/lon).
    unprojected.u = lon_lat[0] * DEG_TO_RAD;
    unprojected.v = lon_lat[1] * DEG_TO_RAD;

    projected = pj_fwd(unprojected, m_proj_context->proj_ptr());

    if ( !m_is_projected ) 
      return Vector2(projected.u * RAD_TO_DEG, projected.v * RAD_TO_DEG);
    else 
      return Vector2(projected.u, projected.v);
  }

}} // vw::cartography
