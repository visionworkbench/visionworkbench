// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file KML.cc
/// 
/// An abstract base class referring to an image on disk.
/// 

// Vision Workbench
#include <iomanip>
#include <vw/FileIO/KML.h>
#include <vw/Math/EulerAngles.h>
#include <vw/Cartography/SimplePointImageManipulation.h>

namespace fs = boost::filesystem;

namespace vw {

  // TAB COUNT ////////////////////////////////////////////////////
  
  std::ostream& operator<<( std::ostream& os, TabCount const& tab) {
    for ( int i = 0; i < tab.count; i++ )
      os << "\t";
    return os;
  }

  // KML FILE /////////////////////////////////////////////////////

  // Constructor / Deconstructor
  KMLFile::KMLFile( std::string filename,
		    std::string name,
		    std::string directory ) : m_filename(filename), m_name(name), m_directory(directory) {
    m_tab.count = 0;
    boost::to_lower( m_name );
    boost::replace_all( m_name, " ", "_" );

    open_kml();
  }

  KMLFile::~KMLFile( void ) {
    close_kml();
  }

  // Lower Level Writing Functions

  void KMLFile::open_bracket( std::string name ) {
    m_bracket_names.push(name);
    m_output_file << m_tab << "<" << name << ">\n";
    m_tab.count++;
  }
  void KMLFile::close_bracket( void ) {
    m_tab.count--;
    m_output_file << m_tab << "</" << m_bracket_names.top() << ">\n";
    m_bracket_names.pop();
  }
  void KMLFile::close_brackets( int i ) {
    while ( i != 0 ) {
      close_bracket();
      i--;
    }
  }
  void KMLFile::close_all_brackets( void ) {
    while ( !m_bracket_names.empty() )
      close_bracket();
  }

  // Low Level Functions
  
  // Enter / Exit Folder 
  // This creates those actual folders in GE that the user can turn on
  // or off.
  void KMLFile::enter_folder( std::string name, 
			      std::string desc ) {
    open_bracket("Folder");
    if ( name != "" )
      m_output_file << m_tab << "<name>"<< name <<"</name>\n";
    if ( desc != "" )
      m_output_file << m_tab << "<description>"<< desc <<"</description>\n";
  }

  void KMLFile::exit_folder(void) {
    close_bracket();
  }

  // Append Objects
  
  // Placemark: This creates the simple generic pushpin marker
  void KMLFile::append_placemark( double lon, double lat,
				  std::string name,
				  std::string description,
				  std::string style,
				  double altitude,
				  bool extrude ) {
    open_bracket("Placemark");
    if ( name != "" )
      m_output_file << m_tab << "<name>"<< name <<"</name>\n";
    if ( description != "" )
      m_output_file << m_tab << "<description>" 
		    << description << "</description>\n";
    if ( style != "") 
      m_output_file << m_tab << "<styleUrl>#"<<style<<"</styleUrl>\n";
    open_bracket("Point");
    if ( extrude )
      m_output_file << m_tab << "<extrude>1</extrude>\n";
    m_output_file << m_tab << "<altitudeMode>absolute</altitudeMode>\n";
    m_output_file << m_tab << "<coordinates>"<< std::setw(10) 
		  << lon <<","<< lat <<"," << altitude
		  << "</coordinates>\n";
    close_brackets(2);
  }

  // Cordinate: This is a specialized display model. This just draws a
  // coordinate frame model at a position with pose. Coordinate vector
  // is in ECEF.
  void KMLFile::append_coordinate( vw::Vector3 position,
				   vw::Quaternion<double> pose,
				   std::string name,
				   std::string description,
				   float scale ) {
    cartography::XYZtoLonLatRadFunctor func;
    Vector3 lon_lat_alt = func(position);

    // Converts from GE's default rotation which is oriented over the
    // site frame to a standard planetocentric rotation.
    Matrix3x3 correction_rot = vw::math::euler_to_rotation_matrix((90-lon_lat_alt(1))*M_PI/180, (90+lon_lat_alt(0))*M_PI/180, 0, "xzy");

    Vector3 angles = rotation_matrix_to_euler_zxy(pose.rotation_matrix()*correction_rot);
    double heading = angles(0)*180/M_PI, tilt = angles(1)*180/M_PI, roll = angles(2)*180/M_PI;

    open_bracket("Placemark");
    if ( name != "" )
      m_output_file << m_tab << "<name>"<< name <<"</name>\n";
    if ( description != "" )
      m_output_file << m_tab << "<description>" 
		    << description << "</description>\n";
    
    open_bracket("LookAt");
    m_output_file << m_tab << "<longitude>" << lon_lat_alt(0) << "</longitude>\n";
    m_output_file << m_tab << "<latitude> " << lon_lat_alt(1) << "</latitude>\n";
    m_output_file << m_tab << "<altitude> " << lon_lat_alt(2)-3396200 << "</altitude>\n";
    m_output_file << m_tab << "<range> " << 1e6 << "</range>\n";
    m_output_file << m_tab << "<tilt>" << 0 << "</tilt>\n";
    m_output_file << m_tab << "<heading>" << 0 << "</heading>\n";
    close_bracket();

    m_output_file << m_tab << "<Model id=\"model_" << name << "\">\n";
    m_tab.count++;
    m_output_file << m_tab << "<altitudeMode>absolute</altitudeMode>\n";
    open_bracket("Location");
    m_output_file << m_tab << "<longitude>" << lon_lat_alt(0) << "</longitude>\n";
    m_output_file << m_tab << "<latitude> " << lon_lat_alt(1) << "</latitude>\n";
    m_output_file << m_tab << "<altitude> " << lon_lat_alt(2)-3396200 << "</altitude>\n";
    close_bracket();
    open_bracket("Orientation");
    m_output_file << m_tab << "<heading>" << heading << "</heading>\n";
    m_output_file << m_tab << "<tilt>" << tilt << "</tilt>\n";
    m_output_file << m_tab << "<roll>" << roll << "</roll>\n";
    close_bracket();
    open_bracket("Scale");
    m_output_file << m_tab << "<x>" << 3000*scale << "</x>\n";
    m_output_file << m_tab << "<y>" << 3000*scale << "</y>\n";
    m_output_file << m_tab << "<z>" << 3000*scale << "</z>\n";
    close_bracket();
    open_bracket("Link");
    m_output_file << m_tab << "<href>/Users/mbroxton/projects/axis.dae</href>\n";
    close_bracket();
    m_tab.count--;
    m_output_file << m_tab << "</Model>\n";

    close_bracket();
  }

  // LatLonAltBox: This is a bounding box, that only displays contents
  // when viewer is inside box.
  void KMLFile::append_latlonaltbox( float north, float south,
				     float east, float west ) {
    open_bracket("LatLonAltBox");
    m_output_file << m_tab << "<north>"<<north<<"</north>\n";
    m_output_file << m_tab << "<south>"<<south<<"</south>\n";
    m_output_file << m_tab << "<east>"<<east<<"</east>\n";
    m_output_file << m_tab << "<west>"<<west<<"</west>\n";
    close_bracket();
  }

  // Lod: This sets the min max pixel viewing range for an object
  void KMLFile::append_lod( float min, float max ) {
    open_bracket("Lod");
    m_output_file << m_tab << "<minLodPixels>"<<min<<"</minLodPixels>\n";
    m_output_file << m_tab << "<maxLodPixels>"<<max<<"</maxLodPixels>\n";
    close_bracket();
  }

  // Style: Defines an Icon to use later
  void KMLFile::append_style( std::string id, std::string color_hex,
			      float scale, std::string image_url ) {
    m_output_file << m_tab << "<Style id=\"" << id << "\">\n";
    m_tab.count++;
    open_bracket("IconStyle");
    if (color_hex != "" )
      m_output_file << m_tab << "<color>" << color_hex << "</color>\n";
    m_output_file << m_tab << "<scale>" << scale << "</scale>\n";
    open_bracket("Icon");
    m_output_file << m_tab << "\t<href>" << image_url << "</href>\n";
    close_brackets(2);
    m_tab.count--;
    m_output_file << m_tab << "</Style>\n";
  }

  // StyleMap: Maps two styles together to create a bipolar icon
  void KMLFile::append_stylemap( std::string id, 
				 std::string style_normal,
				 std::string style_highlight ) {
    m_output_file << m_tab << "<StyleMap id=\"" << id << "\">\n";
    m_tab.count++;
    open_bracket("Pair");
    m_output_file << m_tab << "<key>normal</key>\n";
    m_output_file << m_tab << "<styleUrl>#" << style_normal << "</styleUrl>\n";
    close_bracket();
    open_bracket("Pair");
    m_output_file << m_tab << "<key>highlight</key>\n";
    m_output_file << m_tab << "<styleUrl>#" << style_highlight << "</styleUrl>\n";
    close_bracket();
    m_tab.count--;
    m_output_file << m_tab << "</StyleMap>\n";
  }

  // NetworkLink: Links to another KML file, this also contains an LOD
  // that conditionals the opening of the link.
  void KMLFile::append_network( std::string link,
				double north, double south,
				double east, double west ) {
    open_bracket("NetworkLink");
    open_bracket("Region");
    append_latlonaltbox( north, south,
			 east, west );
    append_lod( 512, -1 );
    close_bracket();
    open_bracket("Link");
    m_output_file << m_tab << "<href>" << link << "</href><viewRefreshMode>onRegion</viewRefreshMode>\n";
    close_brackets(2);
  }

  // Open / Close Stuff
  void KMLFile::open_kml( void ) {
    std::ostringstream path;
    if ( m_directory != "" )
      path << m_directory << "/";
    fs::path kml_path( path.str(), fs::native );
    fs::create_directories( kml_path );
    path << m_filename;
    kml_path = path.str();
    m_output_file.open( kml_path, std::ios::out);

    if (!m_output_file.good())
      vw_throw(IOErr() <<  "An error occured while trying to write KML file.");
    
    m_output_file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    m_output_file << "<kml xmlns=\"http://www.opengis.net/kml/2.2\" xmlns:gx=\"http://www.google.com/kml/ext/2.2\" xmlns:kml=\"http://www.opengis.net/kml/2.2\" xmlns:atom=\"http://www.w3.org/2005/Atom\">\n";
    m_output_file << "<Document>\n";
    
    m_tab.count++;
    
    m_output_file << m_tab << "<name>" << m_name << "</name>\n";
  }

  void KMLFile::close_kml( void ) {
    if (m_output_file.is_open()) {
      if (!m_bracket_names.empty())
	vw_throw(IOErr() << "Error on close out, there seems to be an open bracket somewhere left in the kml.");

      m_tab.count--;
      m_output_file << m_tab << "</Document>\n";
      m_output_file << m_tab << "</kml>\n";

      m_output_file.close();
    }
  }
}
