// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/Camera/BundleAdjustReport.h>

namespace vw {
namespace camera {

  KMLPlaceMark::KMLPlaceMark( std::string filename,
			      std::string name ) {
    m_output_file.open( filename.c_str(), std::ios::out);
    m_tab.count = 0;
    
    if (!m_output_file.good())
      vw_throw(IOErr() <<  "An error occured while trying to write KML file.");
    
    m_output_file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    m_output_file << "<kml xmlns=\"http://www.opengis.net/kml/2.2\" xmlns:gx=\"http://www.google.com/kml/ext/2.2\" xmlns:kml=\"http://www.opengis.net/kml/2.2\" xmlns:atom=\"http://www.w3.org/2005/Atom\">\n";
    m_output_file << "<Document>\n";
    
    m_tab.count++;
    
    m_output_file << m_tab << "<name>" << name << "</name>\n";
    
    // GCP Circle
    m_output_file << m_tab << "<Style id=\"gcp_circle\">\n";
    m_output_file << m_tab << "\t<IconStyle>\n";
    m_output_file << m_tab << "\t\t<scale>1.0</scale>\n";
    m_output_file << m_tab << "\t\t<Icon>\n";
    m_output_file << m_tab << "\t\t\t<href>http://maps.google.com/mapfiles/kml/shapes/placemark_circle.png</href>\n";
    m_output_file << m_tab << "\t\t</Icon>\n";
    m_output_file << m_tab << "\t</IconStyle>\n";
    m_output_file << m_tab << "</Style>\n";
    
    // GCP Circle Highlight
    m_output_file << m_tab << "<Style id=\"gcp_circle_highlight\">\n";
    m_output_file << m_tab << "\t<IconStyle>\n";
    m_output_file << m_tab << "\t\t<scale>1.1</scale>\n";
    m_output_file << m_tab << "\t\t<Icon>\n";
    m_output_file << m_tab << "\t\t\t<href>http://maps.google.com/mapfiles/kml/shapes/placemark_circle_highlight.png</href>\n";
    m_output_file << m_tab << "\t\t</Icon>\n";
    m_output_file << m_tab << "\t</IconStyle>\n";
    m_output_file << m_tab << "</Style>\n";
  
    // GCP Placemark
    m_output_file << m_tab << "<StyleMap id=\"gcp_placemark\">\n";
    m_output_file << m_tab << "\t<Pair>\n";
    m_output_file << m_tab << "\t\t<key>normal</key>\n";
    m_output_file << m_tab << "\t\t<styleUrl>#gcp_circle</styleUrl>\n";
    m_output_file << m_tab << "\t</Pair>\n";
    m_output_file << m_tab << "\t<Pair>\n";
    m_output_file << m_tab << "\t\t<key>highlight</key>\n";
    m_output_file << m_tab << "\t\t<styleUrl>#gcp_circle_highlight</styleUrl>\n";
    m_output_file << m_tab << "\t</Pair>\n";
    m_output_file << m_tab << "</StyleMap>\n";
    
    // Est Circle
    m_output_file << m_tab << "<Style id=\"est_circle\">\n";
    m_output_file << m_tab << "\t<IconStyle>\n";
    m_output_file << m_tab << "\t\t<scale>0.3</scale>\n";
    m_output_file << m_tab << "\t\t<Icon>\n";
    m_output_file << m_tab << "\t\t\t<href>http://maps.google.com/mapfiles/kml/shapes/placemark_circle.png</href>\n";
    m_output_file << m_tab << "\t\t</Icon>\n";
    m_output_file << m_tab << "\t</IconStyle>\n";
    m_output_file << m_tab << "</Style>\n";

    // est Circle Highlight
    m_output_file << m_tab << "<Style id=\"est_circle_highlight\">\n";
    m_output_file << m_tab << "\t<IconStyle>\n";
    m_output_file << m_tab << "\t\t<scale>0.4</scale>\n";
    m_output_file << m_tab << "\t\t<Icon>\n";
    m_output_file << m_tab << "\t\t\t<href>http://maps.google.com/mapfiles/kml/shapes/placemark_circle_highlight.png</href>\n";
    m_output_file << m_tab << "\t\t</Icon>\n";
    m_output_file << m_tab << "\t</IconStyle>\n";
    m_output_file << m_tab << "</Style>\n";
    
    // est Placemark
    m_output_file << m_tab << "<StyleMap id=\"est_placemark\">\n";
    m_output_file << m_tab << "\t<Pair>\n";
    m_output_file << m_tab << "\t\t<key>normal</key>\n";
    m_output_file << m_tab << "\t\t<styleUrl>#est_circle</styleUrl>\n";
    m_output_file << m_tab << "\t</Pair>\n";
    m_output_file << m_tab << "\t<Pair>\n";
    m_output_file << m_tab << "\t\t<key>highlight</key>\n";
    m_output_file << m_tab << "\t\t<styleUrl>#est_circle_highlight</styleUrl>\n";
    m_output_file << m_tab << "\t</Pair>\n";
    m_output_file << m_tab << "</StyleMap>\n";
  }
  
  KMLPlaceMark::~KMLPlaceMark( void ) {
    m_tab.count--;
    m_output_file << m_tab << "</Document>\n";
    m_output_file << m_tab << "</kml>\n";

    m_output_file.close();
  }

  void KMLPlaceMark::write_gcps( ControlNetwork const& cnet ) {
    enter_folder( "Ground Control Points",
		  "Used for Bundle Adjustment in VW" );
    
    unsigned count = 0;
    for ( ControlNetwork::const_iterator iter = cnet.begin();
	  iter != cnet.end(); ++iter ) {
      if ( (*iter).type() == ControlPoint::GroundControlPoint ) {
	count++;
	Vector3 llr = cartography::xyz_to_lon_lat_radius( (*iter).position() );

	std::ostringstream desc;
	desc << "&lt;h1&gt;Viewed by:&lt;/h1&gt;&lt;ol&gt;";
	for ( ControlPoint::const_iterator measure = (*iter).begin();
	      measure != (*iter).end(); ++measure ) {
	  desc << "&lt;li&gt;" << (*measure).image_id() << " " << (*measure).serial() << "&lt;/li&gt;";
	}
	desc << "&lt;/ol&gt;";

	if ( (*iter).id() != "Null" ) {
	  placemark( llr.x(), llr.y(), (*iter).id(), desc.str(), true );
	} else {
	  std::ostringstream gcp_name;
	  gcp_name << "GCP " << count;
	  placemark( llr.x(), llr.y(), gcp_name.str(), desc.str(), true );
	}
      }
    }
    
    exit_folder();
  }
  
  void KMLPlaceMark::write_3d_est( ControlNetwork const& cnet ) {
    enter_folder( "3D Point estimates",
		  "Used for Bundle Adjustment in VW" );

    std::list<Vector2> points; // Lon and Lat only

    unsigned count = 0;
    for ( ControlNetwork::const_iterator iter = cnet.begin();
	  iter != cnet.end(); ++iter ) {
      count++;
      if ( (*iter).type() == ControlPoint::TiePoint ) {
	Vector3 llr = cartography::xyz_to_lon_lat_radius( (*iter).position() );
	points.push_back( Vector2(llr.x(),llr.y()) );
      }
    }
    
    // Grow a bounding box
    BBox2 total;
    for ( std::list<Vector2>::iterator it = points.begin();
	  it != points.end(); ++it )
      total.grow( (*it) );

    // Building tiles of smaller bounding boxes
    Vector2 lower_corner = total.min();
    Vector2 upper_corner = total.max();
    lower_corner.x() = floor( lower_corner.x() );
    lower_corner.y() = floor( lower_corner.y() );
    upper_corner.x() = ceil( upper_corner.x() );
    upper_corner.y() = ceil( upper_corner.y() );
   
    for ( double lon = lower_corner.x(); lon < upper_corner.x(); lon+=9 ) {
      for ( double lat = lower_corner.y(); lat < upper_corner.y(); lat+=9 ) {
	std::ostringstream section_name;
	section_name << lon << " " << lat;
	enter_folder( section_name.str(), "" );
	
	// Writing region
	m_output_file << m_tab << "<Region><LatLonAltBox><north>"<<lat+9<<"</north><south>"<<lat<<"</south><east>"
		      <<lon+9<<"</east><west>"<<lon<<"</west></LatLonAltBox><Lod><minLodPixels>256</minLodPixels>"
		      <<"<maxLodPixels>-1</maxLodPixels></Lod></Region>\n";

	for ( std::list<Vector2>::iterator it = points.begin(); it != points.end(); ++it ) {
	  if ( (*it).x() >= lon && (*it).x() <= (lon+9) && 
	       (*it).y() >= lat && (*it).y() <= (lat+9) ) {
	    placemark( (*it).x(), (*it).y(), "","", false );
	    it = points.erase(it);
	    it--;
	  }
	}

	exit_folder();
      }
    }

    exit_folder();
  }
  
  void KMLPlaceMark::enter_folder( std::string name="", std::string desc="") {
    m_output_file << m_tab << "<Folder>\n";
    m_tab.count++;
    m_output_file << m_tab << "<name>"<< name <<"</name>\n";
    m_output_file << m_tab << "<description>"<< desc <<"</description>\n";
  }

  void KMLPlaceMark::exit_folder(void) {
    m_tab.count--;
    m_output_file << m_tab << "</Folder>\n";
  }
  
  void KMLPlaceMark::placemark( double lon, double lat,
				std::string name="",
				std::string desc="", bool gcp=false ) {
    m_output_file << m_tab << "<Placemark>\n";
    m_tab.count++;
    m_output_file << m_tab << "<name>"<< name <<"</name>\n";
    m_output_file << m_tab << "<description>" << desc << "</description>\n";
    if (gcp) {
      m_output_file << m_tab << "<styleUrl>#gcp_placemark</styleUrl>\n";
    } else {
      m_output_file << m_tab << "<styleUrl>#est_placemark</styleUrl>\n";
    }
    m_output_file << m_tab << "<Point>\n";
    m_output_file << m_tab << "\t<coordinates>"<< std::setw(10) << lon <<","<< lat 
		  <<",0</coordinates>\n";
    m_output_file << m_tab << "</Point>\n";
    m_tab.count--;
    m_output_file << m_tab << "</Placemark>\n";
  }
  
  std::ostream& operator<<( std::ostream& os, TabCount const& tab) {
    for ( int i = 0; i < tab.count; i++ )
      os << "\t";
    return os;
  }

}}
