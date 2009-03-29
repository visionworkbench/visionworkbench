// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/Camera/BundleAdjustReport.h>
namespace fs = boost::filesystem;

namespace vw {
namespace camera {

  KMLPlaceMark::KMLPlaceMark( std::string filename,
			      std::string name,
			      std::string directory="" ) {
    
    m_tab.count = 0;
    m_name = name;
    m_directory = directory;
    boost::to_lower( m_name );
    boost::replace_all( m_name, " ", "_" );

    std::ostringstream path;
    if ( directory != "" )
      path << directory << "/";
    fs::path kml_path( path.str(), fs::native );
    fs::create_directories( kml_path );
    path << filename;
    kml_path = path.str();
    m_output_file.open( kml_path, std::ios::out);

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
    close_kml();
  }

  void KMLPlaceMark::close_kml( void ) {
    if (m_output_file.is_open()) {
      m_tab.count--;
      m_output_file << m_tab << "</Document>\n";
      m_output_file << m_tab << "</kml>\n";

      m_output_file.close();
    }
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
	// GCP data
	desc << "&lt;h2&gt;Ground Control Point&lt;/h2&gt;";
	desc << "&lt;b&gt;Lon:&lt;/b&gt; " << llr.x() << " deg&lt;br&gt;";
	desc << "&lt;b&gt;Lat:&lt;/b&gt; " << llr.y() << " deg&lt;br&gt;";
	desc << "&lt;b&gt;Rad:&lt;/b&gt; " << llr.z() << " m&lt;br&gt;";
	
	// Images viewing
	desc << "&lt;h3&gt;Viewed by:&lt;/h3&gt;&lt;ol&gt;";
	for ( ControlPoint::const_iterator measure = (*iter).begin();
	      measure != (*iter).end(); ++measure ) {
	  desc << "&lt;li&gt;" << "[" << (*measure).image_id() << "] " 
	       << (*measure).serial() << "&lt;/li&gt;";
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

    std::list<Vector3> points; // Lon, Lat, Interest

    unsigned count = 0;
    for ( ControlNetwork::const_iterator iter = cnet.begin();
	  iter != cnet.end(); ++iter ) {
      count++;
      if ( (*iter).type() == ControlPoint::TiePoint ) {
	Vector3 llr = cartography::xyz_to_lon_lat_radius( (*iter).position() );
	points.push_back( Vector3(llr.x(),llr.y(), 1.0/sum((*iter).sigma())) );
      }
    }
    
    // Grow a bounding box
    BBox2 total;
    for ( std::list<Vector3>::iterator it = points.begin();
	  it != points.end(); ++it )
      total.grow( Vector2((*it).x(),(*it).y()) );

    // Building tiles of smaller bounding boxes
    Vector2f lower_corner = total.min();
    Vector2f upper_corner = total.max();
    lower_corner.x() = floor( lower_corner.x() );
    lower_corner.y() = floor( lower_corner.y() );
    upper_corner.x() = ceil( upper_corner.x() );
    upper_corner.y() = ceil( upper_corner.y() );
   
    recursive_placemark_building( points, m_name,
				  upper_corner.y(), lower_corner.y(),
				  upper_corner.x(), lower_corner.x(), 0 );

    exit_folder();
  }
  
  void KMLPlaceMark::recursive_placemark_building( std::list<Vector3>& list,
						   std::string const& name,
						   float& north, float& south,
						   float& east, float& west,
						   int recursive_lvl) {
    // Sub divides
    std::vector<float> north_dv;
    north_dv.push_back(north); north_dv.push_back(north);
    north_dv.push_back(south+(north-south)/2);
    north_dv.push_back( north_dv[2] );
    std::vector<float> south_dv;
    south_dv.push_back( north_dv[3] ); south_dv.push_back( north_dv[3] );
    south_dv.push_back( south ); south_dv.push_back( south );
    std::vector<float> east_dv;
    east_dv.push_back( east ); east_dv.push_back( west + (east-west)/2 );
    east_dv.push_back( east_dv[0] ); east_dv.push_back( east_dv[1] );
    std::vector<float> west_dv;
    west_dv.push_back( east_dv[1] ); west_dv.push_back( west );
    west_dv.push_back( west_dv[0] ); west_dv.push_back( west_dv[1] );

    // Checking list
    int count = 0;
    for ( std::list<Vector3>::iterator it = list.begin();
	  it != list.end(); ++it )
      count++;
    if ( count <= 500 ) {
      // Write a termination file
      enter_folder("","");
      
      // Regioning
      m_output_file << m_tab << "<Region>\n";
      m_tab.count++;
      latlonaltbox( north, south, east, west );
      lod( 512, -1 );
      m_tab.count--;
      m_output_file << m_tab << "</Region>\n";

      // Placemarks
      for ( std::list<Vector3>::iterator it = list.begin();
	    it != list.end(); ++it ) 
	placemark( (*it).x(), (*it).y(), "", "", false );

      exit_folder();
      close_kml();
      return;
    } else {
      // Write a branching file
      list.sort( vector_sorting );

      enter_folder("","");

      for ( char i = 0; i < 4; i++ ) {
	std::ostringstream link;
	if ( recursive_lvl == 0 )
	  link << "data/";
	link << name << int(i) << ".kml";
	network_link( link.str(),
		      north_dv[i], south_dv[i],
		      east_dv[i], west_dv[i] );
      }
	
      enter_folder("","");

      // Regioning
      m_output_file << m_tab << "<Region>\n";
      m_tab.count++;
      latlonaltbox( north, south, east, west );
      lod( 512, -1 );
      m_tab.count--;
      m_output_file << m_tab << "</Region>\n";

      // Placemarks
      count = 500;
      std::list<Vector3>::iterator it = list.begin();
      while ( count > 0 ) {
	placemark( (*it).x(), (*it).y(), "", "", false );
	it = list.erase( it );
	count--;
      }

      exit_folder();
      exit_folder();
    }

    // Making calls to make lower levels
    for ( char i = 0; i < 4; i++ ) {
      std::list<Vector3> temp;
      for ( std::list<Vector3>::iterator it = list.begin();
	    it != list.end(); ++it )
	if ( (*it).y() < north_dv[i] && (*it).y() > south_dv[i] &&
	     (*it).x() < east_dv[i] && (*it).x() > west_dv[i] ) {
	  temp.push_back( *it );
	  it = list.erase( it );
	  it--;
	}
      if ( !temp.empty() ) {
	std::ostringstream new_name;
	new_name << name << int(i);
	std::ostringstream dir;
	if ( m_directory != "" )
	  dir << m_directory << "/";
	if ( recursive_lvl == 0 )
	  dir << "data/";
	KMLPlaceMark subsection( new_name.str() + ".kml", new_name.str(), dir.str() );
	subsection.recursive_placemark_building( temp, new_name.str(),
						 north_dv[i], south_dv[i],
						 east_dv[i], west_dv[i],
						 recursive_lvl+1 );
      }
    }

    if ( !list.empty() )
      std::cout << "Error! Vector not empty!\n";
 
  }

  void KMLPlaceMark::network_link( std::string link, 
				   double north, double south, 
				   double east, double west ) {
    m_output_file << m_tab << "<NetworkLink>\n";
    m_tab.count++;
    m_output_file << m_tab << "<Region>\n";
    m_tab.count++;
    latlonaltbox( north, south,
		  east, west );
    lod( 512, -1 );
    m_tab.count--;
    m_output_file << m_tab << "</Region>\n";
    m_output_file << m_tab << "<Link>\n";
    m_tab.count++;
    m_output_file << m_tab << "<href>" << link << "</href><viewRefreshMode>onRegion</viewRefreshMode>\n";
    m_tab.count--;
    m_output_file << m_tab << "</Link>\n";
    m_tab.count--;
    m_output_file << m_tab << "</NetworkLink>\n";
  }

  void KMLPlaceMark::enter_folder( std::string name="", std::string desc="") {
    m_output_file << m_tab << "<Folder>\n";
    m_tab.count++;
    if ( name != "" )
      m_output_file << m_tab << "<name>"<< name <<"</name>\n";
    if ( desc != "" )
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
  
  void KMLPlaceMark::latlonaltbox( float north, float south,
				   float east, float west ) {
    m_output_file << m_tab << "<LatLonAltBox>\n";
    m_tab.count++;
    m_output_file << m_tab << "<north>"<<north<<"</north>\n";
    m_output_file << m_tab << "<south>"<<south<<"</south>\n";
    m_output_file << m_tab << "<east>"<<east<<"</east>\n";
    m_output_file << m_tab << "<west>"<<west<<"</west>\n";
    m_tab.count--;
    m_output_file << m_tab << "</LatLonAltBox>\n";
  }

  void KMLPlaceMark::lod( float min, float max ) {
    m_output_file << m_tab << "<Lod>\n";
    m_tab.count++;
    m_output_file << m_tab << "<minLodPixels>"<<min<<"</minLodPixels>\n";
    m_output_file << m_tab << "<maxLodPixels>"<<max<<"</maxLodPixels>\n";
    m_tab.count--;
    m_output_file << m_tab << "</Lod>\n";
  }

  std::ostream& operator<<( std::ostream& os, TabCount const& tab) {
    for ( int i = 0; i < tab.count; i++ )
      os << "\t";
    return os;
  }
  
  bool vector_sorting( Vector3 i, Vector3 j) {
    return (i.z() > j.z());
  }

}}
