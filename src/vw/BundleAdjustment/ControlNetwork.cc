// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file ControlNetwork.cc
///

#include <vw/BundleAdjustment/ControlNetwork.h>
#include <boost/algorithm/string.hpp>

// Time Headers
#include <boost/thread/xtime.hpp>
#include <ctime>

#ifdef VW_HAVE_UNISTD_H
#include <unistd.h>
#endif

 ///////////////////////////////////
// Global only to Control Network //
///////////////////////////////////

inline std::string current_posix_time_string() {
  char time_string[2048];
  time_t t = time(0);
  struct tm* time_struct = localtime(&t);
  strftime(time_string, 2048, "%F %T", time_struct);
  return std::string(time_string);
}

inline std::string isis_style_time_string() {
  std::string time = current_posix_time_string();
  boost::erase_all( time, "\n" );
  boost::trim( time );
  boost::replace_all( time, " ", "T" );
  return time;
}

namespace vw {
namespace ba {

  ////////////////////////////
  // Control Measure        //
  ////////////////////////////

  /// Constructor
  ControlMeasure::ControlMeasure( float col, float row,
                                  float col_sigma, float row_sigma,
                                  size_t image_id,
                                  ControlMeasureType type ) :
    m_col(col), m_row(row), m_col_sigma(col_sigma), m_row_sigma(row_sigma), m_image_id(image_id), m_type(type) {

    // Recording time
    m_date_time = isis_style_time_string();

    m_serialNumber = "Null";
    m_description = "Null";
    m_ignore = false;
    m_pixels_dominant = true;

    m_ephemeris_time = m_focalplane_x = m_focalplane_y = m_diameter = 0;
  }

  ControlMeasure::ControlMeasure( ControlMeasureType type ) :
    m_col(0), m_row(0), m_col_sigma(0), m_row_sigma(0), m_image_id(0), m_type(type) {

    // Recording time
    m_date_time = isis_style_time_string();

    m_serialNumber = "Null";
    m_description = "Null";
    m_ignore = false;
    m_pixels_dominant = true;

    m_ephemeris_time = m_focalplane_x = m_focalplane_y = m_diameter = 0;
  }

  /// Write a compressed binary style of measure
  void ControlMeasure::write_binary( std::ostream &f ) {
    // Writing out all the strings first
    f << m_serialNumber << char(0) << m_date_time << char(0)
      << m_description << char(0) << m_chooserName << char(0);
    // Writing the binary data
    f.write((char*)&(m_col), sizeof(m_col));
    f.write((char*)&(m_row), sizeof(m_row));
    f.write((char*)&(m_col_sigma), sizeof(m_col_sigma));
    f.write((char*)&(m_row_sigma), sizeof(m_row_sigma));
    f.write((char*)&(m_diameter), sizeof(m_diameter));
    f.write((char*)&(m_focalplane_x), sizeof(m_focalplane_x));
    f.write((char*)&(m_focalplane_y), sizeof(m_focalplane_y));
    f.write((char*)&(m_ephemeris_time), sizeof(m_ephemeris_time));
    f.write((char*)&(m_image_id), sizeof(m_image_id));
    f.write((char*)&(m_ignore), sizeof(m_ignore));
    f.write((char*)&(m_pixels_dominant), sizeof(m_pixels_dominant));
    f.write((char*)&(m_type), sizeof(m_type));
  }

  /// Reading a compressed binary style of measure
  void ControlMeasure::read_binary( std::istream &f ) {
    // Reading in all the strings
    std::getline( f, m_serialNumber, '\0' );
    std::getline( f, m_date_time, '\0' );
    std::getline( f, m_description, '\0' );
    std::getline( f, m_chooserName, '\0' );
    // Reading the binary data
    f.read((char*)&(m_col),             sizeof(m_col));
    f.read((char*)&(m_row),             sizeof(m_row));
    f.read((char*)&(m_col_sigma),       sizeof(m_col_sigma));
    f.read((char*)&(m_row_sigma),       sizeof(m_row_sigma));
    f.read((char*)&(m_diameter),        sizeof(m_diameter));
    f.read((char*)&(m_focalplane_x),    sizeof(m_focalplane_x));
    f.read((char*)&(m_focalplane_y),    sizeof(m_focalplane_y));
    f.read((char*)&(m_ephemeris_time),  sizeof(m_ephemeris_time));
    f.read((char*)&(m_image_id),        sizeof(m_image_id));
    f.read((char*)&(m_ignore),          sizeof(m_ignore));
    f.read((char*)&(m_pixels_dominant), sizeof(m_pixels_dominant));
    f.read((char*)&(m_type),            sizeof(m_type));
  }

  /// Write an isis style measure
  void ControlMeasure::write_isis( std::ostream &f ) {
    f << "    Group = ControlMeasure\n";
    f << "      SerialNumber   = " << m_serialNumber << std::endl;
    f << "      MeasureType    = ";
    if ( m_type == ControlMeasure::Unmeasured ) {
      f << "Unmeasured\n";
    } else if ( m_type == ControlMeasure::Manual ) {
      f << "Manual\n";
    } else if ( m_type == ControlMeasure::Estimated ) {
      f << "Estimated\n";
    } else if ( m_type == ControlMeasure::Automatic ) {
      f << "Automatic\n";
    } else if ( m_type == ControlMeasure::ValidatedManual ) {
      f << "ValidatedManual\n";
    } else if ( m_type == ControlMeasure::ValidatedAutomatic ) {
      f << "ValidatedAutomatic\n";
    } else {
      vw_throw( vw::NoImplErr() << "Invalid Control Measure type." );
    }
    if ( m_type == ControlMeasure::Unmeasured ) {
      // Why do you even exist?
      f << "      Sample         = Null\n";
      f << "      Line           = Null\n";
    } else {
      f << "      Sample         = " << m_col << "\n";
      f << "      Line           = " << m_row << "\n";
      f << "      ErrorLine      = " << m_col_sigma << "\n";
      f << "      ErrorSample    = " << m_row_sigma << "\n";
      f << "      ErrorMagnitude = " << sigma_magnitude() << "\n";
      f << "      FocalPlaneX    = " << m_focalplane_x << "\n";
      f << "      FocalPlaneY    = " << m_focalplane_y << "\n";
    }
    if ( m_ephemeris_time != 0 )
      f << "      EphemerisTime  = " << m_ephemeris_time << "\n";
    if ( m_diameter > 0 )
      f << "      Diameter       = " << m_diameter << "\n";
    if ( m_date_time != "" )
      f << "      DateTime       = " << m_date_time << "\n";
    if ( m_chooserName != "" )
      f << "      ChooserName    = " << m_chooserName << "\n";
    if ( m_ignore )
      f << "      Ignore         = True\n";
    f << "      Reference      = False\n";    // What is reference?
    if ( m_pixels_dominant )
      f << "      PixelsDominant = True\n";
    else
      f << "      PixelsDominant = False\n";
    f << "    End_Group\n";
  }

  void ControlMeasure::read_isis( std::istream &f ) {

    std::vector<std::string> tokens;
    std::ostringstream ostr;
    std::istringstream converter;
    std::string str;

    // Setting defaults
    m_diameter = 0;
    m_date_time = "";
    m_chooserName = "";
    m_ignore = false;
    m_pixels_dominant = true;

    while (1) {
      if ( f.eof() )
        vw_throw( vw::IOErr() << "Error reading Control Measure, unexpectly hit end of file" );

      // Reading file
      str = "";
      std::getline( f, str );
      boost::split( tokens, str, boost::is_any_of(" =\n") );

      // Cleaning out any tokens that are just ""
      for(std::vector<std::string>::iterator iter = tokens.begin();
          iter != tokens.end(); ++iter ) {
        if ( (*iter) == "" ) {
          iter = tokens.erase(iter);
          iter--;
        }
      }

      // Processing statement
      if ( tokens.empty() )
        continue;

      if ( tokens[0] == "End_Group" )     // End of Control Group
        break;
      else if ( tokens[0] == "SerialNumber" ) {
        read_pvl_property( ostr, tokens );
        m_serialNumber = ostr.str();
      } else if ( tokens[0] == "MeasureType" ) {
        read_pvl_property( ostr, tokens );
        if ( ostr.str() == "Unmeasured" )
          m_type = ControlMeasure::Unmeasured;
        else if ( ostr.str() == "Manual" )
          m_type = ControlMeasure::Manual;
        else if ( ostr.str() == "Estimated" )
          m_type = ControlMeasure::Estimated;
        else if ( ostr.str() == "Automatic" )
          m_type = ControlMeasure::Automatic;
        else if ( ostr.str() == "ValidatedManual" )
          m_type = ControlMeasure::ValidatedManual;
        else if ( ostr.str() == "ValidatedAutomatic" )
          m_type = ControlMeasure::ValidatedAutomatic;
        else
          vw_throw( vw::IOErr() << "Invalid Control Measure type, \""
                    << ostr.str() << "." );
      } else if ( tokens[0] == "Sample" ) {
        read_pvl_property( ostr, tokens );
        if ( ostr.str() == "Null" )
          m_col = 0;
        else {
          converter.str( ostr.str() );
          converter.clear();
          converter >> m_col;
        }
      } else if ( tokens[0] == "Line" ) {
        read_pvl_property( ostr, tokens );
        if ( ostr.str() == "Null" )
          m_row = 0;
        else {
          converter.str( ostr.str() );
          converter.clear();
          converter >> m_row;
        }
      } else if ( tokens[0] == "ErrorLine" ) {
        read_pvl_property( ostr, tokens );
        converter.str( ostr.str() );
        converter.clear();
        converter >> m_col_sigma;
      } else if ( tokens[0] == "ErrorSample" ) {
        read_pvl_property( ostr, tokens );
        converter.str( ostr.str() );
        converter.clear();
        converter >> m_row_sigma;
      } else if ( tokens[0] == "FocalPlaneX" ) {
        read_pvl_property( ostr, tokens );
        converter.str( ostr.str() );
        converter.clear();
        converter >> m_focalplane_x;
      } else if ( tokens[0] == "FocalPlaneY" ) {
        read_pvl_property( ostr, tokens );
        converter.str( ostr.str() );
        converter.clear();
        converter >> m_focalplane_y;
      } else if ( tokens[0] == "EphemerisTime" ) {
        read_pvl_property( ostr, tokens );
        converter.str( ostr.str() );
        converter.clear();
        converter >> m_ephemeris_time;
      } else if ( tokens[0] == "Diameter" ) {
        read_pvl_property( ostr, tokens );
        converter.str( ostr.str() );
        converter.clear();
        converter >> m_diameter;
      } else if ( tokens[0] == "DateTime" ) {
        read_pvl_property( ostr, tokens );
        m_date_time = ostr.str();
      } else if ( tokens[0] == "ChooserName" ) {
        read_pvl_property( ostr, tokens );
        m_chooserName = ostr.str();
      } else if ( tokens[0] == "Ignore" ) {
        m_ignore = true;
      } else if ( tokens[0] == "PixelsDominant" ) {
        boost::to_lower(tokens[1]);
        if ( tokens[1] == "false" )
          m_pixels_dominant = 0;
        else
          m_pixels_dominant = 1;
      }
      if ( m_col_sigma == 0 )
        m_col_sigma = 1;
      if ( m_row_sigma == 0 )
        m_row_sigma = 1;
    }
  }

  ////////////////////////////
  // Control Point          //
  ////////////////////////////

  /// Constructor
  ControlPoint::ControlPoint(ControlPointType type) : m_type(type) {
    m_ignore = false;
    m_id = "Null";
  }

  /// Add a measure
  void ControlPoint::add_measure(ControlMeasure const& measure) {
    m_measures.push_back(measure);
  }

  /// Add multiple measures
  void ControlPoint::add_measures(std::vector<ControlMeasure> const& measures) {
    m_measures.insert(m_measures.end(), measures.begin(), measures.end());
  }

  /// Remove the control point at the specified index.
  void ControlPoint::delete_measure(size_t index) {
    if (index >= this->size())
      vw_throw(LogicErr() << "ControlPoint::delete_control_point -- index " << index << " exceeds control point dimensions.");

    iterator iter = this->begin();
    iter += index;

    m_measures.erase(iter);
  }

  /// Find Control Measure
  size_t ControlPoint::find(ControlMeasure const& query) {
    for (size_t i = 0; i < m_measures.size(); ++i)
      if (m_measures[i] == query)
        return i;
    // If no matches are found, return m_measures.size() (the last index + 1)
    return m_measures.size();
  }

  /// Write a compressed binary style of point
  void ControlPoint::write_binary( std::ostream &f ) {
    // Writing out the string first
    f << m_id << char(0);
    // Writing the binary data
    f.write((char*)&(m_ignore), sizeof(m_ignore));
    f.write((char*)&(m_position[0]), sizeof(m_position[0]));
    f.write((char*)&(m_position[1]), sizeof(m_position[1]));
    f.write((char*)&(m_position[2]), sizeof(m_position[2]));
    f.write((char*)&(m_sigma[0]), sizeof(m_sigma[0]));
    f.write((char*)&(m_sigma[1]), sizeof(m_sigma[1]));
    f.write((char*)&(m_sigma[2]), sizeof(m_sigma[2]));
    f.write((char*)&(m_type), sizeof(m_type));
    int size = m_measures.size();
    f.write((char*)&(size), sizeof(size));
    // Rolling through the measures
    for ( int m = 0; m < size; m++ )
      m_measures[m].write_binary( f );
  }

  /// Reading a compressed binary style of point
  void ControlPoint::read_binary( std::istream &f ) {
    // Reading in the string first
    std::getline( f, m_id, '\0' );
    // Reading in the binary data
    f.read((char*)&(m_ignore),      sizeof (m_ignore));
    f.read((char*)&(m_position[0]), sizeof (m_position[0]));
    f.read((char*)&(m_position[1]), sizeof (m_position[1]));
    f.read((char*)&(m_position[2]), sizeof (m_position[2]));
    f.read((char*)&(m_sigma[0]),    sizeof (m_sigma[0]));
    f.read((char*)&(m_sigma[1]),    sizeof (m_sigma[1]));
    f.read((char*)&(m_sigma[2]),    sizeof (m_sigma[2]));
    f.read((char*)&(m_type),        sizeof (m_type));
    int size;
    f.read((char*)&(size),             sizeof (size));
    m_measures.clear();
    // Reading in all the measures
    for ( int m = 0; m < size; m++ ) {
      m_measures.push_back( ControlMeasure(f, FmtBinary) );
    }
  }

  /// Write an isis style point
  void ControlPoint::write_isis( std::ostream &f ) {
    f << "  Object = ControlPoint\n";
    f << "    PointType = ";
    if ( m_type == ControlPoint::GroundControlPoint ) {
      f << "Ground\n";
    } else if ( m_type == ControlPoint::TiePoint ) {
      f << "Tie\n";
    } else {
      vw_throw( vw::NoImplErr() << "Invalid Control Point type." );
    }
    f << "    PointId   = " << m_id << "\n";
    f << "    Latitude  = " << m_position[1] << "\n";
    f << "    Longitude = " << m_position[0] << "\n";
    f << "    Radius    = " << m_position[2] << "\n";
    if ( m_ignore )
      f << "    Ignore    = True\n";

    // Rolling through measures
    for ( unsigned m = 0; m < size(); m++ ) {
      f << std::endl;
      m_measures[m].write_isis( f );
    }

    f << "  End_Object\n";
  }

  /// Read an isis style point
  void ControlPoint::read_isis( std::istream &f ) {

    std::vector<std::string> tokens;
    std::ostringstream ostr;
    std::istringstream converter;
    std::string str;

    // Setting defaults
    m_ignore = false;

    m_measures.clear();

    while (1) {
      if ( f.eof() )
        vw_throw( vw::IOErr() << "Error reading Control Point, unexpectly hit end of file" );

      // Reading file
      str = "";
      std::getline( f, str );
      boost::split( tokens, str, boost::is_any_of(" =\n") );

      // Cleaning out any tokens that are just ""
      for(std::vector<std::string>::iterator iter = tokens.begin();
          iter != tokens.end(); ++iter ) {
        if ( (*iter) == "" ) {
          iter = tokens.erase(iter);
          iter--;
        }
      }

      // Processing statement
      if ( tokens.empty() )
        continue;

      if ( tokens[0] == "End_Object" )     // End of Control Point
        break;
      else if ( tokens[0] == "PointType" ) {
        read_pvl_property( ostr, tokens );
        if ( ostr.str() == "Ground" )
          m_type = ControlPoint::GroundControlPoint;
        else if ( ostr.str() == "Tie" )
          m_type = ControlPoint::TiePoint;
        else
          vw_throw( vw::IOErr() << "Invalid Control Point type, \""
                    << ostr.str() << "." );
      } else if ( tokens[0] == "PointId" ) {
        read_pvl_property( ostr, tokens );
        m_id = ostr.str();
      } else if ( tokens[0] == "Latitude" ) {
        read_pvl_property( ostr, tokens );
        converter.str( ostr.str() );
        converter.clear();
        converter >> m_position[1];
      } else if ( tokens[0] == "Longitude" ) {
        read_pvl_property( ostr, tokens );
        converter.str( ostr.str() );
        converter.clear();
        converter >> m_position[0];
      } else if ( tokens[0] == "Radius" ) {
        read_pvl_property( ostr, tokens );
        converter.str( ostr.str() );
        converter.clear();
        converter >> m_position[2];
      } else if ( tokens[0] == "Ignore" ) {
        m_ignore = tokens[1] == "True";
      } else if ( tokens[0] == "Group" ) {
        if ( tokens.size() == 1 ) {
          vw_throw( IOErr() << "Failed to read Control Point. Contains incorrect syntax, unlabelled Group" );
        } else if ( tokens[1] == "ControlMeasure" ) {
          m_measures.push_back( ControlMeasure( f, FmtIsisPvl ) );
        } else {
          vw_throw( IOErr() << "Failed to read Control Point. Unkown group \"" << tokens[1] << "\" found." );
        }
      }
    }
  }

  ////////////////////////////
  // Control Network        //
  ////////////////////////////

  /// Constructor
  ControlNetwork::ControlNetwork(std::string id,
                                 ControlNetworkType type,
                                 std::string target_name,
                                 std::string descrip, std::string user_name ) :
    m_targetName(target_name), m_networkId(id), m_description(descrip), m_userName(user_name), m_type(type) {
    // Recording time
    m_created = isis_style_time_string();
  }

  /// Add a single Control Point
  void ControlNetwork::add_control_point(ControlPoint const& point) {
    // Checking for GCPs (if this network supposedly doesn't contain
    // control networks)
    if ( m_type != ControlNetwork::ImageToGround )
      if ( point.type() == ControlPoint::GroundControlPoint )
        m_type = ControlNetwork::ImageToGround;

    m_control_points.push_back(point);
  }

  /// Add a vector of Control Points
  void ControlNetwork::add_control_points(std::vector<ControlPoint> const& points) {
    // Checking for GCPs (if this network supposedly doesn't contain
    // control networks)
    if ( m_type != ControlNetwork::ImageToGround )
      for ( unsigned p = 0; p < points.size(); p++ )
        if ( points[p].type() == ControlPoint::GroundControlPoint )
          m_type = ControlNetwork::ImageToGround;

    m_control_points.insert(m_control_points.end(), points.begin(), points.end());
  }

  // Delete control point
  void ControlNetwork::delete_control_point(size_t index) {
    if (index >= this->size())
      vw_throw(LogicErr() << "ControlNetwork::delete_control_point -- index " << index << " exceeds control network dimensions.");

    iterator iter = this->begin();
    iter += index;

    m_control_points.erase(iter);
  }

  // Find measure
  size_t ControlNetwork::find_measure(ControlMeasure const& query) {
    for (size_t i = 0; i < m_control_points.size(); ++i)
      if (m_control_points[i].find(query) != m_control_points[i].size())
        return i;
    // Otherwise...
    return m_control_points.size();
  }

  /// Write a compressed binary style control network
  void ControlNetwork::write_binary( std::string filename ) {

    // Recording the modified time
    m_modified = isis_style_time_string();

    // Forcing file extension type
    filename = filename.substr(0,filename.rfind("."));
    filename += ".cnet";

    // Opening file
    std::ofstream f( filename.c_str() );

    // Writing out the strings first
    f << m_targetName << char(0) << m_networkId << char(0)
      << m_created << char(0) << m_modified << char(0)
      << m_description << char(0) << m_userName << char(0);
    // Writing the binary data
    f.write((char*)&(m_type), sizeof(m_type));
    int size = m_control_points.size();
    f.write((char*)&(size), sizeof(size));
    // Rolling through the control points
    for ( int p = 0; p < size; p++ )
      m_control_points[p].write_binary( f );

    f.close();
  }

  /// Reading a compressed binary style control network
  void ControlNetwork::read_binary( std::string const& filename ) {

    // Opening file
    std::ifstream f( filename.c_str() );
    if ( !f.is_open() )
      vw_throw( IOErr() << "Failed to open \"" << filename << "\" as a Control Network." );

    // Reading in the strings first
    std::getline( f, m_targetName, '\0' );
    std::getline( f, m_networkId, '\0' );
    std::getline( f, m_created, '\0' );
    std::getline( f, m_modified, '\0' );
    std::getline( f, m_description, '\0' );
    std::getline( f, m_userName, '\0' );

    // Reading in the binary data
    f.read((char*)&(m_type), sizeof(m_type));
    int size;
    f.read((char*)&(size), sizeof(size));

    // Clearing anything left in this control network
    m_control_points.clear();

    // Reading in all the control points
    for ( int p = 0; p < size; p++ ) {
      m_control_points.push_back( ControlPoint( f, FmtBinary ) );
    }

    f.close();
  }

  /// Write an isis style control network
  void ControlNetwork::write_isis( std::string filename ) {
    // Recording the modified time
    m_modified = isis_style_time_string();

    // Forcing file extension type
    filename = filename.substr(0,filename.rfind("."));
    filename += ".net";

    // Making sure all the control points have unique IDs
    for ( unsigned p = 0; p < m_control_points.size(); p++ ) {
      if (m_control_points[p].id() == "Null" ||
          m_control_points[p].id() == "" ) {
        std::ostringstream ostr;
        ostr << std::setw(9) << std::setfill('0') << p;
        m_control_points[p].set_id( ostr.str() );
      }
    }

    // Opening file
    std::ofstream f( filename.c_str() );
    f << std::setprecision( 15 );

    f << "Object = ControlNetwork\n";
    f << "  NetworkId    = " << m_networkId << "\n";
    f << "  NetworkType  = ";
    if ( m_type == ControlNetwork::Singleton ) {
      f << "Singleton\n";
    } else if ( m_type == ControlNetwork::ImageToImage ) {
      f << "ImageToImage\n";
    } else if ( m_type == ControlNetwork::ImageToGround ) {
      f << "ImageToGround\n";
    } else {
      vw_throw( vw::NoImplErr() << "Invalid Control Network type." );
    }
    f << "  TargetName   = " << m_targetName << "\n";
    f << "  UserName     = " << m_userName << "\n";
    f << "  Created      = " << m_created << "\n";
    f << "  LastModified = " << m_modified << "\n";
    f << "  Description  = " << m_description << "\n";

    // Rolling through the control points
    for ( unsigned m = 0; m < size(); m++ ) {
      f << std::endl;
      m_control_points[m].write_isis( f );
    }

    f << "End_Object\nEnd\n";
    f.close();
  }

  /// Read an isis style control network
  void ControlNetwork::read_isis( std::string const& filename ) {

    // Opening file
    std::ifstream f( filename.c_str() );
    if ( !f.is_open() )
      vw_throw( IOErr() << "Failed to open \"" << filename << "\" as a ISIS Control Network." );

    // Clearing anything left in this control network
    m_control_points.clear();

    // Reading file
    std::vector<std::string> tokens;
    std::ostringstream ostr;
    std::string str;
    while ( !f.eof() ) {
      // Reading file
      str = "";
      std::getline( f, str);
      boost::split( tokens, str, boost::is_any_of(" =\n") );

      // Cleaning out any tokens that are just ""
      for(std::vector<std::string>::iterator iter = tokens.begin();
          iter != tokens.end(); ++iter ) {
        if ( (*iter) == "" ) {
          iter = tokens.erase(iter);
          iter--;
        }
      }

      // Processing statement
      if ( tokens.empty() )
        continue;

      if ( tokens[0] == "End" )        // End of file
        break;
      else if ( tokens[0] == "NetworkId" ) {
        read_pvl_property( ostr, tokens );
        m_networkId = ostr.str();
      }
      else if ( tokens[0] == "NetworkType" ) {
        read_pvl_property( ostr, tokens );
        if ( ostr.str() == "Singleton" )
          m_type = ControlNetwork::Singleton;
        else if ( ostr.str() == "ImageToImage" )
          m_type = ControlNetwork::ImageToImage;
        else if ( ostr.str() == "ImageToGround" )
          m_type = ControlNetwork::ImageToGround;
        else
          vw_throw( vw::IOErr() << "Invalid Control Network type, \""
                    << ostr.str() << "." );
      }
      else if ( tokens[0] == "TargetName" ) {
        read_pvl_property( ostr, tokens );
        m_targetName = ostr.str();
      }
      else if ( tokens[0] == "UserName" ) {
        read_pvl_property( ostr, tokens );
        m_userName = ostr.str();
      }
      else if ( tokens[0] == "Created" ) {
        read_pvl_property( ostr, tokens );
        m_created = ostr.str();
      }
      else if ( tokens[0] == "LastModified" ) {
        read_pvl_property( ostr, tokens );
        m_modified = ostr.str();
      }
      else if ( tokens[0] == "Description" ) {
        read_pvl_property( ostr, tokens );
        m_description = ostr.str();
      }
      else if ( tokens[0] == "Object" ) {
        if ( tokens.size() == 1 ) {
          vw_throw( IOErr() << "Failed to open \"" << filename
                    << "\". Contains incorrect syntax, unlabelled Object" );
        } else if ( tokens[1] == "ControlNetwork" ) {
          // Start of file
          continue;
        } else if ( tokens[1] == "ControlPoint" ) {
          m_control_points.push_back( ControlPoint( f, FmtIsisPvl  ) );
        } else if ( tokens[1] == "ControlMeasure" ) {
          vw_throw( IOErr() << "Failed to open \"" << filename
                    << "\". Control Measure found out of order." );
        } else {
          vw_throw( IOErr() << "Failed to open \"" << filename
                    << "\". Unknown Object \"" << tokens[1] << "\" found." );
        }
      }
    }
    f.close();
  }

  ////////////////////////////
  // Generic Ostream        //
  ////////////////////////////

  std::ostream& operator<<( std::ostream& os, ControlMeasure const& measure) {
    os << measure.image_id() << ":" << measure.position();
    return os;
  }

  std::ostream& operator<<( std::ostream& os, ControlPoint const& point) {
    os << "[Control Point: " << point.position() << "] ";
    for (unsigned m=0; m<point.size(); ++m)
      os << point[m] << " ";
    os << "\n";
    return os;
  }

  std::ostream& operator<<( std::ostream& os, ControlNetwork const& cnet) {
    os << "Control Network: " << cnet.size() << " points.\n";
    for (unsigned p=0; p<cnet.size(); ++p)
      os << "\t" << cnet[p];
    os << "\n";
    return os;
  }

  // Read a single pvl propert
  void read_pvl_property( std::ostringstream& ostr,
                          std::vector<std::string>& tokens ) {
    ostr.str("");
    for ( unsigned i = 1; i < tokens.size(); i++ )
      if ( i > 1 )
        ostr << " " << tokens[i];
      else
        ostr << tokens[i];
  }


}} // namespace vw::camera
