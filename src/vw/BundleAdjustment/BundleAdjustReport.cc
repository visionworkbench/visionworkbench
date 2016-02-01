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


#include <vw/BundleAdjustment/BundleAdjustReport.h>
namespace fs = boost::filesystem;

namespace vw {
namespace ba {

  void write_kml_styles( KMLFile& kml ) {
    // GCP Placemark Style
    kml.append_style( "gcp_circle", "", 1.2,
                      "http://maps.google.com/mapfiles/kml/shapes/placemark_circle.png" );
    kml.append_style( "gcp_circle_highlight", "", 1.4,
                      "http://maps.google.com/mapfiles/kml/shapes/placemark_circle_highlight.png" );
    kml.append_stylemap( "gcp_placemark", "gcp_circle",
                         "gcp_circle_highlight" );

    // Est Circle Lvl 1 (Green) Placemark
    kml.append_style( "est_circle_1", "ff00ff00", 0.8,
                      "http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png" );
    kml.append_style( "est_circle_1_highlight", "ff00ff00", 0.9,
                      "http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png" );
    kml.append_stylemap( "est_placemark_1", "est_circle_1",
                         "est_circle_1_highlight" );

    // Est Circle Lvl 2 (Green-Yellow) Placemark
    kml.append_style( "est_circle_2", "ff00ff80", 0.8,
                      "http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png" );
    kml.append_style( "est_circle_2_highlight", "ff00ff80", 0.9,
                      "http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png" );
    kml.append_stylemap( "est_placemark_2", "est_circle_2",
                         "est_circle_2_highlight" );

    // Est Circle Lvl 3 (Yellow) Placemark
    kml.append_style( "est_circle_3", "ff00ffff", 0.8,
                      "http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png" );
    kml.append_style( "est_circle_3_highlight", "ff00ffff", 0.9,
                      "http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png" );
    kml.append_stylemap( "est_placemark_3", "est_circle_3",
                         "est_circle_3_highlight" );

    // Est Circle Lvl 4 (Red-Yellow) Placemark
    kml.append_style( "est_circle_4", "ff0080ff", 0.8,
                      "http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png" );
    kml.append_style( "est_circle_4_highlight", "ff0080ff", 0.9,
                      "http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png" );
    kml.append_stylemap( "est_placemark_4", "est_circle_4",
                         "est_circle_4_highlight" );

    // Est Circle Lvl 5 (Red) Placemark
    kml.append_style( "est_circle_5", "ff0000ff", 0.8,
                      "http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png" );
    kml.append_style( "est_circle_5_highlight", "ff0000ff", 0.9,
                      "http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png" );
    kml.append_stylemap( "est_placemark_5", "est_circle_5",
                         "est_circle_5_highlight" );
  }

  void write_gcps_kml( KMLFile& kml,
                       ControlNetwork const& cnet ) {
    kml.enter_folder( "Ground Control Points",
                      "Used for Bundle Adjustment in VW" );

    unsigned count = 0;
    for ( ControlNetwork::const_iterator iter = cnet.begin();
          iter != cnet.end(); ++iter ) {
      if ( (*iter).type() == ControlPoint::GroundControlPoint ) {
        count++;
        Vector3 llr = cartography::xyz_to_lon_lat_radius_estimate( (*iter).position() ); // TODO: Not accurate!

        std::ostringstream desc;
        // GCP data
        desc << "&lt;h2&gt;Ground Control Point&lt;/h2&gt;";
        desc << "&lt;b&gt;Lon:&lt;/b&gt; " << llr.x() << " deg&lt;br&gt;";
        desc << "&lt;b&gt;Lat:&lt;/b&gt; " << llr.y() << " deg&lt;br&gt;";
        desc << "&lt;b&gt;Rad:&lt;/b&gt; " << std::setprecision(12)
             << llr.z() << " m&lt;br&gt;";

        // Images viewing
        desc << "&lt;h3&gt;Viewed by:&lt;/h3&gt;&lt;ol&gt;";
        for ( ControlPoint::const_iterator measure = (*iter).begin();
              measure != (*iter).end(); ++measure ) {
          desc << "&lt;li&gt;" << "[" << (*measure).image_id() << "] "
               << (*measure).serial() << " ( " << measure->position()[0] << ", "
               << measure->position()[1] <<  " ) " << "&lt;/li&gt;";
        }
        desc << "&lt;/ol&gt;";

        if ( (*iter).id() != "Null" ) {
          kml.append_placemark( llr.x(), llr.y(),
                                (*iter).id(), desc.str(),
                                "gcp_placemark" );
        } else {
          std::ostringstream gcp_name;
          gcp_name << "GCP " << count;
          kml.append_placemark( llr.x(), llr.y(),
                                gcp_name.str(), desc.str(),
                                "gcp_placemark" );
        }
      }
    }

    kml.exit_folder();
  }

  void write_3d_est_kml( KMLFile& kml,
                         ControlNetwork const& cnet,
                         std::vector<double>& image_errors ) {
    kml.enter_folder( "3D Point estimates",
                      "Used for Bundle Adjustment in VW" );

    std::list<Vector3> points; // Lon, Lat, Interest

    unsigned index = 0;
    for ( ControlNetwork::const_iterator iter = cnet.begin();
          iter != cnet.end(); ++iter ) {
      Vector3 llr = cartography::xyz_to_lon_lat_radius_estimate( (*iter).position() ); // TODO: Not accurate!
      double mean_image_error = 0;
      int count_measures = 0;
      for ( ControlPoint::const_iterator m_iter = (*iter).begin();
            m_iter != (*iter).end(); ++m_iter ) {
        mean_image_error += image_errors[index];
        index++;
        count_measures++;
      }
      mean_image_error /= count_measures;
      if ( (*iter).type() == ControlPoint::TiePoint )
        points.push_back( Vector3(llr.x(),llr.y(), mean_image_error) );
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

    // Finding the maximium and minimium error
    double min = 1e20, max = -1;
    for ( std::list<Vector3>::iterator it = points.begin();
          it != points.end(); ++it ) {
      if ( (*it).z() < min )
        min = (*it).z();
      if ( (*it).z() > max )
        max = (*it).z();
    }

    recursive_kml_placemark( kml, points, kml.name(), min, max,
                             upper_corner.y(), lower_corner.y(),
                             upper_corner.x(), lower_corner.x(), 0 );

    kml.exit_folder();
  }

  void recursive_kml_placemark( KMLFile& kml,
                                std::list<Vector3>& list,
                                std::string const& name,
                                double& min, double& max,
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
    double diff = max - min;

    // Checking list
    int count = 0;
    for ( std::list<Vector3>::iterator it = list.begin();
          it != list.end(); ++it )
      count++;
    if ( count <= 500 ) {
      // Write a termination file
      kml.enter_folder("","");

      // Regioning
      kml.open_bracket("Region");
      kml.append_latlonaltbox( north, south, east, west );
      kml.append_lod( 512, -1 );
      kml.close_bracket();

      // Placemarks
      for ( std::list<Vector3>::iterator it = list.begin();
            it != list.end(); ++it ) {
        std::ostringstream desc;
        desc << "Image error: " << (*it).z();
        if ( (*it).z() > 4*diff/5 + min )
          kml.append_placemark( (*it).x(), (*it).y(),
                                "", desc.str(),
                                "est_placemark_5" );
        else if ( (*it).z() > 3*diff/5 + min )
          kml.append_placemark( (*it).x(), (*it).y(),
                                "", desc.str(),
                                "est_placemark_4" );
        else if ( (*it).z() > 2*diff/5 + min )
          kml.append_placemark( (*it).x(), (*it).y(),
                                "", desc.str(),
                                "est_placemark_3" );
        else if ( (*it).z() > 1*diff/5 + min )
          kml.append_placemark( (*it).x(), (*it).y(),
                                "", desc.str(),
                                "est_placemark_2" );
        else
          kml.append_placemark( (*it).x(), (*it).y(),
                                "", desc.str(),
                                "est_placemark_1" );
      }

      kml.exit_folder();
      return;
    } else {
      // Write a branching file
      list.sort( vector_sorting );

      kml.enter_folder("","");

      for ( char i = 0; i < 4; i++ ) {
        std::ostringstream link;
        if ( recursive_lvl == 0 )
          link << "data/";
        link << name << int(i) << ".kml";
        kml.append_network( link.str(),
                            north_dv[i], south_dv[i],
                            east_dv[i], west_dv[i] );
      }

      kml.enter_folder("","");

      // Regioning
      kml.open_bracket("Region");
      kml.append_latlonaltbox( north, south, east, west );
      kml.append_lod( 512, -1 );
      kml.close_bracket();

      // Placemarks
      count = 500;
      std::list<Vector3>::iterator it = list.begin();
      while ( count > 0 ) {
        std::ostringstream desc;
        desc << "Image error: " << (*it).z();
        if ( (*it).z() > 4*diff/5 + min )
          kml.append_placemark( (*it).x(), (*it).y(),
                                "", desc.str(),
                                "est_placemark_5" );
        else if ( (*it).z() > 3*diff/5 + min )
          kml.append_placemark( (*it).x(), (*it).y(),
                                "", desc.str(),
                                "est_placemark_4" );
        else if ( (*it).z() > 2*diff/5 + min )
          kml.append_placemark( (*it).x(), (*it).y(),
                                "", desc.str(),
                                "est_placemark_3" );
        else if ( (*it).z() > 1*diff/5 + min )
          kml.append_placemark( (*it).x(), (*it).y(),
                                "", desc.str(),
                                "est_placemark_2" );
        else
          kml.append_placemark( (*it).x(), (*it).y(),
                                "", desc.str(),
                                "est_placemark_1" );
        it = list.erase( it );
        count--;
      }

      kml.exit_folder();
      kml.exit_folder();
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
        if ( kml.directory() != "" )
          dir << kml.directory() << "/";
        if ( recursive_lvl == 0 )
          dir << "data/";
        KMLFile subsection( new_name.str() + ".kml",
                            new_name.str(), dir.str() );
        write_kml_styles( subsection );
        recursive_kml_placemark( subsection,
                                 temp, new_name.str(),
                                 min, max,
                                 north_dv[i], south_dv[i],
                                 east_dv[i], west_dv[i],
                                 recursive_lvl+1 );
      }
    }

    if ( !list.empty() )
      std::cout << "Error! Vector not empty!\n";

  }

  bool vector_sorting( Vector3 i, Vector3 j) {
    return (i.z() > j.z());
  }

}}
