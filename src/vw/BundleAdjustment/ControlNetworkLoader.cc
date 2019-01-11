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


#include <vw/BundleAdjustment/ControlNetworkLoader.h>
#include <vw/Stereo/StereoModel.h>
#include <vw/InterestPoint/Matcher.h>

using namespace vw;
using namespace vw::ba;

#include <boost/filesystem/fstream.hpp>

namespace fs = boost::filesystem;

struct ContainsEqualIP {
  ip::InterestPoint& m_compare;
  ContainsEqualIP( ip::InterestPoint& comp ) : m_compare(comp) {}

  bool operator()( boost::shared_ptr<IPFeature> in ) {
    if (  m_compare.x == in->m_ip.x &&
          m_compare.y == in->m_ip.y )
      return true;
    return false;
  }
};

// Utility for checking that the point is BA safe
void safe_measurement( ip::InterestPoint& ip ) {
  if ( ip.scale <= 0 ) ip.scale = 10;
}

double vw::ba::triangulate_control_point( ControlPoint& cp,
                                          std::vector<boost::shared_ptr<camera::CameraModel> >
                                          const& camera_models,
                                          double min_angle_radians,
                                          double forced_triangulation_distance) {

  Vector3 position_sum;
  double error = 0, error_sum = 0;
  size_t count = 0;

  // 4.1.) Building a listing of triangulation
  for ( size_t j = 0, k = 1; k < cp.size(); j++, k++ ) {
    // Make sure camera centers are not equal
    size_t j_cam_id = cp[j].image_id();
    size_t k_cam_id = cp[k].image_id();
    if ( norm_2( camera_models[j_cam_id]->camera_center( cp[j].position() ) -
                 camera_models[k_cam_id]->camera_center( cp[k].position() ) ) > 1e-6 ) {
      try {

        double angle_tol = stereo::StereoModel::robust_1_minus_cos(min_angle_radians);

        bool least_squares = false;
        stereo::StereoModel sm( camera_models[ j_cam_id ].get(),
                                camera_models[ k_cam_id ].get(), least_squares,
                                angle_tol );
        
        Vector3 pt = sm( cp[j].position(), cp[k].position(), error );
        // TODO: When forced_triangulation_distance > 0, one can check
        // if the triangulated point is behind the camera, and if yes,
        // to replace it with an artificial point in front of the camera.
        // This will need a good test.
        if (pt != Vector3() ){
          count++;
          position_sum += pt;
          error_sum += error;
        }else if (forced_triangulation_distance <= 0){
          vw_out(WarningMessage,"ba") << "\nCould not triangulate point. If too many such errors, "
                                      << "perhaps your baseline is too small, "
                                      << "or consider decreasing --min-triangulation-angle "
                                      << "or using --forced-triangulation-distance.\n";
        }
      } catch ( std::exception const& e) {
        /* Just let it go */
        vw_out(WarningMessage,"ba") << "\nFailure in triangulation: " << e.what();
      }
    }
  }

  // 4.2.) Summing, averaging, and storing
  if ( count == 0) {
    // vw_out(WarningMessage,"ba") << "\nUnable to triangulate point!\n";
    // At the very least we can provide a point that is some
    // distance out from the camera center and is in the 'general' area.
    size_t j = cp[0].image_id();
    double dist = 10.0; // for backward compatibility
    if (forced_triangulation_distance > 0) 
      dist = forced_triangulation_distance;
    try {
      cp.set_position( camera_models[j]->camera_center(cp[0].position()) +
                       camera_models[j]->pixel_to_vector(cp[0].position())*dist );
    } catch ( const camera::PixelToRayErr& ) {
      cp.set_position( camera_models[j]->camera_center(cp[0].position()) +
                       camera_models[j]->camera_pose(cp[0].position()).rotate(Vector3(0,0,dist)) );
    }
    if (forced_triangulation_distance > 0)
      return 1; // mark triangulation as successful
    return 0;
  } else {
    error_sum /= double(count);
    cp.set_position( position_sum / double(count) );

    return error_sum;
  }
}

bool vw::ba::build_control_network( bool triangulate_control_points,
                                    ba::ControlNetwork& cnet,
                                    std::vector<boost::shared_ptr<camera::CameraModel> >
                                    const& camera_models,
                                    std::vector<std::string> const& image_files,
                                    std::map< std::pair<int, int>, std::string> const& match_files,
                                    size_t min_matches,
                                    double min_angle_radians,
                                    double forced_triangulation_distance) {

  // Note that this statement does not clear the network fully.
  // TODO: Clear all items here. 
  cnet.clear();
  
  // We can't guarantee that image_files is sorted, so we make a
  // std::map to give ourselves a sorted list and access to a binary search.
  std::map<std::string,size_t> image_prefix_map;
  size_t count = 0;
  ba::CameraRelationNetwork<ba::IPFeature> crn;
  BOOST_FOREACH( std::string const& file, image_files ) {
    fs::path file_path(file);
    image_prefix_map[file_path.replace_extension().string()] = count;
    crn.add_node( ba::CameraNode<ba::IPFeature>( count,
                                                 file_path.stem().string() ) );
    cnet.add_image_name(file);
    count++;
  }

  // Look for match files starting with given prefix.
  std::vector<std::string> match_files_vec;
  std::vector<size_t> index1_vec, index2_vec;

  // Searching through the directories available to us.
  typedef std::map<std::string,size_t>::iterator MapIterator;
  int num_images = image_files.size();
  for (int i = 0; i < num_images; i++){ // Loop through all image pairs
    for (int j = i+1; j < num_images; j++){
      std::string image1 = image_files[i];
      std::string image2 = image_files[j];

      // Find the match file for this pair of images
      std::pair<int, int> pair_ind(i, j);
      std::map< std::pair<int, int>, std::string>::const_iterator pair_it
        = match_files.find(pair_ind);
      if (pair_it == match_files.end())
        continue; // This match file was not passed in, that is ok.
      std::string match_file = pair_it->second;

      // 
      std::string prefix1 = fs::path(image1).replace_extension().string();
      std::string prefix2 = fs::path(image2).replace_extension().string();
      MapIterator it1 = image_prefix_map.find( prefix1 );
      MapIterator it2 = image_prefix_map.find( prefix2 );
      if ( it1 == image_prefix_map.end() ||
           it2 == image_prefix_map.end() ) continue;
           
      // Verify that the match file exists
      if (!fs::exists(match_file)) {
        vw_out(WarningMessage) << "Missing match file: " << match_file << std::endl;
        continue;
      }
      match_files_vec.push_back(match_file);
      index1_vec.push_back(it1->second);
      index2_vec.push_back(it2->second);
    }
  }

  // Loop through the match files...
  size_t num_load_rejected = 0, num_loaded = 0;
  for (size_t file_iter = 0; file_iter < match_files_vec.size(); file_iter++){
    std::string match_file = match_files_vec[file_iter];
    size_t index1 = index1_vec[file_iter];
    size_t index2 = index2_vec[file_iter];

    // Actually read in the file as it seems we've found something correct
    std::vector<ip::InterestPoint> ip1, ip2;
    vw_out(DebugMessage,"ba") << "Loading: " << match_file << std::endl;
    ip::read_binary_match_file( match_file, ip1, ip2 );
    if ( ip1.size() < min_matches ) {
      vw_out(DebugMessage,"ba") << "\t" << match_file << "    "
                                << ip1.size() << " matches. [rejected]\n";
      num_load_rejected += ip1.size();
      continue;
    }
    vw_out(DebugMessage,"ba") << "\t" << match_file << "    "
                              << ip1.size() << " matches.\n";
    num_loaded += ip1.size();

    // Remove descriptors from interest points and correct scale
    std::for_each( ip1.begin(), ip1.end(), ip::remove_descriptor );
    std::for_each( ip2.begin(), ip2.end(), ip::remove_descriptor );
    std::for_each( ip1.begin(), ip1.end(), safe_measurement );
    std::for_each( ip2.begin(), ip2.end(), safe_measurement );

    typedef boost::shared_ptr< ba::IPFeature > f_ptr;
    typedef std::list< f_ptr >::iterator f_itr;

    // Checking to see if features already exist, adding if they
    // don't, then linking them.
    vw_out() << "Building the control network for " << match_file <<".\n";
    TerminalProgressCallback progress("ba", "Building: ");
    progress.report_progress(0);
    double inc_prog = 1.0/double(ip1.size());
    for ( size_t k = 0; k < ip1.size(); k++ ) { // Loop through all existing IP
      // Check if either IP is already in the lists
      f_itr ipfeature1 = std::find_if( crn[index1].begin(),
                                       crn[index1].end(),
                                       ContainsEqualIP( ip1[k] ) );
      f_itr ipfeature2 = std::find_if( crn[index2].begin(),
                                       crn[index2].end(),
                                       ContainsEqualIP( ip2[k] ) );
      // If the IP are new, add them to the list
      if ( ipfeature1 == crn[index1].end() ) {
        crn[index1].relations.push_front( f_ptr( new ba::IPFeature( ip1[k], index1 ) ) );
        ipfeature1 = crn[index1].begin();
      }
      if ( ipfeature2 == crn[index2].end() ) {
        crn[index2].relations.push_front( f_ptr( new ba::IPFeature( ip2[k], index2 ) ) );
        ipfeature2 = crn[index2].begin();
      }

      // Doubly linking
      (*ipfeature1)->connection( *ipfeature2, false );
      (*ipfeature2)->connection( *ipfeature1, false );
      progress.report_incremental_progress(inc_prog );
    } // End loop through ip1
    progress.report_finished();
  } // End loop through match files

  if ( num_load_rejected != 0 ) {
    vw_out(WarningMessage,"ba") << "\tDidn't load " << num_load_rejected
                                << " matches due to inadequacy. Decrease the"
                                << " --min-matches parameter to load smaller "
                                << "sets of matches.\n";
    vw_out(WarningMessage,"ba") << "\tLoaded " << num_loaded << " matches.\n";
  }

  // Building control network
  bool success = crn.write_controlnetwork( cnet );

  // Triangulating Positions
  if (triangulate_control_points){
    TerminalProgressCallback progress("ba", "Triangulating: ");
    progress.report_progress(0);
    double inc_prog = 1.0/double(cnet.size());
    BOOST_FOREACH( ba::ControlPoint& cpoint, cnet ) {
      progress.report_incremental_progress(inc_prog );
      ba::triangulate_control_point( cpoint, camera_models, min_angle_radians,
                                     forced_triangulation_distance);
    }
    progress.report_finished();
  }
  return success;
}

void vw::ba::add_ground_control_points(vw::ba::ControlNetwork& cnet,
                                       std::vector<std::string> const& gcp_files,
                                       cartography::Datum const& datum){
  
  namespace fs = boost::filesystem;
  
  std::vector<std::string> const& image_files = cnet.get_image_list();

  // Creating a version of image_files that doesn't contain the path
  typedef std::map<std::string,size_t> LookupType;
  LookupType image_lookup;
  for (size_t i = 0; i < image_files.size(); i++ ) {
    image_lookup[image_files[i]] = i;
    image_lookup[fs::path(image_files[i]).filename().string()] = i;
  }

  std::vector<std::string>::const_iterator gcp_iter = gcp_files.begin();
  std::vector<std::string>::const_iterator gcp_end  = gcp_files.end();

  while ( gcp_iter != gcp_end ) {

    vw_out() << "Loading GCP file: " << *gcp_iter << std::endl;

    if ( !fs::exists( *gcp_iter ) ) {
      vw_throw( ArgumentErr() << "GCP file " << *gcp_iter << " does not exist!" );
    }

    std::ifstream ifile( (*gcp_iter).c_str() );
    std::string line;
    while ( getline(ifile, line, '\n') ){
      // Skip empty lines or lines starting with comments
      if (line.size() == 0) continue;
      if (line.size() > 0 && line[0] == '#') continue;

      boost::replace_all(line, ",", " ");

      // Data to be loaded
      std::vector<Vector4> measure_locations;
      std::vector<std::string> measure_cameras;
      int point_id;
      Vector3 world_location, world_sigma;

      std::istringstream is(line);

      // First elements in the line are the point id, location in
      // the world, and its sigmas
      if (!(is >> point_id >> world_location[0] >> world_location[1]
            >> world_location[2] >> world_sigma[0]
            >> world_sigma[1] >> world_sigma[2])){
        vw_out(WarningMessage) << "Could not parse a ground control point "
                               << "from line: " << line << std::endl;
        continue;
      }

      // Other elements in the line define the position in images
      while(1){
        std::string temp_name;
        Vector4 temp_loc;
        if (is >> temp_name >> temp_loc[0] >> temp_loc[1]
            >> temp_loc[2] >> temp_loc[3]){
          if (temp_loc[2] <= 0 || temp_loc[3] <= 0) {
            vw_throw( ArgumentErr() << "Standard deviations must be positive "
                                    << "when loading ground control points." );
          }
          measure_locations.push_back( temp_loc );
          measure_cameras.push_back( temp_name );
        }else{
                break;
        }
      }

      if (world_sigma[0] <= 0 || world_sigma[1] <= 0 || world_sigma[2] <= 0)
        vw_throw( ArgumentErr() << "Standard deviations must be positive "
                                            << "when loading ground control points." );

      // Make lat,lon into lon,lat
      std::swap(world_location[0], world_location[1]);

      // Building Control Point
      Vector3 xyz = datum.geodetic_to_cartesian(world_location);

      vw_out(VerboseDebugMessage,"ba") << "\t\tLocation: " << xyz << std::endl;
      ControlPoint cpoint(ControlPoint::GroundControlPoint);
      cpoint.set_position(xyz[0],         xyz[1],         xyz[2]        );
      cpoint.set_sigma   (world_sigma[0], world_sigma[1], world_sigma[2]);

      // Adding measures
      std::vector<Vector4    >::iterator m_iter_loc  = measure_locations.begin();
      std::vector<std::string>::iterator m_iter_name = measure_cameras.begin();
      while ( m_iter_loc != measure_locations.end() ) {
        LookupType::iterator it = image_lookup.find(*m_iter_name);
        if ( it != image_lookup.end() ) {
          vw_out(DebugMessage,"ba") << "\t\tAdded Measure: " << *m_iter_name
                                    << " #" << it->second << std::endl;
          ControlMeasure cm( (*m_iter_loc)[0], (*m_iter_loc)[1],
                             (*m_iter_loc)[2], (*m_iter_loc)[3], it->second );
          cpoint.add_measure( cm );
        } else {
          vw_out(WarningMessage,"ba") << "No input image found matching "
                                      << *m_iter_name << std::endl;
        }
        m_iter_loc++;
        m_iter_name++;
      }

      // Append the GCP
      if (cpoint.size() > 0) 
        cnet.add_control_point(cpoint);
      
    } // End line loop
    ifile.close();

    gcp_iter++;
  } // End file loop
}
