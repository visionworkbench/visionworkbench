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
#include <vw/BundleAdjustment/openMVG_tracks.h>
#include <vw/Stereo/StereoModel.h>
#include <vw/InterestPoint/Matcher.h>
#include <vw/Math/RandomSet.h>
#include <vw/Core/Stopwatch.h>

#include <mutex>
#include <boost/filesystem/fstream.hpp>

using namespace vw;
using namespace vw::ba;
namespace fs = boost::filesystem;

const int MAX_TRI_FAILURE_WARNINGS = 100;

struct ContainsEqualIP {
  ip::InterestPoint& m_compare;
  ContainsEqualIP(ip::InterestPoint& comp) : m_compare(comp) {}

  bool operator()(boost::shared_ptr<IPFeature> in) {
    if ( m_compare.x == in->m_ip.x &&
          m_compare.y == in->m_ip.y)
      return true;
    return false;
  }
};

// Utility for checking that the point is BA safe
void safe_measurement(ip::InterestPoint& ip) {
  if (ip.scale <= 0) ip.scale = 10;
}

double vw::ba::triangulate_control_point(ControlPoint& cp,
                                         std::vector<boost::shared_ptr<camera::CameraModel>>
                                         const& camera_models,
                                         double min_angle_radians,
                                         double forced_triangulation_distance) {
  
  Vector3 position_sum(0.0, 0.0, 0.0);
  double error = 0, error_sum = 0;
  size_t count = 0;
  
  double angle_tol = stereo::StereoModel::robust_1_minus_cos(min_angle_radians);

  // Ensure this is initialized
  cp.set_position((Vector3(0, 0, 0)));
  
  // Do pairwise triangulation then average the results. Note that as
  // long as at least two of the rays meet with a triangulation angle
  // no less than min_angle_radians, triangulation will succeed, and
  // these successful triangulation points will be averaged. Hence,
  // while some of the rays can meet at a very small angle, their
  // intersection, being unreliable, will not be added to the
  // mix. This allows, for example, for a triplet of rays, R1, R2, R3,
  // all corresponding to matching features and converging to the same
  // triangulated point, so that R1 and R2 have a very small angle,
  // but R1 and R3, then R2 and R3 have a very solid angle. This
  // triangulated point should be robust enough if the min
  // triangulation angle is, say, no less than 5 degrees, and ideally
  // 15-20 degrees or more.
  for (size_t j = 0, k = 1; k < cp.size(); j++, k++) {
    // Make sure camera centers are not equal
    size_t j_cam_id = cp[j].image_id();
    size_t k_cam_id = cp[k].image_id();
    if (norm_2(camera_models[j_cam_id]->camera_center(cp[j].position()) -
                 camera_models[k_cam_id]->camera_center(cp[k].position())) > 1e-6) {
      try {

        bool least_squares = false;
        stereo::StereoModel sm(camera_models[j_cam_id].get(), camera_models[k_cam_id].get(),
                               least_squares, angle_tol); // use angle_tol0
        
        Vector3 pt = sm(cp[j].position(), cp[k].position(), error);
        // TODO: When forced_triangulation_distance > 0, one can check
        // if the triangulated point is behind the camera, and if yes,
        // to replace it with an artificial point in front of the camera.
        // This will need a good test.
        if (pt != Vector3()) {
          count++;
          position_sum += pt;
          error_sum += error;
        }

      } catch (std::exception const& e) {
        // Just let it go
      }
    }
  }

  // Summing, averaging, and storing
  if (count == 0) {
    // It failed to triangulate. At the very least we can provide a point that is some
    // distance out from the camera center and is in the 'general' area.
    size_t j = cp[0].image_id();
    double dist = 10.0; // for backward compatibility
    if (forced_triangulation_distance > 0) 
      dist = forced_triangulation_distance;
    try {
      cp.set_position(camera_models[j]->camera_center(cp[0].position()) +
                       camera_models[j]->pixel_to_vector(cp[0].position())*dist);
    } catch (...) {
      try {
        cp.set_position(camera_models[j]->camera_center(cp[0].position()) +
                      camera_models[j]->camera_pose(cp[0].position()).rotate(Vector3(0,0,dist)));
      } catch (...) {
        return -1; // nothing can be done
      }
    }
    
    if (forced_triangulation_distance > 0)
      return 1; // mark triangulation as successful
    
    return -1;
  }
  
  error_sum /= double(count);
  Vector3 position = position_sum / double(count);
  cp.set_position(position);

  return error_sum;
}

bool vw::ba::build_control_network(bool triangulate_control_points,
                                   ba::ControlNetwork& cnet,
                                   std::vector<boost::shared_ptr<camera::CameraModel>>
                                   const& camera_models,
                                   std::vector<std::string> const& image_files,
                                   std::map<std::pair<int, int>, std::string> const& match_files,
                                   size_t min_matches,
                                   double min_angle_radians,
                                   double forced_triangulation_distance,
                                   int max_pairwise_matches) {

  typedef std::tuple<double, double, double> ipTriplet; // ip.x, ip.y, ip.scale
  typedef std::pair<std::vector<ip::InterestPoint>, std::vector<ip::InterestPoint>> MATCH_PAIR;
  typedef std::map<std::pair<int, int>, MATCH_PAIR> MATCH_MAP;
  MATCH_MAP ip_subset;

  // Note that this statement does not clear the network fully.
  // TODO: Clear all items here. 
  cnet.clear();
  
  int num_images = image_files.size();

  size_t count = 0;
  BOOST_FOREACH(std::string const& file, image_files) {
    fs::path file_path(file);
    cnet.add_image_name(file);
    count++;
  }

  // Iterate through the match files passed in and record the ones which exist.
  // TODO(oalexan1): Likely this loop and the one below can be
  // integrated and there is no need for match_files_vec, index1_vec, index2_vec.
  std::vector<std::string> match_files_vec;
  std::vector<size_t> index1_vec, index2_vec;
  typedef std::map<std::string,size_t>::iterator MapIterator;
  for (int i = 0; i < num_images; i++){ // Loop through all image pairs
    for (int j = i+1; j < num_images; j++){

      // Find the match file for this pair of images
      std::pair<int, int> pair_ind(i, j);
      auto pair_it = match_files.find(pair_ind);
      if (pair_it == match_files.end())
        continue; // This match file was not passed in
      std::string match_file = pair_it->second;

      // Verify that the match file exists
      if (!fs::exists(match_file)) {
        vw_out(WarningMessage) << "Missing match file: " << match_file << std::endl;
        continue;
      }

      // Also see what happens when GCP is added, when some similar logic is used.
      match_files_vec.push_back(match_file);
      index1_vec.push_back(i);
      index2_vec.push_back(j);
    }
  }

  // Give all interest points in a given image a unique id, and put
  // them in a vector with the id corresponding to the interest point
  std::vector<std::map<ipTriplet, int>> keypoint_map(num_images);

  // Loop through the match files
  size_t num_load_rejected = 0, num_loaded = 0;
  for (size_t file_iter = 0; file_iter < match_files_vec.size(); file_iter++) {
    std::string match_file = match_files_vec[file_iter];
    size_t index1 = index1_vec[file_iter];
    size_t index2 = index2_vec[file_iter];

    // Actually read in the file as it seems we've found something correct
    std::vector<ip::InterestPoint> ip1, ip2;
    // vw_out() << "Loading: " << match_file << std::endl;
    ip::read_binary_match_file(match_file, ip1, ip2);
    if (ip1.size() < min_matches) {
      vw_out(DebugMessage,"ba") << "\t" << match_file << "    "
                                << ip1.size() << " matches. [rejected]\n";
      num_load_rejected += ip1.size();
      continue;
    }
    vw_out() << "Match file " << match_file << " has " << ip1.size() << " matches.\n";

    // Remove descriptors from interest points and correct scale
    std::for_each(ip1.begin(), ip1.end(), ip::remove_descriptor);
    std::for_each(ip2.begin(), ip2.end(), ip::remove_descriptor);
    std::for_each(ip1.begin(), ip1.end(), safe_measurement);
    std::for_each(ip2.begin(), ip2.end(), safe_measurement);

    if (max_pairwise_matches >= 0 && (int)ip1.size() > max_pairwise_matches) {
      vw_out() << "Reducing the number of matches to: " << max_pairwise_matches << ".\n";

      std::vector<int> subset;
      vw::math::pick_random_indices_in_range(ip1.size(), max_pairwise_matches, subset);
      std::sort(subset.begin(), subset.end()); // sort the indices; not strictly necessary
      
      std::vector<ip::InterestPoint> ip1_full, ip2_full;
      ip1_full.swap(ip1);
      ip2_full.swap(ip2);
      
      ip1.resize(max_pairwise_matches);
      ip2.resize(max_pairwise_matches);
      for (size_t it = 0; it < subset.size(); it++) {
        ip1[it] = ip1_full[subset[it]];
        ip2[it] = ip2_full[subset[it]];
      }
    }

    num_loaded += ip1.size();

    for (size_t ip_it = 0; ip_it < ip1.size(); ip_it++) {
      auto dist_left_ip  = ipTriplet(ip1[ip_it].x, ip1[ip_it].y, ip1[ip_it].scale);
      auto dist_right_ip = ipTriplet(ip2[ip_it].x, ip2[ip_it].y, ip2[ip_it].scale);
      // Initialize key keypoint map to zero. Will populate the entities later.
      keypoint_map[index1][dist_left_ip] = 0;
      keypoint_map[index2][dist_right_ip] = 0;
    }

    // Save the matches after getting a subset
    ip_subset[std::make_pair<int, int>(index1, index2)] = std::make_pair(ip1, ip2);
   } // End loop through match files

  if (num_load_rejected != 0)
    vw_out(WarningMessage,"ba")
      << "Did not load " << num_load_rejected << " matches "
      << "due to inadequacy. Decrease the --min-matches parameter "
      << "to load smaller sets of matches.\n";
  vw_out() << "Loaded " << num_loaded << " matches.\n";

  Stopwatch watch;
  watch.start();

  std::vector<std::vector<ipTriplet>> keypoint_vec; // for direct access later
  keypoint_vec.resize(num_images);
  for (size_t cid = 0; cid < num_images; cid++) {
    keypoint_vec[cid].resize(keypoint_map[cid].size());
    int fid = 0;
    for (auto ip_it = keypoint_map[cid].begin(); ip_it != keypoint_map[cid].end();
         ip_it++) {
      auto& dist_ip = ip_it->first;  // alias
      keypoint_map[cid][dist_ip] = fid;
      keypoint_vec[cid][fid] = dist_ip;
      fid++;
    }
  }

  // If feature A in image I matches feather B in image J, which
  // matches feature C in image K, then (A, B, C) belong together in
  // a track, and will have a single triangulated xyz. Build such a track.
  VwOpenMVG::matching::PairWiseMatches match_map;
  for (auto it = ip_subset.begin(); it != ip_subset.end(); it++) {
    std::pair<int, int> const& cid_pair = it->first;     // alias
    
    int left_cid = cid_pair.first;
    int right_cid = cid_pair.second;
    
    MATCH_PAIR const& match_pair = it->second;  // alias
    std::vector<ip::InterestPoint> const& left_ip_vec = match_pair.first;
    std::vector<ip::InterestPoint> const& right_ip_vec = match_pair.second;
    std::vector<VwOpenMVG::matching::IndMatch> mvg_matches;
    for (size_t ip_it = 0; ip_it < left_ip_vec.size(); ip_it++) {
      auto dist_left_ip  = ipTriplet(left_ip_vec[ip_it].x, left_ip_vec[ip_it].y,
                                     left_ip_vec[ip_it].scale);
      auto dist_right_ip = ipTriplet(right_ip_vec[ip_it].x, right_ip_vec[ip_it].y,
                                     right_ip_vec[ip_it].scale);
      
      int left_fid = keypoint_map[left_cid][dist_left_ip];
      int right_fid = keypoint_map[right_cid][dist_right_ip];
      mvg_matches.push_back(VwOpenMVG::matching::IndMatch(left_fid, right_fid));
    }
    match_map[cid_pair] = mvg_matches;
  }

  // Deallocate data that is not needed anymore
  ip_subset.clear(); ip_subset = MATCH_MAP();
  keypoint_map.clear(); keypoint_map.shrink_to_fit();

  {
    // Build tracks
    // De-allocate these as soon as not needed to save memory
    VwOpenMVG::tracks::TracksBuilder trackBuilder;
    trackBuilder.Build(match_map);  // Build:  Efficient fusion of correspondences
    trackBuilder.Filter();          // Filter: Remove tracks that have conflict
    // Export tracks as a map (each entry is a sequence of imageId and featureIndex):
    //  {TrackIndex => {(imageIndex, featureIndex), ... ,(imageIndex, featureIndex)}
    VwOpenMVG::tracks::STLMAPTracks map_tracks;
    trackBuilder.ExportToSTL(map_tracks);
    match_map = VwOpenMVG::matching::PairWiseMatches();  // wipe this, no longer needed
    trackBuilder = VwOpenMVG::tracks::TracksBuilder();   // wipe it
    if (map_tracks.empty())
      return false; 

    // Populate the filtered tracks
    size_t num_elems = map_tracks.size();
    for (auto pid = map_tracks.begin(); pid != map_tracks.end(); pid++) {
      ControlPoint cpoint(ControlPoint::TiePoint);
      for (auto cid_fid = (pid->second).begin(); cid_fid != (pid->second).end(); cid_fid++) {
        int cid = cid_fid->first;
        int fid = cid_fid->second;
        auto const& dist_ip = keypoint_vec.at(cid).at(fid);
        cpoint.add_measure(ControlMeasure(std::get<0>(dist_ip), // x position
                                          std::get<1>(dist_ip), // y position
                                          std::get<2>(dist_ip), // col sigma
                                          std::get<2>(dist_ip), // row sigma
                                          cid)); // image index
      }
      
      if (cpoint.size() > 0)
        cnet.add_control_point(cpoint);
    }
  }
  
  watch.stop();
  vw_out() << "Building the control network took " << watch.elapsed_seconds() << " seconds.\n";
  
  // Triangulating positions
  std::int64_t num_total_points = 0, num_failed_points = 0;
  if (triangulate_control_points) {
    TerminalProgressCallback progress("ba", "Triangulating: ");
    progress.report_progress(0);
    double inc_prog = 1.0/double(cnet.size());
    BOOST_FOREACH(ba::ControlPoint& cpoint, cnet) {
      progress.report_incremental_progress(inc_prog);
      int ans = ba::triangulate_control_point(cpoint, camera_models, min_angle_radians,
                                              forced_triangulation_distance);
      num_total_points++;
      if (ans < 0) 
        num_failed_points++;
    }
    progress.report_finished();
  }

  if (num_failed_points > 0)
    vw_out() << "\n" << "Triangulated successfully "
             << num_total_points - num_failed_points << " out of " << num_total_points
             << " points (ratio: " << 1.0 - num_failed_points / double(num_total_points)
             << "). If too many failures, perhaps your baseline/convergence angle is too small. "
             << "Or consider deleting your run directory and rerunning with more match points, "
             << "decreasing --min-triangulation-angle, or using "
             << "--forced-triangulation-distance.\n";
          
  return true;
}

int vw::ba::add_ground_control_points(vw::ba::ControlNetwork& cnet,
                                       std::vector<std::string> const& gcp_files,
                                       cartography::Datum const& datum){
  
  namespace fs = boost::filesystem;
  
  std::vector<std::string> const& image_files = cnet.get_image_list();

  // Creating a version of image_files that doesn't contain the path.
  typedef std::map<std::string,size_t> LookupType;
  LookupType image_lookup;
  for (size_t i = 0; i < image_files.size(); i++) {
    // TODO(oalexan1): This is fragile logic. What if the same
    // image exists in two directories?
    // TODO(oalexan1): Wipe this. Test if it makes a difference.
    image_lookup[image_files[i]] = i;
    image_lookup[fs::path(image_files[i]).filename().string()] = i;
  }

  std::vector<std::string>::const_iterator gcp_iter = gcp_files.begin();
  std::vector<std::string>::const_iterator gcp_end  = gcp_files.end();
  int num_gcp = 0;
  while (gcp_iter != gcp_end) {

    vw_out() << "Loading GCP file: " << *gcp_iter << std::endl;

    if (!fs::exists(*gcp_iter)) {
      vw_throw(ArgumentErr() << "GCP file " << *gcp_iter << " does not exist!");
    }

    std::ifstream ifile((*gcp_iter).c_str());
    std::string line;
    while (getline(ifile, line, '\n')){
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
            vw_throw(ArgumentErr() << "Standard deviations must be positive "
                                    << "when loading ground control points.");
          }
          measure_locations.push_back(temp_loc);
          measure_cameras.push_back(temp_name);
        }else{
                break;
        }
      }

      if (world_sigma[0] <= 0 || world_sigma[1] <= 0 || world_sigma[2] <= 0)
        vw_throw(ArgumentErr() << "Standard deviations must be positive "
                                            << "when loading ground control points.");

      // Make lat,lon into lon,lat
      std::swap(world_location[0], world_location[1]);

      // Building Control Point
      Vector3 xyz = datum.geodetic_to_cartesian(world_location);

      vw_out(VerboseDebugMessage,"ba") << "\t\tLocation: " << xyz << std::endl;
      ControlPoint cpoint(ControlPoint::GroundControlPoint);
      cpoint.set_position(xyz[0],         xyz[1],         xyz[2]       );
      cpoint.set_sigma   (world_sigma[0], world_sigma[1], world_sigma[2]);

      // Adding measures
      std::vector<Vector4    >::iterator m_iter_loc  = measure_locations.begin();
      std::vector<std::string>::iterator m_iter_name = measure_cameras.begin();
      while (m_iter_loc != measure_locations.end()) {
        LookupType::iterator it = image_lookup.find(*m_iter_name);
        if (it != image_lookup.end()) {
          vw_out(DebugMessage,"ba") << "\t\tAdded Measure: " << *m_iter_name
                                    << " #" << it->second << std::endl;
          ControlMeasure cm((*m_iter_loc)[0], (*m_iter_loc)[1],
                             (*m_iter_loc)[2], (*m_iter_loc)[3], it->second);
          cpoint.add_measure(cm);
        } else {
          vw_out(WarningMessage,"ba") << "No input image found matching "
                                      << *m_iter_name << std::endl;
        }
        m_iter_loc++;
        m_iter_name++;
      }

      // Append the GCP
      if (cpoint.size() > 0) {
        cnet.add_control_point(cpoint);
        num_gcp++;
      }
      
    } // End line loop
    ifile.close();

    gcp_iter++;
  } // End file loop

  return num_gcp;
}
