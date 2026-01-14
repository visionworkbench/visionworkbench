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
#include <vw/InterestPoint/MatcherIO.h>
#include <vw/InterestPoint/InterestPointUtils.h>
#include <vw/Math/RandomSet.h>
#include <vw/Core/Stopwatch.h>

#include <mutex>
#include <boost/filesystem/fstream.hpp>

using namespace vw;
using namespace vw::ba;
namespace fs = boost::filesystem;

namespace {
  // Simple structure to hold a pair of feature indices
  struct IndMatch {
    int m_left;
    int m_right;
    IndMatch(int i=0, int j=0) : m_left(i), m_right(j) {}
  };

  // Map from (cid1, cid2) -> vector of matches
  typedef std::map<std::pair<int, int>, std::vector<IndMatch>> PairWiseMatches;

  // Map from TrackID -> (ImageID -> FeatureID)
  // This replaces VwOpenMVG::tracks::STLMAPTracks
  typedef std::map<int, std::map<int, int>> TrackMap; 
}

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
  
  // Do pairwise triangulation then average the results. For speed, consider
  // only each ray and next 10 rays. Should be enough. Note that as long as at
  // least two of the rays meet with a triangulation angle no less than
  // min_angle_radians, triangulation will succeed, and these successful
  // triangulation points will be averaged. Hence, while some of the rays can
  // meet at a very small angle, their intersection, being unreliable, will not
  // be added to the mix. This allows, for example, for a triplet of rays, R1,
  // R2, R3, all corresponding to matching features and converging to the same
  // triangulated point, so that R1 and R2 have a very small angle, but R1 and
  // R3, then R2 and R3 have a very solid angle. This triangulated point should
  // be robust enough if the min triangulation angle is, say, no less than 5
  // degrees, and ideally 15-20 degrees or more.
  for (size_t j = 0; j < cp.size(); j++) {
    for (size_t k = j + 1; k < std::min(j + 11, cp.size()); k++) {
      size_t j_cam_id = cp[j].image_id();
      size_t k_cam_id = cp[k].image_id();
      
      // Make sure camera centers are not equal
      try {
        // This trips up the CSM frame camera model
        auto const& c1 = camera_models[j_cam_id]->camera_center(cp[j].position());
        auto const& c2 = camera_models[k_cam_id]->camera_center(cp[k].position());
        if (norm_2(c1 - c2) <= 1e-6)
          continue;
      } catch (...) {
        continue;
      }
        
      try {
        stereo::StereoModel sm(camera_models[j_cam_id].get(), camera_models[k_cam_id].get(),
                               angle_tol); // use angle_tol
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
        continue;
      }
    }
  }
  
  // Summing, averaging, and storing
  if (count == 0) {
    // It failed to triangulate. At the very least we can provide a point that is some
    // distance out from the camera center and is in the 'general' area.
    size_t j = cp[0].image_id();
    // Set the point to be a certain distance from the camera center
    // TODO(oalexan1): Later the flag that the triangulation failed is not checked.
    // It is not clear if that should be done. This may end up in too many points
    // thrown out.
    double dist = 1000000.0; // 1000 km
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
        cp.set_ignore(true);
        return -1; // nothing can be done
      }
    }
    
    if (forced_triangulation_distance > 0) 
      return 1; // mark triangulation as successful
    
    cp.set_ignore(true);
    return -1;
  }
  
  error_sum /= double(count);
  Vector3 position = position_sum / double(count);
  cp.set_position(position);

  return error_sum;
}

// Triangulate the points in a control network. Do not triangulate
// GCP or points constrained to a DEM. 
void vw::ba::triangulate_control_network(vw::ba::ControlNetwork& cnet,
                                         std::vector<vw::CamPtr> const& camera_models,
                                         double min_angle_radians,
                                         double forced_triangulation_distance) {

  std::int64_t num_total_points = 0, num_failed_points = 0;
  TerminalProgressCallback progress("ba", "Triangulating: ");
  progress.report_progress(0);
  double inc_prog = 1.0/double(cnet.size());
  BOOST_FOREACH(ba::ControlPoint& cpoint, cnet) {
    progress.report_incremental_progress(inc_prog);
    
    if (cpoint.type() == ControlPoint::GroundControlPoint ||
        cpoint.type() == ControlPoint::PointFromDem)
      continue; // Skip GCPs and points from a DEM
    
    int ans = ba::triangulate_control_point(cpoint, camera_models, min_angle_radians,
                                            forced_triangulation_distance);
    num_total_points++;
    if (ans < 0) 
      num_failed_points++;
  }
  progress.report_finished();

  // This is a verbose message, but it helps when users are confused as to why
  // their bundle adjustment failed.
  double failed_ratio = num_failed_points / double(num_total_points);
  if (num_failed_points > 0 || num_total_points == 0) {
    vw_out() << "Triangulated successfully "
             << num_total_points - num_failed_points << " out of " << num_total_points
             << " control points (ratio: " 
             << 1.0 - failed_ratio << ").\n";
    if (num_total_points < 30 || failed_ratio > 0.5)
      vw_out() << "If too many failures, perhaps your baseline/convergence angle "
              << "is too small. Consider deleting your run directory and rerunning "
              << "with more match points, decreasing --min-triangulation-angle, or using "
              << "--forced-triangulation-distance.\n";
  }
  
  return;
}

// Notation
typedef std::pair<std::vector<ip::InterestPoint>, std::vector<ip::InterestPoint>> MATCH_PAIR;
typedef std::map<std::pair<int, int>, MATCH_PAIR> MATCH_MAP;
typedef std::tuple<double, double, double> ipTriplet; // (ip.x, ip.y, ip.scale)

// Convert the matches to the MVG format, so their track-building code can be used.
void matchesToMvg(MATCH_MAP const& match_map, 
                  std::vector<std::map<ipTriplet, int>> const& keypoint_map,
                  std::vector<std::vector<ipTriplet>> const& keypoint_vec,
                  std::map<std::pair<int, int>, std::string> const& match_files,
                  PairWiseMatches& mvg_match_map) {

  // Wipe the output
  mvg_match_map.clear();
  
  for (auto it = match_map.begin(); it != match_map.end(); it++) {
    std::pair<int, int> const& cid_pair = it->first;     // alias
    
    int left_cid = cid_pair.first;
    int right_cid = cid_pair.second;
    
    // Must not have left_cid == right_cid
    if (left_cid == right_cid) 
      vw::vw_throw(ArgumentErr() 
                   << "Bookkeeping error: Cannot have matches from an image to itself.\n");

    // This potential swap ensures that the pair is always ordered
    // in the same way. It results in almost the same results regardless
    // of the order of images vs order of matches.
    bool swap = (left_cid > right_cid);
    std::pair<int, int> out_pair = cid_pair;
    if (swap) 
      out_pair = std::make_pair(right_cid, left_cid);
      
    MATCH_PAIR const& match_pair = it->second;  // alias
    std::vector<ip::InterestPoint> const& left_ip_vec = match_pair.first;
    std::vector<ip::InterestPoint> const& right_ip_vec = match_pair.second;
    std::vector<IndMatch> mvg_matches;
    for (size_t ip_it = 0; ip_it < left_ip_vec.size(); ip_it++) {
      auto dist_left_ip  = ipTriplet(left_ip_vec[ip_it].x, left_ip_vec[ip_it].y,
                                      left_ip_vec[ip_it].scale);
      auto dist_right_ip = ipTriplet(right_ip_vec[ip_it].x, right_ip_vec[ip_it].y,
                                      right_ip_vec[ip_it].scale);

      auto left_it = keypoint_map[left_cid].find(dist_left_ip);
      if (left_it == keypoint_map[left_cid].end()) 
        vw::vw_throw(ArgumentErr() << "Bookkeeping error in matchesToMvg().\n");
      auto right_it = keypoint_map[right_cid].find(dist_right_ip);
      if (right_it == keypoint_map[right_cid].end()) 
        vw::vw_throw(ArgumentErr() << "Bookkeeping error in matchesToMvg().\n");
      int left_fid = left_it->second;
      int right_fid = right_it->second;

      if (!swap)
        mvg_matches.push_back(IndMatch(left_fid, right_fid));
      else
        mvg_matches.push_back(IndMatch(right_fid, left_fid));
    }

    // Append to vector mvg_match_map[out_pair] the vector mvg_matches.
    // As such, can have both left-to-right and right-to-left matches.
    mvg_match_map[out_pair].insert(mvg_match_map[out_pair].end(), 
                                   mvg_matches.begin(), mvg_matches.end());
  }
  
}

// Build tracks from pairwise matches. This logic was shown to be superior to the 
// OpenMVG :TracksBuilder. It can merge pairs (i, k), (j, k) into (i, j, k), which
// that one could not. It produces more long tracks in general.
// The terminology is as follows. Each track has the form: 
//  {TrackIndex => {(imageIndex, featureIndex), ..., (imageIndex, featureIndex)}
// Below, pid is the track id, cid is camera index, fid is feature index.
void buildTracks(PairWiseMatches const& mvg_match_map,
                   TrackMap & pid_cid_fid) {

  // Wipe the output
  pid_cid_fid.clear();

  // This alternative bookkeeping helps lookup a track based on cid.
  std::map<int, std::map<int, int>> cid_fid_pid;
  
  // In mvg_match_map we have cid pairs as keys. For each such pair, we have fid pairs.
  for (auto it = mvg_match_map.begin(); it != mvg_match_map.end(); it++) {
    std::pair<int, int> const& cid_pair = it->first;     // alias
    std::vector<IndMatch> const& mvg_matches = it->second;

    int left_cid = cid_pair.first;
    int right_cid = cid_pair.second;
    for (size_t ip_it = 0; ip_it < mvg_matches.size(); ip_it++) {
      int left_fid = mvg_matches[ip_it].m_left;
      int right_fid = mvg_matches[ip_it].m_right;
      
      // See if the left feature is already in a track
      auto left_cid_it = cid_fid_pid.find(left_cid);
      int left_pid = -1;
      if (left_cid_it != cid_fid_pid.end()) {
        auto left_fid_it = left_cid_it->second.find(left_fid);
        if (left_fid_it != left_cid_it->second.end())
          left_pid = left_fid_it->second;
      }
      
      // See if the right feature is already in a track
      auto right_cid_it = cid_fid_pid.find(right_cid);
      int right_pid = -1;
      if (right_cid_it != cid_fid_pid.end()) {
        auto right_fid_it = right_cid_it->second.find(right_fid);
        if (right_fid_it != right_cid_it->second.end())
          right_pid = right_fid_it->second;
      }
      
      // If both left and right features are in some track, two options exist.
      // Either they are in the same track, then nothing to do. Or they are in
      // different tracks. It is not clear if merging the tracks is a good
      // thing, as maybe this pair is wrong match. Just do nothing. Works well
      // enough.
      if (left_pid >= 0 && right_pid >= 0)
        continue;
      
      // If the left feature is in a track, but the right one is not, add it to
      // the left feature track.
      if (left_pid >= 0 && right_pid < 0) {
        cid_fid_pid[right_cid][right_fid] = left_pid;
        pid_cid_fid[left_pid][right_cid] = right_fid;
        continue;
      }  
    
      // If the right feature is in a track, but the left one is not, add it to
      // the right feature track.
      if (left_pid < 0 && right_pid >= 0) {
        cid_fid_pid[left_cid][left_fid] = right_pid;
        pid_cid_fid[right_pid][left_cid] = left_fid;
        continue;
      }
      
      // If neither feature is in a track, create a new track.
      if (left_pid < 0 && right_pid < 0) {
        int pid = pid_cid_fid.size(); // one past the last existing pid
        cid_fid_pid[left_cid][left_fid] = pid;
        cid_fid_pid[right_cid][right_fid] = pid;
        pid_cid_fid[pid][left_cid] = left_fid;
        pid_cid_fid[pid][right_cid] = right_fid;
      }
    }
  }
  
  return;
}

// Create tracks and a cnet from matches. Return false if the cnet is empty.
bool matchMapToCnet(std::vector<std::string> const& image_files,
                    std::vector<std::vector<ipTriplet>> const& keypoint_vec,
                    PairWiseMatches const& mvg_match_map,
                    vw::ba::ControlNetwork& cnet) {

  // Wipe fully the network. This does not allow passing a name to it, but that
  // looks unnecessary.
  cnet = ba::ControlNetwork("ASP_control_network");

  // Add image names  
  int num_images = image_files.size();
  for (int i = 0; i < num_images; i++)
    cnet.add_image_name(image_files[i]);

  // If feature A in image I matches feather B in image J, which matches feature
  // C in image K, then (A, B, C) belong together in a track, and will have a
  // single triangulated xyz. Build such tracks.
  TrackMap pid_cid_fid;
  buildTracks(mvg_match_map, pid_cid_fid);

  if (pid_cid_fid.empty())
    return false; 
  
  // Convert to the control network format. Will triangulate later.
  for (auto pid = pid_cid_fid.begin(); pid != pid_cid_fid.end(); pid++) {
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

  return true;
}

bool vw::ba::build_control_network(bool triangulate_control_points,
                                   ba::ControlNetwork& cnet,
                                   std::vector<vw::CamPtr> const& camera_models,
                                   std::vector<std::string> const& image_files,
                                   std::map<std::pair<int, int>, std::string> const& match_files,
                                   size_t min_matches,
                                   double min_angle_radians,
                                   double forced_triangulation_distance,
                                   int max_pairwise_matches,
                                   std::map<std::pair<int, int>, double> const& 
                                   match_sigmas) {

  // TODO(oalexan1): Must be able to handle the case when the matches
  // are from an image later in the list to an image earlier in the list.

  MATCH_MAP match_map;

  // Give all interest points in a given image a unique id, and put
  // them in a vector with the id corresponding to the interest point
  int num_images = image_files.size();
  std::vector<std::map<ipTriplet, int>> keypoint_map(num_images);
  size_t num_load_rejected = 0, num_loaded = 0;
  for (auto it = match_files.begin(); it != match_files.end(); it++) {
    
    std::pair<int, int> pair_ind = it->first;
    std::string const& match_file = it->second; // alias
    int index1 = pair_ind.first; 
    int index2 = pair_ind.second; 
    // Actually read in the file as it seems we've found something correct
    std::vector<ip::InterestPoint> ip1, ip2;
    // vw_out() << "Loading: " << match_file << std::endl;
    ip::read_binary_match_file(match_file, ip1, ip2);
    if (ip1.size() < min_matches) {
      vw_out(DebugMessage,"ba") << "\tRejecting " << match_file << " with "
                                << ip1.size() << " matches.\n";
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
      pick_pair_subset(ip1, ip2, max_pairwise_matches);
    }
    
    num_loaded += ip1.size();

    // The sigma value is optional. It can be used to give more weight (less sigma)
    // to matches from a certain file. These are typically files with few user-selected 
    // matches. 
    auto sigma_it = match_sigmas.find(pair_ind);
    if (sigma_it != match_sigmas.end()) {
      double sigma = sigma_it->second; 
      for (size_t ip_it = 0; ip_it < ip1.size(); ip_it++) {
        ip1[ip_it].scale *= sigma;
        ip2[ip_it].scale *= sigma;
      }
    }
  
    // Initialize the keypoint map to zero. Will populate the entities later.
    for (size_t ip_it = 0; ip_it < ip1.size(); ip_it++) {
      auto dist_left_ip  = ipTriplet(ip1[ip_it].x, ip1[ip_it].y, ip1[ip_it].scale);
      auto dist_right_ip = ipTriplet(ip2[ip_it].x, ip2[ip_it].y, ip2[ip_it].scale);
      keypoint_map[index1][dist_left_ip] = 0;
      keypoint_map[index2][dist_right_ip] = 0;
    }

    // Save the matches after getting a subset
    match_map[pair_ind] = std::make_pair(ip1, ip2);
  } // End loop through match files

  if (num_load_rejected != 0)
    vw_out(WarningMessage,"ba")
      << "Did not load a total of " << num_load_rejected << " matches due to "
      << "pairwise application of the --min-matches parameter. Decrease this "
      << "to load smaller sets of matches.\n";
  vw_out() << "Loaded " << num_loaded << " matches from all files.\n";

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

  // Convert the matches to the MVG format, so their track-building code can be used.
  PairWiseMatches mvg_match_map;
  matchesToMvg(match_map, keypoint_map, keypoint_vec, match_files, mvg_match_map);
  
  // Deallocate data that is not needed anymore
  match_map = MATCH_MAP();
  keypoint_map = std::vector<std::map<ipTriplet, int>>();

  // Build the tracks and the control network using MVG
  bool ans = matchMapToCnet(image_files, keypoint_vec, mvg_match_map, cnet);
  if (!ans) 
   return false;
 
  // Wipe this, no longer needed
  mvg_match_map.clear();
  
  watch.stop();
  vw_out() << "Building the control network took " << watch.elapsed_seconds() 
    << " seconds.\n";
  
  // Triangulate the points in a control network
  if (triangulate_control_points)
    vw::ba::triangulate_control_network(cnet, camera_models, min_angle_radians,
                                        forced_triangulation_distance);

  return true;
}

// A little function to parse the WKT
bool parseDatum(std::string const& wkt, vw::cartography::Datum & datum) {
  
  std::string buff, val;
  std::istringstream is(wkt);

  while (is >> val) {
    // Ignore the pound sign and WKT: 
    if (val.find("#") != std::string::npos) continue;
    if (val.find("WKT:") != std::string::npos) continue;
    buff += val + " ";
  }
  
  // Empty buff means that the string was empty
  if (buff.empty()) 
    return false;
  
  datum.set_wkt(buff);
  
  return true;
}

int vw::ba::add_ground_control_points(vw::ba::ControlNetwork& cnet,
                                      std::vector<std::string> const& gcp_files,
                                      cartography::Datum const& datum,
                                      bool skip_datum_check) {
  
  namespace fs = boost::filesystem;
  
  std::vector<std::string> const& image_files = cnet.get_image_list();

  // Initalize this as the datum passed from outside. 
  vw::cartography::Datum local_datum = datum;

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

    if (!fs::exists(*gcp_iter))
      vw_throw(ArgumentErr() << "GCP file " << *gcp_iter << " does not exist!");

    std::ifstream ifile((*gcp_iter).c_str());
    std::string line;
    int count = 0;
    while (getline(ifile, line, '\n')) {
      
      // Skip empty lines
      if (line.size() == 0) continue;
      
      // The first line can be the WKT string. Then read and validate it.
      if (line.size() > 0 && line[0] == '#' && line.find("WKT:") != std::string::npos &&
          count == 0) { 
        count++; // for next time
        cartography::Datum gcp_datum;
        if (!parseDatum(line, gcp_datum)) 
          continue;
        
        // Use the datum from the GCP file
        local_datum = gcp_datum;
          
        if (skip_datum_check) 
          continue;
        double tol = 1e-6;
        if (std::abs(gcp_datum.semi_major_axis() - datum.semi_major_axis()) > tol ||
            std::abs(gcp_datum.semi_minor_axis() - datum.semi_minor_axis()) > tol ||
            std::abs(gcp_datum.meridian_offset() - datum.meridian_offset()) > tol) {
          vw::vw_out(vw::WarningMessage)
            << "The datum of the GCP file " << *gcp_iter
            << " does not match the datum passed on input. Will use the datum "
            << "from the GCP file.\n";
        } 
        continue;
      }
      count++;
      
      // Skip lines starting with comments
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
            >> world_sigma[1] >> world_sigma[2])) {
        vw_out(WarningMessage) << "Could not parse a ground control point "
                               << "from line: " << line << std::endl;
        continue;
      }

      // Other elements in the line define the position in images
      while (1) {
        std::string image_name;
        Vector4 pix_std;
        if (is >> image_name >> pix_std[0] >> pix_std[1] >> pix_std[2] >> pix_std[3]){
          if (pix_std[2] <= 0 || pix_std[3] <= 0)
            vw_throw(ArgumentErr() << "Standard deviations must be positive "
                                    << "when loading ground control points.");

          measure_locations.push_back(pix_std);
          measure_cameras.push_back(image_name);
        } else {
          break;
        }
      }

      if (world_sigma[0] <= 0 || world_sigma[1] <= 0 || world_sigma[2] <= 0)
        vw_throw(ArgumentErr() << "Standard deviations must be positive "
                                            << "when loading ground control points.");

      // Make lat,lon into lon,lat
      std::swap(world_location[0], world_location[1]);

      // Convert GCP from lon,lat,height to ECEF
      // TODO(oalexan1): Support here projected coordinates. Use
      // the input georef to convert them to ECEF.
      Vector3 xyz = local_datum.geodetic_to_cartesian(world_location);

      vw_out(VerboseDebugMessage,"ba") << "\t\tLocation: " << xyz << std::endl;
      ControlPoint cpoint(ControlPoint::GroundControlPoint);
      cpoint.set_position(xyz[0],         xyz[1],         xyz[2]       );
      cpoint.set_sigma   (world_sigma[0], world_sigma[1], world_sigma[2]);

      // Adding measures
      auto iter_loc  = measure_locations.begin();
      auto iter_name = measure_cameras.begin();
      while (iter_loc != measure_locations.end()) {
        LookupType::iterator it = image_lookup.find(*iter_name);
        if (it != image_lookup.end()) {
          vw_out(DebugMessage,"ba") << "\t\tAdded measure: " << *iter_name
                                    << " #" << it->second << std::endl;
          ControlMeasure cm((*iter_loc)[0], (*iter_loc)[1],
                             (*iter_loc)[2], (*iter_loc)[3], it->second);
          cpoint.add_measure(cm);
        } else {
          vw_out(WarningMessage,"ba") << "No input image found matching "
                                      << *iter_name << std::endl;
        }
        iter_loc++;
        iter_name++;
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
