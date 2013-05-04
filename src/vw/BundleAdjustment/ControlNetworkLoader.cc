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

void vw::ba::triangulate_control_point( ControlPoint& cp,
                                        std::vector<boost::shared_ptr<camera::CameraModel> > const& camera_models,
                                        double const& minimum_angle ) {
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
        stereo::StereoModel sm( camera_models[ j_cam_id ].get(),
                                camera_models[ k_cam_id ].get() );

        if ( sm.convergence_angle( cp[j].position(),
                                   cp[k].position() ) >
             minimum_angle ) {
          count++;
          position_sum += sm( cp[j].position(), cp[k].position(), error );
          error_sum += error;
        }
      } catch ( const camera::PixelToRayErr& ) { /* Just let it go */ }
    }
  }

  // 4.2.) Summing, Averaging, and Storing
  if ( !count ) {
    vw_out(WarningMessage,"ba") << "Unable to triangulation position for point!\n";
    // At the very least we can provide a point that is some
    // distance out from the camera center and is in the 'general'
    // area.
    size_t j = cp[0].image_id();
    try {
      cp.set_position( camera_models[j]->camera_center(cp[0].position()) +
                       camera_models[j]->pixel_to_vector(cp[0].position())*10 );
    } catch ( const camera::PixelToRayErr& ) {
      cp.set_position( camera_models[j]->camera_center(cp[0].position()) +
                       camera_models[j]->camera_pose(cp[0].position()).rotate(Vector3(0,0,10)) );
    }
  } else {
    error_sum /= double(count);
    cp.set_position( position_sum / double(count) );
  }
}

void vw::ba::build_control_network( ba::ControlNetwork& cnet,
                                     std::vector<boost::shared_ptr<camera::CameraModel> > const& camera_models,
                                     std::vector<std::string> const& image_files,
                                     size_t min_matches,
                                     std::vector<std::string> const& directories ) {
  cnet.clear();

  // We can't guarantee that image_files is sorted, so we make a
  // std::map to give ourselves a sorted list and access to a binary
  // search.
  std::map<std::string,size_t> image_prefix_map;
  size_t count = 0;
  ba::CameraRelationNetwork<ba::IPFeature> crn;
  BOOST_FOREACH( std::string const& file, image_files ) {
    fs::path file_path(file);
    image_prefix_map[fs::path(file_path).replace_extension().string()] = count;
    crn.add_node( ba::CameraNode<ba::IPFeature>( count,
                                                 file_path.stem().string() ) );
    count++;
  }

  // Searching through the directories available to us.
  size_t num_load_rejected = 0, num_loaded = 0;
  BOOST_FOREACH( std::string const& directory, directories ) {
    vw_out(VerboseDebugMessage,"ba") << "\tOpening directory \""
                                     << directory << "\".\n";
    fs::directory_iterator end_itr;
    for ( fs::directory_iterator obj( directory ); obj != end_itr; ++obj ) {

      // Skip if not a file
      if ( !fs::is_regular_file( obj->status() ) ) continue;
      // Skip if not a match file
      if ( obj->path().extension() != ".match" ) continue;

      // Pull out the prefixes that made up that match file
      std::string match_base = obj->path().stem().string();
      size_t split_pt = match_base.find("__");
      if ( split_pt == std::string::npos ) continue;
      std::string prefix1 = match_base.substr(0,split_pt);
      std::string prefix2 = match_base.substr(split_pt+2,match_base.size()-split_pt-2);

      // Extract the image indices that correspond to image prefixes.
      typedef std::map<std::string,size_t>::iterator MapIterator;
      MapIterator it1 = image_prefix_map.find( prefix1 );
      MapIterator it2 = image_prefix_map.find( prefix2 );
      if ( it1 == image_prefix_map.end() ||
           it2 == image_prefix_map.end() ) continue;

      // Actually read in the file as it seems we've found something correct
      std::vector<ip::InterestPoint> ip1, ip2;
      ip::read_binary_match_file( obj->path().string(), ip1, ip2 );
      if ( ip1.size() < min_matches ) {
        vw_out(VerboseDebugMessage,"ba") << "\t" << obj->path().string() << "    "
                                         << it1->second << " <-> " << it2->second << " : "
                                         << ip1.size() << " matches. [rejected]\n";
        num_load_rejected += ip1.size();
      } else {
        vw_out(VerboseDebugMessage,"ba") << "\t" << obj->path().string() << "    "
                                         << it1->second << " <-> " << it2->second << " : "
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
        for ( size_t k = 0; k < ip1.size(); k++ ) {
          f_itr ipfeature1 = std::find_if( crn[it1->second].begin(),
                                           crn[it1->second].end(),
                                           ContainsEqualIP( ip1[k] ) );
          f_itr ipfeature2 = std::find_if( crn[it2->second].begin(),
                                           crn[it2->second].end(),
                                           ContainsEqualIP( ip2[k] ) );
          if ( ipfeature1 == crn[it1->second].end() ) {
            crn[it1->second].relations.push_front( f_ptr( new ba::IPFeature( ip1[k], it1->second ) ) );
            ipfeature1 = crn[it1->second].begin();
          }
          if ( ipfeature2 == crn[it2->second].end() ) {
            crn[it2->second].relations.push_front( f_ptr( new ba::IPFeature( ip2[k], it2->second ) ) );
            ipfeature2 = crn[it2->second].begin();
          }

          // Doubly linking
          (*ipfeature1)->connection( *ipfeature2, false );
          (*ipfeature2)->connection( *ipfeature1, false );
        }
      }
    }
  } // end search through directories
  if ( num_load_rejected != 0 ) {
    vw_out(WarningMessage,"ba") << "\tDidn't load " << num_load_rejected
                                << " matches due to inadequacy.\n";
    vw_out(WarningMessage,"ba") << "\tLoaded " << num_loaded << " matches.\n";
  }

  // Building control network
  crn.write_controlnetwork( cnet );

  // Triangulating Positions
  {
    TerminalProgressCallback progress("ba", "Triangulating:");
    progress.report_progress(0);
    double inc_prog = 1.0/double(cnet.size());
    double min_angle = 5.0*M_PI/180.0;
    BOOST_FOREACH( ba::ControlPoint& cpoint, cnet ) {
      progress.report_incremental_progress(inc_prog );
      ba::triangulate_control_point( cpoint, camera_models, min_angle );
    }
    progress.report_finished();
  }
}
