// __BEGIN_LICENSE__
// __END_LICENSE__

#include <vw/BundleAdjustment/ControlNetworkLoader.h>
#include <vw/Stereo/StereoModel.h>
#include <vw/Cartography/SimplePointImageManipulation.h>

using namespace vw;
using namespace vw::ba;

#include <boost/filesystem/operations.hpp>
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

void vw::ba::build_control_network( ControlNetwork& cnet,
                                    std::vector<boost::shared_ptr<camera::CameraModel> > const& camera_models,
                                    std::vector<std::string> const& image_files,
                                    int min_matches ) {
  cnet.clear();

  // 1.) Build CRN (pop on cameras)
  CameraRelationNetwork<IPFeature> crn;
  for ( unsigned i = 0; i < image_files.size(); i++ ) {
    fs::path image_path( image_files[i] );
    crn.add_node( CameraNode<IPFeature>( i,
                                         image_path.stem() ) );
  }

  // 2.) Load up matches into CRN
  {
    TerminalProgressCallback progress("ba","Match Files: ");
    progress.report_progress(0);
    int32 num_load_rejected = 0, num_loaded = 0;
    for ( unsigned i = 0; i < image_files.size(); ++i ) {
      progress.report_progress(float(i)/float(image_files.size()));
      for ( unsigned j = i+1; j < image_files.size(); ++j ) {
        std::string match_filename =
          fs::path( image_files[i] ).replace_extension().string() + "__" +
          fs::path( image_files[j] ).stem() + ".match";

        if ( !fs::exists( match_filename ) )
          continue;

        std::vector<ip::InterestPoint> ip1, ip2;
        ip::read_binary_match_file( match_filename, ip1, ip2 );
        if ( int( ip1.size() ) < min_matches ) {
          vw_out(VerboseDebugMessage,"ba") << "\t" << match_filename << "    "
                                           << i << " <-> " << j << " : "
                                           << ip1.size() << " matches. [rejected]\n";
          num_load_rejected += ip1.size();
        } else {
          vw_out(VerboseDebugMessage,"ba") << "\t" << match_filename << "    "
                                           << i << " <-> " << j << " : "
                                           << ip1.size() << " matches.\n";
          num_loaded += ip1.size();

          typedef boost::shared_ptr< IPFeature > f_ptr;
          typedef std::list< f_ptr >::iterator f_itr;

          // Checking to see if features already exist, adding if they
          // don't, then linking them.
          for ( unsigned k = 0; k < ip1.size(); k++ ) {
            f_itr ipfeature1 = std::find_if( crn[i].begin(),
                                             crn[i].end(),
                                             ContainsEqualIP( ip1[k] ) );
            f_itr ipfeature2 = std::find_if( crn[j].begin(),
                                             crn[j].end(),
                                             ContainsEqualIP( ip2[k] ) );
            if ( ipfeature1 == crn[i].end() ) {
              crn[i].relations.push_front( f_ptr( new IPFeature( ip1[k], i ) ) );
              ipfeature1 = crn[i].begin();
            }
            if ( ipfeature2 == crn[j].end() ) {
              crn[j].relations.push_front( f_ptr( new IPFeature( ip2[k], j ) ) );
              ipfeature2 = crn[j].begin();
            }

            // Doubly linking
            (*ipfeature1)->connection( *ipfeature2, false );
            (*ipfeature2)->connection( *ipfeature1, false );
          }

        }
      } // end j for loop
    }   // end i for loop
    progress.report_finished();
    if ( num_load_rejected != 0 ) {
      vw_out(WarningMessage,"ba") << "\tDidn't load " << num_load_rejected << " matches due to inadequacy.\n";
      vw_out(WarningMessage,"ba") << "\tLoaded " << num_loaded << " matches.\n";
    }
  }

  // 3.) Building Control Network
  crn.write_controlnetwork( cnet );

  // 4.) Triangulating Positions
  {
    TerminalProgressCallback progress("ba", "Triangulating:");
    progress.report_progress(0);
    double inc_prog = 1.0/double(cnet.size());
    for ( ControlNetwork::iterator cpoint = cnet.begin();
          cpoint != cnet.end(); cpoint++ ) {
      progress.report_incremental_progress( inc_prog );

      std::vector< Vector3 > positions;
      double error = 0, error_sum = 0;

      const double min_convergence_angle = 5.0*M_PI/180.0;

      // 4.1.) Building a listing of triangulation
      for ( unsigned j = 0; j < cpoint->size(); j++ ) {
        for ( unsigned k = j+1; k < cpoint->size(); k++ ) {
          // Make sure camera centers are not equal
          int j_cam_id = (*cpoint)[j].image_id();
          int k_cam_id = (*cpoint)[k].image_id();
          if ( norm_2( camera_models[j_cam_id]->camera_center( (*cpoint)[j].position() ) -
                       camera_models[k_cam_id]->camera_center( (*cpoint)[k].position() ) ) > 1e-6 ) {

            stereo::StereoModel sm( camera_models[ j_cam_id ].get(),
                                    camera_models[ k_cam_id ].get() );

            if ( sm.convergence_angle( (*cpoint)[j].position(),
                                       (*cpoint)[k].position() ) >
                 min_convergence_angle ) {
              positions.push_back( sm( (*cpoint)[j].position(),
                                       (*cpoint)[k].position(),
                                       error ) );
              error_sum += error;
            }
          }
        }
      }

      // 4.2.) Summing, Averaging, and Storing
      if ( positions.empty() ) {
        vw_out(WarningMessage,"ba") << "Unable to triangulation position for point!\n";
        // At the very least we can provide a point that is some
        // distance out from the camera center and is in the 'general'
        // area.
        int j = (*cpoint)[0].image_id();
        cpoint->set_position( camera_models[j]->camera_center((*cpoint)[j].position()) +
                              camera_models[j]->pixel_to_vector((*cpoint)[j].position())*10 );
      } else {
        error_sum /= positions.size();
        Vector3 position_avg;
        for ( unsigned j = 0; j < positions.size(); j++ )
          position_avg += positions[j]/positions.size();
        cpoint->set_position( position_avg );
      }

    }
    progress.report_finished();
  }

}

void vw::ba::add_ground_control_points( ControlNetwork& cnet,
                                        std::vector<std::string> const& image_files,
                                        std::vector<std::string> const& gcp_files ) {
  // Creating a version of image_files that doesn't contain the path
  std::vector<std::string> pathless_image_files;
  for ( unsigned i = 0; i < image_files.size(); i++ )
    pathless_image_files.push_back(fs::path(image_files[i]).filename());

  TerminalProgressCallback progress("ba", "Grnd Control:");
  progress.report_progress(0);
  float inc_progress = 1.0/gcp_files.size();
  for ( std::vector<std::string>::const_iterator gcp_name = gcp_files.begin();
        gcp_name != gcp_files.end(); gcp_name++ ) {

    if ( !fs::exists( *gcp_name ) )
      continue;

    progress.report_incremental_progress(inc_progress);

    // Data to be loaded
    std::vector<Vector2> measure_locations;
    std::vector<std::string> measure_cameras;
    Vector3 world_location, world_sigma;

    vw_out(VerboseDebugMessage,"ba") << "\tLoading \"" << *gcp_name
                                     << "\".\n";
    int count = 0;
    std::ifstream ifile( (*gcp_name).c_str() );
    while (!ifile.eof()) {
      if ( count == 0 ) {
        // First line defines position in the world
        ifile >> world_location[0] >> world_location[1]
              >> world_location[2] >> world_sigma[0]
              >> world_sigma[1] >> world_sigma[2];
      } else {
        // Other lines define position in images
        std::string temp_name;
        Vector2 temp_loc;
        ifile >> temp_name >> temp_loc[0] >> temp_loc[1];
        measure_locations.push_back( temp_loc );
        measure_cameras.push_back( temp_name );
      }
      count++;
    }
    ifile.close();

    // Building Control Point
    Vector3 xyz = cartography::lon_lat_radius_to_xyz(world_location);
    vw_out(VerboseDebugMessage,"ba") << "\t\tLocation: "
                                     << xyz << std::endl;
    ControlPoint cpoint(ControlPoint::GroundControlPoint);
    cpoint.set_position(xyz[0],xyz[1],xyz[2]);
    cpoint.set_sigma(world_sigma[0],world_sigma[1],world_sigma[2]);

    // Adding measures
    std::vector<Vector2>::iterator m_iter_loc = measure_locations.begin();
    std::vector<std::string>::iterator m_iter_name = measure_cameras.begin();
    while ( m_iter_loc != measure_locations.end() ) {
      unsigned camera_index;
      for (camera_index = 0; camera_index < image_files.size(); camera_index++ ) {
        if ( *m_iter_name == image_files[camera_index] )
          break;
        else if ( *m_iter_name == pathless_image_files[camera_index])
          break;
      }
      if ( camera_index == image_files.size() ) {
        vw_out(WarningMessage,"ba") << "\t\tWarning: no image found matching "
                                    << *m_iter_name << std::endl;
      } else {
        vw_out(DebugMessage,"ba") << "\t\tAdded Measure: " << *m_iter_name
                                  << " #" << camera_index << std::endl;
        ControlMeasure cm( (*m_iter_loc).x(), (*m_iter_loc).y(),
                           1.0, 1.0, camera_index );
        cpoint.add_measure( cm );
      }
      m_iter_loc++;
      m_iter_name++;
    }

    // Appended GCP
    cnet.add_control_point(cpoint);
  }
  progress.report_finished();

}
