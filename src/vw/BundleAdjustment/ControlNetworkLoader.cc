// __BEGIN_LICENSE__
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

          // Remove descriptors from interest points
          std::for_each( ip1.begin(), ip1.end(), ip::remove_descriptor );
          std::for_each( ip2.begin(), ip2.end(), ip::remove_descriptor );

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
      vw_out(WarningMessage,"ba") << "\tDidn't load " << num_load_rejected
                                  << " matches due to inadequacy.\n";
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
    BOOST_FOREACH( ControlPoint& cpoint, cnet ) {
      progress.report_incremental_progress( inc_prog );

      std::vector< Vector3 > positions;
      double error = 0, error_sum = 0;

      const double min_convergence_angle = 5.0*M_PI/180.0;

      // 4.1.) Building a listing of triangulation
      for ( unsigned j = 0; j < cpoint.size(); j++ ) {
        for ( unsigned k = j+1; k < cpoint.size(); k++ ) {
          // Make sure camera centers are not equal
          int j_cam_id = cpoint[j].image_id();
          int k_cam_id = cpoint[k].image_id();
          if ( norm_2( camera_models[j_cam_id]->camera_center( cpoint[j].position() ) -
                       camera_models[k_cam_id]->camera_center( cpoint[k].position() ) ) > 1e-6 ) {

            try {
              stereo::StereoModel sm( camera_models[ j_cam_id ].get(),
                                      camera_models[ k_cam_id ].get() );

              if ( sm.convergence_angle( cpoint[j].position(),
                                         cpoint[k].position() ) >
                   min_convergence_angle ) {
                positions.push_back( sm( cpoint[j].position(),
                                         cpoint[k].position(),
                                       error ) );
                error_sum += error;
              }
            } catch ( camera::PixelToRayErr e ) { /* Just let it go */ }
          }
        }
      }

      // 4.2.) Summing, Averaging, and Storing
      if ( positions.empty() ) {
        vw_out(WarningMessage,"ba") << "Unable to triangulation position for point!\n";
        // At the very least we can provide a point that is some
        // distance out from the camera center and is in the 'general'
        // area.
        int j = cpoint[0].image_id();
        try {
          cpoint.set_position( camera_models[j]->camera_center(cpoint[0].position()) +
                               camera_models[j]->pixel_to_vector(cpoint[0].position())*10 );
        } catch ( camera::PixelToRayErr e ) {
          cpoint.set_position( camera_models[j]->camera_center(cpoint[0].position()) +
                               camera_models[j]->camera_pose(cpoint[0].position()).rotate(Vector3(0,0,10)) );
        }
      } else {
        error_sum /= positions.size();
        Vector3 position_avg;
        for ( unsigned j = 0; j < positions.size(); j++ )
          position_avg += positions[j]/positions.size();
        cpoint.set_position( position_avg );
      }

    }
    progress.report_finished();
  }
}
