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


#ifndef __VW_BUNDLEADJUSTMENT_CONTROL_NETWORK_LOADER_H__
#define __VW_BUNDLEADJUSTMENT_CONTROL_NETWORK_LOADER_H__

#include <vw/BundleAdjustment/ControlNetwork.h>
#include <vw/BundleAdjustment/CameraRelation.h>
#include <vw/Camera/CameraModel.h>
#include <vw/Cartography/SimplePointImageManipulation.h>

#include <boost/filesystem/operations.hpp>
#include <boost/foreach.hpp>

/// \file ControlNetworkLoader.h Functions for generating control networks

namespace vw {
namespace ba {

  // Builds a control network using given camera models and original
  // image names. This function uses Boost::FS to then find match files
  // that would have been created by 'ipmatch' by searching the entire
  // permutation of the image_files vector.
  void build_control_network( ControlNetwork& cnet,
                              std::vector<boost::shared_ptr<camera::CameraModel> > const& camera_models,
                              std::vector<std::string> const& image_files,
                              size_t min_matches = 30,
                              std::vector<std::string> const& directories = std::vector<std::string>(1,"."),
                              std::string const& prefix = ""
                              );

  void triangulate_control_point( ControlPoint& cp,
                                  std::vector<boost::shared_ptr<camera::CameraModel> > const& camera_models,
                                  double const& minimum_angle );

  // Adds ground control points from individual GCP files to an
  // already built Control Network. The vector image_files serves as a look
  // up chart for relating image names in GCP files to CNET's internal
  // indexing.
  template <class IterT>
  void add_ground_control_points( ControlNetwork& cnet,
                                  std::vector<std::string> const& image_files,
                                  IterT gcp_start, IterT gcp_end ) {
    namespace fs = boost::filesystem;

    // Creating a version of image_files that doesn't contain the path
    typedef std::map<std::string,size_t> LookupType;
    LookupType image_lookup;
    for (size_t i = 0; i < image_files.size(); i++ ) {
      image_lookup[image_files[i]] = i;
      image_lookup[fs::path(image_files[i]).filename().string()] = i;
    }

    while ( gcp_start != gcp_end ) {
      if ( !fs::exists( *gcp_start ) ) {
        gcp_start++;
        continue;
      }

      // Data to be loaded
      std::vector<Vector4> measure_locations;
      std::vector<std::string> measure_cameras;
      Vector3 world_location, world_sigma;

      vw_out(VerboseDebugMessage,"ba") << "\tLoading \"" << *gcp_start
                                       << "\".\n";
      int count = 0;
      std::ifstream ifile( (*gcp_start).c_str() );
      while (!ifile.eof()) {
        if ( count == 0 ) {
          // First line defines position in the world
          ifile >> world_location[0] >> world_location[1]
                >> world_location[2] >> world_sigma[0]
                >> world_sigma[1] >> world_sigma[2];
        } else {
          // Other lines define position in images
          std::string temp_name;
          Vector4 temp_loc;
          if (ifile >> temp_name >> temp_loc[0] >> temp_loc[1]
              >> temp_loc[2] >> temp_loc[3]){
            measure_locations.push_back( temp_loc );
            measure_cameras.push_back( temp_name );
          }
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
      std::vector<Vector4>::iterator m_iter_loc = measure_locations.begin();
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
          vw_out(WarningMessage,"ba") << "\t\tWarning: no image found matching "
                                      << *m_iter_name << std::endl;
        }
        m_iter_loc++;
        m_iter_name++;
      }

      // Appended GCP
      cnet.add_control_point(cpoint);
      gcp_start++;
    }
  }

  template <class IterT>
  void add_ground_control_cnets( ControlNetwork& cnet,
                                 std::vector<std::string> const& image_files,
                                 IterT gcpcnet_start, IterT gcpcnet_end ) {
    namespace fs = boost::filesystem;

    // Creating a version of image_files that doesn't contain the path
    typedef std::map<std::string,size_t> LookupType;
    LookupType image_lookup;
    for (size_t i = 0; i < image_files.size(); i++ ) {
      image_lookup[image_files[i]] = i;
      image_lookup[fs::path(image_files[i]).filename().string()] = i;
    }

    while ( gcpcnet_start != gcpcnet_end ) {
      if ( !fs::exists( *gcpcnet_start ) ) {
        gcpcnet_start++;
        continue;
      }

      vw_out(VerboseDebugMessage,"ba") << "\tLoading \"" << *gcpcnet_start
                                       << "\".\n";

      ControlNetwork gcpcnet("");
      gcpcnet.read_binary(*gcpcnet_start);

      BOOST_FOREACH( ControlPoint & cp, gcpcnet ) {
        bool failed_to_index = false;
        // Fixing indexing
        BOOST_FOREACH( ControlMeasure & cm, cp ) {
          LookupType::iterator it = image_lookup.find(cm.serial());
          if ( it != image_lookup.end() )
            cm.set_image_id(it->second);
          else
            vw_out(WarningMessage,"ba") << "\t\tWarning: no image found matching "
                                        << cm.serial() << std::endl;
        }

        if ( failed_to_index )
          continue;
        cp.set_type( ControlPoint::GroundControlPoint );
        cnet.add_control_point(cp);
        vw_out(DebugMessage,"ba") << "\t\tAdded GCP: " << cp.position() << "\n";
      }

      gcpcnet_start++;
    }
  }

}} //end namespace vw::ba

#endif//__VW_BUNDLEADJUSTMENT_CONTROL_NETWORK_LOADER_H__
