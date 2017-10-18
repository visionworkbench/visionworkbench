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
#include <vw/Cartography/Datum.h>
#include <boost/filesystem/operations.hpp>
#include <boost/foreach.hpp>

/// \file ControlNetworkLoader.h Functions for generating control networks

namespace vw {
namespace ba {

  /// Builds a control network using given camera models and original
  /// image names. This function uses Boost::FS to then find match files
  /// that would have been created by 'ipmatch' by searching the entire
  /// permutation of the image_files vector.
  VW_API bool build_control_network(bool triangulate_points,
                             ControlNetwork& cnet,
                             std::vector<boost::shared_ptr<camera::CameraModel> >
			     const& camera_models,
                             std::vector<std::string> const& image_files,
                             std::map< std::pair<int, int>, std::string> const& match_files,
                             size_t min_matches,
                             double min_angle_radians);
  
  /// Recomputes the world location of a point based on camera observations.
  /// - Returns the mean triangulation error.
  VW_API double triangulate_control_point(ControlPoint& cp,
				   std::vector<boost::shared_ptr<camera::CameraModel> >
				   const& camera_models,
				   double const& min_angle_radians );

  /// Adds ground control points from GCP files to an already built
  /// Control Network. The vector image_files serves as a look up chart
  /// for relating image names in GCP files to CNET's internal
  /// indexing.  Each GCP is a line in the file, containing the point
  /// id, 3D point (as lat,lon,height_above_datum), its sigmas, then,
  /// for each image, the image file name, pixel measurements, and their sigmas.
  VW_API void add_ground_control_points(ControlNetwork& cnet,
				 std::vector<std::string> const& image_files,
				 std::vector<std::string> const& gcp_files,
				 cartography::Datum const& datum);
    
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
