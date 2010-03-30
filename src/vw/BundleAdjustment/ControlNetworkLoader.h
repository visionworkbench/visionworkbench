#ifndef __VW_BUNDLEADJUSTMENT_CONTROL_NETWORK_LOADER_H__
#define __VW_BUNDLEADJUSTMENT_CONTROL_NETWORK_LOADER_H__

#include <vw/BundleAdjustment/ControlNetwork.h>
#include <vw/BundleAdjustment/CameraRelation.h>
#include <vw/Camera/CameraModel.h>

namespace vw {
namespace ba {

  void build_control_network( ControlNetwork& cnet,
                              std::vector<boost::shared_ptr<camera::CameraModel> > const& camera_models,
                              std::vector<std::string> const& image_files,
                              int min_matches = 30 );

  void add_ground_control_points( ControlNetwork& cnet,
                                  std::vector<std::string> const& image_files,
                                  std::vector<std::string> const& gcp_files );

}}

#endif//__VW_BUNDLEADJUSTMENT_CONTROL_NETWORK_LOADER_H__
