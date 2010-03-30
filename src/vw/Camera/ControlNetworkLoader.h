#ifndef __VW_CAMERA_CONTROL_NETWORK_LOADER_H__
#define __VW_CAMERA_CONTROL_NETWORK_LOADER_H__

#include <vw/Camera/ControlNetwork.h>
#include <vw/Camera/CameraRelation.h>
#include <vw/Camera/CameraModel.h>

namespace vw {
namespace camera {

void build_control_network( ControlNetwork& cnet,
                            std::vector<boost::shared_ptr<vw::camera::CameraModel> > const& camera_models,
                            std::vector<std::string> const& image_files,
                            int min_matches = 30 );

}}

#endif//__VW_CAMERA_CONTROL_NETWORK_LOADER_H__
