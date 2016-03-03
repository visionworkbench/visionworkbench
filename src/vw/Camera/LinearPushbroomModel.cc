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

#include <vw/Camera/LinearPushbroomModel.h>

using namespace vw;
using namespace camera;

LinearPushbroomModel::LinearPushbroomModel( double  scan_duration,
                                            int     number_of_lines,
                                            int     samples_per_line,
                                            int     sample_offset,
                                            double  focal_length,
                                            double  along_scan_pixel_size,
                                            double  across_scan_pixel_size,
                                            Vector3 pointing_vec,
                                            Vector3 u_vec,
                                            Quaternion<double> const& camera_pose,
                                            Vector3 const& initial_position,
                                            Vector3 const& velocity_vector) :
LinescanModel<LinearPositionInterpolation,
              ConstantPoseInterpolation>::LinescanModel(number_of_lines,
                                                        samples_per_line,
                                                        sample_offset,
                                                        focal_length,
                                                        along_scan_pixel_size,
                                                        across_scan_pixel_size,
                                                        scan_duration / number_of_lines,
                                                        pointing_vec, u_vec,
                                                        LinearPositionInterpolation(initial_position, velocity_vector),
                                                        ConstantPoseInterpolation(camera_pose)) {}

LinearPushbroomModel::~LinearPushbroomModel() {}
std::string LinearPushbroomModel::type() const { return "LinearPushbroom"; }
