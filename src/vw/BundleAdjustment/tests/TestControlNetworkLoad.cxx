// __BEGIN_LICENSE__
//  Copyright (c) 2006-2026, United States Government as represented by the
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


#include <gtest/gtest_VW.h>

#include <sstream>
#include <vw/BundleAdjustment/ControlNetworkLoader.h>

#include <test/Helpers.h>

using namespace vw;
using namespace vw::test;

// Build a control network with a ground control point and verify the point
// count, its type, and its measures. The old GCP-network file loader this once
// exercised was removed. GCP file loading itself is covered by GCPLoad below.
TEST(ControlNetworkLoad, GroundControlPoint) {

  ControlNetwork net("test");
  net.add_image_name("image1");
  net.add_image_name("image2");

  ControlPoint cp(ControlPoint::GroundControlPoint);
  cp.set_position(Vector3(1,2,3));
  cp.add_measure(ControlMeasure(4,5,1,1,0));
  cp.add_measure(ControlMeasure(7,1,2,2,1));
  net.add_control_point(cp);

  ASSERT_EQ(1u, net.size());
  EXPECT_EQ(ControlPoint::GroundControlPoint, net[0].type());
  EXPECT_EQ(2u, net[0].size()); // two measures
}

TEST(GCPLoad, LoadingGCP) {

  ControlNetwork cnet("testing");
  std::vector<std::string> gcp_files;
  vw::Datum datum;

  cnet.add_image_name("170lo.tif");
  cnet.add_image_name("171lo.tif");
  cnet.add_image_name("172lo.tif");
  cnet.add_image_name("173lo.tif");
  cnet.add_image_name("174lo.tif");
  gcp_files.push_back("sample.gcp");

  vw::add_ground_control_points(cnet, gcp_files, datum);
  EXPECT_EQ(27, cnet.size());
}
