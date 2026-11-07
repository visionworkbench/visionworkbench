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

#include <iostream>
#include <sstream>

#include <vw/BundleAdjustment/CameraRelation.h>
#include <vw/BundleAdjustment/ControlNetwork.h>

#include <test/Helpers.h>

using namespace vw;

class CircleTest : public ::testing::Test {
protected:
  CircleTest() : cnet("TestCNET") {}

  virtual void SetUp() {
    cnet = ControlNetwork("TestCNET");

    for (uint32 i = 0; i < 4; i++) {
      ControlPoint cpoint;

      for (uint32 j = 0; j < i+1; j++) {
        ControlMeasure cm(100*j+i, 100*j+i, 1, 1, j);
        cpoint.add_measure(cm);
      }
      cnet.add_control_point(cpoint);
    }
  }

  ControlNetwork cnet;
};

TEST(CameraRelation, Construction) {
  CameraRelationNetwork<JFeature> crn;
  crn.add_node(CameraNode<JFeature>(0,"monkey"));
  crn.add_node(CameraNode<JFeature>(1,"dog"));
  ASSERT_EQ(crn.size(), 2u);

  {
    // JFeature(point_id, camera_id); set the pixel location separately.
    boost::shared_ptr<JFeature> first(new JFeature(0, 0));
    first->m_location = Vector2f(5, 10);
    boost::shared_ptr<JFeature> second(new JFeature(1, 1));
    second->m_location = Vector2f(10, 5);

    first->connection(second);
    second->connection(first);
    crn[0].relations.push_back(first);
    crn[1].relations.push_back(second);
  }

  std::list<boost::weak_ptr<JFeature>> list;
  (*crn[0].begin())->list_connections(list);
  for (std::list<boost::weak_ptr<JFeature>>::iterator it = list.begin();
       it != list.end(); it++) {
    std::ostringstream ostr;
    ostr << *(*it).lock();
    EXPECT_GT(ostr.str().size(), 20u);
  }
  EXPECT_EQ(list.size(), 2u);
  EXPECT_EQ(list.front().lock()->m_location.x(), 10);
  EXPECT_EQ(list.back().lock()->m_location.x(), 5);
}

TEST_F(CircleTest, IPFeature) {

  CameraRelationNetwork<JFeature> crn;
  crn.from_cnet(cnet);
  EXPECT_EQ(crn.size(), 4u);
  for (uint32 i = 0; i < crn.size(); i++)
    EXPECT_EQ(crn[i].id, i);
  EXPECT_EQ(crn[0].relations.size(), 4u);
  EXPECT_EQ(crn[1].relations.size(), 3u);
  EXPECT_EQ(crn[2].relations.size(), 2u);
  EXPECT_EQ(crn[3].relations.size(), 1u);
  std::list<boost::weak_ptr<JFeature>> list;
  crn[3].relations.front()->list_connections(list);
  EXPECT_EQ(list.size(), 4u);
}

TEST_F(CircleTest, JFeature) {

  // Convert to Camera Relation
  CameraRelationNetwork<JFeature> crn;
  crn.from_cnet(cnet);
  EXPECT_EQ(crn.size(), 4u);
  for (uint32 i = 0; i < crn.size(); i++)
    EXPECT_EQ(crn[i].id, i);
  EXPECT_EQ(crn[0].relations.size(), 4u);
  EXPECT_EQ(crn[1].relations.size(), 3u);
  EXPECT_EQ(crn[2].relations.size(), 2u);
  EXPECT_EQ(crn[3].relations.size(), 1u);
  std::list<boost::weak_ptr<JFeature>> list;
  crn[3].relations.front()->list_connections(list);
  EXPECT_EQ(list.size(), 4u);

}
