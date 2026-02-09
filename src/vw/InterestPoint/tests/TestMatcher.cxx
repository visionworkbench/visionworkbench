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


// TestMatcher.h
#include <gtest/gtest_VW.h>

#include <vw/InterestPoint/Matcher.h>
#include <vw/InterestPoint/InterestPoint.h>
#include <test/Helpers.h>

#include <algorithm>

using namespace vw;
using namespace vw::ip;

TEST( Matcher, IPComparison ) {
  std::list<InterestPoint> ip_list;
  ip_list.push_back( InterestPoint(0,0,0,3) );
  ip_list.push_back( InterestPoint(0,0,0,5) );
  ip_list.push_back( InterestPoint(0,0,0,1) );
  ip_list.push_back( InterestPoint(0,0,0,7) );
  ip_list.push_back( InterestPoint(0,0,0,9) );
  ip_list.push_back( InterestPoint(0,0,0,5.3) );

  ip_list.sort();

  EXPECT_EQ( ip_list.front().interest, 9 );
  EXPECT_EQ( ip_list.back().interest,  1 );
}

TEST( Matcher, DistanceMetric ) {
  InterestPoint ip1(0,0);
  InterestPoint ip2(5,0);
  ip1.descriptor = Vector3(0,1,0);
  ip2.descriptor = Vector3(0,5,0);

  L2NormMetric metric1;
  EXPECT_EQ(metric1(ip1, ip2), 16);

  RelativeEntropyMetric metric2;
  EXPECT_NEAR(metric2(ip1, ip2), -2.321, 1e-3);

  HammingMetric metric3;
  ip1.descriptor = Vector3(1,2,3); // Bit distance is 1, 1, 2 --> Total = 4
  ip2.descriptor = Vector3(0,0,0);
  EXPECT_NEAR(metric3(ip1, ip2), 4, 1e-3);
  
  Vector<float> data, zeros;
  data.set_size(15);
  zeros.set_size(15);
  size_t k=0;
  for (int i=0; i<5; ++i) {
    data[k++] = 1;
    data[k++] = 2;
    data[k++] = 3;
  } 
  zeros.set_all(0);
  ip1.descriptor = data; // Bit distance 4*5 = 20
  ip2.descriptor = zeros;
  EXPECT_NEAR(metric3(ip1, ip2), 20, 1e-3);
}

TEST( Matcher, Constraints ) {
  InterestPoint ip1(0,0,1.0,1.0,0.0);
  InterestPoint ip2(20,0,2.0,2.0,M_PI);
  ip1.descriptor = Vector3(0,1,0);
  ip2.descriptor = Vector3(0,5,0);

  NullConstraint null;
  ScaleOrientationConstraint scale_ori1;
  ScaleOrientationConstraint scale_ori2(0.4, 2.1, -1.2*M_PI, 1.2*M_PI);

  PositionConstraint pos1;
  PositionConstraint pos2(-20,20,-20,20);

  EXPECT_TRUE(  null(ip1,ip2) );
  EXPECT_FALSE( scale_ori1(ip1,ip2) );
  EXPECT_TRUE(  scale_ori2(ip1,ip2) );
  EXPECT_FALSE( pos1(ip1,ip2) );
  EXPECT_TRUE(  pos2(ip1,ip2) );
}

TEST( Matcher, Matcher ) {
  InterestPoint ip1a(0,0,1.0,1.0,0.0);
  InterestPoint ip2a(20,0,2.0,2.0,M_PI);
  InterestPoint ip2b(20,0,2.0,2.0,M_PI);
  InterestPoint ip2c(20,0,2.0,2.0,M_PI);
  InterestPoint ip2d(20,0,2.0,2.0,M_PI);
  InterestPoint ip2e(20,0,2.0,2.0,M_PI);
  ip1a.descriptor = Vector3(0,7.7,0);
  ip2a.descriptor = Vector3(0,5,0);
  ip2b.descriptor = Vector3(0,6,0);
  ip2c.descriptor = Vector3(0,7,0);
  ip2d.descriptor = Vector3(0,8,0);
  ip2e.descriptor = Vector3(0,9,0);

  std::vector<InterestPoint> ip1_list, ip2_list;
  std::vector<InterestPoint> matched_ip1, matched_ip2;
  std::vector<size_t> matched_indexes;

  ip1_list.push_back(ip1a);
  ip2_list.push_back(ip2a);
  ip2_list.push_back(ip2b);
  ip2_list.push_back(ip2c);
  ip2_list.push_back(ip2d);
  ip2_list.push_back(ip2e);

  std::string flann_method = "kmeans";
  InterestPointMatcher<L2NormMetric,NullConstraint> matcher(flann_method);
  matcher(ip1_list, ip2_list, matched_ip1, matched_ip2);
  matcher(ip1_list, ip2_list, matched_indexes );

  EXPECT_EQ( int(matched_ip1.size()), 1);
  EXPECT_EQ( int(matched_ip2.size()), 1);
  EXPECT_EQ( matched_indexes.size(), 1);

  EXPECT_VECTOR_NEAR( matched_ip1.begin()->descriptor,
                      Vector3( 0, 7.7, 0 ), 1e-6 );
  EXPECT_VECTOR_EQ( matched_ip2.begin()->descriptor,
                    Vector3(0,8,0) );

  EXPECT_EQ( matched_indexes[0], 3 );
}


