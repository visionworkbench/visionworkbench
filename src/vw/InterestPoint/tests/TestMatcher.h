// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <cxxtest/TestSuite.h>
#include <vw/InterestPoint/Matcher.h>
#include <vw/InterestPoint/InterestData.h>

using namespace vw;
using namespace vw::ip;

class TestMatcher : public CxxTest::TestSuite
{
  public:

  void test_distance_metric() {
    InterestPoint ip1(0,0);
    InterestPoint ip2(5,0);
    ip1.descriptor = Vector3(0,1,0);
    ip2.descriptor = Vector3(0,5,0);

    L2NormMetric metric1;
    TS_ASSERT_EQUALS(metric1(ip1, ip2), 16);

    RelativeEntropyMetric metric2;
    TS_ASSERT_DELTA(metric2(ip1, ip2), -2.321, 0.001);
  }

  void test_constraints() {
    InterestPoint ip1(0,0,1.0,1.0,0.0);
    InterestPoint ip2(20,0,2.0,2.0,M_PI);
    ip1.descriptor = Vector3(0,1,0);
    ip2.descriptor = Vector3(0,5,0);

    NullConstraint null;
    ScaleOrientationConstraint scale_ori1;
    ScaleOrientationConstraint scale_ori2(0.4, 2.1, -1.2*M_PI, 1.2*M_PI);

    PositionConstraint pos1;
    PositionConstraint pos2(-20,20,-20,20);

    TS_ASSERT_EQUALS(null(ip1,ip2), true);
    TS_ASSERT_EQUALS(scale_ori1(ip1,ip2), false);
    TS_ASSERT_EQUALS(scale_ori2(ip1,ip2), true);
    TS_ASSERT_EQUALS(pos1(ip1,ip2), false);
    TS_ASSERT_EQUALS(pos2(ip1,ip2), true);
  }


  void test_matcher() {
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

    ip1_list.push_back(ip1a);
    ip2_list.push_back(ip2a);
    ip2_list.push_back(ip2b);
    ip2_list.push_back(ip2c);
    ip2_list.push_back(ip2d);
    ip2_list.push_back(ip2e);

    InterestPointMatcher<L2NormMetric,NullConstraint> matcher;
    matcher(ip1_list, ip2_list, matched_ip1, matched_ip2);

    TS_ASSERT_EQUALS(int(matched_ip1.size()), 1);
    TS_ASSERT_EQUALS(int(matched_ip2.size()), 1);

    TS_ASSERT_EQUALS( matched_ip1.begin()->descriptor[0], 0);
    TS_ASSERT_DELTA( matched_ip1.begin()->descriptor[1], 7.7,0.01);
    TS_ASSERT_EQUALS( matched_ip1.begin()->descriptor[2], 0);

    TS_ASSERT_EQUALS( matched_ip2.begin()->descriptor[0], 0);
    TS_ASSERT_EQUALS( matched_ip2.begin()->descriptor[1], 8);
    TS_ASSERT_EQUALS( matched_ip2.begin()->descriptor[2], 0);
  }
};
