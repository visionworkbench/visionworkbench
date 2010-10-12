// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <gtest/gtest.h>

#include <sstream>
#include <vw/BundleAdjustment/ControlNetworkLoader.h>

#include <test/Helpers.h>

using namespace vw;
using namespace vw::ba;
using namespace vw::test;

TEST( ControlNetworkLoad, LoadingGCPNET ) {
  std::vector<std::string> image_names;
  image_names.push_back("image1");
  image_names.push_back("image2");

  // Creating GCP net
  UnlinkName gcpnet("test.cnet");
  std::list<std::string> gcpnets;
  gcpnets.push_back(gcpnet);
  {
    ControlNetwork net("test");
    ControlPoint cp;
    cp.set_position( Vector3(1,2,3) );
    {
      ControlMeasure cm(4,5,1,1,0);
      cm.set_serial("image1");
      cp.add_measure(cm);
    }
    {
      ControlMeasure cm(7,1,2,2,1);
      cm.set_serial("image2");
      cp.add_measure(cm);
    }
    net.add_control_point(cp);
    net.write_binary(gcpnet);
  }

  ControlNetwork net("destination");
  {
    ControlPoint cp;
    cp.set_position( Vector3(9,3,4) );
    {
      ControlMeasure cm(5,5,1,1,0);
      cp.add_measure(cm);
    }
    {
      ControlMeasure cm(5,7,1,1,1);
      cp.add_measure(cm);
    }
    net.add_control_point(cp);
  }
  EXPECT_TRUE( boost::filesystem::exists( gcpnets.front() ) );
  EXPECT_EQ( 1u, net.size() );
  EXPECT_EQ( 1u, gcpnets.size() );

  add_ground_control_cnets( net, image_names,
                            gcpnets.begin(), gcpnets.end() );

  ASSERT_EQ( 2u, net.size() );
  EXPECT_EQ( ControlPoint::GroundControlPoint, net[1].type() );
}
