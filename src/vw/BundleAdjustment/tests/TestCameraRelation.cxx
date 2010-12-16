// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <gtest/gtest.h>

#include <iostream>
#include <sstream>

#include <vw/BundleAdjustment/CameraRelation.h>
#include <vw/BundleAdjustment/ControlNetwork.h>

#include <test/Helpers.h>

using namespace vw;
using namespace vw::ip;
using namespace vw::ba;

class CircleTest : public ::testing::Test {
protected:
  CircleTest() : cnet("TestCNET") {}

  virtual void SetUp() {
    cnet = ControlNetwork( "TestCNET" );

    for ( uint32 i = 0; i < 4; i++ ) {
      ControlPoint cpoint;

      for ( uint32 j = 0; j < i+1; j++ ) {
        ControlMeasure cm( 100*j+i, 100*j+i, 1, 1, j );
        cpoint.add_measure( cm );
      }
      cnet.add_control_point( cpoint );
    }
  }

  ControlNetwork cnet;
};

TEST( CameraRelation, Construction ) {
  CameraRelationNetwork<IPFeature> crn;
  crn.add_node( CameraNode<IPFeature>(0,"monkey") );
  crn.add_node( CameraNode<IPFeature>(1,"dog") );
  ASSERT_EQ( crn.size(), 2u );

  {
    boost::shared_ptr<IPFeature> first( new IPFeature( InterestPoint(5,10),
                                                       0 ) );
    boost::shared_ptr<IPFeature> second( new IPFeature( InterestPoint(10,5),
                                                        1) );

    first->connection( second );
    second->connection( first );
    crn[0].relations.push_back( first );
    crn[1].relations.push_back( second );
  }

  std::list<boost::weak_ptr<IPFeature> > list;
  (*crn[0].begin())->list_connections( list );
  for ( std::list<boost::weak_ptr<IPFeature> >::iterator it = list.begin();
        it != list.end(); it++ ) {
    std::ostringstream ostr;
    ostr << *(*it).lock();
    EXPECT_GT( ostr.str().size(), 20u );
  }
  EXPECT_EQ( list.size(), 2u );
  EXPECT_EQ( list.front().lock()->m_ip.x, 10 );
  EXPECT_EQ( list.back().lock()->m_ip.x, 5 );
}

TEST_F( CircleTest, IPFeature ) {

  // Convert to Camera Relation
  CameraRelationNetwork<IPFeature> crn;
  crn.read_controlnetwork( cnet );
  EXPECT_EQ( crn.size(), 4u );
  for ( uint32 i = 0; i < crn.size(); i++ )
    EXPECT_EQ( crn[i].id, i );
  EXPECT_EQ( crn[0].relations.size(), 4u );
  EXPECT_EQ( crn[1].relations.size(), 3u );
  EXPECT_EQ( crn[2].relations.size(), 2u );
  EXPECT_EQ( crn[3].relations.size(), 1u );
  std::list<boost::weak_ptr<IPFeature> > list;
  crn[3].relations.front()->list_connections( list );
  EXPECT_EQ( list.size(), 4u );

  // Converting back to Control Network
  ControlNetwork cnet2("Copy");
  crn.write_controlnetwork( cnet2 );
  ASSERT_EQ( cnet2.size(), cnet.size() );
  for ( unsigned i = 0; i < cnet2.size(); i++ ) {
    ASSERT_EQ( cnet2[i].size(), cnet[i].size() );
    for ( unsigned j = 0; j < cnet2[i].size(); j++ ) {
      EXPECT_EQ( cnet2[i][j], cnet[i][j] );
    }
  }

}

TEST_F( CircleTest, JFeature ) {

  // Convert to Camera Relation
  CameraRelationNetwork<JFeature> crn;
  crn.read_controlnetwork( cnet );
  EXPECT_EQ( crn.size(), 4u );
  for ( uint32 i = 0; i < crn.size(); i++ )
    EXPECT_EQ( crn[i].id, i );
  EXPECT_EQ( crn[0].relations.size(), 4u );
  EXPECT_EQ( crn[1].relations.size(), 3u );
  EXPECT_EQ( crn[2].relations.size(), 2u );
  EXPECT_EQ( crn[3].relations.size(), 1u );
  std::list<boost::weak_ptr<JFeature> > list;
  crn[3].relations.front()->list_connections( list );
  EXPECT_EQ( list.size(), 4u );

  // Converting back to Control Network
  ControlNetwork cnet2("Copy");
  crn.write_controlnetwork( cnet2 );
  ASSERT_EQ( cnet2.size(), cnet.size() );
  for ( unsigned i = 0; i < cnet2.size(); i++ ) {
    ASSERT_EQ( cnet2[i].size(), cnet[i].size() );
    for ( unsigned j = 0; j < cnet2[i].size(); j++ ) {
      EXPECT_EQ( cnet2[i][j], cnet[i][j] );
    }
  }

}


