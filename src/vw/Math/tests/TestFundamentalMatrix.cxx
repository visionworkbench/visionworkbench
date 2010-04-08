// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <gtest/gtest.h>
#include <test/Helpers.h>
#include <vw/Math/Vector.h>
#include <vw/Math/Matrix.h>
#include <vw/Math/FundamentalMatrix.h>

using namespace vw;
using namespace vw::math;

TEST(Fundamental, 7pFittingFunctor) {
  std::vector<Vector3> ip1, ip2;
  ip1.push_back(Vector3(793,796,1)); // These come from Z's potato
  ip2.push_back(Vector3(187,678,1)); // chip example
  ip1.push_back(Vector3(1647,984,1));
  ip2.push_back(Vector3(974,1016,1));
  ip1.push_back(Vector3(1452,1109,1)); // 3
  ip2.push_back(Vector3(753,1100,1));
  ip1.push_back(Vector3(942,396,1)); // 4
  ip2.push_back(Vector3(394,322,1));
  ip1.push_back(Vector3(681,463,1)); // 5
  ip2.push_back(Vector3(139,345,1));
  ip1.push_back(Vector3(1252,166,1)); // 6
  ip2.push_back(Vector3(755,154,1));
  ip1.push_back(Vector3(977,516,1)); // 7
  ip2.push_back(Vector3(401,441,1));

  Fundamental7FittingFunctor cow;
  Matrix<double> fundie = cow(ip1,ip2);
  EXPECT_EQ(fundie.cols(),3);
  EXPECT_EQ(fundie.rows(),3);
  EXPECT_NEAR( det(fundie), 0, 1e-10 );
  EXPECT_EQ( rank(fundie), 2 );
  for ( size_t i = 0; i < cow.num_solutions(); i++ ) {
    fundie = cow.fundamental_matrix( i );
    EXPECT_NEAR( det(fundie), 0, 1e-10 );
  }
}

TEST(Fundamental, 8pFittingFunctor) {
  std::vector<Vector3> ip1, ip2;
  ip1.push_back(Vector3(971, 605,1)); // 1
  ip1.push_back(Vector3(765, 632,1)); // 2
  ip1.push_back(Vector3(953, 617,1)); // 3
  ip1.push_back(Vector3(984, 606,1)); // 4
  ip1.push_back(Vector3(941, 595,1)); // 5
  ip1.push_back(Vector3(1005, 609,1)); // 6
  ip1.push_back(Vector3(866, 644,1)); // 7
  ip1.push_back(Vector3(728, 635,1)); // 8

  ip2.push_back(Vector3(1318, 661,1));
  ip2.push_back(Vector3(1108, 685,1));
  ip2.push_back(Vector3(1300, 673,1));
  ip2.push_back(Vector3(1331, 662,1));
  ip2.push_back(Vector3(1290, 651,1));
  ip2.push_back(Vector3(1352, 666,1));
  ip2.push_back(Vector3(1213, 698));
  ip2.push_back(Vector3(1076, 685,1));

  Fundamental8FittingFunctor cow;
  Matrix3x3 F = cow(ip1,ip2);
  EXPECT_EQ(F.cols(),3);
  EXPECT_EQ(F.rows(),3);
  EXPECT_NEAR( det(F), 0, 1e-10 );
  EXPECT_EQ( rank(F), 2 );
}
