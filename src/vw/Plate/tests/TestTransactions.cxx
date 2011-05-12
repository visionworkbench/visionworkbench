// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <gtest/gtest.h>
#include <test/Helpers.h>
#include <vw/Plate/FundamentalTypes.h>

using namespace std;
using namespace vw;
using namespace vw::platefile;
using namespace vw::test;

TEST(Transactions, Basic) {
  EXPECT_THROW(Transaction(-1), ArgumentErr);
  EXPECT_NO_THROW(TransactionOrNeg(-1));

  Transaction a1(0), b1(1), c1(2);
  TransactionOrNeg a2(0), b2(1), d2(-1);

  EXPECT_EQ(Transaction(0), a1);
  EXPECT_EQ(Transaction(1), b1);
  EXPECT_EQ(Transaction(2), c1);

  EXPECT_EQ(TransactionOrNeg(0),  a2);
  EXPECT_EQ(TransactionOrNeg(1),  b2);
  EXPECT_EQ(TransactionOrNeg(-1), d2);
  EXPECT_EQ(TransactionOrNeg(-2), d2);

  EXPECT_NO_THROW(EXPECT_EQ(a2.promote(), a1));
  EXPECT_NO_THROW(EXPECT_EQ(b2.promote(), b1));
  EXPECT_THROW(d2.promote(), LogicErr);

  EXPECT_EQ(a1, a2);
  EXPECT_EQ(b1, b2);
  EXPECT_NE(c1, d2);

  EXPECT_EQ(a2, a1);
  EXPECT_EQ(b2, b1);
  EXPECT_NE(d2, c1);
}

TEST(Transactions, Range) {
  EXPECT_THROW(TransactionRange(-1, 0), ArgumentErr);
  EXPECT_NO_THROW(TransactionRange(-1,-1));
}
