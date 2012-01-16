// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// TestDatum.h
#include <gtest/gtest_VW.h>
#include <test/Helpers.h>

#include <vw/Cartography/Datum.h>

using namespace vw;
using namespace vw::cartography;
using namespace vw::test;


#if defined(VW_HAVE_PKG_PROTOBUF) && VW_HAVE_PKG_PROTOBUF==1
TEST( Datum, DatumDesc ) {
  Datum datum("NAD27");

  DatumDesc desc = datum.build_desc();

  Datum datum2(desc);

  EXPECT_EQ(datum.build_desc().DebugString(), datum2.build_desc().DebugString());
}
#endif
