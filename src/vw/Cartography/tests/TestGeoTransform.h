// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// TestGeoTransform.h
#include <cxxtest/TestSuite.h>

#include <vw/Cartography/GeoReference.h>
#include <vw/Cartography/GeoTransform.h>

using namespace std;
using namespace vw;
using namespace vw::cartography;

class TestGeoTransform : public CxxTest::TestSuite
{
public:

  void test_basic_transform()
  {
    GeoReference src_georef;
    src_georef.set_well_known_geogcs("WGS84");

    GeoReference dst_georef;
    dst_georef.set_well_known_geogcs("WGS84");

    GeoTransform geotx(src_georef,dst_georef);

//     TS_TRACE(stringify(geotx.forward(Vector2(0,0))));
//     TS_TRACE(stringify(geotx.reverse(Vector2(0,0))));
   }


};
