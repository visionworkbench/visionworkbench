// __BEGIN_LICENSE__
//
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
//
// Copyright 2006 Carnegie Mellon University. All rights reserved.
//
// This software is distributed under the NASA Open Source Agreement
// (NOSA), version 1.3.  The NOSA has been approved by the Open Source
// Initiative.  See the file COPYING at the top of the distribution
// directory tree for the complete NOSA document.
//
// THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
// KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
// LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
// SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
// A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
// THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
// DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
//
// __END_LICENSE__

// TestGeoReference.h
#include <cxxtest/TestSuite.h>

#include <vw/Cartography/PointImageManipulation.h>

using namespace std;
using namespace vw;
using namespace vw::cartography;

class TestPointImageManipulation : public CxxTest::TestSuite
{
public:

  void test_xyz_to_lon_lat_conversions()
  {
    // Test a full forward and reverse transformation
    Vector3 xyz(-2197110.000000, 1741355.875000, 1898886.875000);
    Vector3 lon_lat_alt = xyz_to_lon_lat_radius(xyz);
    Vector3 xyz2 = lon_lat_radius_to_xyz(lon_lat_alt);

    TS_ASSERT_DELTA(xyz(0), xyz2(0), 1e-2);
    TS_ASSERT_DELTA(xyz(1), xyz2(1), 1e-2);
    TS_ASSERT_DELTA(xyz(2), xyz2(2), 1e-2);

    // Test to see if things still work for West positive coordinate
    // systems.
    lon_lat_alt = xyz_to_lon_lat_radius(xyz,false);
    xyz2 = lon_lat_radius_to_xyz(lon_lat_alt,false);

    TS_ASSERT_DELTA(xyz(0), xyz2(0), 1e-2);
    TS_ASSERT_DELTA(xyz(1), xyz2(1), 1e-2);
    TS_ASSERT_DELTA(xyz(2), xyz2(2), 1e-2);
  }

};

