// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// TestGeoReference.h
#include <gtest/gtest_VW.h>
#include <test/Helpers.h>

#include <vw/Cartography/PointImageManipulation.h>

using namespace vw;
using namespace vw::cartography;
using namespace vw::test;

TEST( PointImageManip, XYZ_to_LonLat ) {
  // Test a full forward and reverse transformation
  Vector3 xyz(-2197110.000000, 1741355.875000, 1898886.875000);
  Vector3 lon_lat_alt = xyz_to_lon_lat_radius(xyz);
  Vector3 xyz2 = lon_lat_radius_to_xyz(lon_lat_alt);

  EXPECT_VECTOR_NEAR( xyz, xyz2, 1e-2 );

  // Test to see if things still work for West positive coordinate
  // systems.
  lon_lat_alt = xyz_to_lon_lat_radius(xyz,false);
  xyz2 = lon_lat_radius_to_xyz(lon_lat_alt,false);

  EXPECT_VECTOR_NEAR( xyz, xyz2, 1e-2 );

  // See if it still works if using 0-360 range
  xyz[1] = -xyz[1];
  lon_lat_alt = xyz_to_lon_lat_radius(xyz,true,false);
  xyz2 = lon_lat_radius_to_xyz(lon_lat_alt);

  EXPECT_VECTOR_NEAR( xyz, xyz2, 1e-2 );
}

