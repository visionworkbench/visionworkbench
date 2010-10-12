// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <gtest/gtest.h>
#include <vw/FileIO.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/ImageIO.h>
#include <test/Helpers.h>

using namespace vw;
using namespace vw::test;

#if defined(VW_HAVE_PKG_GDAL) && VW_HAVE_PKG_GDAL==1

#include <vw/FileIO/DiskImageResourceGDAL.h>

TEST( GDALFeatures, NoDataValue ) {

  ImageView<float> monkey(2,2);
  monkey(0,0) = 0;
  monkey(0,1) = 0.5;
  monkey(1,0) = 0.75;
  monkey(1,1) = 0.8;

  {
    DiskImageResourceGDAL w_rsrc( "nodata.tif", monkey.format() );
    w_rsrc.set_nodata_value(-1);
    EXPECT_NO_THROW( w_rsrc.set_nodata_value( -1 ) );
    write_image( w_rsrc, monkey );
  }

  // Read it back in and confirm
  ImageView<float> monkey_return;
  DiskImageResourceGDAL r_rsrc( "nodata.tif" );
  read_image( monkey_return, r_rsrc );
  ASSERT_TRUE( r_rsrc.has_nodata_value() );
  EXPECT_EQ(  r_rsrc.nodata_value(), -1 );
}

#endif
