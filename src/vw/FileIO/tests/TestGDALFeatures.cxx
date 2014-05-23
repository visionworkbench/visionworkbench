// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NASA Vision Workbench is licensed under the Apache License,
//  Version 2.0 (the "License"); you may not use this file except in
//  compliance with the License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__


#include <gtest/gtest_VW.h>
#include <vw/Image/ImageIO.h>
#include <vw/Image/ImageView.h>
#include <vw/FileIO/DiskImageResourceGDAL.h>
#include <test/Helpers.h>
#include <vw/config.h>

using namespace vw;
using namespace vw::test;

#if defined(VW_HAVE_PKG_GDAL) && VW_HAVE_PKG_GDAL==1

#include <vw/FileIO/DiskImageResourceGDAL.h>


TEST( GDALFeatures, NoDataValue ) {
  UnlinkName nodata("nodata.tif");

  ImageView<float> monkey(2,2);
  monkey(0,0) = 0;
  monkey(0,1) = 0.5;
  monkey(1,0) = 0.75;
  monkey(1,1) = 0.8;

  {
    DiskImageResourceGDAL w_rsrc( nodata, monkey.format() );
    ASSERT_TRUE( w_rsrc.has_nodata_write() );
    EXPECT_NO_THROW( w_rsrc.set_nodata_write(-1) );
    write_image( w_rsrc, monkey );
  }

  // Read it back in and confirm
  ImageView<float> monkey_return;
  DiskImageResourceGDAL r_rsrc( nodata );
  read_image( monkey_return, r_rsrc );
  ASSERT_TRUE( r_rsrc.has_nodata_read() );
  EXPECT_EQ( -1, r_rsrc.nodata_read() );
}

#endif
