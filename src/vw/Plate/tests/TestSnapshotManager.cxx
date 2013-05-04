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
#include <test/Helpers.h>
#include <boost/lexical_cast.hpp>
#include <vw/Image/ImageMath.h>
#include <vw/Plate/PlateFile.h>
#include <vw/Plate/SnapshotManager.h>

using namespace vw;
using namespace vw::platefile;
using namespace vw::test;

class SnapshotTest : public ::testing::Test {
protected:
  typedef PixelGrayA<uint8> PixelT;

  virtual void SetUp() {
    platename = UnlinkName("test.plate");
    platefile =
      boost::shared_ptr<PlateFile>( new PlateFile( Url(platename),
                                                   "equi", "", 256, "png",
                                                   VW_PIXEL_GRAYA,
                                                   VW_CHANNEL_UINT8) );

    // Generate synthetic data .. this is going to make 6 layers of a
    // 1024 px by 1024 px image.
    level = 3;
    bbox = BBox2i(0,0,4,4);
    for ( size_t ch = 0; ch < 6; ch++ ) {
      ImageView<PixelT> input(256,256);
      fill(input, PixelT());
      fill(crop(input,ch,0,2,256), PixelT(ch));

      platefile->transaction_begin("Writing ch:" +
                                   boost::lexical_cast<std::string>(ch),-1);
      platefile->write_request();
      for ( int32 i = bbox.min().x(); i < bbox.max().x(); i++ ) {
        for ( int32 j = bbox.min().y(); j < bbox.max().y(); j++ ) {
          if ( ch == 0 )
            fill(crop(input,ch,0,2,256), PixelT(i + j*bbox.width()));
          platefile->write_update(input,i,j,level);
        }
      }
      platefile->write_complete();
      platefile->transaction_end(true);
    }

    objective.set_size(256,256);
    fill(objective, PixelT() );
    inc_objective.set_size(256,256);
    fill(inc_objective, PixelT() );
    fill(crop(inc_objective,0,0,1,256), PixelT(1,1) );
    fill(crop(objective,1,0,1,256), PixelT(1) );
    fill(crop(objective,2,0,1,256), PixelT(2) );
    fill(crop(objective,3,0,1,256), PixelT(3) );
    fill(crop(objective,4,0,1,256), PixelT(4) );
    fill(crop(objective,5,0,2,256), PixelT(5) );
  }

  ImageView<PixelT> objective;
  ImageView<PixelT> inc_objective;
  UnlinkName platename;
  boost::shared_ptr<PlateFile> platefile;
  uint32 level;
  BBox2i bbox;
};

TEST_F( SnapshotTest, FastOperation ) {
  // Test the test framwork, make sure I understand.
  EXPECT_EQ( 4, platefile->num_levels() );

  std::list<TileHeader> query =
    platefile->search_by_location(0,0,level,TransactionRange(0,100));
  EXPECT_EQ( 6, query.size() );
  query = platefile->search_by_location(3,3,level,TransactionRange(0,100));
  EXPECT_EQ( 6, query.size() );

  // Perform actual test
  SnapshotManager<PixelT> sm(platefile, platefile, false);

  platefile->transaction_begin("Write snapshot",-1);
  platefile->write_request();
  sm.snapshot(level,bbox,TransactionRange(0,6));
  platefile->write_complete();
  platefile->transaction_end(false);

  // Check the results
  query =
    platefile->search_by_location(0,0,level,TransactionRange(0,100));
  EXPECT_EQ( 7, query.size() );
  for ( uint32 i = bbox.min().x(); i < bbox.max().x(); i++ ) {
    for ( uint32 j = bbox.min().y(); j < bbox.max().y(); j++ ) {
      ImageView<PixelT> tile;
      platefile->read(tile,i,j,level,7,true);
      EXPECT_EQ( objective + inc_objective * PixelT(i + j*bbox.width()), tile );
    }
  }
}

TEST_F( SnapshotTest, SmallCache ) {
  uint64 cache_size_before = vw_settings().system_cache_size();

  // Setting the cache size to 7 tiles
  vw_settings().set_system_cache_size( 7 * 256 * 256 * uint32(PixelNumBytes<PixelT>::value));

  SnapshotManager<PixelT> sm(platefile, platefile, false);

  platefile->transaction_begin("Write snapshot",-1);
  platefile->write_request();
  sm.snapshot(level,bbox,TransactionRange(0,6));
  platefile->write_complete();
  platefile->transaction_end(false);

  // Check the results
  std::list<TileHeader> query =
    platefile->search_by_location(0,0,level,TransactionRange(0,100));
  EXPECT_EQ( 7, query.size() );
  for ( uint32 i = bbox.min().x(); i < bbox.max().x(); i++ ) {
    for ( uint32 j = bbox.min().y(); j < bbox.max().y(); j++ ) {
      ImageView<PixelT> tile;
      platefile->read(tile,i,j,level,7,true);
      EXPECT_EQ( objective + inc_objective * PixelT(i + j*bbox.width()), tile );
    }
  }

  vw_settings().set_system_cache_size(cache_size_before);
}

TEST_F( SnapshotTest, TooSmallCache ) {
  uint64 cache_size_before = vw_settings().system_cache_size();

  // Setting the cache size to 3 tiles
  vw_settings().set_system_cache_size( 3 * 256 * 256 * uint32(PixelNumBytes<PixelT>::value));

  SnapshotManager<PixelT> sm(platefile, platefile, false);

  platefile->transaction_begin("Write snapshot",-1);
  platefile->write_request();
  EXPECT_THROW( sm.snapshot(level,bbox,TransactionRange(0,6)),
                LogicErr );
  platefile->write_complete();
  platefile->transaction_end(false);

  vw_settings().set_system_cache_size(cache_size_before);
}
