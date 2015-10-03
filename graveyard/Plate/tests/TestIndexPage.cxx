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
#include <vw/Plate/detail/LocalIndex.h>
#include <vw/Plate/Exception.h>

using namespace std;
using namespace vw;
using namespace vw::platefile;
using namespace vw::platefile::detail;
using namespace vw::test;

class IndexPageTest : public ::testing::Test {
  protected:
  virtual void SetUp() {
    page_path = UnlinkName("IndexPage");
    page.reset(new LocalIndexPage(page_path,0,0,0,1024,1024));
  }
  virtual void TearDown() {
    page.reset();
  }
  UnlinkName page_path;
  boost::shared_ptr<LocalIndexPage> page;
};

TEST_F(IndexPageTest, Empty) {
  EXPECT_EQ(0, page->sparse_size());
}

TEST_F(IndexPageTest, ThrowsOnEmpty) {
  // Try accessing an empty entry.
  EXPECT_THROW( IndexRecord r = page->get(0,0,0),    TileNotFoundErr );
  EXPECT_THROW( IndexRecord r = page->get(1024,0,0), TileNotFoundErr );
  EXPECT_THROW( IndexRecord r = page->get(0,1024,0), TileNotFoundErr );
  EXPECT_THROW( IndexRecord r = page->get(-1,0,0),   TileNotFoundErr );
  EXPECT_THROW( IndexRecord r = page->get(0,-1,0),   TileNotFoundErr );
  ASSERT_EQ(0, page->sparse_size());
}

TEST_F(IndexPageTest, Basic) {

  TileHeader hdr[3];
  IndexRecord rec[4];

  // headers ordered by transaction_id
  hdr[0].set_col(3);
  hdr[0].set_row(5);
  hdr[0].set_transaction_id(1001);

  hdr[1] = hdr[0];
  hdr[1].set_transaction_id(1003);

  hdr[2] = hdr[0];
  hdr[2].set_transaction_id(2001);

  rec[0].set_blob_id(1001);
  rec[0].set_blob_offset(1002);

  rec[1].set_blob_id(1101);
  rec[1].set_blob_offset(1102);

  rec[2].set_blob_id(1201);
  rec[2].set_blob_offset(1202);

  rec[3].set_blob_id(1301);
  rec[3].set_blob_offset(1302);

  IndexRecord out_rec;

  // We're going out of order with the headers intentionally
  page->set(hdr[1], rec[0]);
  out_rec = page->get(hdr[1].col(), hdr[1].row(), hdr[1].transaction_id());

  EXPECT_EQ( rec[0].blob_id(),     out_rec.blob_id() );
  EXPECT_EQ( rec[0].blob_offset(), out_rec.blob_offset() );

  // Add a different record
  page->set(hdr[1], rec[1]);
  out_rec = page->get(hdr[1].col(), hdr[1].row(), hdr[1].transaction_id());
  EXPECT_EQ( rec[1].blob_id(),     out_rec.blob_id() );
  EXPECT_EQ( rec[1].blob_offset(), out_rec.blob_offset() );

  // Try different transaction ids, going forward and back in "time"
  page->set(hdr[0], rec[2]);
  out_rec = page->get(hdr[0].col(), hdr[0].row(), hdr[0].transaction_id());
  EXPECT_EQ( rec[2].blob_id(),     out_rec.blob_id() );
  EXPECT_EQ( rec[2].blob_offset(), out_rec.blob_offset() );

  // Get the id just before hdr[2] (which should be between hdr[1] and hdr[2]
  out_rec = page->get(hdr[2].col(), hdr[2].row(), hdr[2].transaction_id()-1);
  EXPECT_EQ( rec[1].blob_id(),     out_rec.blob_id() );
  EXPECT_EQ( rec[1].blob_offset(), out_rec.blob_offset() );

  out_rec = page->get(hdr[1].col(), hdr[1].row(), -1);
  EXPECT_EQ( rec[1].blob_id(),     out_rec.blob_id() );
  EXPECT_EQ( rec[1].blob_offset(), out_rec.blob_offset() );

  // Try and one that occurs later in the list
  page->set(hdr[2], rec[3]);
  out_rec = page->get(hdr[2].col(), hdr[2].row(), -1);
  EXPECT_EQ( rec[3].blob_id(),     out_rec.blob_id() );
  EXPECT_EQ( rec[3].blob_offset(), out_rec.blob_offset() );

  // Make sure we still only have one entry at this point.
  EXPECT_EQ(1, page->sparse_size());
  EXPECT_THROW( IndexRecord r = page->get(0,0,0), TileNotFoundErr );
}

TEST_F(IndexPageTest, Serialization) {
  TileHeader hdr[3];
  IndexRecord rec[3];

  hdr[0].set_col(3);
  hdr[0].set_row(5);
  hdr[0].set_transaction_id(1);
  hdr[1] = hdr[0];
  hdr[1].set_transaction_id(2);
  hdr[2].set_col(1023);
  hdr[2].set_row(5);
  hdr[2].set_transaction_id(2);

  // Set a few different index records
  rec[0].set_blob_id(1001);
  rec[0].set_blob_offset(1002);
  rec[1].set_blob_id(2001);
  rec[1].set_blob_offset(2002);
  rec[2].set_blob_id(3001);
  rec[2].set_blob_offset(3002);

  page->set(hdr[0],rec[0]);
  page->set(hdr[1],rec[1]);
  page->set(hdr[2],rec[2]);

  EXPECT_EQ(2, page->sparse_size());

  // Now, save the file to disk.
  page->sync();

  boost::shared_ptr<LocalIndexPage> page2(new LocalIndexPage(page_path,0,0,0,1024,1024));

  IndexRecord out_rec;

  out_rec = page2->get(hdr[0].col(), hdr[0].row(), hdr[0].transaction_id());
  EXPECT_EQ( rec[0].blob_id(),     out_rec.blob_id() );
  EXPECT_EQ( rec[0].blob_offset(), out_rec.blob_offset() );

  out_rec = page2->get(hdr[1].col(), hdr[1].row(), hdr[1].transaction_id());
  EXPECT_EQ( rec[1].blob_id(),     out_rec.blob_id() );
  EXPECT_EQ( rec[1].blob_offset(), out_rec.blob_offset() );

  out_rec = page2->get(hdr[2].col(), hdr[2].row(), hdr[2].transaction_id());
  EXPECT_EQ( rec[2].blob_id(),     out_rec.blob_id() );
  EXPECT_EQ( rec[2].blob_offset(), out_rec.blob_offset() );
}
