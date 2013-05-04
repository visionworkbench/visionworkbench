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
#include <vw/Plate/BlobManager.h>
#include <vw/Plate/Exception.h>
#include <boost/filesystem/convenience.hpp>
#include <fstream>

using namespace vw;
using namespace vw::platefile;
using namespace vw::test;
namespace fs = boost::filesystem;

class BlobManagerTest : public ::testing::Test {
  protected:
    UnlinkName blob_dir;

    virtual void SetUp() {
      blob_dir = UnlinkName("BlobDir");
      fs::create_directory(blob_dir);
      bm.reset(new BlobManager(blob_dir));
    }

    void write_to_blob(uint32 blob_id, const char* text, size_t len) {
      std::string fn(bm->name_from_id(blob_id));
      std::ofstream f(fn.c_str(), std::ios::binary|std::ios::trunc|std::ios::out);
      ASSERT_TRUE(f.is_open()) << "Failed to open " << fn;
      ASSERT_TRUE(f.good());
      f.write(text, len);
      ASSERT_TRUE(f.good());
      f.close();
      ASSERT_TRUE(f.good());
    }

    boost::shared_ptr<BlobManager> bm;
};


TEST_F(BlobManagerTest, Basic) {
  ASSERT_EQ( 0, bm->num_blobs() );

  // Test basic lock with one blob_id
  uint32 blob_id = bm->request_lock();
  ASSERT_EQ( 0, blob_id );
  bm->release_lock(blob_id);

  // Test basic lock after release
  blob_id = bm->request_lock();
  ASSERT_EQ( 0, blob_id );
  bm->release_lock(blob_id);
}

TEST_F(BlobManagerTest, TwoRequests) {
  // Test basic lock with two outstanding requests
  uint32 blob_id1 = bm->request_lock();
  uint32 blob_id2 = bm->request_lock();
  ASSERT_EQ( 0, blob_id1 );
  ASSERT_EQ( 1, blob_id2 );
  bm->release_lock(blob_id1);
  bm->release_lock(blob_id2);
}

TEST_F(BlobManagerTest, Sizes) {
  ASSERT_EQ( 0, bm->num_blobs() );
  ASSERT_THROW(bm->blob_size(0), ArgumentErr);

  uint32 blob_id = bm->request_lock();
  ASSERT_EQ(0, bm->blob_size(blob_id));
  bm->release_lock(blob_id);

  ASSERT_EQ(0, bm->blob_size(blob_id));

  blob_id = bm->request_lock();
  write_to_blob(blob_id, "rawr!", 5);
  ASSERT_NO_FATAL_FAILURE();

  ASSERT_EQ(0, bm->blob_size(blob_id));
  bm->release_lock(blob_id);
  ASSERT_EQ(5, bm->blob_size(blob_id));
}

TEST_F(BlobManagerTest, Scan) {
  ASSERT_EQ( 0, bm->num_blobs() );

  uint32 id1 = bm->request_lock(),
         id2 = bm->request_lock(),
         id3 = bm->request_lock();

  write_to_blob(id1, "abcdefg", 2);
  ASSERT_NO_FATAL_FAILURE();
  write_to_blob(id2, "abcdefg", 4);
  ASSERT_NO_FATAL_FAILURE();
  write_to_blob(id3, "abcdefg", 6);
  ASSERT_NO_FATAL_FAILURE();

  bm->release_lock(id1);
  bm->release_lock(id2);
  bm->release_lock(id3);

  ASSERT_EQ(3, bm->num_blobs());

  EXPECT_EQ(2, bm->blob_size(id1));
  EXPECT_EQ(4, bm->blob_size(id2));
  EXPECT_EQ(6, bm->blob_size(id3));

  // Make sure we scanned the directory properly
  bm.reset(new BlobManager(blob_dir));
  ASSERT_EQ(3, bm->num_blobs());

  EXPECT_EQ(2, bm->blob_size(id1));
  EXPECT_EQ(4, bm->blob_size(id2));
  EXPECT_EQ(6, bm->blob_size(id3));
}
