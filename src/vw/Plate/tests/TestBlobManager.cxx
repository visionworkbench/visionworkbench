// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <gtest/gtest.h>
#include <vw/Plate/Exception.h>
#include <vw/Plate/BlobManager.h>

using namespace vw;
using namespace vw::platefile;

class BlobManagerTest : public ::testing::Test {
  protected:
  virtual void SetUp() {
    // Initial blobs = 1, max_blobs = 4;
    bm.reset(new BlobManager(1024, 1, 4));
  }
  boost::shared_ptr<BlobManager> bm;
};


TEST_F(BlobManagerTest, Basic) {
  EXPECT_EQ( 1u, bm->num_blobs() );

  uint64 offset;
  // Test basic lock with one blob_id
  int blob_id = bm->request_lock(offset);
  EXPECT_EQ( 0, blob_id );
  bm->release_lock(blob_id, offset);

  // Test basic lock after release
  blob_id = bm->request_lock(offset);
  EXPECT_EQ( 0, blob_id );
  bm->release_lock(blob_id, offset);
}

TEST_F(BlobManagerTest, TwoRequests) {
  uint64 offset1, offset2;
  // Test basic lock with two outstanding requests
  int blob_id1 = bm->request_lock(offset1);
  int blob_id2 = bm->request_lock(offset2);
  EXPECT_EQ( 0, blob_id1 );
  EXPECT_EQ( 1, blob_id2 );
  bm->release_lock(blob_id1, offset1);
  bm->release_lock(blob_id2, offset2);
}

TEST_F(BlobManagerTest, MaxRequests) {
  uint64 offset;
  // Max out the number of blobs.
  (void)bm->request_lock(offset);
  (void)bm->request_lock(offset);
  (void)bm->request_lock(offset);
  (void)bm->request_lock(offset);

  EXPECT_THROW((void)bm->request_lock(offset), BlobLimitErr);
}

TEST_F(BlobManagerTest, AddToOffset) {
  uint64 offset1, offset2;
  int blob1, blob2;

  blob1 = bm->request_lock(offset1);
  offset1+=1024;
  bm->release_lock(blob1, offset1);

  blob2 = bm->request_lock(offset2);
  ASSERT_EQ(blob1, blob2);
  EXPECT_EQ(offset1, offset2);
  bm->release_lock(blob2, offset2);
}
