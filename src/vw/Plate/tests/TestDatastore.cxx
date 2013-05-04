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
#include <vw/Plate/Datastore.h>
#include <vw/FileIO/TemporaryFile.h>
#include <boost/filesystem/operations.hpp>
namespace fs = boost::filesystem;

using namespace std;
using namespace vw;
using namespace vw::platefile;
using namespace vw::test;

struct IDatastore : public ::testing::TestWithParam<std::string> {
  boost::scoped_ptr<TemporaryDir> m_tmpdir;
  Url m_url;
  IndexHeader m_hdr;

  void SetUp() {
    m_hdr.set_tile_size(256);
    m_hdr.set_tile_filetype("jpg");
    m_hdr.set_pixel_format(VW_PIXEL_RGBA);
    m_hdr.set_channel_type(VW_CHANNEL_UINT8);
    m_hdr.set_type("test");

    m_tmpdir.reset(new TemporaryDir(TEST_OBJDIR));
    m_url.scheme(GetParam());
    m_url.path(m_tmpdir->filename() + "/test.plate");
  }

  void TearDown() {
    fs::remove_all(m_url.path());
  }
};

namespace {
  typedef uint32 val_t;

  static const char TYPE1[] = "test";
  static const char TYPE2[] = "test2";
  static const val_t vA = 6765, vB = 10946, vC = 17711, vD = 28657, vE = 46368, vF = 75025, vG = 121393;

}

TEST_P(IDatastore, Creation) {
  boost::scoped_ptr<Datastore> store;
  ASSERT_THROW(store.reset(Datastore::open(m_url)), ArgumentErr);
  ASSERT_NO_THROW(store.reset(Datastore::open(m_url, m_hdr)));
  store.reset();
  ASSERT_NO_THROW(store.reset(Datastore::open(m_url)));

  IndexHeader hdr = store->index_header();
  EXPECT_EQ(m_hdr.tile_size(),     hdr.tile_size());
  EXPECT_EQ(m_hdr.tile_filetype(), hdr.tile_filetype());
  EXPECT_EQ(m_hdr.pixel_format(),  hdr.pixel_format());
  EXPECT_EQ(m_hdr.channel_type(),  hdr.channel_type());
  EXPECT_EQ(m_hdr.type(),          hdr.type());
  EXPECT_EQ(m_hdr.description(),   hdr.description());

  EXPECT_EQ(0, hdr.num_levels());
  EXPECT_EQ(0, hdr.transaction_read_cursor());
  EXPECT_EQ(1, hdr.transaction_write_cursor());
}

bool SortTilesByTidASC(const Tile& a, const Tile& b) {
  return a.hdr.transaction_id() < b.hdr.transaction_id();
}

TEST_P(IDatastore, TwoInsert) {
  boost::scoped_ptr<Datastore> store;
  ASSERT_NO_THROW(store.reset(Datastore::open(m_url, m_hdr)));

  Datastore::TileSearch r;
  store->head(r, 0, 0, 0, TransactionRange(-1));
  EXPECT_EQ(0, r.size());

  boost::scoped_ptr<WriteState> state1, state2;
  Transaction id1(0), id2(0);

  // ideally this would throw, but not sure how to detect if a transaction has never been started.
  // ASSERT_THROW(state1.reset(store->write_request(1)), LogicErr);
  ASSERT_NO_THROW(id1 = store->transaction_begin("insert test1"));
  ASSERT_NO_THROW(id2 = store->transaction_begin("insert test2"));
  EXPECT_GT(id2, id1);

  ASSERT_NO_THROW(state1.reset(store->write_request(id1)));
  ASSERT_NO_THROW(state2.reset(store->write_request(id2)));

  store->write_update(*state1, 0, 0, 0, TYPE1, reinterpret_cast<const uint8*>(&vA), sizeof(val_t));
  store->write_update(*state2, 0, 0, 0, TYPE1, reinterpret_cast<const uint8*>(&vB), sizeof(val_t));

  store->head(r, 0, 0, 0, TransactionRange(-1));       EXPECT_EQ(1, r.size());
  store->head(r, 0, 0, 0, TransactionRange(id1));      EXPECT_EQ(1, r.size());
  store->head(r, 0, 0, 0, TransactionRange(id2));      EXPECT_EQ(1, r.size());
  store->head(r, 0, 0, 0, TransactionRange(id2+1));    EXPECT_EQ(0, r.size());
  store->head(r, 0, 0, 0, TransactionRange(id1, id2)); EXPECT_EQ(2, r.size());

  // Transactions should be returned in DESC order
  ASSERT_GT(r[0].hdr.transaction_id(), r[1].hdr.transaction_id());

  store->populate(r);
  EXPECT_EQ(2, r.size());
  vector<Tile> ret(r.begin(), r.end());
  sort(ret.begin(), ret.end(), SortTilesByTidASC);

  EXPECT_EQ(vA, *reinterpret_cast<val_t*>(&ret[0].data->operator[](0)));
  EXPECT_EQ(vB, *reinterpret_cast<val_t*>(&ret[1].data->operator[](0)));

  store->write_update(*state1, 0, 0, 0, TYPE2, reinterpret_cast<const uint8*>(&vC), sizeof(val_t));

  store->get(r, 0, 0, 0, TransactionRange(id1));
  EXPECT_EQ(1, r.size());
  EXPECT_EQ(string(TYPE2), r[0].hdr.filetype());
  EXPECT_EQ(vC, *reinterpret_cast<val_t*>(&r[0].data->operator[](0)));

  store->write_complete(*state1);
  store->write_complete(*state2);

  store->transaction_end(id1, true);
  store->transaction_end(id2, true);
}

TEST_P(IDatastore, InsertRegion) {
  boost::scoped_ptr<Datastore> store;
  ASSERT_NO_THROW(store.reset(Datastore::open(m_url, m_hdr)));

  Datastore::TileSearch r;
  store->head(r, 0, 0, 0, TransactionRange(-1));
  EXPECT_EQ(0, r.size());

  boost::scoped_ptr<WriteState> state;
  Transaction id(0);

  ASSERT_NO_THROW(id = store->transaction_begin("insert test3"));
  ASSERT_NO_THROW(state.reset(store->write_request(id)));

  // |A| | | |   (level 2) (row,col)
  // | |B|C| |   (0,0) (1,1) (1,2)
  // | |D| |E|   (2,1) (2,3)
  // | | |F| |   (3,2)
  store->write_update(*state, 2, 0, 0, TYPE1, reinterpret_cast<const uint8*>(&vA), sizeof(val_t));
  store->write_update(*state, 2, 1, 1, TYPE1, reinterpret_cast<const uint8*>(&vB), sizeof(val_t));
  store->write_update(*state, 2, 1, 2, TYPE1, reinterpret_cast<const uint8*>(&vC), sizeof(val_t));
  store->write_update(*state, 2, 2, 1, TYPE1, reinterpret_cast<const uint8*>(&vD), sizeof(val_t));
  store->write_update(*state, 2, 2, 3, TYPE1, reinterpret_cast<const uint8*>(&vE), sizeof(val_t));
  store->write_update(*state, 2, 3, 2, TYPE1, reinterpret_cast<const uint8*>(&vF), sizeof(val_t));

  store->head(r, 0, 0, 0, TransactionRange(-1)); EXPECT_EQ(0, r.size());

  #define CHECK(lvl, zone, count)\
    store->head(r, lvl, zone, TransactionRange(-1));    EXPECT_EQ(count, r.size());\
    store->head(r, lvl, zone, TransactionRange(0, id)); EXPECT_EQ(count, r.size());

  CHECK(2, BBox2u(0,0,1,1), 1);
  CHECK(2, BBox2u(1,1,1,1), 1);
  CHECK(2, BBox2u(1,2,1,1), 1);
  CHECK(2, BBox2u(2,1,1,1), 1);
  CHECK(2, BBox2u(2,3,1,1), 1);
  CHECK(2, BBox2u(3,2,1,1), 1);
  CHECK(2, BBox2u(0,1,1,1), 0);
  CHECK(2, BBox2u(3,3,1,1), 0);

  CHECK(2, BBox2u(0,0,2,2), 2);
  CHECK(2, BBox2u(0,1,2,2), 2);
  CHECK(2, BBox2u(0,2,2,2), 1);
  CHECK(2, BBox2u(1,0,2,2), 2);
  CHECK(2, BBox2u(1,1,2,2), 3);
  CHECK(2, BBox2u(1,2,2,2), 2);
  CHECK(2, BBox2u(2,0,2,2), 1);
  CHECK(2, BBox2u(2,1,2,2), 2);
  CHECK(2, BBox2u(2,2,2,2), 2);

  CHECK(2, BBox2u(0,0,3,3), 4);
  CHECK(2, BBox2u(0,1,3,3), 4);
  CHECK(2, BBox2u(1,0,3,3), 4);
  CHECK(2, BBox2u(1,1,3,3), 5);

  store->write_complete(*state);
  store->transaction_end(id, true);
}

TEST_P(IDatastore, Logging) {
  boost::scoped_ptr<Datastore> store;
  ASSERT_NO_THROW(store.reset(Datastore::open(m_url, m_hdr)));

  store->audit_log()() << "I can't really check that this worked, but I can at least call it" << std::endl;
  store->error_log()() << "This should produce output" << std::endl;
}

std::vector<string> test_urls() {
  std::vector<string> v;
  v.push_back("file");
  v.push_back("dir");
  return v;
}

INSTANTIATE_TEST_CASE_P(URLs, IDatastore, ::testing::ValuesIn(test_urls()));
