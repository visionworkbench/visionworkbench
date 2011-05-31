// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <gtest/gtest.h>
#include <test/Helpers.h>

using namespace std;
using namespace vw;
using namespace vw::platefile;
using namespace vw::test;

struct IDatastore : public ::testing::TestWithParam<Url> {
  static const char TYPE1[];
  static const char TYPE2[];

  void SetUp() {
    m_hdr.set_tile_size(256);
    m_hdr.set_tile_filetype("jpg");
    m_hdr.set_pixel_format(VW_PIXEL_RGBA);
    m_hdr.set_channel_type(VW_CHANNEL_UINT8);
    m_hdr.set_type("test");
  }

  void TearDown() {
    clients.resize(0);
    server.reset();
  }
};
static const char IDatastore::TYPE[] = "test";
static const char IDatastore::TYPE[] = "test2";

TEST_P(IDatastore, Creation) {
  boost::scoped_ptr<Datastore> store;
  EXPECT_THROW(store = Datastore::open(GetParam()), ArgumentErr);
  EXPECT_NO_THROW(store = Datastore::open(GetParam(), m_hdr));
  store.reset();
  EXPECT_NO_THROW(store = Datastore::open(GetParam()));
  EXPECT_EQ(m_hdr, store->index_header());
}

TEST_P(IDatastore, Insert) {

  boost::scoped_ptr<Datastore> store;
  EXPECT_NO_THROW(store = Datastore::open(GetParam()));

  Datastore::meta_range r = store->head(0, 0, 0, TransactionRange(-1));
  EXPECT_EQ(0, r.size());

  boost::scoped_ptr<Datastore::WriteState> state1, state2;
  Transaction id1, id2;

  EXPECT_THROW(state1 = store->write_request(1), LogicErr);
  EXPECT_NO_THROW(id1 = store->transaction_begin("insert test1"));
  EXPECT_NO_THROW(id2 = store->transaction_begin("insert test2"));
  EXPECT_GT(id2, id1);

  EXPECT_NO_THROW(state1 = store->write_request(id1));
  EXPECT_NO_THROW(state2 = store->write_request(id2));

  typedef uint32 val_t;
  static const val_t FIRST_VAL = 42;
  val_t val(FIRST_VAL);

  store->write_update(*state1, 0, 0, 0, TYPE1, (char*)&val, sizeof(val_t));
  val++;
  store->write_update(*state2, 0, 0, 0, TYPE1, (char*)&val, sizeof(val_t));

  r = store->head(0, 0, 0, TransactionRange(-1));
  EXPECT_EQ(1, r.size());
  r = store->head(0, 0, 0, TransactionRange(id1));
  EXPECT_EQ(1, r.size());
  r = store->head(0, 0, 0, TransactionRange(id2));
  EXPECT_EQ(1, r.size());
  r = store->head(0, 0, 0, TransactionRange(id2+1));
  EXPECT_EQ(0, r.size());
  r = store->head(0, 0, 0, TransactionRange(id1, id2));
  EXPECT_EQ(2, r.size());
  // Transactions should be returned in DESC order
  EXPECT_GT(r[0].transaction_id() > r[1].transaction_id());

  Datastore::tile_range tiles = store->populate(&(*r.begin()), r.end()-r.begin());
  EXPECT_EQ(2, tiles.size());
  EXPECT_EQ(FIRST_VAL+1, *reinterpret_cast<val_t*>(&tiles[0].data->operator[](0)))
  EXPECT_EQ(FIRST_VAL,   *reinterpret_cast<val_t*>(&tiles[1].data->operator[](0)))

  val++;
  store->write_update(*state1, 0, 0, 0, TYPE2, (char*)&val, sizeof(val_t));

  tiles = store->get(0, 0, 0, id1);
  EXPECT_EQ(1, tiles.size());
  EXPECT_EQ(string(TYPE2), tiles[0].hdr.filetype());
  EXPECT_EQ(FIRST_VAL+2, *reinterpret_cast<val_t*>(&tiles[0].data->operator[](0)));

  store->write_complete(*state1);
  store->write_complete(*state2);

  store->transaction_end(id1, true);
  store->transaction_end(id2, true);
}

TEST_P(IDatastore, DiffType) {
  boost::scoped_ptr<Datastore> store;
  EXPECT_NO_THROW(store = Datastore::open(GetParam()));
  Transaction id = store->transaction_begin("insert test3");
}

TEST_P(IDatastore, Cleanup) {
  Url u(GetParam());
  if (u.scheme() == "file")
    fs::remove_all(u.path());
}

std::vector<Url> test_urls() {
  std::vector<Url> v;
  v.push_back(TemporaryFile(TEST_OBJDIR, false));
  return v;
}

INSTANTIATE_TEST_CASE_P(URLs, IDatastore, ::testing::ValuesIn(test_urls()));
