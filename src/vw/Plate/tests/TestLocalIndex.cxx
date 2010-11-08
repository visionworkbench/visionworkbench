// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <gtest/gtest.h>
#include <test/Helpers.h>

#include <vw/Plate/LocalIndex.h>
#include <vw/Plate/Exception.h>
#include <vw/Plate/Blob.h>

#include <boost/filesystem/convenience.hpp>
namespace fs = boost::filesystem;

using namespace std;
using namespace vw;
using namespace vw::platefile;
using namespace vw::test;

class LocalIndexTest : public ::testing::Test {
  protected:

  virtual void SetUp() {
    // Set up a sane index default
    index_hdr.set_tile_size(256);
    index_hdr.set_tile_filetype("tiff");
    index_hdr.set_pixel_format(VW_PIXEL_RGB);
    index_hdr.set_channel_type(VW_CHANNEL_UINT8);

    // And a sample tile header to clone
    tile_hdr.set_filetype("tiff");
    tile_hdr.set_col(1);
    tile_hdr.set_row(2);
    tile_hdr.set_level(3);
    tile_hdr.set_transaction_id(4);

    // Create a platefile to play in
    plate_path = UnlinkName("LocalIndex");

    fs::create_directories(plate_path);

    blob_path  = plate_path + "/0.blob";

    // Now construct the index for it
    index.reset(new LocalIndex(plate_path, index_hdr));
    blob.reset(new Blob(blob_path));
  }

  virtual void TearDown() {
    // Don't care about the rest (destructor will get them), but make sure the
    // index and blob go away before the name does.
    blob.reset();
    index.reset();
  }

  IndexHeader index_hdr;
  TileHeader tile_hdr;

  UnlinkName plate_path;

  std::string blob_path;

  boost::shared_ptr<LocalIndex> index;
  boost::shared_ptr<Blob> blob;
};

class LocalIndexTiles : public LocalIndexTest {
  protected:
  virtual void SetUp() {
    LocalIndexTest::SetUp();
    hdrs.reset(new TileHeader[6]);

    TileHeader hdr2 = tile_hdr;
    hdr2.set_col(0);
    hdr2.set_row(0);

    std::fill(hdrs.get(), hdrs.get()+6, hdr2);

    hdrs[0].set_level(0);

    hdrs[1].set_level(1);

    hdrs[2].set_col(1);
    hdrs[2].set_level(1);

    hdrs[3].set_row(1);
    hdrs[3].set_level(1);

    hdrs[4].set_col(1);
    hdrs[4].set_row(1);
    hdrs[4].set_level(1);

    hdrs[5] = hdrs[4];

    test_size = 20;

    test_data.reset(new uint8[test_size]);
    for (size_t i = 0; i < test_size; ++i)
      test_data[i] = i;
  }

  void index_write(const TileHeader &hdr, IndexRecord &rec) {
    uint64 old_offset;
    rec.set_blob_id( index->write_request(old_offset) );
    rec.set_blob_offset(blob->write(hdr, test_data, test_size));
    index->write_update(hdr, rec);
    index->write_complete(rec.blob_id(), rec.blob_offset());
  }

#define check_tile_hdr(expected, actual) do {\
  SCOPED_TRACE("");\
  check_tile_hdr_(expected, actual);\
} while(0)

  void check_tile_hdr_(const TileHeader& expected, const TileHeader& actual) {
    EXPECT_EQ( expected.col(),      actual.col() );
    EXPECT_EQ( expected.row(),      actual.row() );
    EXPECT_EQ( expected.level(),    actual.level() );
    EXPECT_EQ( expected.filetype(), actual.filetype() );
  }


  boost::shared_array<TileHeader> hdrs;

  size_t test_size;
  boost::shared_array<uint8> test_data;
};

TEST(LocalIndex, BasicAccess) {
  UnlinkName file("index");
  std::vector<boost::shared_ptr<IndexLevel> > levels;

  for (int i = 0; i < 5; ++i) {
    boost::shared_ptr<LocalPageGeneratorFactory> page_gen_factory;
    page_gen_factory.reset( new LocalPageGeneratorFactory(file) );
    boost::shared_ptr<IndexLevel> level( new IndexLevel(page_gen_factory, i, 256, 256, 100) );
    levels.push_back(level);
  }
}

TEST(LocalIndex, IndexRecord) {

  UnlinkName name("foo.bar");

  // Test serialization/deserialization
  TileHeader write_hdr;
  write_hdr.set_filetype("tiff");
  write_hdr.set_col(1);
  write_hdr.set_row(1);
  write_hdr.set_level(2);

  std::ofstream ostr(name.c_str(), std::ios::binary);
  write_hdr.SerializeToOstream(&ostr);
  ostr.close();

  std::ifstream istr(name.c_str(), std::ios::binary);
  TileHeader read_hdr;
  read_hdr.ParseFromIstream(&istr);
  istr.close();

  EXPECT_EQ(write_hdr.col(), read_hdr.col());
  EXPECT_EQ(write_hdr.row(), read_hdr.row());
  EXPECT_EQ(write_hdr.level(), read_hdr.level());
  EXPECT_EQ(write_hdr.filetype(), read_hdr.filetype());
}

TEST_F(LocalIndexTest, WriteRead) {

  // Now read the index back in
  LocalIndex index2(plate_path);

  EXPECT_EQ(index->version(),       index2.version());
  EXPECT_EQ(index->tile_size(),     index2.tile_size());
  EXPECT_EQ(index->tile_filetype(), index2.tile_filetype());
}

TEST_F(LocalIndexTest, Transactions) {
  // Store the default transaction id
  int tx0 = index->transaction_cursor();

  // 0 means nothing has happened
  EXPECT_EQ(0, tx0);

  // Start a bunch of transactions
  int tx1 = index->transaction_request("Test transaction #1", -1);
  int tx2 = index->transaction_request("Test transaction #2", -1);
  int tx3 = index->transaction_request("Test transaction #3", -1);
  int tx4 = index->transaction_request("Test transaction #4", -1);
  int tx5 = index->transaction_request("Test transaction #5", -1);
  int tx6 = index->transaction_request("Test transaction #6", -1);

  // Make sure we didn't move the transaction cursor
  EXPECT_EQ( tx0, index->transaction_cursor() );

  // According to the current scheme, each tx id is ++
  EXPECT_EQ(   1, tx1 );
  EXPECT_EQ(   2, tx2 );
  EXPECT_EQ(   3, tx3 );
  EXPECT_EQ(   4, tx4 );
  EXPECT_EQ(   5, tx5 );
  EXPECT_EQ(   6, tx6 );

  // Now complete a transaction, but not the next one in the series.
  index->transaction_complete(tx6, false);
  EXPECT_EQ( tx0, index->transaction_cursor() );

  // Now we complete the next one. the cursor should move forward by one.
  index->transaction_complete(tx1, true);
  EXPECT_EQ( tx1, index->transaction_cursor() );

  // Complete more transactions in reverse order
  index->transaction_complete(tx3, true);
  EXPECT_EQ( tx3, index->transaction_cursor() );
  index->transaction_complete(tx2, true);
  EXPECT_EQ( tx3, index->transaction_cursor() );

  // And the rest
  index->transaction_complete(tx4, false);
  EXPECT_EQ( tx3, index->transaction_cursor() );
  index->transaction_complete(tx5, true);
  EXPECT_EQ( tx5, index->transaction_cursor() );
}

TEST_F(LocalIndexTest, SkipTransactions) {
  // Store the default transaction id
  int tx0 = index->transaction_cursor();

  // 0 means nothing has happened
  EXPECT_EQ(0, tx0);

  // Start a bunch of transactions, including one with an override ID.
  int tx1 = index->transaction_request("Test transaction #1", -1);
  int tx2 = index->transaction_request("Test transaction #2", 13);
  int tx3 = index->transaction_request("Test transaction #2", -1);

  EXPECT_EQ(   1, tx1 );
  EXPECT_EQ(   13, tx2 );
  EXPECT_EQ(   14, tx3 );
}

TEST_F(LocalIndexTiles, BasicReadWrite) {

  IndexRecord in, out;

  index_write(tile_hdr, in);
  out = index->read_request(tile_hdr.col(), tile_hdr.row(), tile_hdr.level(), tile_hdr.transaction_id(), true);

  EXPECT_EQ(in.blob_id()    , out.blob_id());
  EXPECT_EQ(in.blob_offset(), out.blob_offset());
}


TEST_F(LocalIndexTiles, Simple) {

  IndexRecord rec;
  index_write(hdrs[0], rec);
  EXPECT_LE(test_size, rec.blob_offset());

  index_write(hdrs[1], rec);
  index_write(hdrs[2], rec);
  index_write(hdrs[3], rec);
  index_write(hdrs[4], rec);

  index_write(hdrs[1], rec);

  EXPECT_LE(test_size * 6, rec.blob_offset());

  IndexRecord result = index->read_request(hdrs[1].col(), hdrs[1].row(), hdrs[1].level(), -1);
  EXPECT_EQ( rec.blob_id(),     result.blob_id() );
  EXPECT_EQ( rec.blob_offset(), result.blob_offset() );

  // Now let's try some invalid reads/writes
  hdrs[5].set_col(10);
  EXPECT_THROW(index_write(hdrs[5], rec), TileNotFoundErr);

  hdrs[5].set_col(0);
  hdrs[5].set_row(2);
  hdrs[5].set_level(1);
  EXPECT_THROW(index_write(hdrs[5], rec), TileNotFoundErr);

  EXPECT_THROW(index->read_request(0, 0, 2, -1), TileNotFoundErr);
}

TEST_F(LocalIndexTiles, ReadWrite) {

  {
    // Write some data to the Index.
    IndexRecord rec;
    index_write(hdrs[0], rec);
    index_write(hdrs[1], rec);
    index_write(hdrs[2], rec);
    index_write(hdrs[3], rec);
    index_write(hdrs[4], rec);
  }

  // Now, let's save the data to disk, and then read it back.
  {
    index.reset(new LocalIndex(plate_path));
    blob.reset(new Blob(blob_path));

    // Read the data back from the index
    IndexRecord out_rec;
    TileHeader  out_hdr;
    boost::shared_array<uint8> out_data;

    uint64 read_size;

    out_rec  = index->read_request(0, 0, 0, -1);
    out_hdr  = blob->read_header<TileHeader>(out_rec.blob_offset());
    out_data = blob->read_data(out_rec.blob_offset(), read_size);
    ASSERT_RANGE_EQ(&test_data[0], &test_data[test_size], &out_data[0], &out_data[read_size]);
    check_tile_hdr(hdrs[0], out_hdr);

    out_rec  = index->read_request(0, 0, 1, -1);
    out_hdr  = blob->read_header<TileHeader>(out_rec.blob_offset());
    out_data = blob->read_data(out_rec.blob_offset(), read_size);
    ASSERT_RANGE_EQ(&test_data[0], &test_data[test_size], &out_data[0], &out_data[read_size]);
    check_tile_hdr(hdrs[1], out_hdr);

    out_rec = index->read_request(1, 0, 1, -1);
    out_hdr = blob->read_header<TileHeader>(out_rec.blob_offset());
    out_data = blob->read_data(out_rec.blob_offset(), read_size);
    ASSERT_RANGE_EQ(&test_data[0], &test_data[test_size], &out_data[0], &out_data[read_size]);
    check_tile_hdr(hdrs[2], out_hdr);

    out_rec = index->read_request(0, 1, 1, -1);
    out_hdr = blob->read_header<TileHeader>(out_rec.blob_offset());
    out_data = blob->read_data(out_rec.blob_offset(), read_size);
    ASSERT_RANGE_EQ(&test_data[0], &test_data[test_size], &out_data[0], &out_data[read_size]);
    check_tile_hdr(hdrs[3], out_hdr);

    out_rec = index->read_request(1, 1, 1, -1);
    out_hdr = blob->read_header<TileHeader>(out_rec.blob_offset());
    out_data = blob->read_data(out_rec.blob_offset(), read_size);
    ASSERT_RANGE_EQ(&test_data[0], &test_data[test_size], &out_data[0], &out_data[read_size]);
    check_tile_hdr(hdrs[4], out_hdr);
  }
}

TEST_F(LocalIndexTiles, ValidTiles) {

  IndexRecord rec;
  index_write(hdrs[0], rec);
  index_write(hdrs[1], rec);
  index_write(hdrs[2], rec);
  index_write(hdrs[3], rec);
  index_write(hdrs[4], rec);

  MinMaxAccumulator<int64> tid;
  tid(hdrs[0].transaction_id());
  tid(hdrs[1].transaction_id());
  tid(hdrs[2].transaction_id());
  tid(hdrs[3].transaction_id());
  tid(hdrs[4].transaction_id());


  // Test valid_tiles() call at level 0;
  std::list<TileHeader> tiles = index->search_by_region(0, BBox2i(0,0,1,1), tid.minimum(), tid.maximum(), 1, false);
  EXPECT_EQ(1, tiles.size());

  tiles = index->search_by_region(1, BBox2i(0,0,2,2), tid.minimum(), tid.maximum(), 1, false);
  EXPECT_EQ(4, tiles.size());
}
