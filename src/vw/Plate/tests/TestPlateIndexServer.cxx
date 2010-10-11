// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// IMPORTANT NOTE ABOUT THIS UNIT TEST:
//
// This unit test is customized to test against a plateindex server
// that has the TrueMarble.16km.2700x1350_toast.plate platefile.
// Since this platefile is not distributed with VW, these tests may
// not be useful to you.  Rather, they are mostly for internal
// development and rough performance testing, so most of these tests
// will be disabled by default.  Re-enable them if you want to do some
// plateindex_server testing. -mbroxton


#include <vw/Plate/Index.h>

using namespace std;
using namespace vw;
using namespace vw::platefile;

class TestPlateIndexServer : public CxxTest::TestSuite {

public:

  void test_basic_open() {
    std::string url = "pf://index/TrueMarble.16km.2700x1350_toast.plate";
    boost::shared_ptr<Index> idx = Index::construct_open(url);

    // Verify that index header information is correct.
    IndexHeader hdr = idx->index_header();
    TS_ASSERT_EQUALS(hdr.type(), "toast");
    TS_ASSERT_EQUALS(hdr.version(), 2);
    TS_ASSERT_EQUALS(hdr.tile_size(), 256u);
    TS_ASSERT_EQUALS(hdr.tile_filetype(), "png");
    TS_ASSERT_EQUALS(hdr.pixel_format(), VW_PIXEL_RGBA);
    TS_ASSERT_EQUALS(hdr.channel_type(), VW_CHANNEL_UINT8);

    // Test out the max depth RPC.  This is a simple one that ought to work.
    int max_depth = idx->max_depth();
    TS_ASSERT_EQUALS(max_depth, 3);

    // Test basic transaction logic.  This is a little bit more complex.
    std::vector<TileHeader> empty_tileheader_list;

    int start_cursor = idx->transaction_cursor();
    int transaction_id = idx->transaction_request("Test transaction", empty_tileheader_list);
    idx->transaction_complete(transaction_id);
    int finish_cursor = idx->transaction_cursor();
    TS_ASSERT_EQUALS(start_cursor + 1, finish_cursor);
  }

  void test_read() {
    std::string url = "pf://index/TrueMarble.16km.2700x1350_toast.plate";
    boost::shared_ptr<Index> idx = Index::construct_open(url);

    // Read the top level tile information & make sure it's valid.
    IndexRecord rec;
    for (int i = 0; i < 10; ++i) {
      rec = idx->read_request(0,0,0,-1);
      TS_ASSERT_DIFFERS(rec.blob_id(), -1);
    }

    // Read a non-existant tile to make sure it throws an error.
    TS_ASSERT_THROWS( rec = idx->read_request(100,0,0,-1), TileNotFoundErr);

    for (int i = 0; i < 100; ++i) {
      int depth = idx->max_depth();
      TS_ASSERT_EQUALS(depth, 3);
      rec = idx->read_request(0,0,0,-1);
      TS_ASSERT_DIFFERS(rec.blob_id(), -1);
    }
  }

  // void test_create() {
  //   IndexHeader hdr;
  //   hdr.set_tile_size(256);
  //   hdr.set_tile_filetype("tif");
  //   hdr.set_pixel_format(VW_PIXEL_RGBA);
  //   hdr.set_channel_type(VW_CHANNEL_UINT8);
  //   hdr.set_type("toast");
  //   hdr.set_description("");

  //   std::string bad_url = "pf://index/TrueMarble.16km.2700x1350_toast.plate";
  //   boost::shared_ptr<Index> bad_idx;
  //   TS_ASSERT_THROWS( bad_idx = Index::construct_create(bad_url, hdr), PlatefileCreationErr);

  //   std::string url = "pf://index/test.plate";
  //   boost::shared_ptr<Index> idx = Index::construct_create(url, hdr);

  //   // Verify that index header information is correct.
  //   IndexHeader hdr2 = idx->index_header();
  //   TS_ASSERT_EQUALS(hdr2.type(), "toast");
  //   TS_ASSERT_EQUALS(hdr2.version(), 2);
  //   TS_ASSERT_EQUALS(hdr2.tile_size(), 256);
  //   TS_ASSERT_EQUALS(hdr2.tile_filetype(), "tif");
  //   TS_ASSERT_EQUALS(hdr2.pixel_format(), VW_PIXEL_RGBA);
  //   TS_ASSERT_EQUALS(hdr2.channel_type(), VW_CHANNEL_UINT8);
  // }


}; // class TestPlateIndexServer
