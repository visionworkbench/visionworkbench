// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <cxxtest/TestSuite.h>

#include <vw/Plate/Index.h>

using namespace std;
using namespace vw;
using namespace vw::platefile;

class TestIndexNode : public CxxTest::TestSuite {
  boost::shared_array<uint8> m_test_data;


public:

  TestIndexNode() {
    m_test_data = boost::shared_array<uint8>(new uint8[20]);
    for (int i = 0; i < 20; ++i) {
      m_test_data[i] = i;
    }
  }
  
  void check_data(boost::shared_array<uint8> a, boost::shared_array<uint8> b) {
    for (int i = 0; i < 20; ++i) 
      TS_ASSERT_EQUALS( a[i], b[i] );
  }


  void test_index_record() {

    // Test serialization/deserialization
    TileHeader write_hdr;
    write_hdr.set_filetype("tiff");
    write_hdr.set_col(1);
    write_hdr.set_row(1);
    write_hdr.set_depth(2);
    
    std::ofstream ostr("/tmp/foo.bar", std::ios::binary);
    write_hdr.SerializeToOstream(&ostr);
    ostr.close();

    std::ifstream istr("/tmp/foo.bar", std::ios::binary);
    TileHeader read_hdr;
    read_hdr.ParseFromIstream(&istr);
    istr.close();

    TS_ASSERT_EQUALS(write_hdr.col(), read_hdr.col());
    TS_ASSERT_EQUALS(write_hdr.row(), read_hdr.row());
    TS_ASSERT_EQUALS(write_hdr.depth(), read_hdr.depth());
    TS_ASSERT_EQUALS(write_hdr.filetype(), read_hdr.filetype());

    // Clean up
    unlink("/tmp/foo.bar");
  }

  void test_index_write_read() {
    unlink("/tmp/foo.plate/plate.index");
    std::string plate_filename = "/tmp/foo.plate";

    // Write the basic data...
    Index idx(plate_filename, 256, "tif", VW_PIXEL_RGB, VW_CHANNEL_UINT8);

    // And read it back in...
    Index idx2(plate_filename);
    
    TS_ASSERT_EQUALS(idx.version(), idx2.version());
    TS_ASSERT_EQUALS(idx.default_tile_size(), idx2.default_tile_size());
    TS_ASSERT_EQUALS(idx.default_tile_filetype(), idx2.default_tile_filetype());
  }

  void test_index_transactions() {
    unlink("/tmp/foo.plate/plate.index");
    std::string plate_filename = "/tmp/foo.plate";

    // Write the basic data...
    Index idx(plate_filename, 256, "tif", VW_PIXEL_RGB, VW_CHANNEL_UINT8);

    TS_ASSERT_EQUALS(idx.transaction_cursor(), 0);

    // Test one tranaction request
    int tx1 = idx.transaction_request("Test transaction #1");
    TS_ASSERT_EQUALS(tx1, 1);
    TS_ASSERT_EQUALS(idx.transaction_cursor(), 0);

    // And now a few more
    int tx2 = idx.transaction_request("Test transaction #2");
    int tx3 = idx.transaction_request("Test transaction #3");
    int tx4 = idx.transaction_request("Test transaction #4");
    int tx5 = idx.transaction_request("Test transaction #5");
    int tx6 = idx.transaction_request("Test transaction #6");
    TS_ASSERT_EQUALS(tx2, 2);
    TS_ASSERT_EQUALS(tx3, 3);
    TS_ASSERT_EQUALS(tx4, 4);
    TS_ASSERT_EQUALS(tx5, 5);
    TS_ASSERT_EQUALS(tx6, 6);
    TS_ASSERT_EQUALS(idx.transaction_cursor(), 0);

    // Now complet a transaction, but not the next one in the series.
    idx.transaction_complete(tx6);
    TS_ASSERT_EQUALS(idx.transaction_cursor(), 0);

    // Now we complete the next one. the cursor should move forward by one.
    idx.transaction_complete(tx1);
    TS_ASSERT_EQUALS(idx.transaction_cursor(), 1);

    // Complete more transactions in reverse order
    idx.transaction_complete(tx3);
    TS_ASSERT_EQUALS(idx.transaction_cursor(), 1);
    idx.transaction_complete(tx2);
    TS_ASSERT_EQUALS(idx.transaction_cursor(), 3);

    // And the rest
    idx.transaction_complete(tx4);
    TS_ASSERT_EQUALS(idx.transaction_cursor(), 4);
    idx.transaction_complete(tx5);
    TS_ASSERT_EQUALS(idx.transaction_cursor(), 6);
  }



  void test_simple_index() {
    std::string plate_filename = "/tmp/foo.plate";
    std::string index_filename = plate_filename + "/plate.index";
    std::string blob_filename = plate_filename + "/plate_0.blob";
    
    unlink(index_filename.c_str());
    unlink(blob_filename.c_str());

    TileHeader dummy_header0;
    dummy_header0.set_filetype("tif");
    dummy_header0.set_col(0);
    dummy_header0.set_row(0);
    dummy_header0.set_depth(0);

    TileHeader dummy_header1;
    dummy_header1.set_filetype("tif");
    dummy_header1.set_col(0);
    dummy_header1.set_row(0);
    dummy_header1.set_depth(1);

    TileHeader dummy_header2;
    dummy_header2.set_filetype("tif");
    dummy_header2.set_col(1);
    dummy_header2.set_row(0);
    dummy_header2.set_depth(1);

    TileHeader dummy_header3;
    dummy_header3.set_filetype("tif");
    dummy_header3.set_col(0);
    dummy_header3.set_row(1);
    dummy_header3.set_depth(1);

    TileHeader dummy_header4;
    dummy_header4.set_filetype("tif");
    dummy_header4.set_col(1);
    dummy_header4.set_row(1);
    dummy_header4.set_depth(1);

    // Write some data to the Index.
    Index idx(plate_filename, 256, "tif", VW_PIXEL_RGB, VW_CHANNEL_UINT8);
    Blob blob(blob_filename);

    IndexRecord rec;
    rec.set_blob_id( idx.write_request(1024) );
    rec.set_blob_offset(blob.write(dummy_header0, m_test_data, 20));
    rec.set_status(INDEX_RECORD_VALID);
    idx.write_complete(dummy_header0, rec);

    rec.set_blob_id( idx.write_request(1024) );
    rec.set_blob_offset(blob.write(dummy_header1, m_test_data, 20));
    idx.write_complete(dummy_header1, rec);

    rec.set_blob_id( idx.write_request(1024) );
    rec.set_blob_offset(blob.write(dummy_header2, m_test_data, 20));
    idx.write_complete(dummy_header2, rec);

    rec.set_blob_id( idx.write_request(1024) );
    rec.set_blob_offset(blob.write(dummy_header3, m_test_data, 20));
    idx.write_complete(dummy_header3, rec);

    rec.set_blob_id( idx.write_request(1024) );
    rec.set_blob_offset(blob.write(dummy_header4, m_test_data, 20));
    idx.write_complete(dummy_header4, rec);

    // Test re-writing (i.e. changing) an entry.
    TileHeader dummy_header5;
    dummy_header5.set_filetype("tif");
    dummy_header5.set_col(1);
    dummy_header5.set_row(1);
    dummy_header5.set_depth(1);

    rec.set_blob_id( idx.write_request(1024) );
    rec.set_blob_offset(blob.write(dummy_header1, m_test_data, 20));
    idx.write_complete(dummy_header1, rec);

    IndexRecord result = idx.read_request(0, 0, 1);
    TS_ASSERT_EQUALS(result.blob_id(), rec.blob_id());
    TS_ASSERT_EQUALS(result.blob_offset(), rec.blob_offset());

    // Now let's try some invalid reads/writes
    dummy_header5.set_col(10);
    TS_ASSERT_THROWS(idx.write_complete(dummy_header5, rec), IndexErr);
    dummy_header5.set_col(0);
    dummy_header5.set_row(2);
    dummy_header5.set_depth(1);
    TS_ASSERT_THROWS(idx.write_complete(dummy_header5, rec), IndexErr);

    TS_ASSERT_THROWS(idx.read_request(0, 0, 2), TileNotFoundErr);
  }

  void test_index_read_write() {
    std::string plate_filename = "/tmp/foo.plate";
    std::string index_filename = plate_filename + "/plate.index";
    std::string blob_filename = plate_filename + "/plate_0.blob";
    
    unlink(index_filename.c_str());
    unlink(blob_filename.c_str());

    TileHeader dummy_header0;
    dummy_header0.set_filetype("tif");
    dummy_header0.set_col(0);
    dummy_header0.set_row(0);
    dummy_header0.set_depth(0);

    TileHeader dummy_header1;
    dummy_header1.set_filetype("tif");
    dummy_header1.set_col(0);
    dummy_header1.set_row(0);
    dummy_header1.set_depth(1);

    TileHeader dummy_header2;
    dummy_header2.set_filetype("tif");
    dummy_header2.set_col(1);
    dummy_header2.set_row(0);
    dummy_header2.set_depth(1);

    TileHeader dummy_header3;
    dummy_header3.set_filetype("tif");
    dummy_header3.set_col(0);
    dummy_header3.set_row(1);
    dummy_header3.set_depth(1);

    TileHeader dummy_header4;
    dummy_header4.set_filetype("tif");
    dummy_header4.set_col(1);
    dummy_header4.set_row(1);
    dummy_header4.set_depth(1);

    // Write some data to the Index. 
    { 
      Index idx(plate_filename, 256, "tif", VW_PIXEL_RGB, VW_CHANNEL_UINT8);
      Blob blob(blob_filename);

      IndexRecord rec;
      rec.set_blob_id( idx.write_request(1024) );
      rec.set_blob_offset(blob.write(dummy_header0, m_test_data, 20));
      rec.set_status(INDEX_RECORD_VALID);
      idx.write_complete(dummy_header0, rec);

      rec.set_blob_id( idx.write_request(1024) );
      rec.set_blob_offset(blob.write(dummy_header1, m_test_data, 20));
      idx.write_complete(dummy_header1, rec);
      
      rec.set_blob_id( idx.write_request(1024) );
      rec.set_blob_offset(blob.write(dummy_header2, m_test_data, 20));
      idx.write_complete(dummy_header2, rec);
      
      rec.set_blob_id( idx.write_request(1024) );
      rec.set_blob_offset(blob.write(dummy_header3, m_test_data, 20));
      idx.write_complete(dummy_header3, rec);
      
      rec.set_blob_id( idx.write_request(1024) );
      rec.set_blob_offset(blob.write(dummy_header4, m_test_data, 20));
      idx.write_complete(dummy_header4, rec);
    }

    // Now, let's save the data to disk, and then read it back.
    {
      Index idx2(plate_filename);
      Blob blob(blob_filename);

      // Read the data back from the index
      IndexRecord result = idx2.read_request(0, 0, 0);
      TileHeader hdr = blob.read_header<TileHeader>(result.blob_offset());
      boost::shared_array<uint8> retrieved_data = blob.read_data(result.blob_offset());
      check_data(m_test_data, retrieved_data);
      TS_ASSERT_EQUALS(hdr.col(), dummy_header0.col());
      TS_ASSERT_EQUALS(hdr.row(), dummy_header0.row());
      TS_ASSERT_EQUALS(hdr.depth(), dummy_header0.depth());
      TS_ASSERT_EQUALS(hdr.filetype(), dummy_header0.filetype());

      result = idx2.read_request(0, 0, 1);
      hdr = blob.read_header<TileHeader>(result.blob_offset());
      retrieved_data = blob.read_data(result.blob_offset());
      check_data(m_test_data, retrieved_data);
      TS_ASSERT_EQUALS(hdr.col(), dummy_header1.col());
      TS_ASSERT_EQUALS(hdr.row(), dummy_header1.row());
      TS_ASSERT_EQUALS(hdr.depth(), dummy_header1.depth());
      TS_ASSERT_EQUALS(hdr.filetype(), dummy_header1.filetype());

      result = idx2.read_request(1, 0, 1);
      hdr = blob.read_header<TileHeader>(result.blob_offset());
      retrieved_data = blob.read_data(result.blob_offset());
      check_data(m_test_data, retrieved_data);
      TS_ASSERT_EQUALS(hdr.col(), dummy_header2.col());
      TS_ASSERT_EQUALS(hdr.row(), dummy_header2.row());
      TS_ASSERT_EQUALS(hdr.depth(), dummy_header2.depth());
      TS_ASSERT_EQUALS(hdr.filetype(), dummy_header2.filetype());

      result = idx2.read_request(0, 1, 1);
      hdr = blob.read_header<TileHeader>(result.blob_offset());
      retrieved_data = blob.read_data(result.blob_offset());
      check_data(m_test_data, retrieved_data);
      TS_ASSERT_EQUALS(hdr.col(), dummy_header3.col());
      TS_ASSERT_EQUALS(hdr.row(), dummy_header3.row());
      TS_ASSERT_EQUALS(hdr.depth(), dummy_header3.depth());
      TS_ASSERT_EQUALS(hdr.filetype(), dummy_header3.filetype());

      result = idx2.read_request(1, 1, 1);
      hdr = blob.read_header<TileHeader>(result.blob_offset());
      retrieved_data = blob.read_data(result.blob_offset());
      check_data(m_test_data, retrieved_data);
      TS_ASSERT_EQUALS(hdr.col(), dummy_header4.col());
      TS_ASSERT_EQUALS(hdr.row(), dummy_header4.row());
      TS_ASSERT_EQUALS(hdr.depth(), dummy_header4.depth());
      TS_ASSERT_EQUALS(hdr.filetype(), dummy_header4.filetype());
    }

    unlink(index_filename.c_str());
    unlink(blob_filename.c_str());
  }

}; // class TestIndex
