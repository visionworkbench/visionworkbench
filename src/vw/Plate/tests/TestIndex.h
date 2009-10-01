// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
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

  void test_index_record() {

    // Test serialization/deserialization
    TileHeader write_hdr;
    write_hdr.mutable_index_record()->set_blob_id(1);
    write_hdr.mutable_index_record()->set_blob_offset(2);
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
    TS_ASSERT_EQUALS(write_hdr.index_record().blob_id(), 
                     read_hdr.index_record().blob_id());
    TS_ASSERT_EQUALS(write_hdr.index_record().blob_offset(), 
                     read_hdr.index_record().blob_offset());

    // Clean up
    unlink("/tmp/foo.bar");
  }

  void test_index_write_read() {
    unlink("/tmp/foo.plate/plate.index");
    std::string plate_filename = "/tmp/foo.plate";
    boost::shared_ptr<BlobManager> mgr( new BlobManager(2048, 3) );

    // Write the basic data...
    Index idx(plate_filename, 256, "tif", mgr);

    // And read it back in...
    Index idx2(plate_filename);
    
    TS_ASSERT_EQUALS(idx.version(), idx2.version());
    TS_ASSERT_EQUALS(idx.default_tile_size(), idx2.default_tile_size());
    TS_ASSERT_EQUALS(idx.default_tile_filetype(), idx2.default_tile_filetype());
  }


  void test_simple_index() {

    IndexRecord dummy_record0;
    dummy_record0.set_blob_id(0);
    dummy_record0.set_blob_offset(9);

    IndexRecord dummy_record1;
    dummy_record1.set_blob_id(0);
    dummy_record1.set_blob_offset(10);

    IndexRecord dummy_record2;
    dummy_record2.set_blob_id(0);
    dummy_record2.set_blob_offset(11);

    IndexRecord dummy_record3;
    dummy_record3.set_blob_id(0);
    dummy_record3.set_blob_offset(12);

    IndexRecord dummy_record4;
    dummy_record4.set_blob_id(0);
    dummy_record4.set_blob_offset(13);


    // Write some data to the Index.
    std::string plate_filename = "/tmp/foo.plate";
    boost::shared_ptr<BlobManager> mgr( new BlobManager(2048, 3) );
    Index idx(plate_filename, 256, "tif", mgr);

    dummy_record0.set_blob_id( idx.write_request(1024) );
    idx.write_complete(0, 0, 0, dummy_record0);

    dummy_record1.set_blob_id( idx.write_request(1024) );
    idx.write_complete(0, 0, 1, dummy_record1);
 
    dummy_record2.set_blob_id( idx.write_request(1024) );
    idx.write_complete(1, 0, 1, dummy_record2);

    dummy_record3.set_blob_id( idx.write_request(1024) );
    idx.write_complete(0, 1, 1, dummy_record3);

    dummy_record4.set_blob_id( idx.write_request(1024) );
    idx.write_complete(1, 1, 1, dummy_record4);

    // Read the data back from the index
    IndexRecord result = idx.read_request(0, 0, 0);
    TS_ASSERT_EQUALS(result.blob_id(), dummy_record0.blob_id());
    TS_ASSERT_EQUALS(result.blob_offset(), dummy_record0.blob_offset());

    result = idx.read_request(0, 0, 1);
    TS_ASSERT_EQUALS(result.blob_id(), dummy_record1.blob_id());
    TS_ASSERT_EQUALS(result.blob_offset(), dummy_record1.blob_offset());

    result = idx.read_request(1, 0, 1);
    TS_ASSERT_EQUALS(result.blob_id(), dummy_record2.blob_id());
    TS_ASSERT_EQUALS(result.blob_offset(), dummy_record2.blob_offset());

    result = idx.read_request(0, 1, 1);
    TS_ASSERT_EQUALS(result.blob_id(), dummy_record3.blob_id());
    TS_ASSERT_EQUALS(result.blob_offset(), dummy_record3.blob_offset());

    result = idx.read_request(1, 1, 1);
    TS_ASSERT_EQUALS(result.blob_id(), dummy_record4.blob_id());
    TS_ASSERT_EQUALS(result.blob_offset(), dummy_record4.blob_offset());

    // Test re-writing (i.e. changing) an entry.
    IndexRecord dummy_record5;
    dummy_record5.set_blob_id(0);
    dummy_record5.set_blob_offset(15);

    dummy_record5.set_blob_id( idx.write_request(1024) );
    idx.write_complete(0, 0, 1, dummy_record5);

    result = idx.read_request(0, 0, 1);
    TS_ASSERT_EQUALS(result.blob_id(), dummy_record5.blob_id());
    TS_ASSERT_EQUALS(result.blob_offset(), dummy_record5.blob_offset());

    // Now let's try some invalid reads/writes
    dummy_record5.set_blob_id( idx.write_request(1024) );
    TS_ASSERT_THROWS(idx.write_complete(10, 0, 0, dummy_record5), IndexErr);
    dummy_record5.set_blob_id( idx.write_request(1024) );
    TS_ASSERT_THROWS(idx.write_complete(0, 2, 1, dummy_record5), IndexErr);

    TS_ASSERT_THROWS(idx.read_request(0, 0, 2), TileNotFoundErr);
  }

  void test_index_read_write() {
    std::string plate_filename = "/tmp/foo.plate";
    std::string index_filename = plate_filename + "/plate.index";
    std::string blob_filename = plate_filename + "/plate_0.blob";
    
    unlink(index_filename.c_str());
    unlink(blob_filename.c_str());

    TileHeader dummy_record0;
    dummy_record0.set_filetype("tif");
    dummy_record0.set_col(0);
    dummy_record0.set_row(0);
    dummy_record0.set_depth(0);
    dummy_record0.mutable_index_record()->set_blob_id(0);
    dummy_record0.mutable_index_record()->set_blob_offset(9);

    TileHeader dummy_record1;
    dummy_record1.set_filetype("tif");
    dummy_record1.set_col(0);
    dummy_record1.set_row(0);
    dummy_record1.set_depth(1);
    dummy_record1.mutable_index_record()->set_blob_id(0);
    dummy_record1.mutable_index_record()->set_blob_offset(10);

    TileHeader dummy_record2;
    dummy_record2.set_filetype("tif");
    dummy_record2.set_col(1);
    dummy_record2.set_row(0);
    dummy_record2.set_depth(1);
    dummy_record2.mutable_index_record()->set_blob_id(0);
    dummy_record2.mutable_index_record()->set_blob_offset(11);

    TileHeader dummy_record3;
    dummy_record3.set_filetype("tif");
    dummy_record3.set_col(0);
    dummy_record3.set_row(1);
    dummy_record3.set_depth(1);
    dummy_record3.mutable_index_record()->set_blob_id(0);
    dummy_record3.mutable_index_record()->set_blob_offset(12);

    TileHeader dummy_record4;
    dummy_record4.set_filetype("tif");
    dummy_record4.set_col(1);
    dummy_record4.set_row(1);
    dummy_record4.set_depth(1);
    dummy_record4.mutable_index_record()->set_blob_id(0);
    dummy_record4.mutable_index_record()->set_blob_offset(13);

    // Write some data to the Index.
    boost::shared_ptr<BlobManager> mgr( new BlobManager(2048, 3) );
    Index idx(plate_filename, 256, "tif", mgr);
    Blob blob(blob_filename);

    dummy_record0.mutable_index_record()->set_blob_id( idx.write_request(1024) );
    blob.write(dummy_record0, m_test_data, 20);
    idx.write_complete(0, 0, 0, dummy_record0.index_record());

    dummy_record1.mutable_index_record()->set_blob_id( idx.write_request(1024) );
    blob.write(dummy_record1, m_test_data, 20);
    idx.write_complete(0, 0, 1, dummy_record1.index_record());
 
    dummy_record2.mutable_index_record()->set_blob_id( idx.write_request(1024) );
    blob.write(dummy_record2, m_test_data, 20);
    idx.write_complete(1, 0, 1, dummy_record2.index_record());

    dummy_record3.mutable_index_record()->set_blob_id( idx.write_request(1024) );
    blob.write(dummy_record3, m_test_data, 20);
    idx.write_complete(0, 1, 1, dummy_record3.index_record());

    dummy_record4.mutable_index_record()->set_blob_id( idx.write_request(1024) );
    blob.write(dummy_record4, m_test_data, 20);
    idx.write_complete(1, 1, 1, dummy_record4.index_record());
    
    // Now, let's save the data to disk, and then read it back.
    Index idx2(plate_filename);

    // Read the data back from the index
    IndexRecord result = idx2.read_request(0, 0, 0);
    TS_ASSERT_EQUALS(result.blob_id(), dummy_record0.index_record().blob_id());
    TS_ASSERT_EQUALS(result.blob_offset(), dummy_record0.index_record().blob_offset());
    TS_ASSERT_EQUALS(result.valid(), dummy_record0.index_record().valid());

    result = idx2.read_request(0, 0, 1);
    TS_ASSERT_EQUALS(result.blob_id(), dummy_record1.index_record().blob_id());
    TS_ASSERT_EQUALS(result.blob_offset(), dummy_record1.index_record().blob_offset());

    result = idx2.read_request(1, 0, 1);
    TS_ASSERT_EQUALS(result.blob_id(), dummy_record2.index_record().blob_id());
    TS_ASSERT_EQUALS(result.blob_offset(), dummy_record2.index_record().blob_offset());

    result = idx2.read_request(0, 1, 1);
    TS_ASSERT_EQUALS(result.blob_id(), dummy_record3.index_record().blob_id());
    TS_ASSERT_EQUALS(result.blob_offset(), dummy_record3.index_record().blob_offset());

    result = idx2.read_request(1, 1, 1);
    TS_ASSERT_EQUALS(result.blob_id(), dummy_record4.index_record().blob_id());
    TS_ASSERT_EQUALS(result.blob_offset(), dummy_record4.index_record().blob_offset());

    unlink(index_filename.c_str());
  }

}; // class TestIndex
