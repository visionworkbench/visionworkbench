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

public:

  void test_index_record() {

    // Test serialization/deserialization
    IndexRecord write_record;
    write_record.set_blob_id(1);
    write_record.set_blob_offset(2);
    write_record.set_block_size(1024);
    write_record.set_block_filetype("tiff");
    std::ofstream ostr("/tmp/foo.bar", std::ios::binary);
    write_record.SerializeToOstream(&ostr);
    ostr.close();

    std::ifstream istr("/tmp/foo.bar", std::ios::binary);
    IndexRecord read_record;
    read_record.ParseFromIstream(&istr);
    istr.close();

    TS_ASSERT_EQUALS(write_record.blob_id(), read_record.blob_id());
    TS_ASSERT_EQUALS(write_record.blob_offset(), read_record.blob_offset());
    TS_ASSERT_EQUALS(write_record.block_size(), read_record.block_size());
    TS_ASSERT_EQUALS(write_record.block_filetype(), read_record.block_filetype());
    TS_ASSERT_EQUALS(write_record.valid(), read_record.valid());
    
    // Clean up
    unlink("/tmp/foo.bar");
  }

  void test_simple_index() {

    IndexRecord dummy_record0;
    dummy_record0.set_blob_id(0);
    dummy_record0.set_blob_offset(9);
    dummy_record0.set_block_size(1024);
    dummy_record0.set_block_filetype("foo");

    IndexRecord dummy_record1;
    dummy_record1.set_blob_id(0);
    dummy_record1.set_blob_offset(10);
    dummy_record1.set_block_size(1024);
    dummy_record1.set_block_filetype("tiff");

    IndexRecord dummy_record2;
    dummy_record2.set_blob_id(0);
    dummy_record2.set_blob_offset(11);
    dummy_record2.set_block_size(1024);
    dummy_record2.set_block_filetype("png");

    IndexRecord dummy_record3;
    dummy_record3.set_blob_id(0);
    dummy_record3.set_blob_offset(12);
    dummy_record3.set_block_size(1024);
    dummy_record3.set_block_filetype("jpg");

    IndexRecord dummy_record4;
    dummy_record4.set_blob_id(0);
    dummy_record4.set_blob_offset(13);
    dummy_record4.set_block_size(1024);
    dummy_record4.set_block_filetype("tif");


    // Write some data to the Index.
    boost::shared_ptr<BlobManager> mgr( new BlobManager(2048, 3) );
    Index idx(mgr, 256, "tif");

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
    TS_ASSERT_EQUALS(result.block_size(), dummy_record0.block_size());
    TS_ASSERT_EQUALS(result.block_filetype(), dummy_record0.block_filetype());

    result = idx.read_request(0, 0, 1);
    TS_ASSERT_EQUALS(result.blob_id(), dummy_record1.blob_id());
    TS_ASSERT_EQUALS(result.blob_offset(), dummy_record1.blob_offset());
    TS_ASSERT_EQUALS(result.block_size(), dummy_record1.block_size());
    TS_ASSERT_EQUALS(result.block_filetype(), dummy_record1.block_filetype());

    result = idx.read_request(1, 0, 1);
    TS_ASSERT_EQUALS(result.blob_id(), dummy_record2.blob_id());
    TS_ASSERT_EQUALS(result.blob_offset(), dummy_record2.blob_offset());
    TS_ASSERT_EQUALS(result.block_size(), dummy_record2.block_size());
    TS_ASSERT_EQUALS(result.block_filetype(), dummy_record2.block_filetype());

    result = idx.read_request(0, 1, 1);
    TS_ASSERT_EQUALS(result.blob_id(), dummy_record3.blob_id());
    TS_ASSERT_EQUALS(result.blob_offset(), dummy_record3.blob_offset());
    TS_ASSERT_EQUALS(result.block_size(), dummy_record3.block_size());
    TS_ASSERT_EQUALS(result.block_filetype(), dummy_record3.block_filetype());

    result = idx.read_request(1, 1, 1);
    TS_ASSERT_EQUALS(result.blob_id(), dummy_record4.blob_id());
    TS_ASSERT_EQUALS(result.blob_offset(), dummy_record4.blob_offset());
    TS_ASSERT_EQUALS(result.block_size(), dummy_record4.block_size());
    TS_ASSERT_EQUALS(result.block_filetype(), dummy_record4.block_filetype());

    // Test re-writing (i.e. changing) an entry.
    IndexRecord dummy_record5;
    dummy_record5.set_blob_id(0);
    dummy_record5.set_blob_offset(15);
    dummy_record5.set_block_size(1024);
    dummy_record5.set_block_filetype("jp2k");

    dummy_record5.set_blob_id( idx.write_request(1024) );
    idx.write_complete(0, 0, 1, dummy_record5);

    result = idx.read_request(0, 0, 1);
    TS_ASSERT_EQUALS(result.blob_id(), dummy_record5.blob_id());
    TS_ASSERT_EQUALS(result.blob_offset(), dummy_record5.blob_offset());
    TS_ASSERT_EQUALS(result.block_size(), dummy_record5.block_size());
    TS_ASSERT_EQUALS(result.block_filetype(), dummy_record5.block_filetype());

    // Now let's try some invalid reads/writes
    dummy_record5.set_blob_id( idx.write_request(1024) );
    TS_ASSERT_THROWS(idx.write_complete(10, 0, 0, dummy_record5), IndexErr);
    dummy_record5.set_blob_id( idx.write_request(1024) );
    TS_ASSERT_THROWS(idx.write_complete(0, 2, 1, dummy_record5), IndexErr);

    TS_ASSERT_THROWS(idx.read_request(0, 0, 2), TileNotFoundErr);
  }

  void test_index_read_write() {

    IndexRecord dummy_record0;
    dummy_record0.set_blob_id(0);
    dummy_record0.set_blob_offset(9);
    dummy_record0.set_block_size(1024);
    dummy_record0.set_block_filetype("foo");

    IndexRecord dummy_record1;
    dummy_record1.set_blob_id(0);
    dummy_record1.set_blob_offset(10);
    dummy_record1.set_block_size(1024);
    dummy_record1.set_block_filetype("tiff");

    IndexRecord dummy_record2;
    dummy_record2.set_blob_id(0);
    dummy_record2.set_blob_offset(11);
    dummy_record2.set_block_size(1024);
    dummy_record2.set_block_filetype("png");

    IndexRecord dummy_record3;
    dummy_record3.set_blob_id(0);
    dummy_record3.set_blob_offset(12);
    dummy_record3.set_block_size(1024);
    dummy_record3.set_block_filetype("jpg");

    IndexRecord dummy_record4;
    dummy_record4.set_blob_id(0);
    dummy_record4.set_blob_offset(13);
    dummy_record4.set_block_size(1024);
    dummy_record4.set_block_filetype("tif");

    // Write some data to the Index.
    boost::shared_ptr<BlobManager> mgr( new BlobManager(2048, 3) );
    Index idx(mgr, 256, "tif");

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

    
    // Now, let's save the data to disk, and then read it back.
    idx.save("/tmp/foo.index");
    Index idx2("/tmp/foo.index");

    // Read the data back from the index
    IndexRecord result = idx2.read_request(0, 0, 0);
    TS_ASSERT_EQUALS(result.blob_id(), dummy_record0.blob_id());
    TS_ASSERT_EQUALS(result.blob_offset(), dummy_record0.blob_offset());
    TS_ASSERT_EQUALS(result.block_size(), dummy_record0.block_size());
    TS_ASSERT_EQUALS(result.block_filetype(), dummy_record0.block_filetype());
    TS_ASSERT_EQUALS(result.valid(), dummy_record0.valid());

    result = idx2.read_request(0, 0, 1);
    TS_ASSERT_EQUALS(result.blob_id(), dummy_record1.blob_id());
    TS_ASSERT_EQUALS(result.blob_offset(), dummy_record1.blob_offset());
    TS_ASSERT_EQUALS(result.block_size(), dummy_record1.block_size());
    TS_ASSERT_EQUALS(result.block_filetype(), dummy_record1.block_filetype());

    result = idx2.read_request(1, 0, 1);
    TS_ASSERT_EQUALS(result.blob_id(), dummy_record2.blob_id());
    TS_ASSERT_EQUALS(result.blob_offset(), dummy_record2.blob_offset());
    TS_ASSERT_EQUALS(result.block_size(), dummy_record2.block_size());
    TS_ASSERT_EQUALS(result.block_filetype(), dummy_record2.block_filetype());

    result = idx2.read_request(0, 1, 1);
    TS_ASSERT_EQUALS(result.blob_id(), dummy_record3.blob_id());
    TS_ASSERT_EQUALS(result.blob_offset(), dummy_record3.blob_offset());
    TS_ASSERT_EQUALS(result.block_size(), dummy_record3.block_size());
    TS_ASSERT_EQUALS(result.block_filetype(), dummy_record3.block_filetype());

    result = idx2.read_request(1, 1, 1);
    TS_ASSERT_EQUALS(result.blob_id(), dummy_record4.blob_id());
    TS_ASSERT_EQUALS(result.blob_offset(), dummy_record4.blob_offset());
    TS_ASSERT_EQUALS(result.block_size(), dummy_record4.block_size());
    TS_ASSERT_EQUALS(result.block_filetype(), dummy_record4.block_filetype());

    //unlink("/tmp/foo.index");
  }

}; // class TestIndex
