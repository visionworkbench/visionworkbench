// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <cxxtest/TestSuite.h>

#include <vw/Plate/Blob.h>

using namespace std;
using namespace vw;
using namespace vw::platefile;

class TestBlobIterator : public CxxTest::TestSuite
{
  boost::shared_array<uint8> m_test_data;
  boost::shared_array<uint8> m_verify_data;
public:

  TestBlobIterator() {
    m_test_data = boost::shared_array<uint8>(new uint8[20]);
    for (int i = 0; i < 20; ++i) {
      m_test_data[i] = i;
    }
  }

  void test_iterator()
  {

    unlink("/tmp/foo.blob");
    boost::shared_ptr<Blob> blob(new Blob("/tmp/foo.blob"));

    // First, create some entries for the blob
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
    
    // Write some dummy data to the file with various index entries
    int64 offset = blob->write(dummy_record0, m_test_data, 20);
    offset = blob->write(dummy_record1, m_test_data, 20);
    offset = blob->write(dummy_record2, m_test_data, 20);
    
    // Create an iterator
    Blob::iterator iter = blob->begin();

    IndexRecord result = *iter;
    TS_ASSERT_EQUALS(result.blob_id(), dummy_record0.blob_id());
    TS_ASSERT_EQUALS(result.blob_offset(), dummy_record0.blob_offset());
    TS_ASSERT_EQUALS(result.block_size(), dummy_record0.block_size());
    TS_ASSERT_EQUALS(result.block_filetype(), dummy_record0.block_filetype());
    TS_ASSERT_DIFFERS(iter, blob->end());

    ++iter;
    result = *iter;
    TS_ASSERT_EQUALS(result.blob_id(), dummy_record1.blob_id());
    TS_ASSERT_EQUALS(result.blob_offset(), dummy_record1.blob_offset());
    TS_ASSERT_EQUALS(result.block_size(), dummy_record1.block_size());
    TS_ASSERT_EQUALS(result.block_filetype(), dummy_record1.block_filetype());
    TS_ASSERT_DIFFERS(iter, blob->end());

    ++iter;
    result = *iter;
    TS_ASSERT_EQUALS(result.blob_id(), dummy_record2.blob_id());
    TS_ASSERT_EQUALS(result.blob_offset(), dummy_record2.blob_offset());
    TS_ASSERT_EQUALS(result.block_size(), dummy_record2.block_size());
    TS_ASSERT_EQUALS(result.block_filetype(), dummy_record2.block_filetype());
    TS_ASSERT_DIFFERS(iter, blob->end());

    ++iter;
    TS_ASSERT_EQUALS(iter, blob->end());
  }

    
}; // class TestBlobIterator
