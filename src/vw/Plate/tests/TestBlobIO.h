// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <cxxtest/TestSuite.h>

#include <vw/Plate/Blob.h>
#include <vw/Plate/BlobManager.h>
#include <vw/Plate/IndexRecord.pb.h>

using namespace std;
using namespace vw;
using namespace vw::platefile;

class TestBlobIO : public CxxTest::TestSuite
{
  boost::shared_array<uint8> m_test_data;
  boost::shared_array<uint8> m_verify_data;
public:

  TestBlobIO() {
    m_test_data = boost::shared_array<uint8>(new uint8[20]);
    for (int i = 0; i < 20; ++i) {
      m_test_data[i] = i;
    }
  }

  void test_write_then_read()
  {
    IndexRecord rec; 
    rec.set_block_size(22);
    unlink("/tmp/foo.blob");

    // First test, creates a new blob file.
    {
      Blob blob("/tmp/foo.blob");
    
      // Write the data to the file.
      int64 offset = blob.write(rec, m_test_data, 20);

      // Read it back in.
      boost::shared_array<uint8> m_verify_data = blob.read_data(offset);

      for (int i = 0; i < 20; ++i) 
        TS_ASSERT_EQUALS( m_test_data[i], m_verify_data[i] );
    }


    // Second test, appends to a blob file.
    {
      Blob blob("/tmp/foo.blob");
    
      // Write the data to the file.
      
      int64 offset = blob.write(rec, m_test_data, 20);
      
      // Read it back in.
      boost::shared_array<uint8> m_verify_data = blob.read_data(offset);

      for (int i = 0; i < 20; ++i) 
        TS_ASSERT_EQUALS( m_test_data[i], m_verify_data[i] );
    }

  }

  void test_file_write_read() {
    IndexRecord rec; 
    rec.set_block_size(22);

    const char* f1 = "/tmp/foo.blob";
    const char* f2 = "/tmp/foo3.blob";
    unlink(f2);

    Blob blob("/tmp/foo2.blob");
    
    // Do one loop through the blob file, placing f1 into the file,
    // and then reading it back out and saving it as f2.
    int64 offset;
    int32 size;
    blob.write_from_file(f1, rec, offset, size);
    blob.read_to_file(f2, offset, size);

    // ----

    // Now, check to make sure it worked!!
    std::ifstream istr1(f1, ios::binary);
    std::ifstream istr2(f2, ios::binary);

    // Check to see if files are the same size
    istr1.seekg (0, std::ios::end);
    istr2.seekg (0, std::ios::end);
    int64 size1 = istr1.tellg();
    int64 size2 = istr2.tellg();
    TS_ASSERT_EQUALS(size1, size2);

    if (size1 == size2) {
      istr1.seekg (0, std::ios::beg);
      istr2.seekg (0, std::ios::beg);
      
      boost::shared_array<char> data1(new char[size1]);
      boost::shared_array<char> data2(new char[size2]);
      
      istr1.read(data1.get(), size1);
      istr2.read(data2.get(), size2);

      for (int i=0; i < size1; ++i) {
        if (data1[i] != data2[i]) {
          TS_FAIL("Write/Read test of blob file failed. Data does not agree between input and output files.");
        }
      }
    }
    
    istr1.close();
    istr2.close();

    unlink(f1);
    unlink(f2);
    unlink("/tmp/foo2.blob");
  }

  void test_blob_manager() {

    BlobManager bm(1024, 4);
    TS_ASSERT_EQUALS(bm.num_blobs(), 4);

    std::cout << "\nTesting blob manager.  If this hangs, then there is a problem!!\n";
    int blob_id = bm.request_lock(1024);

    blob_id = bm.request_lock(1024);
    blob_id = bm.request_lock(1024);
    blob_id = bm.request_lock(1024);
    bm.release_lock(0);
    blob_id = bm.request_lock(1024);
  }
    
}; // class TestBlobIO
