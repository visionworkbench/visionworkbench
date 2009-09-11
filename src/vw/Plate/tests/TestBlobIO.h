// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <cxxtest/TestSuite.h>

#include <vw/Plate/BlobIO.h>

using namespace std;
using namespace vw;
using namespace vw::platefile;

class TestBlobIO : public CxxTest::TestSuite
{
  boost::shared_array<int> m_test_data;
  boost::shared_array<int> m_verify_data;
public:

  TestBlobIO() {
    m_test_data = boost::shared_array<int>(new int[20]);
    for (int i = 1; i < 21; ++i) {
      m_test_data[i-1] = i;
    }
  }

  void test_write_then_read()
  {
    unlink("/tmp/foo.blob");

    // First test, creates a new blob file.
    {
      Blob blob("/tmp/foo.blob");
      TS_ASSERT_EQUALS( blob.version(), 1 );
    
      // Write the data to the file.
      size_t offset = blob.write(m_test_data, 20 * sizeof(int));

      // Read it back in.
      boost::shared_array<int> m_verify_data = blob.read<int>(offset, 20 * sizeof(int));

      for (int i = 0; i < 20; ++i) 
        TS_ASSERT_EQUALS( m_test_data[i], m_verify_data[i] );
    }


    // Second tests, appends to a blob file.
    {
      Blob blob("/tmp/foo.blob");
      TS_ASSERT_EQUALS( blob.version(), 1 );
    
      // Write the data to the file.
      size_t offset = blob.write(m_test_data, 20 * sizeof(int));
      
      // Read it back in.
      boost::shared_array<int> m_verify_data = blob.read<int>(offset, 20 * sizeof(int));

      for (int i = 0; i < 20; ++i) 
        TS_ASSERT_EQUALS( m_test_data[i], m_verify_data[i] );

    }


  }

  void test_file_write_read() {
    const char* f1 = "/tmp/foo.blob";
    const char* f2 = "/tmp/foo3.blob";
    unlink(f2);

    Blob blob("/tmp/foo2.blob");
    
    // Do one loop through the blob file, placing f1 into the file,
    // and then reading it back out and saving it as f2.
    size_t offset, size;
    blob.write_from_file(f1, offset, size);
    blob.read_as_file(f2, offset, size);

    // ----

    // Now, check to make sure it worked!!
    std::ifstream istr1(f1, ios::binary);
    std::ifstream istr2(f2, ios::binary);

    // Check to see if files are the same size
    istr1.seekg (0, std::ios::end);
    istr2.seekg (0, std::ios::end);
    size_t size1 = istr1.tellg();
    size_t size2 = istr2.tellg();
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

    

}; // class TestBlobIO
