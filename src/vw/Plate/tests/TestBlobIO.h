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
      std::cout << "Offset : " << offset << "\n";

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
      std::cout << "Offset : " << offset << "\n";
      
      // Read it back in.
      boost::shared_array<int> m_verify_data = blob.read<int>(offset, 20 * sizeof(int));

      for (int i = 0; i < 20; ++i) 
        TS_ASSERT_EQUALS( m_test_data[i], m_verify_data[i] );

    }


  }

}; // class TestBlobIO
