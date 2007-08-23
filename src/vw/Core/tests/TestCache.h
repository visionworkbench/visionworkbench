// __BEGIN_LICENSE__
// 
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
// 
// Copyright 2006 Carnegie Mellon University. All rights reserved.
// 
// This software is distributed under the NASA Open Source Agreement
// (NOSA), version 1.3.  The NOSA has been approved by the Open Source
// Initiative.  See the file COPYING at the top of the distribution
// directory tree for the complete NOSA document.
// 
// THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
// KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
// LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
// SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
// A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
// THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
// DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
// 
// __END_LICENSE__

// TestFunctors.h
#include <cxxtest/TestSuite.h>

#include <vw/Core/Cache.h>
#include <vw/Core/FundamentalTypes.h>

using namespace std;
using namespace vw;



class TestCache : public CxxTest::TestSuite
{
  // A dummy BlockGenerator that generates blocks of 1-byte data.
  class BlockGenerator {
    int m_dimension;
    vw::uint8 m_fill_value;

  public:
    BlockGenerator(int dimension, vw::uint8 fill_value = 0) :
      m_dimension(dimension), m_fill_value(fill_value) {}
  
    typedef vw::uint8 value_type;
  
  
    size_t size() const {
      return m_dimension * m_dimension * sizeof(vw::uint8);
    }
    
    boost::shared_ptr< value_type > generate() const {
      
      boost::shared_ptr< value_type > ptr(new vw::uint8[m_dimension*m_dimension]);
      (&(*ptr))[0] = m_fill_value;
      return ptr;
    }
  };

public:
  void test_cache_line()
  {

    static int dimension = 1024;
    static int num_actual_blocks = 3;
    static int num_cache_blocks = 2;

    //std::cout << "\n\n";
    //set_debug_level(VerboseDebugMessage+1);

    vw::Cache cache( num_cache_blocks*dimension*dimension );  

    std::vector<Cache::Handle<BlockGenerator> > cache_handles(num_actual_blocks);
    for (int i = 0; i < num_actual_blocks; ++i) {
      cache_handles[i] = cache.insert( BlockGenerator( dimension, i ) );
    }

    for (int i = 0; i < num_actual_blocks; ++i) {
      TS_ASSERT_EQUALS(i,(&(*cache_handles[i]))[0]);
    }
  }
};
