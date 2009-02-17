// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// TestFunctors.h
#include <memory>
#include <cxxtest/TestSuite.h>

#include <vw/Core/Cache.h>
#include <vw/Core/FundamentalTypes.h>

using namespace std;
using namespace vw;

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

class TestCache : public CxxTest::TestSuite
{
  typedef Cache::Handle<BlockGenerator> HandleT;
  typedef std::vector<HandleT> ContainerT;

  auto_ptr<vw::Cache> cache;
  auto_ptr<ContainerT> cache_handles;

  const static int dimension = 1024;
  const static int num_actual_blocks = 10;
  const static int num_cache_blocks = 3;

public:
  void setUp() {
    cache = auto_ptr<vw::Cache>( new vw::Cache(num_cache_blocks*dimension*dimension ) );
    cache_handles = auto_ptr<ContainerT>(new ContainerT(num_actual_blocks));

    for (int i = 0; i < num_actual_blocks; ++i) {
      cache_handles->at(i) = cache->insert( BlockGenerator( dimension, i ) );
    }
  }

  void tearDown() {
    cache_handles.reset();
    cache.reset();
  }

  void test_cache_line_lru()
  {
    // Test assumes cache size less than total data
    TS_ASSERT_LESS_THAN(num_cache_blocks, num_actual_blocks);

    for (int i = 0; i < num_actual_blocks; ++i) {

      Cache::Handle<BlockGenerator> &h = cache_handles->at(i);

      TS_ASSERT_EQUALS(h.size(), dimension*dimension*sizeof(BlockGenerator::value_type));

      TS_ASSERT(!h.valid());
      TS_ASSERT_EQUALS(i, *h);
      TS_ASSERT(h.valid());

      boost::shared_ptr<BlockGenerator::value_type> ptr = h;
      TS_ASSERT_EQUALS(i, *ptr);

      // LRU cache, so the last num_cache_blocks should still be valid
      for (int j = 0; i-j >= 0; ++j)
      {
        if (j < num_cache_blocks) {
          TS_ASSERT(cache_handles->at(i-j).valid());
        } else {
          TS_ASSERT(!cache_handles->at(i-j).valid());
        }
      }
    }
  }

  void test_cache_priority() {
    // Test assumes cache size less than total data
    TS_ASSERT_LESS_THAN(num_cache_blocks, num_actual_blocks);

    // prime the cache
    for (int i = 0; i < num_actual_blocks; ++i) {
      TS_ASSERT_EQUALS(*cache_handles->at(i), i);
    }

    //make sure last element is valid
    TS_ASSERT(cache_handles->at(num_actual_blocks-1).valid());

    // deprioritize the last element, and ensure it's still valid
    cache_handles->at(num_actual_blocks-1).deprioritize();
    TS_ASSERT(cache_handles->at(num_actual_blocks-1).valid());

    // bring the first element back into cache
    TS_ASSERT(!cache_handles->at(0).valid());
    TS_ASSERT_EQUALS(*cache_handles->at(0), 0);
    TS_ASSERT(cache_handles->at(0).valid());

    // make sure that the deprioritized element dropped out of cache
    TS_ASSERT(!cache_handles->at(num_actual_blocks-1).valid())
  }

};
