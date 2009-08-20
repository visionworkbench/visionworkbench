// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <memory>
#include <gtest/gtest.h>

#include <vw/Core/Cache.h>
#include <vw/Core/FundamentalTypes.h>

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

class CacheTest : public ::testing::Test {
  public:
    typedef Cache::Handle<BlockGenerator> HandleT;
    typedef std::vector<HandleT> ContainerT;

    vw::Cache cache;
    ContainerT cache_handles;

    const static int dimension = 1024;
    const static int num_actual_blocks = 10;
    const static int num_cache_blocks = 3;

    CacheTest() :
      cache( num_cache_blocks*dimension*dimension ), cache_handles( num_actual_blocks ) {
        for (int i = 0; i < num_actual_blocks; ++i) {
          cache_handles[i] = cache.insert( BlockGenerator( dimension, i ) );
        }
    }
};

TEST_F(CacheTest, CacheLineLRU) {
  // Test assumes cache size less than total data
  ASSERT_LT(+num_cache_blocks, +num_actual_blocks);

  for (int i = 0; i < num_actual_blocks; ++i) {

    Cache::Handle<BlockGenerator> &h = cache_handles[i];

    ASSERT_EQ(h.size(), dimension*dimension*sizeof(BlockGenerator::value_type));

    ASSERT_FALSE(h.valid());
    EXPECT_EQ(i, *h);
    EXPECT_TRUE(h.valid());

    boost::shared_ptr<BlockGenerator::value_type> ptr = h;
    EXPECT_EQ(i, *ptr);

    // LRU cache, so the last num_cache_blocks should still be valid
    for (int j = 0; i-j >= 0; ++j)
    {
      SCOPED_TRACE(::testing::Message() << "Cache block " << i);
      if (j < num_cache_blocks) {
        EXPECT_TRUE(cache_handles[i-j].valid());
      } else {
        EXPECT_FALSE(cache_handles[i-j].valid());
      }
    }
  }
}

TEST_F(CacheTest, Priority) {
  // Test assumes cache size less than total data
  ASSERT_LT(+num_cache_blocks, +num_actual_blocks);

  // prime the cache
  for (int i = 0; i < num_actual_blocks; ++i) {
    EXPECT_EQ(i, *cache_handles[i]);
  }

  //make sure last element is valid
  ASSERT_TRUE(cache_handles[num_actual_blocks-1].valid());

  // deprioritize the last element, and ensure it's still valid
  cache_handles[num_actual_blocks-1].deprioritize();
  ASSERT_TRUE(cache_handles[num_actual_blocks-1].valid());

  // bring the first element back into cache
  ASSERT_FALSE(cache_handles[0].valid());
  EXPECT_EQ( 0, *cache_handles[0] );
  EXPECT_TRUE(cache_handles[0].valid());

  // make sure that the deprioritized element dropped out of cache
  EXPECT_FALSE(cache_handles[num_actual_blocks-1].valid());
}
