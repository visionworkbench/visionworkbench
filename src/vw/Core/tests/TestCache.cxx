// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NASA Vision Workbench is licensed under the Apache License,
//  Version 2.0 (the "License"); you may not use this file except in
//  compliance with the License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__

#include <numeric>
#include <gtest/gtest_VW.h>

#include <vw/Core/Cache.h>
#include <vw/Core/FundamentalTypes.h>
#include <vw/Core/ThreadPool.h>

using namespace vw;

// A dummy BlockGenerator that generates blocks of 1-byte data.
class BlockGenerator {

protected:
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

    boost::shared_ptr< value_type > ptr(new vw::uint8[m_dimension*m_dimension], boost::checked_array_deleter<value_type>());
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
        for (uint8 i = 0; i < num_actual_blocks; ++i) {
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
    EXPECT_NO_THROW( h.release() );
    EXPECT_TRUE(h.valid());

    boost::shared_ptr<BlockGenerator::value_type> ptr = h;
    EXPECT_EQ(i, *ptr);
    h.release();

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
    EXPECT_NO_THROW( cache_handles[i].release() );
  }

  //make sure last element is valid
  ASSERT_TRUE(cache_handles[num_actual_blocks-1].valid());

  // deprioritize the last element, and ensure it's still valid
  cache_handles[num_actual_blocks-1].deprioritize();
  ASSERT_TRUE(cache_handles[num_actual_blocks-1].valid());

  // bring the first element back into cache
  ASSERT_FALSE(cache_handles[0].valid());
  EXPECT_EQ( 0, *cache_handles[0] );
  EXPECT_NO_THROW( cache_handles[0].release() );
  EXPECT_TRUE(cache_handles[0].valid());

  // make sure that the deprioritized element dropped out of cache
  EXPECT_FALSE(cache_handles[num_actual_blocks-1].valid());
}

// Every copy increases the fill_value by one
class GenGen : public BlockGenerator {
  public:
    GenGen() : BlockGenerator(1, 0) {}

    GenGen(const GenGen& obj) :
      BlockGenerator(obj.m_dimension, boost::numeric_cast<uint8>(obj.m_fill_value+1)) {}
};

TEST(Cache, Types) {
  // This test makes sure cache handles can hold polymorphic types
  // and gives the user control over copying vs. not.
  typedef Cache::Handle<GenGen> obj_t;
  typedef Cache::Handle<boost::shared_ptr<GenGen> > sptr_t;

  vw::Cache cache(sizeof(vw::uint8));

  GenGen value;
  boost::shared_ptr<GenGen> shared(new GenGen());

  obj_t  h1 = cache.insert(value);
  sptr_t h2 = cache.insert(shared);

  // Make sure the value hasn't generated yet
  ASSERT_FALSE(h1.valid());
  ASSERT_FALSE(h2.valid());

  // value should have copied once
  EXPECT_EQ(1, *h1);
  EXPECT_NO_THROW( h1.release() );
  // shared_ptr should not have copied
  EXPECT_EQ(0, *h2);
  EXPECT_NO_THROW( h2.release() );
}

// Here's a more aggressive test that uses many threads plus a good
// chunk of memory (24k).
class ArrayDataGenerator {
public:
  typedef std::vector<char> value_type;
  ArrayDataGenerator() {}
  size_t size() const { return 1024; }
  boost::shared_ptr<value_type> generate() const {
    char number = rand() % 255;
    boost::shared_ptr<value_type> ptr( new value_type(1024,number) );
    return ptr;
  }
};

class TestTask : public vw::Task {
  std::vector<Cache::Handle<ArrayDataGenerator> > m_handles;

public:
  TestTask( std::vector<Cache::Handle<ArrayDataGenerator> > const& handles ) : m_handles( handles ) {}

  virtual ~TestTask() {}
  virtual void operator()() {
    size_t index_a = rand() % m_handles.size(),
      index_b = rand() % m_handles.size();
    std::vector<uint8> a( (const std::vector<uint8>&)(*m_handles[index_a]) );
    m_handles[index_a].release();
    std::vector<uint8> b( (const std::vector<uint8>&)(*m_handles[index_b]) );
    m_handles[index_b].release();
    VW_ASSERT( a.size() != 0, IOErr() << "Size is wrong!" );
    VW_ASSERT( b.size() != 0, IOErr() << "Size is wrong!" );
    std::vector<uint8> c( a.size() );
    for ( size_t i = 0; i < a.size(); i++ ) {
      c[i] = a[i] + b[i];
    }
    // Don't let the compiler optimize this away.
    volatile int result = std::accumulate(c.begin(), c.end(), 0 );
    EXPECT_EQ(result, result); // Rm compiler warning due to unused variable
  }
};

TEST(Cache, StressTest) {
  typedef Cache::Handle<ArrayDataGenerator> handle_t;
  vw::Cache cache( 6*1024 );

  // Creating handles for our data
  std::vector<handle_t> handles;
  for ( size_t i = 0; i < 24; i++ ) {
    handles.push_back( cache.insert( ArrayDataGenerator() ) );
  }

  // Create a FIFO buffer with a 1000 threads ... 12 running at
  // anytime that try to access the cache randomly.
  FifoWorkQueue queue(12);
  for ( size_t i = 0; i < 1000; i++ ) {
    boost::shared_ptr<Task> task( new TestTask(handles) );
    queue.add_task( task );
  }

  // If this test fails ... it will deadlock. Maybe we should have
  // this entire test case in a different thread where I can measure
  // its time?
  EXPECT_NO_THROW( queue.join_all(); );
}
