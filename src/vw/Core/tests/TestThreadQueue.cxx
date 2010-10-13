// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <gtest/gtest.h>
#include <vw/Core/ThreadQueue.h>

using namespace vw;

TEST(ThreadQueue, Basic) {
  ThreadQueue<uint32> q;

  ASSERT_TRUE(q.empty());
  for (uint32 i = 0; i < 50; ++i)
    q.push(i);
  ASSERT_FALSE(q.empty());

  uint32 pop;
  for (uint32 i = 0; i < 50; ++i) {
    ASSERT_FALSE(q.empty());
    q.wait_pop(pop);
    EXPECT_EQ(i, pop);
  }
}

class PushTask {
    ThreadQueue<uint32>& m_queue;
    unsigned m_count, m_value;
  public:
    PushTask(ThreadQueue<uint32>& q, uint32 count, uint32 value) : m_queue(q), m_count(count), m_value(value) {}
    void operator()() {
      for (uint32 i = 0; i < m_count; ++i) {
        m_queue.push(m_value);
      }
    }
};

TEST(ThreadQueue, Threaded) {
  typedef boost::shared_ptr<PushTask>    TheTask;
  typedef boost::shared_ptr<vw::Thread>  TheThread;

  ThreadQueue<uint32> q;
  ASSERT_TRUE(q.empty());

  std::vector<std::pair<TheTask, TheThread> > threads(20);

  for (size_t i = 0; i < threads.size(); ++i) {
    TheTask task(new PushTask(q, 10, uint32(i)));
    TheThread thread( new Thread(task) );
    threads[i] = std::make_pair(task, thread);
  }

  for (size_t i = 0; i < threads.size(); ++i) {
    threads[i].second->join();
  }

  ASSERT_FALSE(q.empty());

  std::vector<uint32> ret(threads.size());
  uint32 value = uint32(ret.size())+1;

  while (!q.empty()) {
    EXPECT_TRUE(q.timed_wait_pop(value, 0));
    EXPECT_LT(value, ret.size());
    ret[value]++;
  }
  EXPECT_FALSE(q.timed_wait_pop(value, 0));

  for (size_t i = 0; i < ret.size(); ++i) {
    EXPECT_EQ(10u, ret[i]);
  }
}
