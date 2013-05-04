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


#include <gtest/gtest_VW.h>

using namespace std;

#if defined(VW_HAVE_APACHE) && (VW_HAVE_APACHE==1)
#  define HAS_APACHE(x) x
#else
#  define HAS_APACHE(x) DISABLED_ ## x
#endif

TEST(ModPlate, HAS_APACHE(Assumptions)) {
  // I'm storing an intptr_t in an apr_table, which expects a const char*.
  // This test makes sure that's okay...ish

  ASSERT_GE(sizeof(intptr_t), sizeof(const char*));

  volatile const char *buf1;
  volatile intptr_t buf2;

  intptr_t zero = 0,
           big   = numeric_limits<intptr_t>::max(),
           small = numeric_limits<intptr_t>::min();

  buf1 = (const char*)zero;
  buf2 = (intptr_t)buf1;
  EXPECT_EQ(zero, buf2);

  buf1 = (const char*)big;
  buf2 = (intptr_t)buf1;
  EXPECT_EQ(big, buf2);

  buf1 = (const char*)small;
  buf2 = (intptr_t)buf1;
  EXPECT_EQ(small, buf2);

  // There doesn't seem to be any reliable way to LINK against apr/apache. Sigh.
#if 0
  apr_status_t ret;
  apr_pool_t *pool;
  apr_table_t *table;

  ret = apr_pool_create(&pool, NULL);
  ASSERT_EQ(APR_SUCCESS, ret);

  apr_table_t *t = apr_table_make(pool, 4);
  ASSERT_TRUE(t);

  apr_table_setn(table, "zero",  (const char*)zero);
  apr_table_setn(table, "big",   (const char*)big);
  apr_table_setn(table, "small", (const char*)small);

  EXPECT_EQ(zero,  reinterpret_cast<intptr_t>(apr_table_get(table, "zero")));
  EXPECT_EQ(big,   reinterpret_cast<intptr_t>(apr_table_get(table, "big")));
  EXPECT_EQ(small, reinterpret_cast<intptr_t>(apr_table_get(table, "small")));

  apr_pool_destroy(pool);
#endif
}
