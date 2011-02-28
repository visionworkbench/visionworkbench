// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <gtest/gtest.h>
#include <test/Helpers.h>
#include <vw/Core/TypeDeduction.h>

using namespace vw;
using namespace std;

#if 0
These tests generated with this python script (requires >= python 2.7):
  import itertools
  t   = ['%sint%s' % (u,s) for s in 8, 16, 32, 64 for u in '', 'u'] + ['float32', 'float64', 'user_t']
  pad = 1 + max(map(len, t))
  for l,r in itertools.combinations_with_replacement(t, 2):
      print 'TRIAL(%s, %s, %s);' % tuple(i.rjust(pad) for i in (l,r,''))
#endif

#define DEDUCE(a,b,c) (boost::is_same<typename TypeDeductionHelper<a,b>::type, c>::value)

#define TRIAL_(a,b,c) \
  EXPECT_TRUE(DEDUCE(a,b,c)) \
      << "[" << #a << "] and [" << #b << "] produced [" << gi::GetTypeName<typename TypeDeductionHelper<a,b>::type>() \
      << "] expected [" << #c << "]";

#define TRIAL(a,b,c) \
  TRIAL_(a,b,c);      \
  TRIAL_(b,a,c);

struct user_t {};

TEST(TypeDeduction, All) {
  TRIAL(    int8,     int8,     int8);
  TRIAL(    int8,    uint8,    uint8);
  TRIAL(    int8,    int16,    int16);
  TRIAL(    int8,   uint16,   uint16);
  TRIAL(    int8,    int32,    int32);
  TRIAL(    int8,   uint32,   uint32);
  TRIAL(    int8,    int64,    int64);
  TRIAL(    int8,   uint64,   uint64);
  TRIAL(    int8,  float32,  float32);
  TRIAL(    int8,  float64,  float64);
  TRIAL(    int8,   user_t,   user_t);
  TRIAL(   uint8,    uint8,    uint8);
  TRIAL(   uint8,    int16,    int16);
  TRIAL(   uint8,   uint16,   uint16);
  TRIAL(   uint8,    int32,    int32);
  TRIAL(   uint8,   uint32,   uint32);
  TRIAL(   uint8,    int64,    int64);
  TRIAL(   uint8,   uint64,   uint64);
  TRIAL(   uint8,  float32,  float32);
  TRIAL(   uint8,  float64,  float64);
  TRIAL(   uint8,   user_t,   user_t);
  TRIAL(   int16,    int16,    int16);
  TRIAL(   int16,   uint16,   uint16);
  TRIAL(   int16,    int32,    int32);
  TRIAL(   int16,   uint32,   uint32);
  TRIAL(   int16,    int64,    int64);
  TRIAL(   int16,   uint64,   uint64);
  TRIAL(   int16,  float32,  float32);
  TRIAL(   int16,  float64,  float64);
  TRIAL(   int16,   user_t,   user_t);
  TRIAL(  uint16,   uint16,   uint16);
  TRIAL(  uint16,    int32,    int32);
  TRIAL(  uint16,   uint32,   uint32);
  TRIAL(  uint16,    int64,    int64);
  TRIAL(  uint16,   uint64,   uint64);
  TRIAL(  uint16,  float32,  float32);
  TRIAL(  uint16,  float64,  float64);
  TRIAL(  uint16,   user_t,   user_t);
  TRIAL(   int32,    int32,    int32);
  TRIAL(   int32,   uint32,   uint32);
  TRIAL(   int32,    int64,    int64);
  TRIAL(   int32,   uint64,   uint64);
  TRIAL(   int32,  float32,  float32);
  TRIAL(   int32,  float64,  float64);
  TRIAL(   int32,   user_t,   user_t);
  TRIAL(  uint32,   uint32,   uint32);
  TRIAL(  uint32,    int64,    int64);
  TRIAL(  uint32,   uint64,   uint64);
  TRIAL(  uint32,  float32,  float32);
  TRIAL(  uint32,  float64,  float64);
  TRIAL(  uint32,   user_t,   user_t);
  TRIAL(   int64,    int64,    int64);
  TRIAL(   int64,   uint64,   uint64);
  TRIAL(   int64,  float32,  float32);
  TRIAL(   int64,  float64,  float64);
  TRIAL(   int64,   user_t,   user_t);
  TRIAL(  uint64,   uint64,   uint64);
  TRIAL(  uint64,  float32,  float32);
  TRIAL(  uint64,  float64,  float64);
  TRIAL(  uint64,   user_t,   user_t);
  TRIAL( float32,  float32,  float32);
  TRIAL( float32,  float64,  float64);
  TRIAL( float32,   user_t,   user_t);
  TRIAL( float64,  float64,  float64);
  TRIAL( float64,   user_t,   user_t);
  TRIAL(  user_t,   user_t,   user_t);

  // and finally two different user-made types. this should be a compile error, so it's commented out.
  //struct user2_t {};
  //TRIAL(  user_t,   user2_t,   user_t);
}
