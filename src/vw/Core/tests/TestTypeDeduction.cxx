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

// These tests generated with this python script (requires >= python 2.7):
#if 0
  import itertools
  import sys

  assert len(sys.argv) > 1, 'Usage: %s <yes|no>\nArgument determines whether the standard types are used'
  assert sys.argv[1] in ('yes', 'no'), 'Usage: %s <yes|no>\nArgument determines whether the standard types are used'
  USE_STANDARD = True if sys.argv[1] == 'yes' else False

# In increasing order of preference
  if USE_STANDARD:
      typelist = ['char', 'signed char', 'unsigned char', 'short', 'unsigned short', 'int', 'unsigned int', 'long', 'unsigned long', 'long long', 'unsigned long long', 'float', 'double', 'long double']
  else:
      typelist = ['%sint%s' % (u,s) for s in 8, 16, 32, 64 for u in '', 'u'] + ['float32', 'float64', 'user_t']

  def prefer(a,b):
      return a if typelist.index(a) > typelist.index(b) else b

  def note(a,b):
      if USE_STANDARD:
          p = prefer(a,b)
          return ' /*int*/' if p != 'int' and prefer('int', p) == 'int' else ''
      else:
          return ''

  pad = 1 + max(map(len, typelist))
  for l,r in itertools.combinations_with_replacement(typelist, 2):
      print 'TRIAL(%s,%s,%s%s);' % tuple([i.rjust(pad) for i in (l,r,prefer(l,r))] + [note(l,r)])
#endif


#define DEDUCE(a,b,c) (boost::is_same<TypeDeductionHelper<a,b>::type, c>::value)

#define TRIAL_(a,b,c) do {\
  EXPECT_TRUE(DEDUCE(a,b,c)) \
      << "[" << #a << "] and [" << #b << "] produced [" << gi::GetTypeName<TypeDeductionHelper<a,b>::type>() \
      << "] expected [" << #c << "]"; \
} while (0)

#define TRIAL(a,b,c) do {\
  TRIAL_(a,b,c);         \
  TRIAL_(b,a,c);         \
} while (0)

struct user_t {};
struct user2_t {};

// This should be a compile error, so it's commented out (duplicate index)
//namespace vw { namespace core { namespace detail {
//  _VW_INTERNAL_TYPE_DEDUCTION(user_t, TypeDeductionIndex<int8>::value);
//}}}

TEST(TypeDeduction, VW) {
  TRIAL(    int8,    int8,    int8);
  TRIAL(    int8,   uint8,   uint8);
  TRIAL(    int8,   int16,   int16);
  TRIAL(    int8,  uint16,  uint16);
  TRIAL(    int8,   int32,   int32);
  TRIAL(    int8,  uint32,  uint32);
  TRIAL(    int8,   int64,   int64);
  TRIAL(    int8,  uint64,  uint64);
  TRIAL(    int8, float32, float32);
  TRIAL(    int8, float64, float64);
  TRIAL(    int8,  user_t,  user_t);
  TRIAL(   uint8,   uint8,   uint8);
  TRIAL(   uint8,   int16,   int16);
  TRIAL(   uint8,  uint16,  uint16);
  TRIAL(   uint8,   int32,   int32);
  TRIAL(   uint8,  uint32,  uint32);
  TRIAL(   uint8,   int64,   int64);
  TRIAL(   uint8,  uint64,  uint64);
  TRIAL(   uint8, float32, float32);
  TRIAL(   uint8, float64, float64);
  TRIAL(   uint8,  user_t,  user_t);
  TRIAL(   int16,   int16,   int16);
  TRIAL(   int16,  uint16,  uint16);
  TRIAL(   int16,   int32,   int32);
  TRIAL(   int16,  uint32,  uint32);
  TRIAL(   int16,   int64,   int64);
  TRIAL(   int16,  uint64,  uint64);
  TRIAL(   int16, float32, float32);
  TRIAL(   int16, float64, float64);
  TRIAL(   int16,  user_t,  user_t);
  TRIAL(  uint16,  uint16,  uint16);
  TRIAL(  uint16,   int32,   int32);
  TRIAL(  uint16,  uint32,  uint32);
  TRIAL(  uint16,   int64,   int64);
  TRIAL(  uint16,  uint64,  uint64);
  TRIAL(  uint16, float32, float32);
  TRIAL(  uint16, float64, float64);
  TRIAL(  uint16,  user_t,  user_t);
  TRIAL(   int32,   int32,   int32);
  TRIAL(   int32,  uint32,  uint32);
  TRIAL(   int32,   int64,   int64);
  TRIAL(   int32,  uint64,  uint64);
  TRIAL(   int32, float32, float32);
  TRIAL(   int32, float64, float64);
  TRIAL(   int32,  user_t,  user_t);
  TRIAL(  uint32,  uint32,  uint32);
  TRIAL(  uint32,   int64,   int64);
  TRIAL(  uint32,  uint64,  uint64);
  TRIAL(  uint32, float32, float32);
  TRIAL(  uint32, float64, float64);
  TRIAL(  uint32,  user_t,  user_t);
  TRIAL(   int64,   int64,   int64);
  TRIAL(   int64,  uint64,  uint64);
  TRIAL(   int64, float32, float32);
  TRIAL(   int64, float64, float64);
  TRIAL(   int64,  user_t,  user_t);
  TRIAL(  uint64,  uint64,  uint64);
  TRIAL(  uint64, float32, float32);
  TRIAL(  uint64, float64, float64);
  TRIAL(  uint64,  user_t,  user_t);
  TRIAL( float32, float32, float32);
  TRIAL( float32, float64, float64);
  TRIAL( float32,  user_t,  user_t);
  TRIAL( float64, float64, float64);
  TRIAL( float64,  user_t,  user_t);
  TRIAL(  user_t,  user_t,  user_t);

  // and finally two different user-made types. this should be a compile error, so it's commented out.
  //TRIAL(  user_t,   user2_t,   user_t);
}

TEST(TypeDeduction, Standard) {
  // The commented out ints are the places where C++ does automatic int promotion, and we don't
  TRIAL(               char,               char,               char /*int*/);
  TRIAL(               char,        signed char,        signed char /*int*/);
  TRIAL(               char,      unsigned char,      unsigned char /*int*/);
  TRIAL(               char,              short,              short /*int*/);
  TRIAL(               char,     unsigned short,     unsigned short /*int*/);
  TRIAL(               char,                int,                int);
  TRIAL(               char,       unsigned int,       unsigned int);
  TRIAL(               char,               long,               long);
  TRIAL(               char,      unsigned long,      unsigned long);
  TRIAL(               char,              float,              float);
  TRIAL(               char,             double,             double);
  TRIAL(               char,        long double,        long double);
  TRIAL(        signed char,        signed char,        signed char /*int*/);
  TRIAL(        signed char,      unsigned char,      unsigned char /*int*/);
  TRIAL(        signed char,              short,              short /*int*/);
  TRIAL(        signed char,     unsigned short,     unsigned short /*int*/);
  TRIAL(        signed char,                int,                int);
  TRIAL(        signed char,       unsigned int,       unsigned int);
  TRIAL(        signed char,               long,               long);
  TRIAL(        signed char,      unsigned long,      unsigned long);
  TRIAL(        signed char,              float,              float);
  TRIAL(        signed char,             double,             double);
  TRIAL(        signed char,        long double,        long double);
  TRIAL(      unsigned char,      unsigned char,      unsigned char /*int*/);
  TRIAL(      unsigned char,              short,              short /*int*/);
  TRIAL(      unsigned char,     unsigned short,     unsigned short /*int*/);
  TRIAL(      unsigned char,                int,                int);
  TRIAL(      unsigned char,       unsigned int,       unsigned int);
  TRIAL(      unsigned char,               long,               long);
  TRIAL(      unsigned char,      unsigned long,      unsigned long);
  TRIAL(      unsigned char,              float,              float);
  TRIAL(      unsigned char,             double,             double);
  TRIAL(      unsigned char,        long double,        long double);
  TRIAL(              short,              short,              short /*int*/);
  TRIAL(              short,     unsigned short,     unsigned short /*int*/);
  TRIAL(              short,                int,                int);
  TRIAL(              short,       unsigned int,       unsigned int);
  TRIAL(              short,               long,               long);
  TRIAL(              short,      unsigned long,      unsigned long);
  TRIAL(              short,              float,              float);
  TRIAL(              short,             double,             double);
  TRIAL(              short,        long double,        long double);
  TRIAL(     unsigned short,     unsigned short,     unsigned short /*int*/);
  TRIAL(     unsigned short,                int,                int);
  TRIAL(     unsigned short,       unsigned int,       unsigned int);
  TRIAL(     unsigned short,               long,               long);
  TRIAL(     unsigned short,      unsigned long,      unsigned long);
  TRIAL(     unsigned short,              float,              float);
  TRIAL(     unsigned short,             double,             double);
  TRIAL(     unsigned short,        long double,        long double);
  TRIAL(                int,                int,                int);
  TRIAL(                int,       unsigned int,       unsigned int);
  TRIAL(                int,               long,               long);
  TRIAL(                int,      unsigned long,      unsigned long);
  TRIAL(                int,              float,              float);
  TRIAL(                int,             double,             double);
  TRIAL(                int,        long double,        long double);
  TRIAL(       unsigned int,       unsigned int,       unsigned int);
  TRIAL(       unsigned int,               long,               long);
  TRIAL(       unsigned int,      unsigned long,      unsigned long);
  TRIAL(       unsigned int,              float,              float);
  TRIAL(       unsigned int,             double,             double);
  TRIAL(       unsigned int,        long double,        long double);
  TRIAL(               long,               long,               long);
  TRIAL(               long,      unsigned long,      unsigned long);
  TRIAL(               long,              float,              float);
  TRIAL(               long,             double,             double);
  TRIAL(               long,        long double,        long double);
  TRIAL(      unsigned long,      unsigned long,      unsigned long);
  TRIAL(      unsigned long,              float,              float);
  TRIAL(      unsigned long,             double,             double);
  TRIAL(      unsigned long,        long double,        long double);
  TRIAL(              float,              float,              float);
  TRIAL(              float,             double,             double);
  TRIAL(              float,        long double,        long double);
  TRIAL(             double,             double,             double);
  TRIAL(             double,        long double,        long double);
  TRIAL(        long double,        long double,        long double);

#if defined(BOOST_HAS_LONG_LONG)
  TRIAL(               char,          long long,          long long);
  TRIAL(               char, unsigned long long, unsigned long long);
  TRIAL(        signed char,          long long,          long long);
  TRIAL(        signed char, unsigned long long, unsigned long long);
  TRIAL(      unsigned char,          long long,          long long);
  TRIAL(      unsigned char, unsigned long long, unsigned long long);
  TRIAL(              short,          long long,          long long);
  TRIAL(              short, unsigned long long, unsigned long long);
  TRIAL(     unsigned short,          long long,          long long);
  TRIAL(     unsigned short, unsigned long long, unsigned long long);
  TRIAL(                int,          long long,          long long);
  TRIAL(                int, unsigned long long, unsigned long long);
  TRIAL(       unsigned int,          long long,          long long);
  TRIAL(       unsigned int, unsigned long long, unsigned long long);
  TRIAL(               long,          long long,          long long);
  TRIAL(               long, unsigned long long, unsigned long long);
  TRIAL(      unsigned long,          long long,          long long);
  TRIAL(      unsigned long, unsigned long long, unsigned long long);
  TRIAL(          long long,          long long,          long long);
  TRIAL(          long long, unsigned long long, unsigned long long);
  TRIAL(          long long,              float,              float);
  TRIAL(          long long,             double,             double);
  TRIAL(          long long,        long double,        long double);
  TRIAL( unsigned long long, unsigned long long, unsigned long long);
  TRIAL( unsigned long long,              float,              float);
  TRIAL( unsigned long long,             double,             double);
  TRIAL( unsigned long long,        long double,        long double);
#endif

}
