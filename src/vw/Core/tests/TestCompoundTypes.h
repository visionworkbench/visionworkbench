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

// TestCompoundTypes.h
#include <cxxtest/TestSuite.h>

#include <vw/Core/CompoundTypes.h>

using namespace vw;

// A simple compound type.  We only test with this one type here, and thus 
// we do not exercise all of the specialized code paths.  However, those 
// do largely get exercized by the pixel math tests in the Image module.
template <class ChannelT>
class TestCompound {
  ChannelT values[2];
public:
  TestCompound( ChannelT a, ChannelT b ) { values[0]=a; values[1]=b; }
  ChannelT& operator[]( int i ) { return values[i]; }
  ChannelT const& operator[]( int i ) const { return values[i]; }
};

// Simple functions to test the compound_apply logic.
template <class T> T add( T const& a, T const& b ) { return a + b; }
template <class T> void add_in_place( T& a, T const& b ) { a += b; }
template <class T> T add_one( T const& val ) { return val + 1; }
template <class T> void add_one_in_place( T& val ) { val += 1; }

// A dummy type that is neither a scalar nor a compound type.
class DummyType {};

namespace vw {
  template<class ChannelT> struct CompoundChannelType<TestCompound<ChannelT> > { typedef ChannelT type; };
  template<class ChannelT> struct CompoundNumChannels<TestCompound<ChannelT> > { static const int32 value = 2; };
  template<class InT, class OutT> struct CompoundChannelCast<TestCompound<InT>, OutT> { typedef TestCompound<OutT> type; };
}

class TestCompoundTypes : public CxxTest::TestSuite
{
public:

  template <class T1, class T2>
  static bool is_of_type( T2 ) {
    return boost::is_same<T1,T2>::value;
  }

  void test_TestCompound() {
    TS_ASSERT(( boost::is_same<CompoundChannelType<TestCompound<float> >::type, float>::value ));
    TS_ASSERT_EQUALS( CompoundNumChannels<TestCompound<float> >::value, 2 );
    TS_ASSERT(( boost::is_same<CompoundChannelCast<TestCompound<float>,int>::type, TestCompound<int> >::value ));
  }
  
  void test_traits() {
    TS_ASSERT(( !IsCompound<uint8>::value ));
    TS_ASSERT(( !IsCompound<double>::value ));
    TS_ASSERT(( IsCompound<TestCompound<uint8> >::value ));
    TS_ASSERT(( IsCompound<TestCompound<double> >::value ));
    TS_ASSERT(( !IsCompound<DummyType>::value ));
    TS_ASSERT(( !IsCompound<const uint8>::value ));
    TS_ASSERT(( !IsCompound<const double>::value ));
    TS_ASSERT(( IsCompound<TestCompound<const uint8> >::value ));
    TS_ASSERT(( IsCompound<TestCompound<const double> >::value ));
    TS_ASSERT(( !IsCompound<const DummyType>::value ));
    
    TS_ASSERT(( IsScalarOrCompound<uint8>::value ));
    TS_ASSERT(( IsScalarOrCompound<double>::value ));
    TS_ASSERT(( IsScalarOrCompound<TestCompound<uint8> >::value ));
    TS_ASSERT(( IsScalarOrCompound<TestCompound<double> >::value ));
    TS_ASSERT(( !IsScalarOrCompound<DummyType>::value ));
    TS_ASSERT(( IsScalarOrCompound<const uint8>::value ));
    TS_ASSERT(( IsScalarOrCompound<const double>::value ));
    TS_ASSERT(( IsScalarOrCompound<const TestCompound<uint8> >::value ));
    TS_ASSERT(( IsScalarOrCompound<const TestCompound<double> >::value ));
    TS_ASSERT(( !IsScalarOrCompound<const DummyType>::value ));
    
    TS_ASSERT(( CompoundIsCompatible<double,double>::value ));
    TS_ASSERT(( CompoundIsCompatible<uint8,double>::value ));
    TS_ASSERT(( CompoundIsCompatible<TestCompound<double>,TestCompound<double> >::value ));
    TS_ASSERT(( CompoundIsCompatible<TestCompound<uint8>,TestCompound<double> >::value ));
    TS_ASSERT(( !CompoundIsCompatible<TestCompound<double>,double>::value ));
    TS_ASSERT(( !CompoundIsCompatible<double,TestCompound<double> >::value ));
    TS_ASSERT(( CompoundIsCompatible<const double,double>::value ));
    TS_ASSERT(( CompoundIsCompatible<const uint8,double>::value ));
    TS_ASSERT(( CompoundIsCompatible<const TestCompound<double>,TestCompound<double> >::value ));
    TS_ASSERT(( CompoundIsCompatible<const TestCompound<uint8>,TestCompound<double> >::value ));
    TS_ASSERT(( !CompoundIsCompatible<const TestCompound<double>,double>::value ));
    TS_ASSERT(( !CompoundIsCompatible<const double,TestCompound<double> >::value ));
    TS_ASSERT(( CompoundIsCompatible<double,const double>::value ));
    TS_ASSERT(( CompoundIsCompatible<uint8,const double>::value ));
    TS_ASSERT(( CompoundIsCompatible<TestCompound<double>,const TestCompound<double> >::value ));
    TS_ASSERT(( CompoundIsCompatible<TestCompound<uint8>,const TestCompound<double> >::value ));
    TS_ASSERT(( !CompoundIsCompatible<TestCompound<double>,const double>::value ));
    TS_ASSERT(( !CompoundIsCompatible<double,const TestCompound<double> >::value ));
    TS_ASSERT(( CompoundIsCompatible<const double,const double>::value ));
    TS_ASSERT(( CompoundIsCompatible<const uint8,const double>::value ));
    TS_ASSERT(( CompoundIsCompatible<const TestCompound<double>,const TestCompound<double> >::value ));
    TS_ASSERT(( CompoundIsCompatible<const TestCompound<uint8>,const TestCompound<double> >::value ));
    TS_ASSERT(( !CompoundIsCompatible<const TestCompound<double>,const double>::value ));
    TS_ASSERT(( !CompoundIsCompatible<const double,const TestCompound<double> >::value ));

    TS_ASSERT(( boost::is_same<CompoundAccumulatorType<uint8>::type,AccumulatorType<uint8>::type>::value ));
    TS_ASSERT(( boost::is_same<CompoundAccumulatorType<TestCompound<uint8> >::type,TestCompound<AccumulatorType<uint8>::type> >::value ));
    TS_ASSERT(( boost::is_same<CompoundAccumulatorType<const uint8>::type,AccumulatorType<uint8>::type>::value ));
    TS_ASSERT(( boost::is_same<CompoundAccumulatorType<const TestCompound<uint8> >::type,TestCompound<AccumulatorType<uint8>::type> >::value ));
  }

  void test_compound_select_channel() {
    uint8 vali = 3;
    TS_ASSERT_EQUALS( compound_select_channel<const uint8&>(vali,0), 3 );
    compound_select_channel<uint8&>(vali,0) = 5;
    TS_ASSERT_EQUALS( compound_select_channel<const uint8&>(vali,0), 5 );
    float valf = 4.0;
    TS_ASSERT_EQUALS( compound_select_channel<const float&>(valf,0), 4.0 );
    compound_select_channel<float&>(valf,0) = 6;
    TS_ASSERT_EQUALS( compound_select_channel<const float&>(valf,0), 6.0 );
    TestCompound<uint8> valci( 1, 2 );
    TS_ASSERT_EQUALS( compound_select_channel<const uint8&>(valci,0), 1 );
    TS_ASSERT_EQUALS( compound_select_channel<const uint8&>(valci,1), 2 );
    compound_select_channel<uint8&>(valci,0) = 3;
    compound_select_channel<uint8&>(valci,1) = 4;
    TS_ASSERT_EQUALS( compound_select_channel<const uint8&>(valci,0), 3 );
    TS_ASSERT_EQUALS( compound_select_channel<const uint8&>(valci,1), 4 );
    TestCompound<float> valcf( 2.0, 3.0 );
    TS_ASSERT_EQUALS( compound_select_channel<const float&>(valcf,0), 2.0 );
    TS_ASSERT_EQUALS( compound_select_channel<const float&>(valcf,1), 3.0 );
    compound_select_channel<float&>(valcf,0) = 3.0;
    compound_select_channel<float&>(valcf,1) = 4.0;
    TS_ASSERT_EQUALS( compound_select_channel<const float&>(valcf,0), 3.0 );
    TS_ASSERT_EQUALS( compound_select_channel<const float&>(valcf,1), 4.0 );
  }

  void test_binary_compound_apply() {
    uint8 ai=1, bi=2, ci=compound_apply(&add<uint8>, ai, bi);
    TS_ASSERT_EQUALS( ci, 3 );
    float af=1.0, bf=2.0, cf=compound_apply(&add<float>, af, bf);
    TS_ASSERT_EQUALS( cf, 3.0 );
    TestCompound<uint8> aci(1,2), bci(3,4), cci=compound_apply(&add<uint8>, aci, bci);
    TS_ASSERT_EQUALS( cci[0], 4 );
    TS_ASSERT_EQUALS( cci[1], 6 );
    TestCompound<float> acf(1,2), bcf(3,4), ccf=compound_apply(&add<float>, acf, bcf);
    TS_ASSERT_EQUALS( ccf[0], 4 );
    TS_ASSERT_EQUALS( ccf[1], 6 );
  }

  void test_binary_compound_apply_in_place() {
    uint8 ai=1, bi=2;
    compound_apply_in_place(&add_in_place<uint8>, ai, bi);
    TS_ASSERT_EQUALS( ai, 3 );
    float af=1.0, bf=2.0;
    compound_apply_in_place(&add_in_place<float>, af, bf);
    TS_ASSERT_EQUALS( af, 3.0 );
    TestCompound<uint8> aci(1,2), bci(3,4);
    compound_apply_in_place(&add_in_place<uint8>, aci, bci);
    TS_ASSERT_EQUALS( aci[0], 4 );
    TS_ASSERT_EQUALS( aci[1], 6 );
    TestCompound<float> acf(1,2), bcf(3,4);
    compound_apply_in_place(&add_in_place<float>, acf, bcf);
    TS_ASSERT_EQUALS( acf[0], 4 );
    TS_ASSERT_EQUALS( acf[1], 6 );
  }

  void test_unary_compound_apply() {
    uint8 ai=1, bi=compound_apply(&add_one<uint8>, ai);
    TS_ASSERT_EQUALS( bi, 2 );
    uint8 af=1.0, bf=compound_apply(&add_one<float>, af);
    TS_ASSERT_EQUALS( bf, 2.0 );
    TestCompound<uint8> aci(1,2), bci=compound_apply(&add_one<uint8>, aci);
    TS_ASSERT_EQUALS( bci[0], 2 );
    TS_ASSERT_EQUALS( bci[1], 3 );
    TestCompound<float> acf(1,2), bcf=compound_apply(&add_one<float>, acf);
    TS_ASSERT_EQUALS( bcf[0], 2 );
    TS_ASSERT_EQUALS( bcf[1], 3 );
  }

  void test_unary_compound_apply_in_place() {
    uint8 vali = 0;
    compound_apply_in_place( &add_one_in_place<uint8>, vali );
    TS_ASSERT_EQUALS( vali, 1.0 );
    double valf = 0.0;
    compound_apply_in_place( &add_one_in_place<double>, valf );
    TS_ASSERT_EQUALS( valf, 1.0 );
    TestCompound<int> valc( 1, 2 );
    compound_apply_in_place( &add_one_in_place<int>, valc );
    TS_ASSERT_EQUALS( valc[0], 2);
    TS_ASSERT_EQUALS( valc[1], 3);
  }

};
