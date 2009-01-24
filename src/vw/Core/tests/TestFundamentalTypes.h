#include <cxxtest/TestSuite.h>

#include <vw/Core/FundamentalTypes.h>
#include <sstream>

template <typename T>
T smallest() {return vw::ScalarTypeLimits<T>::smallest();}

template <typename T>
T lowest() {return vw::ScalarTypeLimits<T>::lowest();}

template <typename T>
T highest() {return vw::ScalarTypeLimits<T>::highest();}

class TestFundamentalTypes : public CxxTest::TestSuite
{
public:
  void test_Ranges() {
      // TODO: commented out some checks because there's no reliable way to declare 
      // a 64-bit literal. this should probably be extended to check whether
      // -ffast-math (... broken math) is enabled.
      TS_ASSERT_EQUALS(smallest<vw::int8>(),    1);
      TS_ASSERT_EQUALS(smallest<vw::int16>(),   1);
      TS_ASSERT_EQUALS(smallest<vw::int32>(),   1);
      TS_ASSERT_EQUALS(smallest<vw::int64>(),   1);
      TS_ASSERT_EQUALS(smallest<vw::uint8>(),   1);
      TS_ASSERT_EQUALS(smallest<vw::uint16>(),  1);
      TS_ASSERT_EQUALS(smallest<vw::uint32>(),  1);
      TS_ASSERT_EQUALS(smallest<vw::uint64>(),  1);
      TS_ASSERT_DELTA( smallest<vw::float32>(), 0, 0.0000001);
      TS_ASSERT_DELTA( smallest<vw::float64>(), 0, 0.0000001);

      TS_ASSERT_EQUALS(lowest<vw::int8>(),    -128L);
      TS_ASSERT_EQUALS(lowest<vw::int16>(),   -32768L);
      TS_ASSERT_EQUALS(lowest<vw::int32>(),   -2147483648L);
      //TS_ASSERT_EQUALS(lowest<vw::int64>(),   -9223372036854775808LL);
      TS_ASSERT_EQUALS(lowest<vw::uint8>(),   0UL);
      TS_ASSERT_EQUALS(lowest<vw::uint16>(),  0UL);
      TS_ASSERT_EQUALS(lowest<vw::uint32>(),  0UL);
      TS_ASSERT_EQUALS(lowest<vw::uint64>(),  0UL);
      //TS_ASSERT_DELTA( lowest<vw::float32>(), 0, );
      //TS_ASSERT_DELTA( lowest<vw::float64>(), 0, );

      TS_ASSERT_EQUALS(highest<vw::int8>(),    127L);
      TS_ASSERT_EQUALS(highest<vw::int16>(),   32767L);
      TS_ASSERT_EQUALS(highest<vw::int32>(),   2147483647L);
      //TS_ASSERT_EQUALS(highest<vw::int64>(),   9223372036854775807LL);
      TS_ASSERT_EQUALS(highest<vw::uint8>(),   255UL);
      TS_ASSERT_EQUALS(highest<vw::uint16>(),  65535UL);
      TS_ASSERT_EQUALS(highest<vw::uint32>(),  4294967295UL);
      //TS_ASSERT_EQUALS(highest<vw::uint64>(),  18446744073709551615ULL);
      //TS_ASSERT_DELTA( highest<vw::float32>(), 0, );
      //TS_ASSERT_DELTA( highest<vw::float64>(), 0, );
  }

  void test_numeric() {
      std::stringstream s;
      s << vw::_numeric(vw::uint8(42)) << vw::_numeric(vw::int8(-42));
      TS_ASSERT_SAME_DATA(s.str().c_str(), "42-42", 6);
      TS_ASSERT_DELTA(vw::_numeric(4.2), 4.2, 0.000001);
  }

  void test_FundamentalType() {
      vw::FundamentalTypeClass<vw::int32> t, t2(42);
      TS_ASSERT(t ==  0);
      TS_ASSERT(t2 == 42);
  }

};
