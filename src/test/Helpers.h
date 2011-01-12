// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_TESTS_CONFIG_TEST_H__
#define __VW_TESTS_CONFIG_TEST_H__

#include <gtest/gtest.h>
#include <cmath>
#include <complex>
#include <string>
#include <boost/function.hpp>
#include <queue>
#include <cstdlib>

#include <vw/config.h>
#include <vw/Core/Log.h>
#include <vw/Math/Vector.h>
#include <vw/Math/Matrix.h>
#include <vw/Image/PixelTypeInfo.h>
#include <vw/Image/PixelMath.h>

#if defined(VW_ENABLE_EXCEPTIONS) && (VW_ENABLE_EXCEPTIONS==1)
#define HAS_EXCEPTIONS(x) x
#else
#define HAS_EXCEPTIONS(x) DISABLED_ ## x
#endif

#if defined(VW_ENABLE_CONFIG_FILE) && (VW_ENABLE_CONFIG_FILE==1)
#define HAS_CONFIG_FILE(x) x
#else
#define HAS_CONFIG_FILE(x) DISABLED_ ## x
#endif

namespace gi = ::testing::internal;

namespace vw { namespace test { }}
namespace t  = vw::test;

namespace vw {
  namespace test {

using namespace ::testing;

#ifndef TEST_OBJDIR
#error TEST_OBJDIR is not defined! Define it before including this header.
#endif

// Create a temporary filename that is unlinked when constructed and destructed
class UnlinkName : public std::string {
  public:
    UnlinkName() {}
    UnlinkName(const std::string& base, const std::string& directory=TEST_OBJDIR);
    UnlinkName(const char *base,        const std::string& directory=TEST_OBJDIR);
    ~UnlinkName();
};

// A getenv with a default value
std::string getenv2(const char *key, const std::string& Default);

// reduce the damage from using gtest internal bits, and make sure uint8 is
// seen as numeric.
template <typename T>
gi::String format(const T& x) {
  return gi::FormatForFailureMessage(_numeric(x));
}


// A version of std::mismatch that returns the set of differences rather than
// just the first
template <class InputIterator1, class InputIterator2, class Pred>
void mismatch_queue(InputIterator1 first1, InputIterator1 last1, InputIterator2 first2,
                    std::queue<std::pair<InputIterator1, InputIterator2> >& answer,
                    const Pred& p)
{
  while ( first1!=last1 ) {
    if (!p(*first1, *first2))
      answer.push(std::make_pair(first1, first2));
    ++first1; ++first2;
  }
}

template <typename ImplT>
class CmpWorker {
  public:
    inline ImplT& impl() { return static_cast<ImplT&>(*this); }
    inline ImplT const& impl() const { return static_cast<ImplT const&>(*this); }

    template <typename T1, typename T2>
    bool operator()(const T1& e, const T2& a) const {return impl()(e, a);}

    template <typename T1, typename T2>
    Message what(const std::string& ename, const std::string& aname, const T1& e, const T2& a) const {
      return impl().what(ename, aname, e, a);
    }
};

class CmpEqual : public CmpWorker<CmpEqual> {
  public:
    template <typename T1, typename T2>
    bool operator()(const T1& a, const T2& b) const { return a == b; }

    template <typename T1, typename T2>
    Message what(const std::string& ename, const std::string& aname, const T1& e, const T2& a) const {
      return Message(gi::EqFailure(ename.c_str(), aname.c_str(), t::format(e), t::format(a), false).message());
    }
};

class CmpTypeEqual : public CmpWorker<CmpTypeEqual> {
  public:
    template <typename T1, typename T2>
    bool operator()(const T1& a, const T2& b) const {
      return boost::is_same<T1,T2>::value && a == b;
    }

    template <typename T1, typename T2>
    Message what(const std::string& ename, const std::string& aname, const T1& e, const T2& a) const {
      if (!boost::is_same<T1,T2>::value)
        return Message() << ename << " and " << aname << " are not the same type";
      return Message() << gi::EqFailure(ename.c_str(), aname.c_str(), t::format(e), t::format(a), false);
    }
};

template <typename DeltaT>
class CmpNear : public CmpWorker<CmpNear<DeltaT> > {
    const char *dexpr;
    const DeltaT& delta;
  public:
    CmpNear(const char* dexpr, const DeltaT& delta) : dexpr(dexpr), delta(delta) {}

    template <typename T1, typename T2>
    bool operator()(const T1& a, const T2& b) const { return value_diff(a, b) <= delta; }

    template <typename T1, typename T2>
    Message what(const std::string& ename, const std::string& aname, const T1& e, const T2& a) const {
      Message msg;
      msg << "The difference between "
          << ename   << " and " << aname
          << " is "  << t::format(value_diff(e, a))
          << ", which exceeds " << dexpr << ", where\n"
          << ename << " evaluates to " << t::format(e) << ",\n"
          << aname << " evaluates to " << t::format(a)
          << ", and\n" << dexpr << " evaluates to " << t::format(delta) << ".";
      return msg;
    }
};

template <typename DeltaT>
CmpNear<DeltaT> cmp_near(const char* dexpr, const DeltaT& delta) {
  return CmpNear<DeltaT>(dexpr, delta);
}

#if 0
template <typename ExpectT>
struct CmpULP : public CmpWorker<CmpULP<ExpectT> > {
  BOOST_STATIC_ASSERT(boost::is_floating_point<ExpectT>::value);
  public:
    template <typename ActualT>
    bool operator()(const ExpectT& a, const ActualT& b) const {
      const gi::FloatingPoint<ExpectT> lhs(a), rhs(b);
      return lhs.AlmostEquals(rhs);
    }

    template <typename ActualT>
    Message what(const std::string& ename, const std::string& aname, const ExpectT& e, const ActualT& a) const {
      std::ostringstream es, as;
      es << std::setprecision(std::numeric_limits<ExpectT>::digits10 + 2) << e;
      as << std::setprecision(std::numeric_limits<ExpectT>::digits10 + 2) << a;
      return Message() << gi::EqFailure(ename.c_str(), aname.c_str(), es.str(), as.str(), false);
    }
};
#endif

template <typename CmpT>
class _CheckOne {
    const CmpT& cmp;
  public:
    _CheckOne() : cmp(CmpT()) {}
    _CheckOne(const CmpT& cmp) : cmp(cmp) {}

    template <typename ExpectT, typename ActualT>
    AssertionResult operator()(const char* ename, const char* aname, const ExpectT& e, const ActualT& a) const
    {
      if (cmp(e, a))
        return AssertionSuccess();
      return AssertionFailure() << cmp.what(ename, aname, e, a);
    }
};
template <typename CmpT>
_CheckOne<CmpT> check_one(const CmpT& cmp) {
  return _CheckOne<CmpT>(cmp);
}

template <typename CmpT>
class _CheckRange {
    const CmpT& cmp;
  public:
    _CheckRange() : cmp(CmpT()) {}
    _CheckRange(const CmpT& cmp) : cmp(cmp) {}

    template <typename Range1T, typename Range2T>
    AssertionResult operator()(const char* ename, const char* aname,
                               const Range1T& e, const Range2T& a) const
    {
      if (a.size() != e.size())
        return AssertionFailure()
                << "Iterator ranges (" << ename << ") and (" << aname
                << ") represent different-sized ranges.";

      typedef std::pair<typename Range1T::const_iterator, typename Range2T::const_iterator> range_t;

      std::queue<range_t> ret;
      mismatch_queue(e.begin(), e.end(), a.begin(), ret, cmp);

      if (ret.empty())
        return AssertionSuccess();

      Message msg;
      bool need_newline = false;
      while (!ret.empty()) {
        const range_t& r = ret.front();
        const std::string idx = "[" + stringify(std::distance(e.begin(), r.first)) + "]";
        if (need_newline)
          msg << std::endl;
        msg << cmp.what(ename + idx, aname + idx, *(r.first), *(r.second));
        need_newline = true;
        ret.pop();
      }
      return AssertionFailure(msg);
    }

    template <typename Iter1T, typename Iter2T>
    AssertionResult operator()(const char* e0name, const char* /*e1name*/,
                               const char* a0name, const char* /*a1name*/,
                               const Iter1T& e0, const Iter1T& e1,
                               const Iter2T& a0, const Iter2T& a1) const
    {
      return this->operator()(e0name, a0name, boost::make_iterator_range(e0, e1), boost::make_iterator_range(a0, a1));
    }
};

template <typename CmpT>
_CheckRange<CmpT> check_range(const CmpT& cmp) {
  return _CheckRange<CmpT>(cmp);
}

template <typename CmpT>
class _CheckMatrix {
    const CmpT& cmp;
  public:
    _CheckMatrix() : cmp(CmpT()) {}
    _CheckMatrix(const CmpT& cmp) : cmp(cmp) {}

    template <typename T1, typename T2>
    AssertionResult operator()(const char* ename, const char* aname, const MatrixBase<T1>& expect, const MatrixBase<T2>& actual)
    {
      const size_t rows = expect.impl().rows(),
                   cols = expect.impl().cols();

      if ( rows != actual.impl().rows() || cols != actual.impl().cols() )
        return AssertionFailure() << "Cannot compare " << ename << " and " << aname << ": Different size: "
                                  << rows << "x" << cols << " != "
                                  << actual.impl().rows() << "x" << actual.impl().cols();

      Message msg;
      bool failed = false;
      bool need_newline = false;
      for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
          if (need_newline)
            msg << std::endl;
          const std::string idx = "[" + stringify(i) +  ","  + stringify(j) + "]";
          if (!cmp(expect.impl()(i,j), actual.impl()(i,j))) {
            failed = true;
            msg << cmp.what(ename + idx, aname + idx, expect.impl()(i,j), actual.impl()(i,j));
          }
        }
      }

      if (failed)
        return AssertionFailure(msg);
      return AssertionSuccess();
    }
};

template <typename CmpT>
_CheckMatrix<CmpT> check_matrix(const CmpT& cmp) {
  return _CheckMatrix<CmpT>(cmp);
}

#define EXPECT_RANGE_EQ(expect0, expect1, actual0, actual1) \
  EXPECT_PRED_FORMAT4(t::check_range(t::CmpEqual()), expect0, expect1, actual0, actual1)
#define ASSERT_RANGE_EQ(expect0, expect1, actual0, actual1) \
  ASSERT_PRED_FORMAT4(t::check_range(t::CmpEqual()), expect0, expect1, actual0, actual1)
#define EXPECT_RANGE_NEAR(expect0, expect1, actual0, actual1, delta) \
  EXPECT_PRED_FORMAT4(t::check_range(t::cmp_near(#delta, delta)), expect0, expect1, actual0, actual1)
#define ASSERT_RANGE_NEAR(expect0, expect1, actual0, actual1, delta) \
  ASSERT_PRED_FORMAT4(t::check_range(t::cmp_near(#delta, delta)), expect0, expect1, actual0, actual1)

#define EXPECT_MATRIX_EQ(expect, actual, delta)\
  EXPECT_PRED_FORMAT2(t::check_matrix(t::CmpEqual()), expect, actual)
#define ASSERT_MATRIX_EQ(expect, actual, delta)\
  ASSERT_PRED_FORMAT2(t::check_matrix(t::CmpEqual()), expect, actual)
#define EXPECT_MATRIX_NEAR(expect, actual, delta)\
  EXPECT_PRED_FORMAT2(t::check_matrix(t::cmp_near(#delta, delta)), expect, actual)
#define ASSERT_MATRIX_NEAR(expect, actual, delta)\
  ASSERT_PRED_FORMAT2(t::check_matrix(t::cmp_near(#delta, delta)), expect, actual)

#define EXPECT_VECTOR_EQ(expect, actual)\
  EXPECT_PRED_FORMAT2(t::check_range(t::CmpEqual()), expect, actual)
#define ASSERT_VECTOR_EQ(expect, actual)\
  ASSERT_PRED_FORMAT2(t::check_range(t::CmpEqual()), expect, actual)
#define EXPECT_VECTOR_NEAR(expect, actual, delta)\
  EXPECT_PRED_FORMAT2(t::check_range(t::cmp_near(#delta, delta)), expect, actual)
#define ASSERT_VECTOR_NEAR(expect, actual, delta)\
  ASSERT_PRED_FORMAT2(t::check_range(t::cmp_near(#delta, delta)), expect, actual)

#define EXPECT_TYPE_EQ( expect, actual )\
  EXPECT_PRED_FORMAT2(t::check_one(t::CmpTypeEqual()), expect, actual )
#define ASSERT_TYPE_EQ( expect, actual )\
  ASSERT_PRED_FORMAT2(t::check_one(t::CmpTypeEqual()), expect, actual )

#define EXPECT_PIXEL_NEAR(expect, actual, delta)\
  EXPECT_PRED_FORMAT2(t::check_one(t::cmp_near(#delta, delta)), expect, actual)
#define ASSERT_PIXEL_NEAR(expect, actual, delta)\
  ASSERT_PRED_FORMAT2(t::check_one(t::cmp_near(#delta, delta)), expect, actual)
#define EXPECT_PIXEL_EQ(expect, actual)\
  EXPECT_PRED_FORMAT2(t::check_one(t::CmpTypeEqual()), expect, actual)
#define ASSERT_PIXEL_EQ(expect, actual)\
  ASSERT_PRED_FORMAT2(t::check_one(t::CmpTypeEqual()), expect, actual)

#define EXPECT_VW_EQ(expect, actual)\
  EXPECT_PRED_FORMAT2(t::check_one(t::CmpTypeEqual()), expect, actual)
#define ASSERT_VW_EQ(expect, actual)\
  ASSERT_PRED_FORMAT2(t::check_one(t::CmpTypeEqual()), expect, actual)

// DEPRECATED
#define EXPECT_MATRIX_FLOAT_EQ(e, a)  EXPECT_MATRIX_NEAR(e, a, 1e20)
#define EXPECT_MATRIX_DOUBLE_EQ(e, a) EXPECT_MATRIX_NEAR(e, a, 1e45)
#define EXPECT_COMPLEX_MATRIX_NEAR(e, a, d)  EXPECT_MATRIX_NEAR(e, a, d)
#define EXPECT_VECTOR_FLOAT_EQ(e, a)  EXPECT_VECTOR_NEAR(e, a, 1e20)
#define EXPECT_VECTOR_DOUBLE_EQ(e, a) EXPECT_VECTOR_NEAR(e, a, 1e45)

template <typename T1, typename T2>
T1 value_diff_round(const T2& x) {
  if (boost::is_integral<T1>::value)
    return boost::numeric_cast<T1>(::ceil(x));
  else
    return boost::numeric_cast<T1>(x);
}

template <typename ElemT>
typename boost::enable_if<boost::is_integral<ElemT>, ElemT>::type
value_diff_helper(const ElemT& a, const ElemT& b) {
  return boost::numeric_cast<ElemT>(::abs(a - b));
}

template <typename ElemT>
typename boost::enable_if<boost::is_floating_point<ElemT>, ElemT>::type
value_diff_helper(const ElemT& a, const ElemT& b) {
  return boost::numeric_cast<ElemT>(::fabs(a - b));
}

template <typename ElemT, typename Elem2T>
typename PixelChannelType<ElemT>::type value_diff(const vw::PixelMathBase<ElemT>& a, const vw::PixelMathBase<Elem2T>& b) {
  BOOST_STATIC_ASSERT((boost::is_same<ElemT, Elem2T>::value));
  typedef typename PixelChannelType<ElemT>::type channel_type;
  channel_type acc = channel_type();
  ElemT diff = a - b;
  for( size_t c=0; c < PixelNumChannels<ElemT>::value; ++c ) {
    const channel_type& x = compound_select_channel<const channel_type&>(diff,c);
    acc += x * x;
  }
  return value_diff_round<channel_type>(::sqrt(acc));
}

template <typename ElemT>
ElemT value_diff(const std::complex<ElemT>& a, const std::complex<ElemT>& b) {
  return boost::numeric_cast<ElemT>(std::abs(a - b));
}

template <typename T1, typename T2>
struct both_are_arithmetic : boost::mpl::and_<boost::is_arithmetic<T1>, boost::is_arithmetic<T2> > {};

template <typename T1, typename T2>
typename boost::enable_if<both_are_arithmetic<T1,T2>, typename DifferenceType<T1,T2>::type>::type
value_diff(const T1& a, const T2& b) {
  BOOST_STATIC_ASSERT(boost::is_arithmetic<T1>::value);
  BOOST_STATIC_ASSERT(boost::is_arithmetic<T2>::value);
  typedef typename DifferenceType<T1,T2>::type diff_t;
  return value_diff_helper<diff_t>(a, b);
}

}} // namespace t

#endif
