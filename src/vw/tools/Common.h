// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_TOOLS_COMMON_H__
#define __VW_TOOLS_COMMON_H__

#include <vw/Core/Exception.h>

#include <boost/preprocessor.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/algorithm/string/case_conv.hpp>
#include <map>

#define VW_PP_ELEM_TO_ELEM_STRING(s, data, elem) elem, BOOST_PP_STRINGIZE(elem)
#define VW_PP_ELEM_TO_STRING_ELEM(s, data, elem) BOOST_PP_STRINGIZE(elem), elem

namespace vw {
namespace tools {

  VW_DEFINE_EXCEPTION(Usage, Exception);

namespace detail {

template <typename T>
static const typename T::mapped_type& get(const T& m, const typename T::key_type& k) {
  typename T::const_iterator i = m.find(k);
  VW_ASSERT(i != m.end(), Usage() << "Unknown key: " << k);
  return i->second;
}
} // namespace detail

template <typename T>
class Tristate {
    T m_value;
    bool m_req;
    bool m_set;
  public:
    bool set() const {return m_set;}
    Tristate(const T& val = T(), bool required = false, bool is_set = false)
      : m_value(val), m_req(required), m_set(is_set) {}
    friend std::istream& operator>>(std::istream& in, Tristate& val) {
      T v;
      in >> v;
      val = v;
      return in;
    }
    Tristate& operator= (const T& value) {
      m_set = true;
      m_value = value;
      return *this;
    }
    operator const T&() const {
      VW_ASSERT(m_set || !m_req, LogicErr() << "Tried to access an unset required value");
      return m_value;
    }
    const T& value() const {
      VW_ASSERT(m_set || !m_req, LogicErr() << "Tried to access an unset required value");
      return m_value;
    }
};

}} // namespace vw::tools

#define VW_DEFINE_ENUM(name, len, tuple) \
class name {\
  public:\
    enum Value {BOOST_PP_TUPLE_REM(len)tuple};\
    static const Value elems[len];\
    static size_t count() {return len;}\
  private:\
    typedef std::map<Value, std::string> s_t;\
    typedef std::map<std::string, Value> v_t;\
    static const s_t m_sm; \
    static const v_t m_vm; \
    Value m_val; \
  public:\
    name() {} \
    name(Value v) : m_val(v) {} \
    name(const std::string& v) : m_val(vw::tools::detail::get(m_vm, boost::to_upper_copy(v))) {}\
    const std::string& string() const {\
      return vw::tools::detail::get(m_sm, m_val);\
    }\
    operator Value() const {return value();} \
    Value value() const {\
      return m_val; \
    }\
    friend std::istream& operator>>(std::istream& in, name& val) {\
      std::string str;\
      in >> str;\
      val = name(str); \
      return in;\
    }\
    friend std::ostream& operator<<(std::ostream& out, const name& val) {\
      out << val.string();\
      return out;\
    }\
    static std::string list(const std::string& delim = ", ") { \
      std::string ret;\
      BOOST_FOREACH(const Value& v, elems) {\
        if (v != elems[0])\
          ret += delim;\
        ret += name(v).string();\
      }\
      return ret;\
    } \
};\
const name::Value name::elems[len] = {BOOST_PP_TUPLE_REM(len)tuple}; \
const name::s_t name::m_sm = boost::assign::map_list_of BOOST_PP_SEQ_TRANSFORM(VW_PP_ELEM_TO_ELEM_STRING, 0, BOOST_PP_TUPLE_TO_SEQ(len, tuple));\
const name::v_t name::m_vm = boost::assign::map_list_of BOOST_PP_SEQ_TRANSFORM(VW_PP_ELEM_TO_STRING_ELEM, 0, BOOST_PP_TUPLE_TO_SEQ(len, tuple));


#endif
