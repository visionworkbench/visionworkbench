// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_PLATE_FUNDAMENTALTYPES_H__
#define __VW_PLATE_FUNDAMENTALTYPES_H__

#include <vw/Core/FundamentalTypes.h>
#include <vw/Core/Exception.h>
#include <boost/operators.hpp>

namespace google { namespace protobuf {
  class Closure;
}}

namespace vw {
namespace platefile {

namespace detail {

  struct RequireCall {
    typedef ::google::protobuf::Closure Func;
    Func *call;
    RequireCall(Func* call);
    ~RequireCall();
  };

  // These numbers carefully chosen to match reality and the proper ordering...
  // "No Transaction" should be -1 and should always compare greater than
  // anything in the valid transaction range. Every invalid transaction should
  // compare equal (which is to say, subclasses must make sure every invalid
  // transaction id is mapped to NO_TRANSACTION)
  static const uint32 NO_TRANSACTION  = static_cast<uint32>(int32(-1));
  static const uint32 MAX_TRANSACTION = 2147483647u; // std::numeric_limits<int32>::max();
  BOOST_STATIC_ASSERT(static_cast<int32>(NO_TRANSACTION) == int32(-1));
  BOOST_STATIC_ASSERT(MAX_TRANSACTION < NO_TRANSACTION);

  template <typename ImplT>
  //class TransactionBase : boost::totally_ordered<ImplT> {
  class TransactionBase {
    private:
      inline ImplT& impl() { return static_cast<ImplT&>(*this); }
      inline ImplT const& impl() const { return static_cast<ImplT const&>(*this); }
    protected:
      // Maintain contract:
      // Must be [0..MAX_TRANSACTION] or NO_TRANSACTION
      uint32 m_id;
      bool invalid() const {return m_id == NO_TRANSACTION;}
    public:
      friend std::istream& operator>>(std::istream& i, TransactionBase& val) {
        uint32 id;
        i >> id;
        val.impl().set(id);
        return i;
      }
      friend std::ostream& operator<<(std::ostream& o, const TransactionBase& val) {
        if (val.m_id == NO_TRANSACTION)
          o << "inf";
        else
          o << val.m_id;
        return o;
      }
  };
}

class Transaction;

class TransactionOrNeg : public detail::TransactionBase<TransactionOrNeg>, boost::totally_ordered<TransactionOrNeg> {
  public:
    TransactionOrNeg() VW_NOTHROW;
    TransactionOrNeg(int32 id) VW_NOTHROW;
    // promote to a concrete Transaction. Will throw if transaction is negative.
    Transaction promote() const;
    // Set the id. All negative ids are considered a request for the newest transaction
    void set(int32 id) VW_NOTHROW;
    // Is this a request for the newest transaction id?
    bool newest() const VW_NOTHROW;
    bool operator<(const TransactionOrNeg& x) const {
      return m_id < x.m_id;
    }
    bool operator==(const TransactionOrNeg& x) const {
      return m_id == x.m_id;
    }
};

class Transaction : public detail::TransactionBase<Transaction> {
  public:
    Transaction(uint32 id);
    // Set the id. Will throw if transaction id is -1.
    void set(uint32 id);
    operator uint32() const VW_NOTHROW;
    operator TransactionOrNeg() const VW_NOTHROW;
};

}} // vw::platefile

#endif
