// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Plate/FundamentalTypes.h>
#include <vw/Plate/IndexData.pb.h>
#include <google/protobuf/service.h>

namespace vw {
namespace platefile {
  namespace detail {
    RequireCall::RequireCall(Func* call) : call(call) {}
    RequireCall::~RequireCall() {call->Run();}
  }

TransactionOrNeg::TransactionOrNeg() VW_NOTHROW {
  set(detail::NO_TRANSACTION);
}
TransactionOrNeg::TransactionOrNeg(int32 id) VW_NOTHROW {
  set(id);
}

Transaction TransactionOrNeg::promote() const {
  // catch the error here, rather than in Transaction's constructor
  VW_ASSERT(!invalid(), LogicErr() << "Cannot promote a negative transaction id!");
  return Transaction(m_id);
}

void TransactionOrNeg::set(int32 id) VW_NOTHROW {
  // map negative ids to NO_TRANSACTION
  uint32 idp = static_cast<uint32>(id);
  if (idp > detail::MAX_TRANSACTION) m_id = detail::NO_TRANSACTION;
  else m_id = idp;
}

bool TransactionOrNeg::newest() const VW_NOTHROW {
  return invalid();
}

#define SAME_OP(Type) \
    bool Type::operator< (const TransactionOrNeg& x) const { return m_id <  x.m_id; } \
    bool Type::operator< (const uint32& x) const           { return m_id <  x; }      \
    bool Type::operator==(const uint32& x) const           { return m_id == x; }      \
    bool Type::operator< (const Transaction& x) const      { return m_id <  x.m_id; } \
    bool Type::operator==(const Transaction& x) const      { return m_id == x.m_id; }

bool Transaction::operator>     ( const TransactionOrNeg& x) const {return m_id > x.m_id; }
bool Transaction::operator>     ( const uint32& x) const           {return m_id > x; }
bool TransactionOrNeg::operator>( const Transaction& x) const      {return m_id > x.m_id; }
bool TransactionOrNeg::operator>( const uint32& x) const           {return m_id > x; }
bool TransactionOrNeg::operator==(const TransactionOrNeg& x) const { return m_id == x.m_id; }

SAME_OP(Transaction);
SAME_OP(TransactionOrNeg);

Transaction::Transaction(uint32 id) {
  set(id);
}

void Transaction::set(uint32 id) {
  VW_ASSERT(id <= detail::MAX_TRANSACTION, ArgumentErr() << "Transaction ID greater than max: " << id);
  m_id = id;
}

Transaction::operator uint32() const VW_NOTHROW {return m_id;}
Transaction::operator TransactionOrNeg() const VW_NOTHROW {return TransactionOrNeg(static_cast<int32>(m_id));}

TransactionRange::TransactionRange(TransactionOrNeg first, TransactionOrNeg last)
  : Super(first, last)
{
  VW_ASSERT(first <= last, ArgumentErr() << "TransactionRange: first must be <= last");
}

TransactionRange::TransactionRange(TransactionOrNeg only)
  : Super(only, only)
{}

TransactionOrNeg TransactionRange::first() const {
  return Super::first;
}
TransactionOrNeg TransactionRange::last() const {
  if (second.newest())
    return second;
  else
    return second.promote();
}

std::ostream& operator<<(std::ostream& o, const vw::platefile::TransactionRange& range) {
  return (o << "range(" << range.first() << "," << range.last() << ")");
}

std::ostream& operator<<(std::ostream& o, const vw::platefile::TileHeader& hdr) {
  return (o << hdr.col() << "," << hdr.row() << "@" << hdr.level() << " (t_id = " << hdr.transaction_id() << ")");
}

bool operator==(const vw::platefile::TileHeader& a, const vw::platefile::TileHeader& b) {
  if (a.col() != b.col())                       return false;
  if (a.row() != b.row())                       return false;
  if (a.level() != b.level())                   return false;
  if (a.transaction_id() != b.transaction_id()) return false;
  if (a.filetype() != b.filetype())             return false;
  return true;
}

bool operator!=(const vw::platefile::TileHeader& a, const vw::platefile::TileHeader& b) {
  return !(a == b);
}

bool OrderHeaderByTidDesc::operator()(const vw::platefile::TileHeader& a, const vw::platefile::TileHeader& b) const {
  return Transaction(a.transaction_id()) > Transaction(b.transaction_id());
}

}} // vw::platefile
