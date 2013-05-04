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


#include <vw/Plate/FundamentalTypes.h>
#include <vw/Plate/IndexData.pb.h>
#include <vw/Plate/Datastore.h>
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

uint32 Transaction::MAX_POSSIBLE() {
  return detail::MAX_TRANSACTION;
}

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
  return Super::second;
}

std::ostream& operator<<(std::ostream& o, const vw::platefile::TransactionRange& range) {
  return (o << "range(" << range.first() << "," << range.last() << ")");
}

std::ostream& operator<<(std::ostream& o, const vw::platefile::TileHeader& hdr) {
  return (o << hdr.col() << "," << hdr.row() << "@" << hdr.level() << " (t_id = " << hdr.transaction_id() << ", type=" << hdr.filetype() << ")");
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

bool OrderTileByTidDesc::operator()(const vw::platefile::Tile& a, const vw::platefile::Tile& b) const {
  return Transaction(a.hdr.transaction_id()) > Transaction(b.hdr.transaction_id());
}

}} // vw::platefile
