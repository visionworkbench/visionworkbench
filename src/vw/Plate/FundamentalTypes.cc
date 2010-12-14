#include <vw/Plate/FundamentalTypes.h>
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

Transaction::Transaction(uint32 id) {
  set(id);
}

void Transaction::set(uint32 id) {
  VW_ASSERT(id <= detail::MAX_TRANSACTION, ArgumentErr() << "Transaction ID greater than max: " << id);
  m_id = id;
}

Transaction::operator uint32() const VW_NOTHROW {return m_id;}
Transaction::operator TransactionOrNeg() const VW_NOTHROW {return TransactionOrNeg(static_cast<int32>(m_id));}

}} // vw::platefile
