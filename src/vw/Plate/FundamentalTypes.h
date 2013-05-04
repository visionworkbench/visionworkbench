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


#ifndef __VW_PLATE_FUNDAMENTALTYPES_H__
#define __VW_PLATE_FUNDAMENTALTYPES_H__

#include <vw/Core/FundamentalTypes.h>
#include <vw/Core/Exception.h>
#include <boost/operators.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/shared_container_iterator.hpp>
#include <boost/range/iterator_range.hpp>

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

class TransactionOrNeg : public detail::TransactionBase<TransactionOrNeg>,
                         boost::totally_ordered<TransactionOrNeg,
                         boost::totally_ordered<TransactionOrNeg, uint32,
                         boost::totally_ordered<TransactionOrNeg, Transaction> > >
{
  protected:
    friend class Transaction;
  public:
    TransactionOrNeg() VW_NOTHROW;
    TransactionOrNeg(int32 id) VW_NOTHROW;
    // promote to a concrete Transaction. Will throw if transaction is negative.
    Transaction promote() const;
    // Set the id. All negative ids are considered a request for the newest transaction
    void set(int32 id) VW_NOTHROW;
    // Is this a request for the newest transaction id?
    bool newest() const VW_NOTHROW;

    bool operator<(const TransactionOrNeg& x) const;
    bool operator==(const TransactionOrNeg& x) const;

    bool operator<(const uint32& x) const;
    bool operator>(const uint32& x) const;
    bool operator==(const uint32& x) const;

    bool operator<(const Transaction& x) const;
    bool operator>(const Transaction& x) const;
    bool operator==(const Transaction& x) const;
};

class Transaction : public detail::TransactionBase<Transaction>,
                    boost::totally_ordered<Transaction,
                    boost::totally_ordered<Transaction, uint32> >
{
  protected:
    friend class TransactionOrNeg;
  public:
    Transaction(uint32 id);
    // Set the id. Will throw if transaction id is -1.
    void set(uint32 id);
    operator uint32() const VW_NOTHROW;
    operator TransactionOrNeg() const VW_NOTHROW;

    bool operator<(const TransactionOrNeg& x) const;
    bool operator>(const TransactionOrNeg& x) const;

    bool operator<(const uint32& x) const;
    bool operator>(const uint32& x) const;
    bool operator==(const uint32& x) const;

    bool operator<(const Transaction& x) const;
    bool operator==(const Transaction& x) const;

    static uint32 MAX_POSSIBLE();
};

class TransactionRange : private std::pair<TransactionOrNeg, TransactionOrNeg> {
  typedef std::pair<TransactionOrNeg, TransactionOrNeg> Super;
  public:
    // INCLUSIVE RANGE
    TransactionRange(TransactionOrNeg first, TransactionOrNeg last);
    TransactionRange(TransactionOrNeg only);
    TransactionOrNeg first() const;
    TransactionOrNeg last() const;
};

std::ostream& operator<<(std::ostream& o, const vw::platefile::TransactionRange& range);

// We can't edit the protobuf-generated code, so this is next best place for this
class TileHeader;
class Tile;
std::ostream& operator<<(std::ostream& o, const vw::platefile::TileHeader& hdr);
bool operator==(const vw::platefile::TileHeader& a, const vw::platefile::TileHeader& b);
bool operator!=(const vw::platefile::TileHeader& a, const vw::platefile::TileHeader& b);

struct OrderHeaderByTidDesc {
  bool operator()(const vw::platefile::TileHeader& a, const vw::platefile::TileHeader& b) const;
};
struct OrderTileByTidDesc {
  bool operator()(const vw::platefile::Tile& a, const vw::platefile::Tile& b) const;
};

// The boost make_shared_iterator_range returns a pair rather than an iterator_range.
// Provide a helper that does the right thing.
template <typename ContainerT>
boost::iterator_range<boost::shared_container_iterator<const ContainerT> >
make_const_shared_range(boost::shared_ptr<ContainerT> c) {
  typedef boost::shared_container_iterator<const ContainerT> iter_t;
  iter_t begin(c->begin(), c);
  iter_t     end(c->end(), c);
  return boost::make_iterator_range(begin, end);
}

}} // vw::platefile

#endif
