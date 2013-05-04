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


#ifndef __VW_PLATE_DETAIL_MIPMAPHELPERS_H__
#define __VW_PLATE_DETAIL_MIPMAPHELPERS_H__

#include <boost/assign/list_of.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <vw/Plate/IndexData.pb.h>
#include <vw/Math/BBox.h>
#include <vw/Core/ProgressCallback.h>

namespace vw { namespace platefile { namespace detail {

using namespace vw;

typedef boost::tuple<uint32, uint32, uint32> rowcoltid_t;
typedef boost::tuple<uint32, uint32>         rowcol_t;
typedef boost::tuple<uint32, rowcoltid_t>    tile_order_t;

inline uint32 therow(const platefile::TileHeader& h) { return h.row(); }
inline uint32 thecol(const platefile::TileHeader& h) { return h.col(); }
inline uint32 thetid(const platefile::TileHeader& h) { return h.transaction_id(); }
inline uint32 therow(const rowcoltid_t& h) { return h.get<0>(); }
inline uint32 thecol(const rowcoltid_t& h) { return h.get<1>(); }
inline uint32 thetid(const rowcoltid_t& h) { return h.get<2>(); }
inline uint32 therow(const rowcol_t& h) { return h.get<0>(); }
inline uint32 thecol(const rowcol_t& h) { return h.get<1>(); }

inline BBox2i move_down(const BBox2i& input, uint32 level_change) {
  return input * (1 << level_change);
}

inline rowcol_t parent_tile(const uint32 row, const uint32 col) {
  return rowcol_t(row/2, col/2);
}

// given a parent tile (1 level up) and a hdr, calculate the composite id;
// [0==UL, 1==UR, 2==LL, 3==LR]
template <typename T1, typename T2>
uint32 calc_composite_id(const T1& parent, const T2& hdr) {
  typedef std::map<rowcol_t, uint32> map_t;
  static const map_t lookup = boost::assign::map_list_of
    (rowcol_t(0,0), 0)
    (rowcol_t(0,1), 1)
    (rowcol_t(1,0), 2)
    (rowcol_t(1,1), 3);

  rowcol_t offset(therow(hdr) - therow(parent) * 2,
                  thecol(hdr) - thecol(parent) * 2);

  map_t::const_iterator i = lookup.find(offset);
  VW_ASSERT(i != lookup.end(), LogicErr() << "Cannot determine composite id for hdr " << thecol(hdr) << "," << therow(hdr) << " and parent "
      << thecol(parent) << "," << therow(parent));
  return i->second;
}

struct SortByTidDesc {
  bool operator()(const tile_order_t& a, const tile_order_t& b)
  {
    return b.get<0>() < a.get<0>();
  }
  template <typename T>
  bool operator()(const T& a, const T& b) {
    return thetid(b) < thetid(a);
  }
};

// Ascending order
struct SortByTid {
  bool operator()(const tile_order_t& a, const tile_order_t& b) {
    return a.get<0>() < b.get<0>();
  }
  template <typename T>
  bool operator()(const T& a, const T& b) {
    return thetid(a) < thetid(b);
  }
};

class RememberCallback : public SubProgressCallback {
    mutable double m_count;
    double m_total;
  public:
    RememberCallback(const ProgressCallback &parent, double percent, double total)
      : SubProgressCallback(parent, parent.progress(), parent.progress() + percent), m_count(0), m_total(total) {}
    void tick(uint32 count = 1) const {
      m_count += count;
      if (m_count > m_total) m_count = m_total;
      this->report_fractional_progress(m_count, m_total);
    }
    ~RememberCallback() {
      // Calling a virtual function from a destructor is bad :/
      if (m_count < m_total)
        SubProgressCallback::report_progress(1);
    }
};
}

inline std::ostream& operator<<(std::ostream& o, const detail::rowcoltid_t& hdr) {
  return (o << detail::thecol(hdr) << "," << detail::therow(hdr) << " (t_id = " << detail::thetid(hdr) << ")");
}

}}

#endif
