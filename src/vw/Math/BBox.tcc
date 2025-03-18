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


/// \file BBox.tcc
///
/// Provides a generic bounding box.

namespace vw {
namespace math {


// *****************************************************************************************
// class BBox
// *****************************************************************************************

// A general n-dimensional axis-aligned bounding box class,
// represented by vectors pointing to the minimal and maximal corners.

//-----------------------
// Constructors 

template <class RealT, size_t DimN>
BBox<RealT, DimN>::BBox() {
  // Make sure we have a type for which we know limits
  BOOST_STATIC_ASSERT(std::numeric_limits<RealT>::is_specialized);
  if (std::numeric_limits<RealT>::is_integer) {
    for (size_t i = 0; i < m_min.size(); ++i) {
      m_min[i] = std::numeric_limits<RealT>::max();
      m_max[i] = std::numeric_limits<RealT>::min();
    }
  }
  else { // Not an integer
    for (size_t i = 0; i < m_min.size(); ++i) {
      m_min[i] = std::numeric_limits<RealT>::max();
      m_max[i] = static_cast<RealT>(-std::numeric_limits<RealT>::max());
    }
  }
}

template <class RealT, size_t DimN>
template <class RealT1, size_t DimN1>
BBox<RealT, DimN>::BBox( BBox<RealT1, DimN1> const& bbox )
  : m_min(bbox.min()), m_max(bbox.max()) {}



template <class RealT, size_t DimN>
BBox<RealT, DimN>::BBox( RealT minx, RealT miny, RealT width, RealT height )
  : m_min(Vector<RealT,2>(minx,      miny       )),
    m_max(Vector<RealT,2>(minx+width,miny+height)) {
  checkLengthEqual<2>();
}

template <class RealT, size_t DimN>
BBox<RealT, DimN>::BBox( RealT minx, RealT miny, RealT minz, RealT width, RealT height, RealT depth )
  : m_min(Vector<RealT,3>(minx,      miny,       minz      )),
    m_max(Vector<RealT,3>(minx+width,miny+height,minz+depth)) {
  checkLengthEqual<3>();
}

template <class RealT, size_t DimN>
template <class RealT1, size_t DimN1>
BBox<RealT, DimN>& BBox<RealT, DimN>::operator=( BBox<RealT1, DimN1> const& bbox ) {
  this->min() = bbox.min();
  this->max() = bbox.max();
  return *this;
}


// End Constructors 
//-----------------------

template <class RealT, size_t DimN>
template <class VectorT>
void BBox<RealT, DimN>::grow(VectorBase<VectorT> const& point) {
  VW_ASSERT(m_min.size() == 0 || point.impl().size() == m_min.size(), 
            ArgumentErr() << "Vector must have dimension " << m_min.size() << ".");
  if (m_min.size() == 0) {
    m_min = point;
    m_max = point;
  }
  else {
    for (size_t i = 0; i < m_min.size(); ++i) {
      if (point.impl()[i] > m_max[i])
        m_max[i] = RealT(point.impl()[i]);
      if (point.impl()[i] < m_min[i])
        m_min[i] = RealT(point.impl()[i]);
    }
  }
}

/// Grows a bounding box to include the given bounding box.
template <class RealT, size_t DimN>
template <class RealT1, size_t DimN1>
void BBox<RealT, DimN>::grow( BBox<RealT1, DimN1> const& bbox ){
  if (bbox.empty())
    return;
  grow(bbox.min()); grow(bbox.max());
}

template <class RealT, size_t DimN>
template <class RealT1, size_t DimN1>
void BBox<RealT, DimN>::crop( BBox<RealT1, DimN1> const& bbox ) {
  VW_ASSERT(m_min.size() == 0 || bbox.min().size() == m_min.size(), 
            ArgumentErr() << "BBox must have dimension " << m_min.size() << ".");
  for( size_t i=0; i<m_min.size(); ++i ) {
    if( m_min[i] < bbox.min()[i] ) 
      m_min[i] = bbox.min()[i];

    if( m_max[i] > bbox.max()[i] ) 
      m_max[i] = bbox.max()[i];
  }
}

template <class RealT, size_t DimN>
template <class VectorT>
bool BBox<RealT, DimN>::contains( const VectorBase<VectorT> &point ) const {
  using namespace vector_containment_comparison;
  VW_ASSERT(m_min.size() == 0 || point.impl().size() == m_min.size(), 
            ArgumentErr() << "Vector must have dimension " << m_min.size() << ".");
  return ((m_min.size() != 0) && (point >= m_min) && (point < m_max));
}

template <class RealT, size_t DimN>
template <class RealT1, size_t DimN1>
bool BBox<RealT, DimN>::contains( const BBox<RealT1, DimN1> &bbox ) const {
  using namespace vector_containment_comparison;
  VW_ASSERT(m_min.size() == 0 || bbox.min().size() == m_min.size(), 
            ArgumentErr() << "BBox must have dimension " << m_min.size() << ".");
  return ((m_min.size() != 0) && (bbox.min() >= m_min) && (bbox.max() <= m_max));
}

template <class RealT, size_t DimN>
template <class RealT1, size_t DimN1>
bool BBox<RealT, DimN>::intersects( const BBox<RealT1, DimN1>& bbox ) const {
  VW_ASSERT(m_min.size() == 0 || bbox.min().size() == m_min.size(), 
            ArgumentErr() << "BBox must have dimension " << m_min.size() << ".");
  for( size_t i=0; i<m_min.size(); ++i ) {
    if( m_min[i] >= bbox.max()[i] ||
        m_max[i] <= bbox.min()[i] )
      return false;
  }
  return (m_min.size() != 0);
}

// TODO(oalexan1): In many places below use fully_empty() rather than empty()
  
template <class RealT, size_t DimN>
bool BBox<RealT, DimN>::empty() const {
  for( size_t i=0; i<m_min.size(); ++i )
    if (m_min[i] >= m_max[i]) return true;
  return (m_min.size() <= 0);
}

template <class RealT, size_t DimN>
bool BBox<RealT, DimN>::fully_empty() const {
  for( size_t i=0; i<m_min.size(); ++i )
    if (m_min[i] > m_max[i]) return true;
  return (m_min.size() <= 0);
}

template <class RealT, size_t DimN>
RealT BBox<RealT, DimN>::width() const {
  checkLengthGte<1>();
  if (empty()) return 0.0; // Bug fix for underflow/overflow
  return this->max()[0] - this->min()[0];
}

template <class RealT, size_t DimN>
RealT BBox<RealT, DimN>::height() const {
  checkLengthGte<2>();
  if (empty()) return 0.0; // Bug fix for underflow/overflow
  return this->max()[1] - this->min()[1];
}

template <class RealT, size_t DimN>
RealT BBox<RealT, DimN>::depth() const {
  checkLengthGte<3>();
  if (empty()) return 0.0; // Bug fix for underflow/overflow
  return this->max()[2] - this->min()[2];
}


template <class RealT, size_t DimN>
RealT BBox<RealT, DimN>::area() const {
  checkLengthGte<2>();
  if (empty()) return 0.0; // Bug fix for underflow/overflow
  return (this->max()[0] - this->min()[0]) *
         (this->max()[1] - this->min()[1]);
}

template <class RealT, size_t DimN>
RealT BBox<RealT, DimN>::volume() const {
  checkLengthGte<3>();
  if (empty()) return 0.0; // Bug fix for underflow/overflow
  return (this->max()[0] - this->min()[0]) *
         (this->max()[1] - this->min()[1]) *
         (this->max()[2] - this->min()[2]);
}

template <class RealT, size_t DimN>
Vector<RealT, DimN> BBox<RealT, DimN>::size() const {
  // To do: This bugfix makes some tests fail, need to investigate.
  // Turning it off for now.
  //if (empty()) return Vector<RealT, DimN>(); // Bug fix for underflow/overflow
  return (m_max - m_min);
}

template <class RealT, size_t DimN>
Vector<RealT, DimN> BBox<RealT, DimN>::center() const {
  if (empty()) return Vector<RealT, DimN>(); // Bug fix for underflow/overflow
  return 0.5 * (m_min + m_max);
}

template <class RealT, size_t DimN>
void BBox<RealT, DimN>::set_size(Vector<RealT, DimN> size) {
  m_max = m_min + size;
}

template <class RealT, size_t DimN>
void BBox<RealT, DimN>::expand( RealT offset ) {
  if (empty()) return; // Bug fix for underflow/overflow
  for( size_t i=0; i<m_min.size(); ++i ) {
    m_min[i] -= offset;
    m_max[i] += offset;
  }
}

template <class RealT, size_t DimN>
template <class VectorT>
void BBox<RealT, DimN>::expand( VectorBase<VectorT> const& vec ) {
  VW_ASSERT(m_min.size() == 0 || vec.impl().size() == m_min.size(), 
            ArgumentErr() << "Vector must have dimension " << m_min.size() << ".");
  if (empty()) return; // Bug fix for underflow/overflow
  for( size_t i=0; i<m_min.size(); ++i ) {
    m_min[i] -= vec.impl()[i];
    m_max[i] += vec.impl()[i];
  }
}

template <class RealT, size_t DimN>
void BBox<RealT, DimN>::contract( RealT offset ) {
  if (empty()) return; // Bug fix for overflow/underflow
  for( size_t i=0; i<m_min.size(); ++i ) {
    m_min[i] += offset;
    m_max[i] -= offset;
  }
}

template <class RealT, size_t DimN>
template <class ScalarT>
BBox<RealT, DimN>& BBox<RealT, DimN>::operator*=( ScalarT s ) {
  VW_ASSERT( s > 0, ArgumentErr()
             << "Cannot scale a box by a non-positive number." );
  if (empty()) return *this; // Bug fix for underflow/overflow
  m_min *= s;
  m_max *= s;
  return *this;
}

template <class RealT, size_t DimN>
template <class ScalarT>
BBox<RealT, DimN>& BBox<RealT, DimN>::operator/=( ScalarT s ) {
  VW_ASSERT( s > 0, ArgumentErr()
             << "Cannot scale a box by a non-positive number." );
  if (empty()) return *this; // Bug fix for underflow/overflow
  m_min /= s;
  m_max /= s;
  return *this;
}

template <class RealT, size_t DimN>
template <class VectorT>
BBox<RealT, DimN>& BBox<RealT, DimN>::operator+=( VectorBase<VectorT> const& v ) {
  if (empty()) return *this; // Bug fix for underflow/overflow
  m_min += v;
  m_max += v;
  return *this;
}

template <class RealT, size_t DimN>
template <class VectorT>
BBox<RealT, DimN>& BBox<RealT, DimN>::operator-=( VectorBase<VectorT> const& v ) {
  if (empty()) return *this; // Bug fix for underflow/overflow
  m_min -= v;
  m_max -= v;
  return *this;
}


// Overloaded functions for operating on BBox
// - Some of these operations should probably have been named class functions.

template <class RealT, size_t DimN, class ScalarT>
inline BBox<RealT, DimN> operator*( BBox<RealT, DimN> const& bbox, ScalarT s ) {
  BBox<RealT, DimN> result = bbox;
  result *= s;
  return result;
}

template <class RealT, size_t DimN, class ScalarT>
inline BBox<RealT, DimN> operator/( BBox<RealT, DimN> const& bbox, ScalarT s ) {
  BBox<RealT, DimN> result = bbox;
  result /= s;
  return result;
}

template <class RealT, size_t DimN, class VectorT>
inline BBox<RealT, DimN> operator+( BBox<RealT, DimN> const& bbox, VectorBase<VectorT> const& v ) {
  BBox<RealT, DimN> result = bbox;
  result += v.impl();
  return result;
}

template <class RealT, size_t DimN, class VectorT>
inline BBox<RealT, DimN> operator-( BBox<RealT, DimN> const& bbox, VectorBase<VectorT> const& v ) {
  BBox<RealT, DimN> result = bbox;
  result -= v.impl();
  return result;
}

template <class RealT1, size_t DimN1, class RealT2, size_t DimN2>
inline bool operator==( BBox<RealT1,DimN1> const& bbox1, BBox<RealT2,DimN2> const& bbox2 ) {
  return bbox1.min()==bbox2.min() && bbox1.max()==bbox2.max();
}

template <class RealT1, size_t DimN1, class RealT2, size_t DimN2>
inline bool operator!=( BBox<RealT1,DimN1> const& bbox1, BBox<RealT2,DimN2> const& bbox2 ) {
  return bbox1.min()!=bbox2.min() || bbox1.max()!=bbox2.max();
}


template <class RealT>
std::ostream& operator<<( std::ostream& os, BBox<RealT,2> const& bbox ) {
  return os << "Min: (" << bbox.min()[0] << ", " << bbox.min()[1] << ") width: " << bbox.width() << " height: " << bbox.height();
}

template <class RealT, size_t DimN, class VectorT>
inline BBox<RealT, DimN> elem_prod( BBox<RealT, DimN> const& bbox, VectorBase<VectorT> const& v ) {
  if (bbox.empty()) 
    return BBox<RealT, DimN>(); // Bug fix for underflow/overflow
  return BBox<RealT, DimN>( elem_prod(bbox.min(),v), elem_prod(bbox.max(),v) );
}

template <class RealT, size_t DimN, class VectorT>
inline BBox<RealT, DimN> elem_quot( BBox<RealT, DimN> const& bbox, VectorBase<VectorT> const& v ) {
  if (bbox.empty()) 
    return BBox<RealT, DimN>(); // Bug fix for underflow/overflow
  return BBox<RealT, DimN>( elem_quot(bbox.min(),v), elem_quot(bbox.max(),v) );
}



} // namespace math


// A helper function to grow a floating-point bounding box
// to the smallest enclosing integer bounding box.
template <class RealT1, size_t DimN1>
inline typename boost::enable_if<boost::is_float<RealT1>,BBox<int32,DimN1> >::type
grow_bbox_to_int( math::BBox<RealT1, DimN1> const& bbox ) {
  BBox<int32,DimN1> result;
  for ( size_t i = 0; i < DimN1; i++ ) {

    // Bug fix: Cannot grow a completely empty box.
    // Note: A box with bbox.min()[i] == bbox.max()[i]
    // is in principle empty, however such boxes are created
    // by VW when growing an empty box with a point.
    // VW is messed up here, but we must use below ">"
    // rather than ">=" to not break existing behavior.
    if (bbox.min()[i] > bbox.max()[i]){
      return  BBox<int32,DimN1>();
    }

    // Note: This is not quite correct, as we grow the box
    // even when it is already integer, but exiting
    // behavior again depends on this.
    result.min()[i] = (int32)floor(bbox.min()[i]);
    result.max()[i] = (int32)floor(bbox.max()[i])+1;
  }
  return result;
}

} // namespace vw

