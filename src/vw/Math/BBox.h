// __BEGIN_LICENSE__
// 
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
// 
// Copyright 2006 Carnegie Mellon University. All rights reserved.
// 
// This software is distributed under the NASA Open Source Agreement
// (NOSA), version 1.3.  The NOSA has been approved by the Open Source
// Initiative.  See the file COPYING at the top of the distribution
// directory tree for the complete NOSA document.
// 
// THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
// KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
// LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
// SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
// A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
// THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
// DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
// 
// __END_LICENSE__

/// \file BBox.h
///
/// Provides a generic bounding-box type based on vectors pointing 
/// to the minimal and maximal (e.g. upper-left and lower-right) corners.
#ifndef __VW_MATH__BBOX_H__
#define __VW_MATH__BBOX_H__

#include <limits>

#include <boost/static_assert.hpp>

#include <vw/Math/Vector.h>

namespace vw {
namespace math {

  /// \cond INTERNAL
  namespace vector_containment_comparison {
    template <class ElemT, int SizeN>
    inline bool operator<( Vector<ElemT,SizeN> const& v1, Vector<ElemT,SizeN> const& v2 ) {
      for( int i=0; i<SizeN; ++i)
        if( ! ( v1[i] < v2[i] ) ) return false;
      return true;
    }

    template <class ElemT, int SizeN>
    inline bool operator<=( Vector<ElemT,SizeN> const& v1, Vector<ElemT,SizeN> const& v2 ) {
      for( int i=0; i<SizeN; ++i)
        if( ! ( v1[i] <= v2[i] ) ) return false;
      return true;
    }

    template <class ElemT, int SizeN>
    inline bool operator>( Vector<ElemT,SizeN> const& v1, Vector<ElemT,SizeN> const& v2 ) {
      for( int i=0; i<SizeN; ++i)
        if( ! ( v1[i] > v2[i] ) ) return false;
      return true;
    }

    template <class ElemT, int SizeN>
    inline bool operator>=( Vector<ElemT,SizeN> const& v1, Vector<ElemT,SizeN> const& v2 ) {
      for( int i=0; i<SizeN; ++i)
        if( ! ( v1[i] >= v2[i] ) ) return false;
      return true;
    }

  } // namespace vector_containment_comparison
  /// \endcond


  // *******************************************************************
  // class BBox
  // *******************************************************************

  /// A general n-dimensional axis-aligned bounding box class,
  /// represented by vectors pointing to the minimal and maximal
  /// corners.
  template <class RealT, int DimN>
  class BBox {
  public:

    /// Default constructor.  Constructs the ultimate empty bounding 
    /// box, whose limits are at the opposite corners of the underlying 
    /// numeric space.  This is a useful starting point if you intend 
    /// to grow your bounding box to fit a collection of items.
    BBox() {
      // Make sure we have a type for which we know limits
      BOOST_STATIC_ASSERT(std::numeric_limits<RealT>::is_specialized);
      if (std::numeric_limits<RealT>::is_integer) {
        for (int i = 0; i < DimN; i++) {
          m_min[i] = std::numeric_limits<RealT>::max();
          m_max[i] = std::numeric_limits<RealT>::min();
        }
      }
      else {
        for (int i = 0; i < DimN; i++) {
          m_min[i] = std::numeric_limits<RealT>::max();
          m_max[i] = -std::numeric_limits<RealT>::max();
        }
      }
    }

    /// Constructs a bounding box with the given minimal and maximal 
    /// points.
    BBox( Vector<RealT, DimN> const& min, Vector<RealT, DimN> const& max ) :
      m_min( min ), m_max( max ) {}

    /// Constructs a 2D bounding box with the given minimal point
    /// coordinates and dimensions.  (Only valid for 2D bouding
    /// boxes.)
    BBox( RealT minx, RealT miny, RealT width, RealT height )
      : m_min( Vector<RealT,2>(minx,miny) ), 
        m_max( Vector<RealT,2>(minx+width,miny+height) )
    {
      BOOST_STATIC_ASSERT( DimN==2 );
    }

    /// Grows a bounding box to include the given point.
    void grow( Vector<RealT, DimN> const& point ) {
      for (int i = 0; i < DimN; i++) {
	if (point[i] > m_max[i])
	  m_max[i] = point[i];
	if (point[i] < m_min[i])
	  m_min[i] = point[i];
      }
    }
    
    /// Grows a bounding box to include the given bounding box.
    void grow( BBox const& bbox ) {
      grow(bbox.min());
      grow(bbox.max());
    }

    /// Crops (intersects) this bounding box to the given bounding box.
    void crop( BBox const& bbox ) {
      for( int i=0; i<DimN; ++i ) {
        if( m_min[i] < bbox.m_min[i] )
          m_min[i] = bbox.m_min[i];
        if( m_max[i] > bbox.m_max[i] )
          m_max[i] = bbox.m_max[i];
      }
    }

    /// Expands this bounding box by the given offset in every direction.
    void expand( RealT offset ) {
      for( int i=0; i<DimN; ++i ) {
        m_min[i] -= offset;
        m_max[i] += offset;
      }
    }

    /// Contracts this bounding box by the given offset in every direction.
    void contract( RealT offset ) {
      for( int i=0; i<DimN; ++i ) {
        m_min[i] += offset;
        m_max[i] -= offset;
      }
    }

    /// Returns true if the given point is contained in the bounding box.
    bool contains( const Vector<RealT, DimN> &point ) const {
      using namespace vector_containment_comparison;
      return ((point >= m_min) && (point < m_max));
    }

    /// Returns true if the given bounding box is entirely contained
    /// in this bounding box.
    bool contains( const BBox &bbox ) const {
      using namespace vector_containment_comparison;
      return ((bbox.m_min >= m_min) && (bbox.m_max <= m_max));
    }

    /// Returns true if the given bounding box intersects this
    /// bounding box.
    bool intersects( const BBox& bbox ) const {
      for( int i=0; i<DimN; ++i ) {
        if( m_min[i] >= bbox.m_max[i] ||
            m_max[i] <= bbox.m_min[i] )
          return false;
      }
      return true;
    }

    /// Returns the size (i.e. the diagonal vector) of the bounding box.
    Vector<RealT, DimN> size() const { return (m_max - m_min); }

    /// Returns the center point of the bounding box.
    Vector<RealT, DimN> center() const { return 0.5 * (m_min + m_max); }

    /// Returns the minimal point of the bounding box.
    Vector<RealT, DimN> const& min() const { return m_min; }

    /// Returns the maximal point of the bounding box.
    Vector<RealT, DimN> const& max() const { return m_max; }

    /// Returns the minimal point of the bounding box (non-const overload).
    Vector<RealT, DimN> &min() { return m_min; }

    /// Returns the maximal point of the bounding box (non-const overload).
    Vector<RealT, DimN> &max() { return m_max; }

    /// Returns the width (i.e. size in the first dimension) of the
    /// bounding box.
    RealT width() const {
      BOOST_STATIC_ASSERT( DimN >= 1 );
      return m_max[0] - m_min[0];
    }

    /// Returns the height (i.e. size in the second dimension) of the
    /// bounding box.  Only valid for bouding boxes of dimension two
    /// or greater.
    RealT height() const {
      BOOST_STATIC_ASSERT( DimN >= 2 );
      return m_max[1] - m_min[1];
    }

    /// Returns true if the bounding box is empty (i.e. degenerate).
    bool empty() const {
      for( int i=0; i<DimN; ++i )
        if( m_min[i] >= m_max[i] ) return true;
      return false;
    }

    /// Scales the bounding box relative to the origin.
    template <class ScalarT>
    BBox& operator*=( ScalarT s ) {
      m_min *= s;
      m_max *= s;
      return *this;
    }

    /// Scales the bounding box relative to the origin.
    template <class ScalarT>
    BBox& operator/=( ScalarT s ) {
      m_min /= s;
      m_max /= s;
      return *this;
    }

    /// Offsets the bounding box by the given vector.
    BBox& operator+=( Vector<RealT,DimN> const& v ) {
      m_min += v;
      m_max += v;
      return *this;
    }

    /// Offsets the bounding box by the negation of the given vector.
    BBox& operator-=( Vector<RealT,DimN> const& v ) {
      m_min -= v;
      m_max -= v;
      return *this;
    }

  private:
    Vector<RealT, DimN> m_min, m_max;
  };

  /// Scales a bounding box relative to the origin.
  template <class RealT, int DimN, class ScalarT>
  inline BBox<RealT,DimN> operator*( BBox<RealT,DimN> const& bbox, ScalarT s ) {
    BBox<RealT,DimN> result = bbox;
    result *= s;
    return result;
  }

  /// Scales a bounding box relative to the origin.
  template <class RealT, int DimN, class ScalarT>
  inline BBox<RealT,DimN> operator/( BBox<RealT,DimN> const& bbox, ScalarT s ) {
    BBox<RealT,DimN> result = bbox;
    result /= s;
    return result;
  }

  /// Scales a bounding box relative to the origin.
  template <class RealT, int DimN, class ScalarT>
  inline BBox<RealT,DimN> operator*( ScalarT s, BBox<RealT,DimN> const& bbox ) {
    return bbox * s;
  }
  
  /// Offsets a bounding box by the given vector.
  template <class RealT, int DimN, class VectorT>
  inline BBox<RealT,DimN> operator+( BBox<RealT,DimN> const& bbox, VectorBase<VectorT> const& v ) {
    BBox<RealT,DimN> result = bbox;
    result += v.impl();
    return result;
  }

  /// Offsets a bounding box by the given vector.
  template <class RealT, int DimN, class VectorT>
  inline BBox<RealT,DimN> operator+( VectorBase<VectorT> const& v, BBox<RealT,DimN> const& bbox ) {
    return bbox + v;
  }

  /// Offsets a bounding box by the negation of the given vector.
  template <class RealT, int DimN, class VectorT>
  inline BBox<RealT,DimN> operator-( BBox<RealT,DimN> const& bbox, VectorBase<VectorT> const& v ) {
    BBox<RealT,DimN> result = bbox;
    result -= v.impl();
    return result;
  }

  /// Equality of two bounding boxes.
  template <class Real1T, class Real2T, int DimN>
  inline bool operator==( BBox<Real1T,DimN> const& bbox1, BBox<Real2T,DimN> const& bbox2 ) {
    return bbox1.min()==bbox2.min() && bbox1.max()==bbox2.max();
  }

  /// Inequality of two bounding boxes.
  template <class Real1T, class Real2T, int DimN>
  inline bool operator!=( BBox<Real1T,DimN> const& bbox1, BBox<Real2T,DimN> const& bbox2 ) {
    return bbox1.min()!=bbox2.min() || bbox1.max()!=bbox2.max();
  }

  /// Writes a bounding box to an ostream.
  template <class RealT, int DimN>
  std::ostream& operator<<( std::ostream& os, BBox<RealT,DimN> const& bbox ) {
    return os << "(" << bbox.min() << "-" << bbox.max() << ")";
  }

} // namespace math

  // Convenience typedefs
  using math::BBox;
  typedef BBox<float64, 2> BBox2;
  typedef BBox<float64, 3> BBox3;
  typedef BBox<float64, 4> BBox4;
  typedef BBox<float32, 2> BBox2f;
  typedef BBox<float32, 3> BBox3f;
  typedef BBox<float32, 4> BBox4f;
  typedef BBox<int32, 2> BBox2i;
  typedef BBox<int32, 3> BBox3i;
  typedef BBox<int32, 4> BBox4i;

} // namespace vw

#endif // __VW_MATH__BBOX_H__
