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
/// Provides generic bounding shapes.
#ifndef __VW_MATH__BBOX_H__
#define __VW_MATH__BBOX_H__

#include <limits>
#include <math.h>

#include <boost/static_assert.hpp>

#include <vw/Math/Vector.h>

namespace vw {
namespace math {

  /// \cond INTERNAL
  namespace vector_containment_comparison {
    template <class VectorT1, class VectorT2>
    inline bool operator<( VectorBase<VectorT1> const& v1, VectorBase<VectorT2> const& v2 ) {
      VW_ASSERT( v1.impl().size() == v2.impl().size(), ArgumentErr() << "Cannot compare vectors of different length." );
      for( int i=0; i<v1.impl().size(); ++i)
        if( ! ( v1.impl()[i] < v2.impl()[i] ) ) return false;
      return true;
    }

    template <class VectorT1, class VectorT2>
    inline bool operator<=( VectorBase<VectorT1> const& v1, VectorBase<VectorT2> const& v2 ) {
      VW_ASSERT( v1.impl().size() == v2.impl().size(), ArgumentErr() << "Cannot compare vectors of different length." );
      for( int i=0; i<v1.impl().size(); ++i)
        if( ! ( v1.impl()[i] <= v2.impl()[i] ) ) return false;
      return true;
    }

    template <class VectorT1, class VectorT2>
    inline bool operator>( VectorBase<VectorT1> const& v1, VectorBase<VectorT2> const& v2 ) {
      VW_ASSERT( v1.impl().size() == v2.impl().size(), ArgumentErr() << "Cannot compare vectors of different length." );
      for( int i=0; i<v1.impl().size(); ++i)
        if( ! ( v1.impl()[i] > v2.impl()[i] ) ) return false;
      return true;
    }

    template <class VectorT1, class VectorT2>
    inline bool operator>=( VectorBase<VectorT1> const& v1, VectorBase<VectorT2> const& v2 ) {
      VW_ASSERT( v1.impl().size() == v2.impl().size(), ArgumentErr() << "Cannot compare vectors of different length." );
      for( int i=0; i<v1.impl().size(); ++i)
        if( ! ( v1.impl()[i] >= v2.impl()[i] ) ) return false;
      return true;
    }

    template <class VectorT1, class VectorT2>
    inline VectorT1 max( VectorBase<VectorT1> const& v1, VectorBase<VectorT2> const& v2 ) {
      VW_ASSERT( v1.impl().size() == v2.impl().size(), ArgumentErr() << "Cannot compute max of vectors of different length." );
      VectorT1 v3;
      v3.set_size(v1.impl().size());
      for( int i=0; i<v1.impl().size(); ++i)
        v3[i] = std::max(v1.impl()[i], v2.impl()[i]);
      return v3;
    }

    template <class VectorT1, class VectorT2>
    inline VectorT1 min( VectorBase<VectorT1> const& v1, VectorBase<VectorT2> const& v2 ) {
      VW_ASSERT( v1.impl().size() == v2.impl().size(), ArgumentErr() << "Cannot compute min of vectors of different length." );
      VectorT1 v3;
      v3.set_size(v1.impl().size());
      for( int i=0; i<v1.impl().size(); ++i)
        v3[i] = std::min(v1.impl()[i], v2.impl()[i]);
      return v3;
    }

  } // namespace vector_containment_comparison
  /// \endcond


  // *******************************************************************
  // class BShapeBase
  // *******************************************************************

  /// A CRTP base class for general n-dimensional bounding shapes.  
  /// Provides a mechanism for restricting function arguments to 
  /// bounding shapes, provides general bounding-shape operations,
  /// and provides the various arithmetic assignment operators.
  template <class BShapeT, class RealT, int DimN>
  class BShapeBase {
  typedef Vector<RealT, DimN> BShapeVectorT;
  public:

    /// Returns the derived implementation type.
    BShapeT& shape_impl() { return *static_cast<BShapeT*>(this); }
    
    /// Returns the derived implementation type.
    BShapeT const& shape_impl() const { return *static_cast<BShapeT const*>(this); }

    /// Grows a bounding shape to include the given point.
    template <class VectorT>
    void grow( VectorBase<VectorT> const& point ) {
      shape_impl().grow(point);
    }

    /// Expands this bounding shape by the given offset in every direction.
    void expand( RealT offset ) {
      shape_impl().expand(offset);
    }

    /// Contracts this bounding shape by the given offset in every direction.
    void contract( RealT offset ) {
      shape_impl().contract(offset);
    }

    /// Returns true if the given point is contained in the bounding shape.
    template <class VectorT>
    bool contains( VectorBase<VectorT> const& point ) const {
      return shape_impl().contains(point);
    }

    /// Returns the center point of the bounding shape.
    BShapeVectorT center() const {
      return shape_impl().center();
    }

    /// Returns true if the bounding shape is empty (i.e. degenerate).
    bool empty() const {
      return shape_impl().empty();
    }

    /// Scales the bounding shape relative to the origin.
    template <class ScalarT>
    BShapeT& operator*=( ScalarT s ) {
      return shape_impl() = shape_impl() * s;
    }

    /// Scales the bounding shape relative to the origin.
    template <class ScalarT>
    BShapeT& operator/=( ScalarT s ) {
      return shape_impl() = shape_impl() / s;
    }

    /// Offsets the bounding shape by the given vector.
    template <class VectorT>
    BShapeT& operator+=( VectorBase<VectorT> const& v ) {
      return shape_impl() = shape_impl() + v;
    }

    /// Offsets the bounding shape by the negation of the given vector.
    template <class VectorT>
    BShapeT& operator-=( VectorBase<VectorT> const& v ) {
      return shape_impl() = shape_impl() - v;
    }
  };
  
  /// Scales a bounding shape relative to the origin.
  template <class BShapeT, class RealT, int DimN, class ScalarT>
  inline BShapeT operator*( BShapeBase<BShapeT, RealT, DimN> const& bshape, ScalarT s ) {
    BShapeT result = bshape.shape_impl();
    result *= s;
    return result;
  }

  /// Scales a bounding shape relative to the origin.
  template <class BShapeT, class RealT, int DimN, class ScalarT>
  inline BShapeT operator/( BShapeBase<BShapeT, RealT, DimN> const& bshape, ScalarT s ) {
    BShapeT result = bshape.shape_impl();
    result /= s;
    return result;
  }

  /// Scales a bounding shape relative to the origin.
  template <class BShapeT, class RealT, int DimN, class ScalarT>
  inline BShapeT operator*( ScalarT s, BShapeBase<BShapeT, RealT, DimN> const& bshape ) {
    return bshape * s;
  }
  
  /// Offsets a bounding shape by the given vector.
  template <class BShapeT, class RealT, int DimN, class VectorT>
  inline BShapeT operator+( BShapeBase<BShapeT, RealT, DimN> const& bshape, VectorBase<VectorT> const& v ) {
    BShapeT result = bshape.shape_impl();
    result += v.impl();
    return result;
  }

  /// Offsets a bounding shape by the given vector.
  template <class BShapeT, class RealT, int DimN, class VectorT>
  inline BShapeT operator+( VectorBase<VectorT> const& v, BShapeBase<BShapeT, RealT, DimN> const& bshape ) {
    return bshape + v;
  }

  /// Offsets a bounding shape by the negation of the given vector.
  template <class BShapeT, class RealT, int DimN, class VectorT>
  inline BShapeT operator-( BShapeBase<BShapeT, RealT, DimN> const& bshape, VectorBase<VectorT> const& v ) {
    BShapeT result = bshape.shape_impl();
    result -= v.impl();
    return result;
  }
  
  /// Equality of two bounding shapes.
  template <class BShapeT, class RealT, int DimN>
  inline bool operator==( BShapeBase<BShapeT, RealT, DimN> const& bshape1, BShapeBase<BShapeT, RealT, DimN> const& bshape2 ) {
    return bshape1.shape_impl() == bshape2.shape_impl();
  }
  
  /// Inequality of two bounding shapes.
  template <class BShapeT, class RealT, int DimN>
  inline bool operator!=( BShapeBase<BShapeT, RealT, DimN> const& bshape1, BShapeBase<BShapeT, RealT, DimN> const& bshape2 ) {
    return bshape1.shape_impl() != bshape2.shape_impl();
  }

  /// Writes a bounding shape to an ostream.
  template <class BShapeT, class RealT, int DimN>
  std::ostream& operator<<( std::ostream& os, BShapeBase<BShapeT, RealT, DimN> const& bshape ) {
    return os << bshape.shape_impl();
  }

  // *******************************************************************
  // class BBoxBase
  // *******************************************************************

  /// A general n-dimensional axis-aligned bounding box class,
  /// represented by vectors pointing to the minimal and maximal
  /// corners.
  template <class BBoxT, class RealT, int DimN>
  class BBoxBase : public BShapeBase<BBoxBase<BBoxT, RealT, DimN>, RealT, DimN> {
  public:

    /// Default constructor.  Constructs the ultimate empty bounding 
    /// box, whose limits are at the opposite corners of the underlying 
    /// numeric space.  This is a useful starting point if you intend 
    /// to grow your bounding box to fit a collection of items.
    BBoxBase() {
      // Make sure we have a type for which we know limits
      BOOST_STATIC_ASSERT(std::numeric_limits<RealT>::is_specialized);
      if (std::numeric_limits<RealT>::is_integer) {
        for (unsigned i = 0; i < m_min.size(); ++i) {
          m_min[i] = std::numeric_limits<RealT>::max();
          m_max[i] = std::numeric_limits<RealT>::min();
        }
      }
      else {
        for (unsigned i = 0; i < m_min.size(); ++i) {
          m_min[i] = std::numeric_limits<RealT>::max();
          m_max[i] = -std::numeric_limits<RealT>::max();
        }
      }
    }

    /// Constructs a bounding box with the given minimal and maximal 
    /// points.
    template <class VectorT1, class VectorT2>
    BBoxBase( VectorBase<VectorT1> const& min, VectorBase<VectorT2> const& max ) :
      m_min( min ), m_max( max ) {}

    /// Returns the derived implementation type.
    BBoxT& impl() { return *static_cast<BBoxT*>(this); }
    
    /// Returns the derived implementation type.
    BBoxT const& impl() const { return *static_cast<BBoxT const*>(this); }

    /// Grows a bounding box to include the given point.
    template <class VectorT>
    void grow( VectorBase<VectorT> const& point ) {
      VW_ASSERT(m_min.size() == 0 || point.impl().size() == m_min.size(), ArgumentErr() << "Vector must have dimension " << m_min.size() << ".");
      if (m_min.size() == 0) {
        m_min = point;
        m_max = point;
      }
      else {
        for (unsigned i = 0; i < m_min.size(); ++i) {
          if (point.impl()[i] > m_max[i])
            m_max[i] = point.impl()[i];
          if (point.impl()[i] < m_min[i])
            m_min[i] = point.impl()[i];
        }
      }
    }
    
    /// Grows a bounding box to include the given bounding box.
    template <class BBoxT1, class RealT1, int DimN1>
    void grow( BBoxBase<BBoxT1, RealT1, DimN1> const& bbox ) {
      grow(bbox.min());
      grow(bbox.max());
    }

    /// Crops (intersects) this bounding box to the given bounding box.
    template <class BBoxT1, class RealT1, int DimN1>
    void crop( BBoxBase<BBoxT1, RealT1, DimN1> const& bbox ) {
      VW_ASSERT(m_min.size() == 0 || bbox.min().size() == m_min.size(), ArgumentErr() << "BBox must have dimension " << m_min.size() << ".");
      for( unsigned i=0; i<m_min.size(); ++i ) {
        if( m_min[i] < bbox.min()[i] )
          if( m_max[i] < bbox.min()[i] )
            m_min[i] = m_max[i]; 
          else
            m_min[i] = bbox.min()[i];

        if( m_max[i] > bbox.max()[i] )
          if ( m_min[i] > bbox.max()[i] )
            m_max[i] = m_min[i];
          else
            m_max[i] = bbox.max()[i];
      }
    }

    /// Expands this bounding box by the given offset in every direction.
    void expand( RealT offset ) {
      for( int i=0; i<m_min.size(); ++i ) {
        m_min[i] -= offset;
        m_max[i] += offset;
      }
    }

    /// Contracts this bounding box by the given offset in every direction.
    void contract( RealT offset ) {
      for( unsigned i=0; i<m_min.size(); ++i ) {
        m_min[i] += offset;
        m_max[i] -= offset;
      }
    }

    /// Returns true if the given point is contained in the bounding box.
    template <class VectorT>
    bool contains( const VectorBase<VectorT> &point ) const {
      using namespace vector_containment_comparison;
      VW_ASSERT(m_min.size() == 0 || point.impl().size() == m_min.size(), ArgumentErr() << "Vector must have dimension " << m_min.size() << ".");
      return ((point >= m_min) && (point < m_max));
    }

    /// Returns true if the given bounding box is entirely contained
    /// in this bounding box.
    template <class BBoxT1, class RealT1, int DimN1>
    bool contains( const BBoxBase<BBoxT1, RealT1, DimN1> &bbox ) const {
      using namespace vector_containment_comparison;
      VW_ASSERT(m_min.size() == 0 || bbox.min().size() == m_min.size(), ArgumentErr() << "BBox must have dimension " << m_min.size() << ".");
      return ((bbox.min() >= m_min) && (bbox.max() <= m_max));
    }

    /// Returns true if the given bounding box intersects this
    /// bounding box.
    template <class BBoxT1, class RealT1, int DimN1>
    bool intersects( const BBoxBase<BBoxT1, RealT1, DimN1>& bbox ) const {
      VW_ASSERT(m_min.size() == 0 || bbox.min().size() == m_min.size(), ArgumentErr() << "BBox must have dimension " << m_min.size() << ".");
      for( unsigned i=0; i<m_min.size(); ++i ) {
        if( m_min[i] >= bbox.max()[i] ||
            m_max[i] <= bbox.min()[i] )
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
      return impl().width();
    }

    /// Returns the height (i.e. size in the second dimension) of the
    /// bounding box.  Only valid for bouding boxes of dimension two
    /// or greater.
    RealT height() const {
      return impl().height();
    }

    /// Returns the depth (i.e. size in the third dimension) of the
    /// bounding box.  Only valid for bouding boxes of dimension three
    /// or greater.
    RealT depth() const {
      return impl().depth();
    }

    /// Returns true if the bounding box is empty (i.e. degenerate).
    bool empty() const {
      for( unsigned i=0; i<m_min.size(); ++i )
        if( m_min[i] >= m_max[i] ) return true;
      return (m_min.size() <= 0);
    }

    /// Scales the bounding box relative to the origin.
    template <class ScalarT>
    BBoxBase& operator*=( ScalarT s ) {
      m_min *= s;
      m_max *= s;
      return *this;
    }

    /// Scales the bounding box relative to the origin.
    template <class ScalarT>
    BBoxBase& operator/=( ScalarT s ) {
      m_min /= s;
      m_max /= s;
      return *this;
    }

    /// Offsets the bounding box by the given vector.
    template <class VectorT>
    BBoxBase& operator+=( VectorBase<VectorT> const& v ) {
      m_min += v;
      m_max += v;
      return *this;
    }

    /// Offsets the bounding box by the negation of the given vector.
    template <class VectorT>
    BBoxBase& operator-=( VectorBase<VectorT> const& v ) {
      m_min -= v;
      m_max -= v;
      return *this;
    }

  protected:
    Vector<RealT, DimN> m_min, m_max;
  };
  
  /// Scales a bounding box relative to the origin.
  template <class BBoxT, class RealT, int DimN, class ScalarT>
  inline BBoxT operator*( BBoxBase<BBoxT, RealT, DimN> const& bbox, ScalarT s ) {
    BBoxT result = bbox.impl();
    result *= s;
    return result;
  }

  /// Scales a bounding box relative to the origin.
  template <class BBoxT, class RealT, int DimN, class ScalarT>
  inline BBoxT operator/( BBoxBase<BBoxT, RealT, DimN> const& bbox, ScalarT s ) {
    BBoxT result = bbox.impl();
    result /= s;
    return result;
  }

  /// Scales a bounding box relative to the origin.
  template <class BBoxT, class RealT, int DimN, class ScalarT>
  inline BBoxT operator*( ScalarT s, BBoxBase<BBoxT, RealT, DimN> const& bbox ) {
    return bbox * s;
  }
  
  /// Offsets a bounding box by the given vector.
  template <class BBoxT, class RealT, int DimN, class VectorT>
  inline BBoxT operator+( BBoxBase<BBoxT, RealT, DimN> const& bbox, VectorBase<VectorT> const& v ) {
    BBoxT result = bbox.impl();
    result += v.impl();
    return result;
  }

  /// Offsets a bounding box by the given vector.
  template <class BBoxT, class RealT, int DimN, class VectorT>
  inline BBoxT operator+( VectorBase<VectorT> const& v, BBoxBase<BBoxT, RealT, DimN> const& bbox ) {
    return bbox + v;
  }

  /// Offsets a bounding box by the negation of the given vector.
  template <class BBoxT, class RealT, int DimN, class VectorT>
  inline BBoxT operator-( BBoxBase<BBoxT, RealT, DimN> const& bbox, VectorBase<VectorT> const& v ) {
    BBoxT result = bbox.impl();
    result -= v.impl();
    return result;
  }

  /// Equality of two bounding boxes.
  template <class BBoxT1, class RealT1, int DimN1, class BBoxT2, class RealT2, int DimN2>
  inline bool operator==( BBoxBase<BBoxT1,RealT1,DimN1> const& bbox1, BBoxBase<BBoxT2,RealT2,DimN2> const& bbox2 ) {
    return bbox1.min()==bbox2.min() && bbox1.max()==bbox2.max();
  }

  /// Inequality of two bounding boxes.
  template <class BBoxT1, class RealT1, int DimN1, class BBoxT2, class RealT2, int DimN2>
  inline bool operator!=( BBoxBase<BBoxT1,RealT1,DimN1> const& bbox1, BBoxBase<BBoxT2,RealT2,DimN2> const& bbox2 ) {
    return bbox1.min()!=bbox2.min() || bbox1.max()!=bbox2.max();
  }
  
  /// Writes a bounding box to an ostream.
  template <class BBoxT, class RealT, int DimN>
  std::ostream& operator<<( std::ostream& os, BBoxBase<BBoxT,RealT,DimN> const& bbox ) {
    return os << "(" << bbox.min() << "-" << bbox.max() << ")";
  }

  // *******************************************************************
  // class BBox
  // *******************************************************************

  /// A general fixed-dimensional axis-aligned bounding box class,
  /// represented by vectors pointing to the minimal and maximal
  /// corners.
  template <class RealT, int DimN = 0>
  class BBox : public BBoxBase<BBox<RealT, DimN>, RealT, DimN> {
  public:

    /// Default constructor.  Constructs the ultimate empty bounding 
    /// box, whose limits are at the opposite corners of the underlying 
    /// numeric space.  This is a useful starting point if you intend 
    /// to grow your bounding box to fit a collection of items.
    BBox() : BBoxBase<BBox<RealT, DimN>, RealT, DimN>() {}

    /// Constructs a bounding box with the given minimal and maximal 
    /// points.
    BBox( Vector<RealT, DimN> const& min, Vector<RealT, DimN> const& max ) :
      BBoxBase<BBox<RealT, DimN>, RealT, DimN>( min, max ) {}

    /// Constructs a 2D bounding box with the given minimal point
    /// coordinates and dimensions.  (Only valid for 2D bouding
    /// boxes.)
    BBox( RealT minx, RealT miny, RealT width, RealT height )
      : BBoxBase<BBox<RealT, DimN>, RealT, DimN>( 
          Vector<RealT,2>(minx,miny), 
          Vector<RealT,2>(minx+width,miny+height) ) {
      BOOST_STATIC_ASSERT( DimN==2 );
    }

    /// Constructs a 3D bounding box with the given minimal point
    /// coordinates and dimensions.  (Only valid for 3D bouding
    /// boxes.)
    BBox( RealT minx, RealT miny, RealT minz, RealT width, RealT height, RealT depth )
      : BBoxBase<BBox<RealT, DimN>, RealT, DimN>( 
          Vector<RealT,3>(minx,miny,minz), 
          Vector<RealT,3>(minx+width,miny+height,minz+depth) ) {
      BOOST_STATIC_ASSERT( DimN==3 );
    }

    /// Copy constructor.
    template <class BBoxT1, class RealT1, int DimN1>
    BBox( BBoxBase<BBoxT1, RealT1, DimN1> const& bbox )
      : BBoxBase<BBox<RealT, DimN>, RealT, DimN>( bbox.min(), bbox.max() ) {}
    
    /// Copy assignment operator.
    template <class BBoxT1, class RealT1, int DimN1>
    BBox& operator=( BBoxBase<BBoxT1, RealT1, DimN1> const& bbox ) {
      this->min() = bbox.min();
      this->max() = bbox.max();
      return *this;
    }

    /// Returns the width (i.e. size in the first dimension) of the
    /// bounding box.
    RealT width() const {
      BOOST_STATIC_ASSERT( DimN >= 1 );
      return this->max()[0] - this->min()[0];
    }

    /// Returns the height (i.e. size in the second dimension) of the
    /// bounding box.  Only valid for bouding boxes of dimension two
    /// or greater.
    RealT height() const {
      BOOST_STATIC_ASSERT( DimN >= 2 );
      return this->max()[1] - this->min()[1];
    }

    /// Returns the depth (i.e. size in the third dimension) of the
    /// bounding box.  Only valid for bouding boxes of dimension three
    /// or greater.
    RealT depth() const {
      BOOST_STATIC_ASSERT( DimN >= 3 );
      return this->max()[2] - this->min()[2];
    }
  };

  /// A general arbitrary-dimensional axis-aligned bounding box class,
  /// represented by vectors pointing to the minimal and maximal
  /// corners.
  template <class RealT>
  class BBox<RealT, 0> : public BBoxBase<BBox<RealT, 0>, RealT, 0> {
  public:

    /// Default constructor.  Constructs the ultimate empty bounding 
    /// box, whose limits are at the opposite corners of the underlying 
    /// numeric space.  This is a useful starting point if you intend 
    /// to grow your bounding box to fit a collection of items.
    BBox() : BBoxBase<BBox<RealT, 0>, RealT, 0>() {}

    /// Constructs a bounding box with the given minimal and maximal 
    /// points.
    BBox( Vector<RealT, 0> const& min, Vector<RealT, 0> const& max ) :
      BBoxBase<BBox<RealT, 0>, RealT, 0>( min, max ) {}

    /// Constructs a 2D bounding box with the given minimal point
    /// coordinates and dimensions.  (Only valid for 2D bouding
    /// boxes.)
    BBox( RealT minx, RealT miny, RealT width, RealT height )
      : BBoxBase<BBox<RealT, 0>, RealT, 0>( 
          Vector<RealT,2>(minx,miny), 
          Vector<RealT,2>(minx+width,miny+height) ) {}

    /// Constructs a 3D bounding box with the given minimal point
    /// coordinates and dimensions.  (Only valid for 3D bouding
    /// boxes.)
    BBox( RealT minx, RealT miny, RealT minz, RealT width, RealT height, RealT depth )
      : BBoxBase<BBox<RealT, 0>, RealT, 0>( 
          Vector<RealT,3>(minx,miny,minz), 
          Vector<RealT,3>(minx+width,miny+height,minz+depth) ) {}

    /// Copy constructor.
    template <class BBoxT1, class RealT1, int DimN1>
    BBox( BBoxBase<BBoxT1, RealT1, DimN1> const& bbox )
      : BBoxBase<BBox<RealT, 0>, RealT, 0>( bbox.min(), bbox.max() ) {}
    
    /// Copy assignment operator.
    template <class BBoxT1, class RealT1, int DimN1>
    BBox& operator=( BBoxBase<BBoxT1, RealT1, DimN1> const& bbox ) {
      this->min() = bbox.min();
      this->max() = bbox.max();
      return *this;
    }

    /// Returns the width (i.e. size in the first dimension) of the
    /// bounding box.
    RealT width() const {
      VW_ASSERT(this->min().size() >= 1, LogicErr() << "BBox must be of dimension >= 1 to get width.");
      return this->max()[0] - this->min()[0];
    }

    /// Returns the height (i.e. size in the second dimension) of the
    /// bounding box.  Only valid for bouding boxes of dimension two
    /// or greater.
    RealT height() const {
      VW_ASSERT(this->min().size() >= 2, LogicErr() << "BBox must be of dimension >= 2 to get height.");
      return this->max()[1] - this->min()[1];
    }

    /// Returns the depth (i.e. size in the third dimension) of the
    /// bounding box.  Only valid for bouding boxes of dimension three
    /// or greater.
    RealT depth() const {
      VW_ASSERT(this->min().size() >= 3, LogicErr() << "BBox must be of dimension >= 3 to get depth.");
      return this->max()[2] - this->min()[2];
    }
  };

  // *******************************************************************
  // class BBallBase
  // *******************************************************************

  /// A general n-dimensional bounding ball class,
  /// represented by a vector pointing to the center, and a radius.
  template <class BBallT, class RealT, int DimN>
  class BBallBase : public BShapeBase<BBallBase<BBallT, RealT, DimN>, RealT, DimN> {
  public:

    /// Default constructor. Constructs an empty bounding ball.
    BBallBase() : m_radius( 0 ) {}

    /// Constructs a bounding ball with the given center and radius.
    template <class VectorT>
    BBallBase( VectorBase<VectorT> const& center, RealT radius ) :
      m_center( center ), m_radius( radius ) {}

    /// Constructs a 2D bounding ball with the given center point
    /// coordinates and radius.  (Only valid for 2D bouding
    /// balls.)
    BBallBase( RealT centerx, RealT centery, RealT radius )
      : m_center( Vector<RealT,2>(centerx,centery) ), 
        m_radius( radius ) {}

    /// Constructs a 3D bounding ball with the given center point
    /// coordinates and radius.  (Only valid for 3D bouding
    /// balls.)
    BBallBase( RealT centerx, RealT centery, RealT centerz, RealT radius )
      : m_center( Vector<RealT,3>(centerx,centery,centerz) ), 
        m_radius( radius ) {}

    /// Returns the derived implementation type.
    BBallT& impl() { return *static_cast<BBallT*>(this); }
    
    /// Returns the derived implementation type.
    BBallT const& impl() const { return *static_cast<BBallT const*>(this); }

    /// Grows a bounding ball to include the given point.
    template <class VectorT>
    void grow( VectorBase<VectorT> const& point ) {
      VW_ASSERT(point.impl().size() == m_center.size(), ArgumentErr() << "Vector must have dimension " << m_center.size() << ".");
      Vector<RealT, DimN> v = point - m_center;
      RealT dist = norm_2(v);
      RealT radius_d;
      if (dist > m_radius) {
        radius_d = (dist - m_radius)/2;
        m_center += v/dist*radius_d;
        m_radius += radius_d;
      }
    }
    
    /// Grows a bounding ball to include the given bounding ball.
    template <class BBallT1, class RealT1, int DimN1>
    void grow( BBallBase<BBallT1, RealT1, DimN1> const& bball ) {
      Vector<RealT, DimN> v = bball.m_center - m_center;
      grow(bball.m_center + normalize(v)*bball.m_radius);
    }

    /// Crops (intersects) this bounding ball to the given bounding ball.
    template <class BBallT1, class RealT1, int DimN1>
    void crop( BBallBase<BBallT1, RealT1, DimN1> const& bball ) {
      Vector<RealT, DimN> v = bball.m_center - m_center;
      RealT dist = norm_2(v);
      if (dist > (m_radius + bball.m_radius)) {
        // intersects(bball) == false
        m_radius = 0;
      }
      else if (dist <= (m_radius - bball.m_radius)) {
        // contains(bball) == true
        m_center = bball.m_center;
        m_radius = bball.m_radius;
      }
      else if (dist <= (bball.m_radius - m_radius)) {
        // bball.contains(*this) == true
        // smallest BBall is *this
        return;
      }
      else if (m_radius >= dist) {
        // smallest BBall is bball
        m_center = bball.m_center;
        m_radius = bball.m_radius;
      }
      else if (bball.m_radius >= dist) {
        // smallest BBall is *this
        return;
      }
      else {
        // normal intersection
        RealT x = (dist + (m_radius*m_radius - bball.m_radius*bball.m_radius)/dist)/2;
        m_center += v/dist*x;
        m_radius = std::sqrt(m_radius*m_radius - x*x);
      }
    }

    /// Expands this bounding ball by the given offset in every direction.
    void expand( RealT offset ) {
      m_radius += offset;
    }

    /// Contracts this bounding ball by the given offset in every direction.
    void contract( RealT offset ) {
      m_radius -= offset;
    }

    /// Returns true if the given point is contained in the bounding ball.
    template <class VectorT>
    bool contains( const VectorBase<VectorT> &point ) const {
      VW_ASSERT(point.impl().size() == m_center.size(), ArgumentErr() << "Vector must have dimension " << m_center.size() << ".");
      return (norm_2(m_center - point) <= m_radius);
    }

    /// Returns true if the given bounding ball is entirely contained
    /// in this bounding ball.
    template <class BBallT1, class RealT1, int DimN1>
    bool contains( const BBallBase<BBallT1, RealT1, DimN1> &bball ) const {
      return (norm_2(m_center - bball.m_center) <= (m_radius - bball.m_radius));
    }

    /// Returns true if the given bounding ball intersects this
    /// bounding ball.
    template <class BBallT1, class RealT1, int DimN1>
    bool intersects( const BBallBase<BBallT1, RealT1, DimN1>& bball ) const {
      return (norm_2(m_center - bball.m_center) <= (m_radius + bball.m_radius));
    }

    /// Returns the size (i.e. the diameter) of the bounding ball.
    RealT size() const { return 2*m_radius; }

    /// Returns the center point of the bounding ball.
    Vector<RealT, DimN> center() const { return m_center; }

    /// Returns the center point of the bounding ball.
    Vector<RealT, DimN> const& center_() const { return m_center; }

    /// Returns the center point of the bounding ball.
    Vector<RealT, DimN> &center_() { return m_center; }

    /// Returns the radius of the bounding ball.
    Vector<RealT, DimN> const& radius() const { return m_radius; }

    /// Returns the radius of the bounding ball.
    Vector<RealT, DimN> &radius() { return m_radius; }

    /// Returns true if the bounding ball is empty (i.e. degenerate).
    bool empty() const {
      return (m_radius <= 0);
    }

    /// Scales the bounding ball relative to the origin.
    template <class ScalarT>
    BBallBase& operator*=( ScalarT s ) {
      m_center *= s;
      m_radius *= s;
      return *this;
    }

    /// Scales the bounding ball relative to the origin.
    template <class ScalarT>
    BBallBase& operator/=( ScalarT s ) {
      m_center /= s;
      m_radius /= s;
      return *this;
    }

    /// Offsets the bounding ball by the given vector.
    template <class VectorT>
    BBallBase& operator+=( VectorBase<VectorT> const& v ) {
      m_center += v;
      return *this;
    }

    /// Offsets the bounding ball by the negation of the given vector.
    template <class VectorT>
    BBallBase& operator-=( VectorBase<VectorT> const& v ) {
      m_center -= v;
      return *this;
    }

  private:
    Vector<RealT, DimN> m_center;
    RealT m_radius;
  };
  
  /// Scales a bounding ball relative to the origin.
  template <class BBallT, class RealT, int DimN, class ScalarT>
  inline BBallT operator*( BBallBase<BBallT, RealT, DimN> const& bball, ScalarT s ) {
    BBallT result = bball.impl();
    result *= s;
    return result;
  }

  /// Scales a bounding ball relative to the origin.
  template <class BBallT, class RealT, int DimN, class ScalarT>
  inline BBallT operator/( BBallBase<BBallT, RealT, DimN> const& bball, ScalarT s ) {
    BBallT result = bball.impl();
    result /= s;
    return result;
  }

  /// Scales a bounding ball relative to the origin.
  template <class BBallT, class RealT, int DimN, class ScalarT>
  inline BBallT operator*( ScalarT s, BBallBase<BBallT, RealT, DimN> const& bball ) {
    return bball * s;
  }
  
  /// Offsets a bounding ball by the given vector.
  template <class BBallT, class RealT, int DimN, class VectorT>
  inline BBallT operator+( BBallBase<BBallT, RealT, DimN> const& bball, VectorBase<VectorT> const& v ) {
    BBallT result = bball.impl();
    result += v.impl();
    return result;
  }

  /// Offsets a bounding ball by the given vector.
  template <class BBallT, class RealT, int DimN, class VectorT>
  inline BBallT operator+( VectorBase<VectorT> const& v, BBallBase<BBallT, RealT, DimN> const& bball ) {
    return bball + v;
  }

  /// Offsets a bounding ball by the negation of the given vector.
  template <class BBallT, class RealT, int DimN, class VectorT>
  inline BBallT operator-( BBallBase<BBallT, RealT, DimN> const& bball, VectorBase<VectorT> const& v ) {
    BBallT result = bball.impl();
    result -= v.impl();
    return result;
  }
  
  /// Equality of two bounding balls.
  template <class BBallT1, class RealT1, int DimN1, class BBallT2, class RealT2, int DimN2>
  inline bool operator==( BBallBase<BBallT1,RealT1,DimN1> const& bball1, BBallBase<BBallT2,RealT2,DimN2> const& bball2 ) {
    return bball1.center_()==bball2.center_() && bball1.radius()==bball2.radius();
  }
  
  /// Inequality of two bounding balls.
  template <class BBallT1, class RealT1, int DimN1, class BBallT2, class RealT2, int DimN2>
  inline bool operator!=( BBallBase<BBallT1,RealT1,DimN1> const& bball1, BBallBase<BBallT2,RealT2,DimN2> const& bball2 ) {
    return bball1.center_()!=bball2.center_() || bball1.radius()!=bball2.radius();
  }
  
  /// Writes a bounding ball to an ostream.
  template <class BBallT, class RealT, int DimN>
  std::ostream& operator<<( std::ostream& os, BBallBase<BBallT,RealT,DimN> const& bball ) {
    return os << "(" << bball.center_() << "-" << bball.radius() << ")";
  }

  // *******************************************************************
  // class BBall
  // *******************************************************************

  /// A general fixed-dimensional bounding ball class,
  /// represented by a vector pointing to the center, and a radius.
  template <class RealT, int DimN = 0>
  class BBall : public BBallBase<BBall<RealT, DimN>, RealT, DimN> {
  public:

    /// Default constructor. Constructs an empty bounding ball.
    BBall() : BBallBase<BBall<RealT, DimN>, RealT, DimN>() {}

    /// Constructs a bounding ball with the given center and radius.
    template <class VectorT>
    BBall( VectorBase<VectorT> const& center, RealT radius ) :
      BBallBase<BBall<RealT, DimN>, RealT, DimN>( center, radius ) {}

    /// Constructs a 2D bounding ball with the given center point
    /// coordinates and radius.  (Only valid for 2D bouding
    /// balls.)
    BBall( RealT centerx, RealT centery, RealT radius )
      : BBallBase<BBall<RealT, DimN>, RealT, DimN>( centerx, centery, radius ) {
      BOOST_STATIC_ASSERT( DimN==2 );
    }

    /// Constructs a 3D bounding ball with the given center point
    /// coordinates and radius.  (Only valid for 3D bouding
    /// balls.)
    BBall( RealT centerx, RealT centery, RealT centerz, RealT radius )
      : BBallBase<BBall<RealT, DimN>, RealT, DimN>( centerx, centery, centerz, radius ) {
      BOOST_STATIC_ASSERT( DimN==3 );
    }

    /// Copy constructor.
    template <class BBallT1, class RealT1, int DimN1>
    BBall( BBallBase<BBallT1, RealT1, DimN1> const& bball )
      : BBallBase<BBall<RealT, DimN>, RealT, DimN>( bball.center_() , bball.radius() ) {}
    
    /// Copy assignment operator.
    template <class BBallT1, class RealT1, int DimN1>
    BBall& operator=( BBallBase<BBallT1, RealT1, DimN1> const& bball ) {
      this->center_() = bball.center_();
      this->radius() = bball.radius();
      return *this;
    }
  };

  /// A general arbitrary-dimensional bounding ball class,
  /// represented by a vector pointing to the center, and a radius.
  template <class RealT>
  class BBall<RealT, 0> : public BBallBase<BBall<RealT, 0>, RealT, 0> {
  public:

    /// Default constructor. Constructs an empty bounding ball.
    BBall() : BBallBase<BBall<RealT, 0>, RealT, 0>() {}

    /// Constructs a bounding ball with the given center and radius.
    template <class VectorT>
    BBall( VectorBase<VectorT> const& center, RealT radius ) :
      BBallBase<BBall<RealT, 0>, RealT, 0>( center, radius ) {}

    /// Constructs a 2D bounding ball with the given center point
    /// coordinates and radius.  (Only valid for 2D bouding
    /// balls.)
    BBall( RealT centerx, RealT centery, RealT radius )
      : BBallBase<BBall<RealT, 0>, RealT, 0>( centerx, centery, radius ) {}

    /// Constructs a 3D bounding ball with the given center point
    /// coordinates and radius.  (Only valid for 3D bouding
    /// balls.)
    BBall( RealT centerx, RealT centery, RealT centerz, RealT radius )
      : BBallBase<BBall<RealT, 0>, RealT, 0>( centerx, centery, centerz, radius ) {}

    /// Copy constructor.
    template <class BBallT1, class RealT1, int DimN1>
    BBall( BBallBase<BBallT1, RealT1, DimN1> const& bball )
      : BBallBase<BBall<RealT, 0>, RealT, 0>( bball.center_() , bball.radius() ) {}
    
    /// Copy assignment operator.
    template <class BBallT1, class RealT1, int DimN1>
    BBall& operator=( BBallBase<BBallT1, RealT1, DimN1> const& bball ) {
      this->center_() = bball.center_();
      this->radius() = bball.radius();
      return *this;
    }
  };

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
  typedef BBox<float64> BBoxN;
  typedef BBox<float32> BBoxNf;
  typedef BBox<int32> BBoxNi;

  using math::BBall;
  typedef BBall<float64, 2> BBall2;
  typedef BBall<float64, 3> BBall3;
  typedef BBall<float64, 4> BBall4;
  typedef BBall<float32, 2> BBall2f;
  typedef BBall<float32, 3> BBall3f;
  typedef BBall<float32, 4> BBall4f;
  typedef BBall<int32, 2> BBall2i;
  typedef BBall<int32, 3> BBall3i;
  typedef BBall<int32, 4> BBall4i;
  typedef BBall<float64> BBallN;
  typedef BBall<float32> BBallNf;
  typedef BBall<int32> BBallNi;
} // namespace vw

#endif // __VW_MATH__BBOX_H__
