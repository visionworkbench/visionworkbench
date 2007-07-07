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

/// \file BShape.h
///
/// Provides a base for generic bounding shapes.
#ifndef __VW_MATH__BSHAPE_H__
#define __VW_MATH__BSHAPE_H__

#include <iostream>

#include <vw/Math/Vector.h>

namespace vw {
namespace math {

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

}} // namespace vw::math

#endif // __VW_MATH__BSHAPE_H__
