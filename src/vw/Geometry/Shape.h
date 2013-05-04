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


/// \file Shape.h
///
/// Provides a base class for generic shapes.
#ifndef __VW_GEOMETRY_SHAPE_H__
#define __VW_GEOMETRY_SHAPE_H__

#include <vw/Math/Vector.h>

namespace vw {
namespace geometry {

  // *******************************************************************
  // class ShapeBase
  // *******************************************************************

  /// A CRTP base class for general n-dimensional bounding shapes.
  /// Provides a mechanism for restricting function arguments to
  /// shapes, provides general shape operations, and provides the
  /// various arithmetic assignment operators.
  template <class ShapeT, class RealT, int DimN>
  class ShapeBase {
  typedef Vector<RealT, DimN> ShapeVectorT;
  public:

    /// Returns the derived implementation type.
    ShapeT& shape_impl() { return *static_cast<ShapeT*>(this); }

    /// Returns the derived implementation type.
    ShapeT const& shape_impl() const { return *static_cast<ShapeT const*>(this); }
/* These introduce ambiguities in Box, and it's not clear that we really want them....
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
    ShapeVectorT center() const {
      return shape_impl().center();
    }

    /// Returns true if the bounding shape is empty (i.e. degenerate).
    bool empty() const {
      return shape_impl().empty();
    }

    /// Scales the bounding shape relative to the origin.
    template <class ScalarT>
    ShapeT& operator*=( ScalarT s ) {
      return shape_impl() = shape_impl() * s;
    }

    /// Scales the bounding shape relative to the origin.
    template <class ScalarT>
    ShapeT& operator/=( ScalarT s ) {
      return shape_impl() = shape_impl() / s;
    }

    /// Offsets the bounding shape by the given vector.
    template <class VectorT>
    ShapeT& operator+=( VectorBase<VectorT> const& v ) {
      return shape_impl() = shape_impl() + v;
    }

    /// Offsets the bounding shape by the negation of the given vector.
    template <class VectorT>
    ShapeT& operator-=( VectorBase<VectorT> const& v ) {
      return shape_impl() = shape_impl() - v;
    }
*/
  };

  /// Scales a bounding shape relative to the origin.
  template <class ShapeT, class RealT, int DimN, class ScalarT>
  inline ShapeT operator*( ShapeBase<ShapeT, RealT, DimN> const& shape, ScalarT s ) {
    ShapeT result = shape.shape_impl();
    result *= s;
    return result;
  }

  /// Scales a bounding shape relative to the origin.
  template <class ShapeT, class RealT, int DimN, class ScalarT>
  inline ShapeT operator/( ShapeBase<ShapeT, RealT, DimN> const& shape, ScalarT s ) {
    ShapeT result = shape.shape_impl();
    result /= s;
    return result;
  }

  /// Scales a bounding shape relative to the origin.
  template <class ShapeT, class RealT, int DimN, class ScalarT>
  inline ShapeT operator*( ScalarT s, ShapeBase<ShapeT, RealT, DimN> const& shape ) {
    return shape * s;
  }

  /// Offsets a bounding shape by the given vector.
  template <class ShapeT, class RealT, int DimN, class VectorT>
  inline ShapeT operator+( ShapeBase<ShapeT, RealT, DimN> const& shape, VectorBase<VectorT> const& v ) {
    ShapeT result = shape.shape_impl();
    result += v.impl();
    return result;
  }

  /// Offsets a bounding shape by the given vector.
  template <class ShapeT, class RealT, int DimN, class VectorT>
  inline ShapeT operator+( VectorBase<VectorT> const& v, ShapeBase<ShapeT, RealT, DimN> const& shape ) {
    return shape + v;
  }

  /// Offsets a bounding shape by the negation of the given vector.
  template <class ShapeT, class RealT, int DimN, class VectorT>
  inline ShapeT operator-( ShapeBase<ShapeT, RealT, DimN> const& shape, VectorBase<VectorT> const& v ) {
    ShapeT result = shape.shape_impl();
    result -= v.impl();
    return result;
  }

  /// Equality of two bounding shapes.
  template <class ShapeT, class RealT, int DimN>
  inline bool operator==( ShapeBase<ShapeT, RealT, DimN> const& shape1, ShapeBase<ShapeT, RealT, DimN> const& shape2 ) {
    return shape1.shape_impl() == shape2.shape_impl();
  }

  /// Inequality of two bounding shapes.
  template <class ShapeT, class RealT, int DimN>
  inline bool operator!=( ShapeBase<ShapeT, RealT, DimN> const& shape1, ShapeBase<ShapeT, RealT, DimN> const& shape2 ) {
    return shape1.shape_impl() != shape2.shape_impl();
  }

}} // namespace vw::geometry

#endif // __VW_GEOMETRY_SHAPE_H__
