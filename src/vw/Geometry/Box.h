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


/// \file Box.h
///
/// A box shape class.
#ifndef __VW_GEOMETRY_BOX_H__
#define __VW_GEOMETRY_BOX_H__

#include <iostream>
#include <fstream>
#include <limits>
#include <vector>
#include <cmath>

#include <boost/static_assert.hpp>

#include <vw/Math/Vector.h>
#include <vw/Math/BBox.h>
#include <vw/Geometry/Shape.h>

namespace vw {
namespace geometry {

  // *******************************************************************
  // class Box
  // *******************************************************************

  /// A box shape class, which derives almost all of its implementation
  /// from the math::BBox superclass.
  template <class RealT, int DimN = 0>
  class Box : public ShapeBase<Box<RealT, DimN>, RealT, DimN>,
              public BBox<RealT,DimN> {
  public:

    /// Default constructor.  Constructs the ultimate empty box,
    /// whose limits are at the opposite corners of the underlying
    /// numeric space.
    Box() {}

    /// Constructs a box with the given minimal and maximal points.
    Box( Vector<RealT, DimN> const& min, Vector<RealT, DimN> const& max )
      : BBox<RealT, DimN>( min, max ) {}

    /// Constructs a 2D box with the given minimal point coordinates
    /// and dimensions.  (Only valid for 2D boxes.)
    Box( RealT minx, RealT miny, RealT width, RealT height )
      : BBox<RealT, DimN>( Vector<RealT,2>(minx,miny),
                           Vector<RealT,2>(minx+width,miny+height) )
    {
      BOOST_STATIC_ASSERT( DimN==2 );
    }

    /// Constructs a 3D box with the given minimal point coordinates
    /// and dimensions.  (Only valid for 3D boxes.)
    Box( RealT minx, RealT miny, RealT minz, RealT width, RealT height, RealT depth )
      : BBox<RealT, DimN>( Vector<RealT,3>(minx,miny,minz),
                           Vector<RealT,3>(minx+width,miny+height,minz+depth) )
    {
      BOOST_STATIC_ASSERT( DimN==3 );
    }

    /// Generalized copy constructor.
    template <class RealT1, int DimN1>
    Box( BBox<RealT1, DimN1> const& bbox )
      : BBox<RealT, DimN>( bbox.min(), bbox.max() ) {}

    /// Generalized copy assignment operator.
    template <class RealT1, int DimN1>
    Box& operator=( BBox<RealT1, DimN1> const& bbox ) {
      this->min() = bbox.min();
      this->max() = bbox.max();
      return *this;
    }

    /// We are our own bounding box, so we just return our superclass.
    BBox<RealT,DimN> const& bbox() const {
      return *this;
    }

    /// Writes the box.
    void write( std::ostream& os = std::cout, bool binary = false ) const {
      std::vector<Vector<RealT, DimN> > points;
      points.push_back(this->min());
      points.push_back(this->max());
      write_point_list(os, points, binary);
    }

    /// Writes the box.
    void write( const char *fn, bool binary = false ) const {
      if (binary) {
        std::ofstream of( fn, std::ofstream::out | std::ofstream::binary );
        write(of, binary);
        of.close();
      }
      else {
        std::ofstream of( fn );
        write(of, binary);
        of.close();
      }
    }

  };

  /// Scales a box relative to the origin.
  template <class RealT, int DimN, class ScalarT>
  inline Box<RealT, DimN> operator*( Box<RealT, DimN> const& box, ScalarT s ) {
    return box.bbox() * s;
  }

  /// Scales a box relative to the origin.
  template <class RealT, int DimN, class ScalarT>
  inline Box<RealT, DimN> operator/( Box<RealT, DimN> const& box, ScalarT s ) {
    return box.bbox() / s;
  }

  /// Scales a box relative to the origin.
  template <class RealT, int DimN, class ScalarT>
  inline Box<RealT, DimN> operator*( ScalarT s, Box<RealT, DimN> const& box ) {
    return s * box.bbox();
  }

  /// Translates a box by the given vector.
  template <class RealT, int DimN, class VectorT>
  inline Box<RealT, DimN> operator+( Box<RealT, DimN> const& box, VectorBase<VectorT> const& v ) {
    return box.bbox() + v;
  }

  /// Translates a box by the given vector.
  template <class RealT, int DimN, class VectorT>
  inline Box<RealT, DimN> operator+( VectorBase<VectorT> const& v, Box<RealT, DimN> const& box ) {
    return v + box.bbox();
  }

  /// Translates a box by the negation of the given vector.
  template <class RealT, int DimN, class VectorT>
  inline Box<RealT, DimN> operator-( Box<RealT, DimN> const& box, VectorBase<VectorT> const& v ) {
    return box.bbox() - v;
  }

  /// Writes a box to a file.
  template <class RealT, int DimN>
  void write_box( std::string const& filename, Box<RealT,DimN> const& box, bool binary = false ) {
    box.write(filename.c_str(), binary);
  }

  /// Reads a box from a file.
  template <class RealT, int DimN>
  void read_box( std::string const& filename, Box<RealT,DimN>& box, bool binary = false ) {
    std::vector<Vector<RealT,DimN> > points;
    read_point_list(filename, points, binary);
    VW_ASSERT( points.size() == 2, IOErr() << "Incorrect number of points in box file!" );
    box = Box<RealT,DimN>(points[0], points[1]);
  }

} // namespace geometry

  // Convenience typedefs
  using geometry::Box;
  typedef Box<float64, 2> Box2;
  typedef Box<float64, 3> Box3;
  typedef Box<float64, 4> Box4;
  typedef Box<float64> BoxN;
  typedef Box<float32, 2> Box2f;
  typedef Box<float32, 3> Box3f;
  typedef Box<float32, 4> Box4f;
  typedef Box<float32> BoxNf;

} // namespace vw

#endif // __VW_GEOMETRY_BOX_H__
