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

/// \file BConvex.h
///
/// Provides a generic convex bounding shape.
#ifndef __VW_MATH__BCONVEX_H__
#define __VW_MATH__BCONVEX_H__

// In addition to the standard incremental grow() method, BConvex provides a
// constructor that takes a std::vector of vw::Vectors. If qhull is available,
// this constructor uses it to construct the initial convex hull, which is much
// better than using the incremental grow() method for a large number of
// points. Otherwise, this constructor is just like calling the incremental
// grow() method for each point yourself. Note that you will probably not
// notice the difference immediately--because of lazy evaluation in the Parma
// Polyhedra Library (PPL), the std::vector constructor will be a bit faster
// without qhull installed, but, if qhull is not installed, the first operation
// that requires the hull can take a very long time while the hull is computed.

#include <iostream>
#include <limits>
#include <math.h>
#include <vector>

#include <boost/static_assert.hpp>

#include <vw/Math/BShape.h>
#include <vw/Math/BBox.h>
#include <vw/Math/Vector.h>
#include <vw/Math/Functors.h>

namespace vw {
namespace math {

  /// \cond INTERNAL
  namespace bconvex_promote {
    struct Promoted {
      //NOTE: we use these types because they are the longest input types supported by GMP
      typedef long int SignedT;
      typedef long unsigned int UnsignedT;
      typedef double FloatT;
      union Val {
        SignedT s;
        UnsignedT u;
        FloatT f;
      };
      Promoted() : is_signed( false ), is_integral( true ) {
        val.u = 0;
      }
      Val val;
      bool is_signed;
      bool is_integral;
    };

    std::ostream& operator<<( std::ostream& os, Promoted const& r );

    template <class ValT, bool IntegralN=true, bool SignedN=true>
    struct PromoteBehavior {
      static inline void set( ValT val, Promoted &r ) {
        r.is_signed = true;
        r.is_integral = true;
        r.val.s = (Promoted::SignedT)val;
      }
      static inline ValT get( Promoted const& r ) {
        return (ValT)(r.val.s);
      }
    };

    template <class ValT>
    struct PromoteBehavior<ValT,true,false> {
      static inline void set( ValT val, Promoted &r ) {
        r.is_signed = false;
        r.is_integral = true;
        r.val.u = (Promoted::UnsignedT)val;
      }
      static inline ValT get( Promoted const& r ) {
        return (ValT)(r.val.u);
      }
    };

    template <class ValT>
    struct PromoteBehavior<ValT,false,true> {
      static inline void set( ValT val, Promoted &r ) {
        r.is_signed = true;
        r.is_integral = false;
        r.val.f = (Promoted::FloatT)val;
      }
      static inline ValT get( Promoted const& r ) {
        return (ValT)(r.val.f);
      }
    };

    template <class ValT>
    struct PromoteBehavior<ValT,false,false> {
      static inline void set( ValT val, Promoted &r ) {
        r.is_signed = false;
        r.is_integral = false;
        r.val.f = (Promoted::FloatT)val;
      }
      static inline ValT get( Promoted const& r ) {
        return (ValT)(r.val.f);
      }
    };

    template <class ValT>
    struct PromoteFuncs {
      static void set( ValT val, Promoted &r ) {
        // Make sure we have a type for which we know limits
        BOOST_STATIC_ASSERT(std::numeric_limits<ValT>::is_specialized);
        PromoteBehavior<ValT,std::numeric_limits<ValT>::is_integer,std::numeric_limits<ValT>::is_signed>::set(val, r);
      }
      static ValT get( Promoted const& r ) {
        // Make sure we have a type for which we know limits
        BOOST_STATIC_ASSERT(std::numeric_limits<ValT>::is_specialized);
        return PromoteBehavior<ValT,std::numeric_limits<ValT>::is_integer,std::numeric_limits<ValT>::is_signed>::get(r);
      }
    };
    
    template <class ValT, bool IntegralN=true>
    struct EpsilonBehavior {
      static ValT epsilon() {
        return 1;
      }
    };
    
    template <class ValT>
    struct EpsilonBehavior<ValT,false> {
      static ValT epsilon() {
        return std::numeric_limits<ValT>::epsilon();
      }
    };
    
    template <class ValT>
    struct EpsilonFunctor {
      static ValT epsilon() {
        // Make sure we have a type for which we know limits
        BOOST_STATIC_ASSERT(std::numeric_limits<ValT>::is_specialized);
        return EpsilonBehavior<ValT,std::numeric_limits<ValT>::is_integer>::epsilon();
      }
    };

  } // namespace bconvex_promote
  /// \endcond
  
  // *******************************************************************
  // class BConvex
  // *******************************************************************
  
  /// A general arbitrary-dimensional convex shape class.
  class BConvex : public BShapeBase<BConvex, double, 0> {
  public:
    /// Default constructor.
    BConvex() {
      m_poly = 0;
    }

    /// Polyhedron-only constructor.
    BConvex( const void *poly ) {
      m_poly = new_poly(poly);
    }
    
    /// Vector-of-points constructor.
    template <class ContainerT>
    BConvex( std::vector<ContainerT> const& points ) {
      unsigned num_points = points.size();
      VW_ASSERT( num_points > 0, ArgumentErr() << "No points provided to BConvex vector-of-points constructor!" );
      unsigned dim = points[0].size();
      if (have_qhull()) {
        double *p = new double[num_points*dim];
        for (unsigned i = 0, k = 0; i < num_points; i++)
          for (unsigned int j = 0; j < dim; j++, k++)
            p[k] = points[i][j];
        init_with_qhull(dim, num_points, p);
        delete[] p;
      }
      else {
        m_poly = new_poly(dim);
        for (unsigned i = 0; i < num_points; i++)
          grow(points[i]);
      }
    }

    /// Destructor.
    ~BConvex() {
      delete_poly(m_poly);
      m_poly = 0;
    }
    
    /// Copy constructor.
    BConvex( BConvex const& bconv ) {
      m_poly = new_poly(bconv.poly());
    }
    
    /// Copy assignment operator.
    BConvex& operator=( BConvex const& bconv ) {
      copy_poly(bconv.poly(), m_poly);
      return *this;
    }
    
    /// Returns the internal polyhedron.
    void *poly() {return m_poly;}
    
    /// Returns the internal polyhedron.
    const void *poly() const {return m_poly;}
    
    /// Grows a convex shape to include the given point.
    template <class VectorT>
    void grow( VectorBase<VectorT> const& point ) {
      Vector<bconvex_promote::Promoted> p;
      new_poly_if_needed(point.impl().size(), m_poly);
      grow_(promote_vector(point, p));
    }
    
    /// Grows a convex shape to include the given convex shape.
    void grow( BConvex const& bconv );
    
    /// Crops (intersects) this convex shape to the given convex shape.
    void crop( BConvex const& bconv );
    
    /// Expands this convex shape by the given offset in every direction.
    void expand( double offset );

    /// Contracts this convex shape by the given offset in every direction.
    void contract( double offset ) {
      expand(-offset);
    }
    
    /// Returns true if the given point is contained in the convex shape.
    template <class VectorT>
    bool contains( VectorBase<VectorT> const& point ) const {
      Vector<bconvex_promote::Promoted> p;
      if (!m_poly)
        return false;
      return contains_(promote_vector(point, p));
    }
    
    /// Returns true if the given convex shape is entirely contained
    /// in this convex shape.
    bool contains( BConvex const& bconv ) const;
    
    /// Returns true if the given convex shape intersects this
    /// convex shape.
    bool intersects( BConvex const& bconv ) const;
    
    /// Returns the size (i.e. twice the largest distance between
    /// the center and a vertex) of the convex shape.
    double size() const;

    /// Returns the center point of the convex shape.
    Vector<double> center() const;

    /// Returns true if the convex shape is empty (i.e. degenerate).
    bool empty() const;
    
    /// Returns the number of facets.
    unsigned num_facets() const;
    
    /// Returns the number of vertices.
    unsigned num_vertices() const;
    
    /// Returns a bounding box that contains this convex shape.
    BBoxN bounding_box() const;

    /// Returns true if the given convex shape is equal to this
    /// convex shape.
    bool equal( BConvex const& bconv ) const;

    /// Prints the convex shape.
    void print( std::ostream& os = std::cout ) const;

    /// Prints the convex shape in vrml.
    void write_vrml( std::ostream& os = std::cout ) const;

    /// Prints the convex shape in vrml.
    void write_vrml( const char *fn ) const;

    /// Prints the convex shape in oogl.
    void write_oogl( std::ostream& os = std::cout ) const;

    /// Prints the convex shape in oogl.
    void write_oogl( const char *fn ) const;

    /// Scales the convex shape relative to the origin.
    template <class ScalarT>
    BConvex& operator*=( ScalarT s ) {
      using namespace bconvex_promote;
      VW_ASSERT( !empty(), LogicErr() << "Cannot multiply an empty polyhedron by a scalar!" );
      Promoted r;
      operator_mult_eq_(promote_scalar(s, r));
      return *this;
    }

    /// Scales the convex shape relative to the origin.
    template <class ScalarT>
    BConvex& operator/=( ScalarT s ) {
      using namespace bconvex_promote;
      VW_ASSERT( !empty(), LogicErr() << "Cannot divide an empty polyhedron by a scalar!" );
      Promoted r;
      operator_div_eq_(promote_scalar(s, r));
      return *this;
    }

    /// Offsets the convex shape by the given vector.
    template <class VectorT>
    BConvex& operator+=( VectorBase<VectorT> const& v ) {
      VW_ASSERT( !empty(), LogicErr() << "Cannot add a vector to an empty polyhedron!" );
      Vector<bconvex_promote::Promoted> p;
      operator_plus_eq_(promote_vector(v, p));
      return *this;
    }

    /// Offsets the convex shape by the negation of the given vector.
    template <class VectorT>
    BConvex& operator-=( VectorBase<VectorT> const& v ) {
      VW_ASSERT( !empty(), LogicErr() << "Cannot subtract a vector from an empty polyhedron!" );
      Vector<bconvex_promote::Promoted> p;
      operator_minus_eq_(promote_vector(v, p));
      return *this;
    }

    /// Converts a Vector to a Promoted Vector.
    template <class VectorT>
    static Vector<bconvex_promote::Promoted> const& promote_vector( VectorBase<VectorT> const& point, Vector<bconvex_promote::Promoted> &p ) {
      using namespace bconvex_promote;
      unsigned dim = point.impl().size();
      p.set_size(dim);
      for (unsigned i = 0; i < dim; i++)
        PromoteFuncs<typename VectorT::value_type>::set(point.impl()(i), p(i));
      return p;
    }
    
    /// Converts a scalar to a Promoted.
    template <class ScalarT>
    static inline bconvex_promote::Promoted const& promote_scalar( ScalarT s, bconvex_promote::Promoted &r ) {
      using namespace bconvex_promote;
      PromoteFuncs<ScalarT>::set(s, r);
      return r;
    }

  protected:
    /// Creates a polyhedron with the given dimension.
    static void *new_poly( unsigned dim );

    /// Creates a polyhedron that is a copy of the given polyhedron.
    static void *new_poly( const void *poly );

    /// Creates a polyhedron with the given dimension, if it has
    /// not already been created.
    static void new_poly_if_needed( unsigned dim, void *&poly );
    
    /// Deletes the given polyhedron.
    static void delete_poly( void *poly );
    
    /// Copies polyhedron from_poly to polyhedron to_poly. 
    static void copy_poly( const void *from_poly, void *&to_poly );
    
    /// Returns whether qhull is available.
    static bool have_qhull();
    
    /// Initialize polyhedron using qhull.
    void init_with_qhull( unsigned dim, unsigned num_points, double *p );

    /// Grows a convex shape to include the given point.
    void grow_( Vector<bconvex_promote::Promoted> const& point );

    /// Returns true if the given point is contained in the convex shape.
    bool contains_( Vector<bconvex_promote::Promoted> const& point ) const;

    /// Scales the convex shape relative to the origin.
    void operator_mult_eq_( bconvex_promote::Promoted const& s );

    /// Scales the convex shape relative to the origin.
    void operator_div_eq_( bconvex_promote::Promoted const& s );

    /// Offsets the convex shape by the given vector.
    void operator_plus_eq_( Vector<bconvex_promote::Promoted> const& v );

    /// Offsets the convex shape by the negation of the given vector.
    void operator_minus_eq_( Vector<bconvex_promote::Promoted> const& v );
    
    void *m_poly;
  };
  
  /// Scales a convex shape relative to the origin.
  template <class ScalarT>
  inline BConvex operator*( BConvex const& bconv, ScalarT s ) {
    BConvex result = bconv;
    result *= s;
    return result;
  }

  /// Scales a convex shape relative to the origin.
  template <class ScalarT>
  inline BConvex operator/( BConvex const& bconv, ScalarT s ) {
    BConvex result = bconv;
    result /= s;
    return result;
  }

  /// Scales a convex shape relative to the origin.
  template <class ScalarT>
  inline BConvex operator*( ScalarT s, BConvex const& bconv ) {
    return bconv * s;
  }
  
  /// Offsets a convex shape by the given vector.
  template <class VectorT>
  inline BConvex operator+( BConvex const& bconv, VectorBase<VectorT> const& v ) {
    BConvex result = bconv;
    result += v.impl();
    return result;
  }

  /// Offsets a convex shape by the given vector.
  template <class VectorT>
  inline BConvex operator+( VectorBase<VectorT> const& v, BConvex const& bconv ) {
    return bconv + v;
  }

  /// Offsets a convex shape by the negation of the given vector.
  template <class VectorT>
  inline BConvex operator-( BConvex const& bconv, VectorBase<VectorT> const& v ) {
    BConvex result = bconv;
    result -= v.impl();
    return result;
  }

  /// Equality of two convex shapes.
  inline bool operator==( BConvex const& bconv1, BConvex const& bconv2 ) {
    return bconv1.equal(bconv2);
  }

  /// Inequality of two convex shapes.
  inline bool operator!=( BConvex const& bconv1, BConvex const& bconv2 ) {
    return !bconv1.equal(bconv2);
  }
  
  /// Writes a convex shape to an ostream.
  std::ostream& operator<<( std::ostream& os, BConvex const& bconv );

} // namespace math

  // Include BConvex in vw namespace.
  using math::BConvex;
} // namespace vw

#endif // __VW_MATH__BCONVEX_H__
