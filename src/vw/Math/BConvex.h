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
  namespace bconvex_rational {
    struct Rational {
      //NOTE: we use these types because they are the longest input types supported by GMP
      typedef long int SignedT;
      typedef long unsigned int UnsignedT;
      Rational() : signed_num( 0 ), unsigned_num( 0 ), den( 1 ), is_signed( false ) {}
      SignedT signed_num;
      UnsignedT unsigned_num;
      UnsignedT den;
      bool is_signed;
    };

    std::ostream& operator<<( std::ostream& os, Rational const& r );

    template <class ValT, bool IntegralN=true, bool SignedN=true>
    struct RationalBehavior {
      static inline void set( ValT val, ValT max_abs, Rational &r ) {
        ArgAbsFunctor abs_func;
        if (abs_func(val) > max_abs)
          std::cout << val << "->" << abs_func(val) << " > " << max_abs << std::endl;
        VW_ASSERT( abs_func(val) <= max_abs, ArgumentErr() << "Trying to convert val to rational where abs(val) > max_abs." );
        r.is_signed = true;
        r.den = 1;
        r.signed_num = (Rational::SignedT)val;
      }
      static inline ValT get( Rational const& r ) {
        return (ValT)(r.signed_num);
      }
    };

    template <class ValT>
    struct RationalBehavior<ValT,true,false> {
      static inline void set( ValT val, ValT max_abs, Rational &r ) {
        VW_ASSERT( val <= max_abs, ArgumentErr() << "Trying to convert val to rational where abs(val) > max_abs." );
        r.is_signed = false;
        r.den = 1;
        r.unsigned_num = (Rational::UnsignedT)val;
      }
      static inline ValT get( Rational const& r ) {
        return (ValT)(r.unsigned_num);
      }
    };

    template <class ValT>
    struct RationalBehavior<ValT,false,true> {
      static inline void set( ValT val, ValT max_abs, Rational &r ) {
        ArgAbsFunctor abs_func;
        VW_ASSERT( abs_func(val) <= max_abs, ArgumentErr() << "Trying to convert val to rational where abs(val) > max_abs." );
        r.is_signed = true;
        ArgCeilFunctor ceil_func;
        Rational::SignedT den = std::numeric_limits<Rational::SignedT>::max() / (Rational::SignedT)ceil_func(max_abs);
        r.den = (Rational::UnsignedT)den;
        r.signed_num = (Rational::SignedT)(val * den);
      }
      static inline ValT get( Rational const& r ) {
        return (ValT)(r.signed_num) / (ValT)(r.den);
      }
    };

    template <class ValT>
    struct RationalBehavior<ValT,false,false> {
      static inline void set( ValT val, ValT max_abs, Rational &r ) {
        VW_ASSERT( val <= max_abs, ArgumentErr() << "Trying to convert val to rational where abs(val) > max_abs." );
        r.is_signed = false;
        ArgCeilFunctor ceil_func;
        r.den = std::numeric_limits<Rational::UnsignedT>::max() / (Rational::UnsignedT)ceil_func(max_abs);
        r.unsigned_num = (Rational::UnsignedT)(val * r.den);
      }
      static inline ValT get( Rational const& r ) {
        return (ValT)(r.unsigned_num) / (ValT)(r.den);
      }
    };

    template <class ValT>
    struct RationalFuncs {
      static void set( ValT val, ValT max_abs, Rational &r ) {
        // Make sure we have a type for which we know limits
        BOOST_STATIC_ASSERT(std::numeric_limits<ValT>::is_specialized);
        RationalBehavior<ValT,std::numeric_limits<ValT>::is_integer,std::numeric_limits<ValT>::is_signed>::set(val, max_abs, r);
      }
      static ValT get( Rational const& r ) {
        // Make sure we have a type for which we know limits
        BOOST_STATIC_ASSERT(std::numeric_limits<ValT>::is_specialized);
        return RationalBehavior<ValT,std::numeric_limits<ValT>::is_integer,std::numeric_limits<ValT>::is_signed>::get(r);
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

  } // namespace bconvex_rational
  /// \endcond
  
  // *******************************************************************
  // class BConvex
  // *******************************************************************
  
  /// A general arbitrary-dimensional convex shape class.
  class BConvex : public BShapeBase<BConvex, double, 0> {
  public:
    /// Dimension-only constructor.
    BConvex( unsigned dim ) {
      m_poly = new_poly(dim);
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
      Vector<bconvex_rational::Rational> p;
      grow_(convert_vector(point, p));
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
      Vector<bconvex_rational::Rational> p;
      return contains_(convert_vector(point, p));
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
    
    /// Returns a bounding box that contains this convex shape.
    BBoxN bounding_box() const;

    /// Returns true if the given convex shape is equal to this
    /// convex shape.
    bool equal( BConvex const& bconv ) const;

    /// Prints the convex shape.
    void print( std::ostream& os = std::cout ) const;

    /// Scales the convex shape relative to the origin.
    template <class ScalarT>
    BConvex& operator*=( ScalarT s ) {
      using namespace bconvex_rational;
      Rational r;
      operator_mult_eq_(convert_scalar(s, r));
      return *this;
    }

    /// Scales the convex shape relative to the origin.
    template <class ScalarT>
    BConvex& operator/=( ScalarT s ) {
      using namespace bconvex_rational;
      Rational r;
      operator_div_eq_(convert_scalar(s, r));
      return *this;
    }

    /// Offsets the convex shape by the given vector.
    template <class VectorT>
    BConvex& operator+=( VectorBase<VectorT> const& v ) {
      Vector<bconvex_rational::Rational> p;
      operator_plus_eq_(convert_vector(v, p));
      return *this;
    }

    /// Offsets the convex shape by the negation of the given vector.
    template <class VectorT>
    BConvex& operator-=( VectorBase<VectorT> const& v ) {
      Vector<bconvex_rational::Rational> p;
      operator_minus_eq_(convert_vector(v, p));
      return *this;
    }

    /// Converts a Vector to a Rational Vector.
    template <class VectorT>
    static Vector<bconvex_rational::Rational> const& convert_vector( VectorBase<VectorT> const& point, Vector<bconvex_rational::Rational> &p ) {
      using namespace bconvex_rational;
      unsigned dim = point.impl().size();
      p.set_size(dim);
      typename VectorT::value_type max_abs = norm_inf(point);
      typename VectorT::value_type eps = EpsilonFunctor<typename VectorT::value_type>::epsilon();
      max_abs = std::max(max_abs, eps);
      for (unsigned i = 0; i < dim; i++)
        RationalFuncs<typename VectorT::value_type>::set(point.impl()(i), max_abs, p(i));
      return p;
    }
    
    /// Converts a scalar to a Rational.
    template <class ScalarT>
    static inline bconvex_rational::Rational const& convert_scalar( ScalarT s, bconvex_rational::Rational &r ) {
      using namespace bconvex_rational;
      RationalFuncs<ScalarT>::set(s, s, r);
      return r;
    }

  protected:
    /// Creates a polyhedron with the given dimension.
    static void *new_poly( unsigned dim );

    /// Creates a polyhedron that is a copy of the given polyhedron.
    static void *new_poly( const void *poly );
    
    /// Deletes the given polyhedron.
    static void delete_poly( void *poly );
    
    /// Copies polyhedron from_poly to polyhedron to_poly. 
    static void copy_poly( const void *from_poly, void *to_poly );
    
    /// Returns whether qhull is available.
    static bool have_qhull();
    
    /// Initialize polyhedron using qhull.
    void init_with_qhull( unsigned dim, unsigned num_points, double *p );

    /// Grows a convex shape to include the given point.
    void grow_( Vector<bconvex_rational::Rational> const& point );

    /// Returns true if the given point is contained in the convex shape.
    bool contains_( Vector<bconvex_rational::Rational> const& point ) const;

    /// Scales the convex shape relative to the origin.
    void operator_mult_eq_( bconvex_rational::Rational const& s );

    /// Scales the convex shape relative to the origin.
    void operator_div_eq_( bconvex_rational::Rational const& s );

    /// Offsets the convex shape by the given vector.
    void operator_plus_eq_( Vector<bconvex_rational::Rational> const& v );

    /// Offsets the convex shape by the negation of the given vector.
    void operator_minus_eq_( Vector<bconvex_rational::Rational> const& v );
    
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
