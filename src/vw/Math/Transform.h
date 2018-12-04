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


// \file Transform.h
//
// Transform classes used other places in VW.
//
#ifndef __VW_MATH_TRANSFORM_H__
#define __VW_MATH_TRANSFORM_H__

#include <vw/Core/Features.h>
#include <vw/Core/Log.h>
#include <vw/Math/Vector.h>
#include <vw/Math/Matrix.h>
#include <vw/Math/BBox.h>
#include <vw/Math/LevenbergMarquardt.h>

namespace vw {

  enum FunctionType {
    ConvexFunction,
    ContinuousFunction,
    DiscontinuousFunction
  };


  // The abstract base class of all transform functors.  You should
  // not subclass Transform directly, but rather indirectly via
  // TransformBase or TransformHelper.  Similarly, you should not use
  // Transform pointers or references directly, but rather indirectly via TransformRef.
  class Transform {
    double m_tolerance;

    class ForwardLMA : public math::LeastSquaresModelBase<ForwardLMA> {
      const Transform* m_tx;
    public:
      ForwardLMA( const Transform* tx ) : m_tx(tx) {}
      typedef Vector2   result_type;
      typedef Vector2   domain_type;
      typedef Matrix2x2 jacobian_type;

      inline result_type operator()( domain_type const& x ) const {
        return m_tx->reverse( x );
      }
    };

    class ReverseLMA : public math::LeastSquaresModelBase<ReverseLMA> {
      const Transform* m_tx;
    public:
      ReverseLMA( const Transform* tx ) : m_tx(tx) {}
      typedef Vector2   result_type;
      typedef Vector2   domain_type;
      typedef Matrix2x2 jacobian_type;

      inline result_type operator()( domain_type const& x ) const {
        return m_tx->forward( x );
      }
    };

  public:
    Transform() : m_tolerance(0.0) {}

    virtual ~Transform() {}

    // 'Forward' defines the transformation from coordinates in the source
    // to coordinates in the destination.  That routine is
    // not actually needed to perform  the transformation itself, but
    // it can be used to determine the appropriate dimensions for the output.
    //
    // 'Reverse' defines the transformation from coordinates in our
    // target back to coordinates in the original.
    //
    // The default implementations here actually performs a least
    // squares search using the reverse method. This means we get
    // access to both methods even when a function only implements
    // one. It may not be accurate however.
    virtual Vector2 forward( Vector2 const& point ) const {
      int status;
      return levenberg_marquardt( ForwardLMA( this ), point, point, status );
    }
    virtual Vector2 reverse( Vector2 const& point ) const {
      int status;
      return levenberg_marquardt( ReverseLMA( this ), point, point, status );
    }

    /// Specifies the properties of the forward mapping function.
    virtual FunctionType forward_type() const { return DiscontinuousFunction; }

    /// Specifies the properties of the reverse mapping function.
    virtual FunctionType reverse_type() const { return DiscontinuousFunction; }

    /// This applies the forward transformation to an entire bounding box of pixels.
    virtual BBox2i forward_bbox( BBox2i const& /*output_bbox*/ ) const { 
      vw_throw( NoImplErr() << "forward_bbox() is not implemented for this transform." ); 
      return BBox2i(); 
    }

    /// This applies the reverse transformation to an entire bounding box of pixels.
    virtual BBox2i reverse_bbox( BBox2i const& /*input_bbox*/ ) const { 
      vw_throw( NoImplErr() << "reverse_bbox() is not implemented for this transform." ); 
      return BBox2i(); 
    }

    // TODO: Do we use this for anything???

    /// Set the tolerance (in pixels) to which this Transform is willing to be approximated.
    virtual void set_tolerance( double tolerance ) { m_tolerance = tolerance; }

    /// Get the tolerance (in pixels) to which this Transform is willing to be approximated.
    virtual double tolerance() const { return m_tolerance; }

  }; // End class Transform


  // The CRTP base class for all transform functors.
  //
  // Specific transform classes should derive from this class (or
  // equivalently TransformHelper, below) rather than the abstract
  // Transform class.  At this point the primary purposes of this
  // class are to aid in disambiguating function overloads and to
  // implement bounding box transformations, taking into account the
  // derived type's forward and reverse transform type traits.
  template <class ImplT>
  class TransformBase : public Transform {
  public:
    inline ImplT      & impl()       { return static_cast<ImplT      &>(*this); }
    inline ImplT const& impl() const { return static_cast<ImplT const&>(*this); }

    virtual BBox2i forward_bbox( BBox2i const& bbox ) const {
      ImplT const& txform = impl();
      BBox2 transformed_bbox;
      if (bbox.empty()) 
        return transformed_bbox; // bugfix
      switch( txform.forward_type() ) {
      case ConvexFunction:
        transformed_bbox.grow( txform.forward( Vector2(bbox.min().x(),  bbox.min().y()  ) ) ); // Top left
        transformed_bbox.grow( txform.forward( Vector2(bbox.max().x()-1,bbox.min().y()  ) ) ); // Top right
        transformed_bbox.grow( txform.forward( Vector2(bbox.min().x(),  bbox.max().y()-1) ) ); // Bottom left
        transformed_bbox.grow( txform.forward( Vector2(bbox.max().x()-1,bbox.max().y()-1) ) ); // Bottom right
        break;
      case ContinuousFunction:
        for( int32 x=bbox.min().x(); x<bbox.max().x(); ++x ) { // Top and bottom
          transformed_bbox.grow( txform.forward( Vector2(x,bbox.min().y()) ) );
          transformed_bbox.grow( txform.forward( Vector2(x,bbox.max().y()-1) ) );
        }
        for( int32 y=bbox.min().y()+1; y<bbox.max().y()-1; ++y ) { // Left and right
          transformed_bbox.grow( txform.forward( Vector2(bbox.min().x(),y) ) );
          transformed_bbox.grow( txform.forward( Vector2(bbox.max().x()-1,y) ) );
        }
        break;
      case DiscontinuousFunction:
        for( int32 y=bbox.min().y(); y<bbox.max().y(); ++y )
          for( int32 x=bbox.min().x(); x<bbox.max().x(); ++x )
            transformed_bbox.grow( txform.forward( Vector2(x,y) ) );
        break;
      }
      return grow_bbox_to_int( transformed_bbox );
    }

    virtual BBox2i reverse_bbox( BBox2i const& bbox ) const {
      ImplT const& txform = impl();
      BBox2 transformed_bbox;
      if (bbox.empty()) 
        return transformed_bbox; // bugfix
      switch( txform.reverse_type() ) {
      case ConvexFunction:
        transformed_bbox.grow( txform.reverse( Vector2(bbox.min().x(),  bbox.min().y()  ) ) ); // Top left
        transformed_bbox.grow( txform.reverse( Vector2(bbox.max().x()-1,bbox.min().y()  ) ) ); // Top right
        transformed_bbox.grow( txform.reverse( Vector2(bbox.min().x(),  bbox.max().y()-1) ) ); // Bottom left
        transformed_bbox.grow( txform.reverse( Vector2(bbox.max().x()-1,bbox.max().y()-1) ) ); // Bottom right
        break;
      case ContinuousFunction:
        for( int32 x=bbox.min().x(); x<bbox.max().x(); ++x ) { // Top and bottom
          transformed_bbox.grow( txform.reverse( Vector2(x,bbox.min().y()) ) );
          transformed_bbox.grow( txform.reverse( Vector2(x,bbox.max().y()-1) ) );
        }
        for( int32 y=bbox.min().y()+1; y<bbox.max().y()-1; ++y ) { // Left and right
          transformed_bbox.grow( txform.reverse( Vector2(bbox.min().x(),y) ) );
          transformed_bbox.grow( txform.reverse( Vector2(bbox.max().x()-1,y) ) );
        }
        break;
      case DiscontinuousFunction:
        for( int32 y=bbox.min().y(); y<bbox.max().y(); ++y )
          for( int32 x=bbox.min().x(); x<bbox.max().x(); ++x )
            transformed_bbox.grow( txform.reverse( Vector2(x,y) ) );
        break;
      }
      return grow_bbox_to_int( transformed_bbox );
    }

    // This function is deprecated, and provided for backwards
    // compatibility only.  Use reverse(BBox2i) instead.
    BBox2i compute_input_bbox( BBox2i const& output_bbox ) const VW_DEPRECATED {
      return reverse_bbox( output_bbox );
    }
  }; // End class TransformBase


  // A helper CRTP base class for transform functors.
  //
  // Deriving from this base class is equivalent to deriving
  // from TransformBase, except it allows you to conveniently
  // specify the forward and reverse function type at the
  // same time.
  template <class ImplT, int ForwardType=DiscontinuousFunction, int ReverseType=DiscontinuousFunction>
  class TransformHelper : public TransformBase<ImplT> {
  public:
    virtual FunctionType forward_type() const { return (FunctionType)ForwardType; }
    virtual FunctionType reverse_type() const { return (FunctionType)ReverseType; }
  };

  // Identity transform functor
  //
  // Does nothing!
  class IdentityTransform : public TransformHelper<IdentityTransform,ConvexFunction,ConvexFunction> {
  public:
    IdentityTransform(){}

    inline Vector2 reverse( Vector2 const& p ) const { return p; }
    inline Vector2 forward( Vector2 const& p ) const { return p; }
  };

  // Resample transform functor
  //
  // Transform points by applying a scaling in x and y.
  class ResampleTransform : public TransformHelper<ResampleTransform,ConvexFunction,ConvexFunction> {
    double m_xfactor, m_yfactor;
  public:
    ResampleTransform( double x_scaling, double y_scaling ) :
      m_xfactor( x_scaling ) , m_yfactor( y_scaling ) {}

    template <class VectorT>
    ResampleTransform( VectorBase<VectorT> const& v ) {
      VW_ASSERT( v.impl().size() == 2,
                 ArgumentErr() << "Vector must have 2 dimensions" );
      m_xfactor = v.impl()[0];
      m_yfactor = v.impl()[1];
    }

    inline Vector2 reverse( Vector2 const& p ) const {
      return Vector2( p(0) / m_xfactor, p(1) / m_yfactor );
    }

    inline Vector2 forward( Vector2 const& p ) const {
      return Vector2( p(0) * m_xfactor, p(1) * m_yfactor );
    }
  };


  // Translate transform functor
  class TranslateTransform : public TransformHelper<TranslateTransform,ConvexFunction,ConvexFunction> {
    double m_xtrans, m_ytrans;
  public:
    TranslateTransform(double x_translation, double y_translation) :
      m_xtrans( x_translation ) , m_ytrans( y_translation ) {}

    template <class VectorT>
    TranslateTransform( VectorBase<VectorT> const& v ) {
      VW_ASSERT( v.impl().size() == 2,
                 ArgumentErr() << "Vector must have 2 dimensions" );
      m_xtrans = v.impl()[0];
      m_ytrans = v.impl()[1];
    }

    inline Vector2 reverse(const Vector2 &p) const {
      return Vector2( p(0) - m_xtrans, p(1) - m_ytrans );
    }

    inline Vector2 forward(const Vector2 &p) const {
      return Vector2( p(0) + m_xtrans, p(1) + m_ytrans );
    }
  };


  // Linear function (i.e. 2x2 matrix) transform functor
  class LinearTransform : public TransformHelper<LinearTransform,ConvexFunction,ConvexFunction> {
  protected:
    Matrix2x2 m_matrix, m_matrix_inverse;
    LinearTransform() {}
  public:
    LinearTransform( Matrix2x2 const& matrix )
      : m_matrix( matrix ), m_matrix_inverse( inverse( matrix ) ) {}

    inline Vector2 forward( Vector2 const& p ) const {
      return m_matrix * p;
    }

    inline Vector2 reverse( Vector2 const& p ) const {
      return m_matrix_inverse * p;
    }
  };

  // Affine function (i.e. linear plus translation) transform functor
  class AffineTransform : public TransformHelper<AffineTransform,ConvexFunction,ConvexFunction> {
  protected:
    double m_a, m_b, m_c, m_d;
    double m_ai, m_bi, m_ci, m_di;
    double m_x, m_y;
    AffineTransform() {}
  public:
    AffineTransform( Matrix2x2 const& matrix, Vector2 const& offset ) :
      m_a( matrix(0,0) ), m_b( matrix(0,1) ),
      m_c( matrix(1,0) ), m_d( matrix(1,1) ),
      m_x( offset(0) ), m_y( offset(1) )
    {
      Matrix2x2 inv = inverse( matrix );
      m_ai = inv(0,0);
      m_bi = inv(0,1);
      m_ci = inv(1,0);
      m_di = inv(1,1);
    }

    inline Vector2 forward( Vector2 const& p ) const {
      return Vector2(m_a*p[0]+m_b*p[1] + m_x,
                     m_c*p[0]+m_d*p[1] + m_y);
    }

    inline Vector2 reverse( Vector2 const& p ) const {
      double px = p[0]-m_x;
      double py = p[1]-m_y;
      return Vector2(m_ai*px+m_bi*py,
                     m_ci*px+m_di*py);
    }
   friend std::ostream& operator<<(std::ostream&, const AffineTransform&); 
  };
  inline std::ostream& operator<<(std::ostream& os, const AffineTransform& trans);

  // Helper functions to pull and push a matrix to an affine transform
  inline Matrix3x3 affine2mat(AffineTransform const& transform);
  inline AffineTransform mat2affine(Matrix3x3 const& T);

  // Rotate transform functor
  class RotateTransform : public AffineTransform {
  public:
    RotateTransform( double theta, Vector2 const& translate ) {
      double c = cos(theta), s=sin(theta);
      Matrix2x2 rotate( c, -s, s, c );
      Vector2 rhs = translate - (rotate*translate);
      m_a    = rotate(0,0); m_b = rotate(0,1);
      m_c    = rotate(1,0); m_d = rotate(1,1);
      m_x    = rhs(0);
      m_y    = rhs(1);
      rotate = inverse(rotate);
      m_ai   = rotate(0,0); m_bi = rotate(0,1);
      m_ci   = rotate(1,0); m_di = rotate(1,1);
    }
  };

  // Homography transform functor
  //
  // Transform points by applying a linear operator
  // in homogeneous coordinates, i.e. a 3x3 homography.
  class HomographyTransform : public TransformHelper<HomographyTransform,ConvexFunction,ConvexFunction> {
  protected:
    Matrix3x3 m_H, m_H_inverse;
    HomographyTransform() {}
  public:

    HomographyTransform(Matrix3x3 H) : m_H(H), m_H_inverse( inverse(H) ) {}

    inline Vector2 reverse(const Vector2 &p) const {
      double w = m_H_inverse(2,0) * p(0) + m_H_inverse(2,1) * p(1) + m_H_inverse(2,2);
      return Vector2( ( m_H_inverse(0,0) * p(0) + m_H_inverse(0,1) * p(1) + m_H_inverse(0,2) ) / w,
                      ( m_H_inverse(1,0) * p(0) + m_H_inverse(1,1) * p(1) + m_H_inverse(1,2) ) / w);
    }

    inline Vector2 forward(const Vector2 &p) const {
      double w = m_H(2,0) * p(0) + m_H(2,1) * p(1) + m_H(2,2);
      return Vector2( ( m_H(0,0) * p(0) + m_H(0,1) * p(1) + m_H(0,2) ) / w,
                      ( m_H(1,0) * p(0) + m_H(1,1) * p(1) + m_H(1,2) ) / w);
    }
    friend std::ostream& operator<<(std::ostream&, const HomographyTransform&);
  };
  inline std::ostream& operator<<(std::ostream& os, const HomographyTransform& trans);


  // InverseTransform inverts a transform functor
  template <class TxT>
  class InverseTransform : public TransformBase<InverseTransform<TxT> >
  {
    TxT tx;
  public:
    InverseTransform( TxT const& tx ) : tx(tx) {}

    inline Vector2 forward( Vector2 const& p ) const {
      return tx.reverse( p );
    }

    inline Vector2 reverse( Vector2 const& p ) const {
      return tx.forward( p );
    }

    virtual FunctionType forward_type() const { return tx.reverse_type(); }
    virtual FunctionType reverse_type() const { return tx.forward_type(); }
  };

  // Inverts a transform functor via an InverseTransform object.
  template <class TxT>
  inline InverseTransform<TxT> inverse( TransformBase<TxT> const& tx ) {
    return InverseTransform<TxT>( tx.impl() );
  }


  // CompositionTransform transform functor adaptor
  //
  // This is a wrapper class that allows you to composite two transform
  // functors.  The arguments to the constructor are in the usual
  // function composition order.  That is, CompositionTransform(tx1,tx2)
  // yields a functor whose forward mapping is tx1.forward(tx2.forward(v)),
  // and thus whose reverse mapping is tx2.reverse(tx1.reverse(v)).
  // The usual way to construct a CompositionTransform is with the
  // compose() free function, below.
  template <class Tx1T, class Tx2T>
  class CompositionTransform : public TransformBase<CompositionTransform<Tx1T,Tx2T> >
  {
    Tx1T tx1;
    Tx2T tx2;
  public:
    CompositionTransform( Tx1T const& tx1, Tx2T const& tx2 ) : tx1(tx1), tx2(tx2) {}

    inline Vector2 forward( Vector2 const& p ) const {
      return tx1.forward( tx2.forward( p ) );
    }

    inline Vector2 reverse( Vector2 const& p ) const {
      return tx2.reverse( tx1.reverse( p ) );
    }

    static FunctionType compose_type( FunctionType f, FunctionType g ) {
      if( f==DiscontinuousFunction || g==DiscontinuousFunction ) return DiscontinuousFunction;
      if( f==ContinuousFunction || g==ContinuousFunction ) return ContinuousFunction;
      return ConvexFunction;
    }

    virtual FunctionType forward_type() const { return compose_type( tx1.forward_type(), tx2.forward_type() ); }
    virtual FunctionType reverse_type() const { return compose_type( tx1.reverse_type(), tx2.reverse_type() ); }
  };

  // Composes two transform functors via a CompositionTransform object.
  template <class Tx1T, class Tx2T>
  CompositionTransform<Tx1T,Tx2T>
  inline compose( TransformBase<Tx1T> const& tx1,
                  TransformBase<Tx2T> const& tx2 ) {
    typedef CompositionTransform<Tx1T,Tx2T> result_type;
    return result_type( tx1.impl(), tx2.impl() );
  }

  // Composes three transform functors via a CompositionTransform object.
  template <class Tx1T, class Tx2T, class Tx3T>
  CompositionTransform<Tx1T,CompositionTransform<Tx2T,Tx3T> >
  inline compose( TransformBase<Tx1T> const& tx1,
                  TransformBase<Tx2T> const& tx2,
                  TransformBase<Tx3T> const& tx3 ) {
    typedef CompositionTransform<Tx1T,CompositionTransform<Tx2T,Tx3T> > result_type;
    return result_type( tx1.impl(), compose( tx2, tx3 ) );
  }

  // Composes four transform functors via a CompositionTransform object.
  template <class Tx1T, class Tx2T, class Tx3T, class Tx4T>
  CompositionTransform<Tx1T,CompositionTransform<Tx2T,CompositionTransform<Tx3T,Tx4T> > >
  inline compose( TransformBase<Tx1T> const& tx1,
                  TransformBase<Tx2T> const& tx2,
                  TransformBase<Tx3T> const& tx3,
                  TransformBase<Tx4T> const& tx4 ) {
    typedef CompositionTransform<Tx1T,CompositionTransform<Tx2T,CompositionTransform<Tx3T,Tx4T> > > result_type;
    return result_type( tx1.impl(), compose( tx2, tx3, tx4 ) );
  }

  // TransformRef virtualized transform functor adaptor
  class TransformRef : public TransformBase<TransformRef> {
    boost::shared_ptr<Transform> m_transform;
  public:
    template <class TransformT> explicit TransformRef( TransformT const& transform ) : m_transform( new TransformT( transform ) ) {}
    explicit TransformRef( boost::shared_ptr<Transform> const& transform ) : m_transform( transform ) {}

    FunctionType forward_type () const { return m_transform->forward_type(); }
    FunctionType reverse_type () const { return m_transform->reverse_type(); }
    double       tolerance    () const { return m_transform->tolerance   (); }

    Vector2      forward      ( Vector2 const& point ) const { return m_transform->forward( point );     }
    Vector2      reverse      ( Vector2 const& point ) const { return m_transform->reverse( point );     }
    BBox2i       forward_bbox ( BBox2i  const& bbox  ) const { return m_transform->forward_bbox( bbox ); }
    BBox2i       reverse_bbox ( BBox2i  const& bbox  ) const { return m_transform->reverse_bbox( bbox ); }
    void         set_tolerance( double  tolerance    )       { m_transform->set_tolerance( tolerance );  }
  };



//==================================================================================
// Function definitions
  
Matrix3x3 affine2mat(AffineTransform const& transform){
  Vector2 O = transform.forward(Vector2(0, 0));
  Vector2 A = transform.forward(Vector2(1, 0));
  Vector2 B = transform.forward(Vector2(0, 1));

  A -= O;
  B -= O;
  Matrix3x3 T;
  T.set_identity();
  T(0, 0) = A[0]; T(1, 0) = A[1];
  T(0, 1) = B[0]; T(1, 1) = B[1];
  T(0, 2) = O[0]; T(1, 2) = O[1];
  return T;
}

AffineTransform mat2affine(Matrix3x3 const& T){
  return AffineTransform(submatrix(T, 0, 0, 2, 2), Vector2(T(0, 2), T(1, 2)));
}

std::ostream& operator<<( std::ostream& os, AffineTransform const& trans ) {
  std::ostringstream oss; // To use custom precision
  oss.precision(10);
  oss << "AffineTransform: " << trans.m_a  << ", " << trans.m_b  << ", " 
                             << trans.m_c  << ", " << trans.m_d  << std::endl;
  oss << "               : " << trans.m_ai << ", " << trans.m_bi << ", " 
                             << trans.m_ci << ", " << trans.m_di << std::endl;
  oss << "               : " << trans.m_x  << ", " << trans.m_y  << std::endl;
  os << oss.str();
  return os;
}

std::ostream& operator<<( std::ostream& os, HomographyTransform const& trans ) {
  std::ostringstream oss; // To use custom precision
  oss.precision(10);
  oss << "HomograpyTransform: " << trans.m_H << std::endl;
  os << oss.str();
  return os;
}

} // namespace vw

#endif // __VW_MATH_TRANSFORM_H__
