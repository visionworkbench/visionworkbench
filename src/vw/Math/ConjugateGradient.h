// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file Math/ConjugateGradient.h
///
/// Conjugate gradient and related descent methods and line
/// search methods.  Note!  No stopping criteria have been
/// implemented besdies iteration count!
///
/// To use these functions, you must first express your
/// problem in a particular form.  That form is a little
/// more convoluted right now than I would like.  You supply
/// a cost functor that meets these requirements:
/// * Defines a result_type that is the type returned by
///   evaluating the functor.  Typically float or double.
/// * Defines a domain_type that is the type of the search
///   space.  Often a Vector<foo>, but can reflect other
///   topologies if needed.
/// * Defines a gradient_type corresponding to the space of
///   tangent vectors.  Typically Vector<foo>.
/// * Defines a method: result_type operator()( domain_type const& x ) const;
///   that evaluates the cost function at the given point.
/// * Defines a method: gradient_type gradient( domain_type const& x ) const;
///   that evaluates the gradient of the cost function at the
///   given point.
/// * Defines a method: unsigned dimension() const; that returns
///   the dimension of the tangent space (i.e. the dimension of
///   the vector returned by the gradient() method).
/// * The domain_type must implement a method: domain_type domain_type::operator+( gradient_type const& g ) const;
///   that adds a tangent vector to a position.  You get this for
///   free if both domain_type and gradient_type are Vector<foo>.
///   This is where you do most of the hard work if domain_type
///   represents some non-trivial topological space.
/// * The gradient_type must implement scalar multiplication on
///   the left.  You get this for free in the usual case when
///   gradient_type is just Vector<foo>.
///
/// Once you've set up your problem in that form, just create
/// a functor object and a domain_type object corresponding to
/// your initial guess and then call the optimization function
/// like this:
///   result = conjugate_gradient( cost_functor, initial_guess, ArmijoStepSize(max_stepsize), numiters);
/// I recommend using the Armijo step-size rule for now.  The
/// stepsize that you pass to the constructor should be an
/// over-estimate of the stepsize, i.e. of the ratio between
/// the distance travelled and the slope of the cost function.
/// In a bundle adjustment problem where the cost function
/// represents the sum of a very large number of square pixel
/// errors the stepsize will often be quite small.  I've been
/// using 1e-6 for my problem with good results.  Play around.
///
/// I also offer a constant step size, which there's almost no
/// good reason to use, as well as a Charalambous step size rule,
/// which may be buggy and certainly is underperforming Armijo
/// for me at the moment.  I also provide a steepest_descent()
/// method for comparison to conjugate_gradient().

#ifndef __VW_MATH_CONJUGATEGRADIENT_H__
#define __VW_MATH_CONJUGATEGRADIENT_H__

#include <vw/Core/Log.h>

#define VW_CONJGRAD_MAX_ITERS_BETWEEN_SPACER_STEPS 20

namespace vw {
namespace math{

  class ConstantStepSize {
    double stepsize;
  public:
    ConstantStepSize( double stepsize ) : stepsize(stepsize) {}
    template <class FuncT>
    typename FuncT::domain_type operator()( FuncT const& /*func*/,
                                            typename FuncT::domain_type const& pos,
                                            typename FuncT::result_type const& /*val*/,
                                            typename FuncT::gradient_type const& /*grad*/,
                                            typename FuncT::gradient_type const& dir ) const {
      return pos + stepsize * dir;
    }
  };


  class ArmijoStepSize {
    double m_initial_stepsize, m_beta, m_sigma;
    static const int m_max_iterations = 1000;
  public:
    ArmijoStepSize( double initial_stepsize=1.0, double beta=0.1, double sigma=0.001 )
      : m_initial_stepsize(initial_stepsize), m_beta(beta), m_sigma(sigma) {}
    template <class FuncT>
    typename FuncT::domain_type operator()( FuncT const& func,
                                            typename FuncT::domain_type const& pos,
                                            typename FuncT::result_type const& val,
                                            typename FuncT::gradient_type const& grad,
                                            typename FuncT::gradient_type const& dir ) const {
      double stepsize = m_initial_stepsize;
      double thresh = m_sigma*m_initial_stepsize*dot_prod(grad,dir);
      int count = 0;
      while( true ) {
        typename FuncT::domain_type new_pos = pos + stepsize * dir;
        typename FuncT::result_type new_val = func(new_pos);
        if( new_val - val <= thresh || count > m_max_iterations ) return new_pos;
        //if( ( stepsize < 1e-20 && val-new_val > 0 ) || stepsize < 1e-40 ) {
        //  vw_out(DebugMessage, "math") << "ArmijoStepSize punting!  (slope=" << dot_prod(grad,dir) << ", delta=" << (new_val-val) << ", thresh=" << thresh << ")" << std::endl;
        //  return new_pos;
        //}
        stepsize *= m_beta;
        thresh *= m_beta;
        count++;
      }
    }
  };


  class CharalambousStepSize {
    double initial_stepsize;
    double mu, sigma, gamma;
  public:
    CharalambousStepSize( double initial_stepsize=1.0, double mu=1e-3, double sigma=0.01, double gamma=1e-2 )
      : initial_stepsize(initial_stepsize), mu(mu), sigma(sigma), gamma(gamma) {}
    template <class FuncT>
    typename FuncT::domain_type operator()( FuncT const& func,
                                            typename FuncT::domain_type const& pos,
                                            typename FuncT::result_type const& val,
                                            typename FuncT::gradient_type const& grad,
                                            typename FuncT::gradient_type const& dir ) const {
      double stepdelta = initial_stepsize;
      double base_slope = dot_prod(grad,dir);
      double cur_stepsize = 0;
      double cur_val = val;
      double cur_slope = base_slope;
      vw_out(DebugMessage, "math") << "*** CharalambousStepSize ***" << std::endl;
      vw_out(DebugMessage, "math") << "Initial val: " << cur_val << std::endl;
      vw_out(DebugMessage, "math") << "Initial slope: " << cur_slope << std::endl;
      while( true ) {
        double new_stepsize = cur_stepsize + stepdelta;
        vw_out(DebugMessage, "math") << "step: " << new_stepsize << std::endl;
        typename FuncT::domain_type new_pos = pos + new_stepsize * dir;
        typename FuncT::result_type new_val = func(new_pos);
        typename FuncT::gradient_type new_grad = func.gradient(new_pos);
        double new_slope = dot_prod(new_grad,dir);
        vw_out(DebugMessage, "math") << "val: " << new_val << std::endl;
        vw_out(DebugMessage, "math") << "slope: " << new_slope << std::endl;
        if( new_val <= val + mu*new_stepsize*base_slope && fabs(new_slope) <= -sigma*base_slope ) {
          vw_out(DebugMessage, "math") << "*** ACCEPTING! ***" << std::endl;
          return new_pos;
        }
        else if( new_slope <= 0 ) {
          if( new_val >= val ) {
            vw_out(DebugMessage, "math") << "Too large!" << std::endl;
            stepdelta *= 0.1;
          }
          else {
            vw_out(DebugMessage, "math") << "Too small! (Retaining.)" << std::endl;
            cur_stepsize = new_stepsize;
            cur_val = new_val;
            cur_slope = new_slope;
            stepdelta *= 10;
          }
        }
        else {
          vw_out(DebugMessage, "math") << "Interpolating!" << std::endl;
          double a = 3*((cur_slope+new_slope)*stepdelta - 2*(new_val-cur_val))/(stepdelta*stepdelta*stepdelta);
          double b = 2*(3*(new_val-cur_val) - (2*cur_slope+new_slope)*stepdelta)/(stepdelta*stepdelta);
          vw_out(DebugMessage, "math") << "a=" << a << " b=" << b << " c=" << cur_slope << std::endl;
          double step = 0;
          if( fabs(a)*stepdelta < 1e-4*fabs(b) ) {
            vw_out(DebugMessage, "math") << "Quadratic" << std::endl;
            step = -cur_slope/b;
          }
          else {
            vw_out(DebugMessage, "math") << "Cubic" << std::endl;
            step = (-b+sqrt(b*b-4*a*cur_slope))/(2*a);
          }
          vw_out(DebugMessage, "math") << "Before: " << step << std::endl;
          if( step < gamma*stepdelta ) step = gamma*stepdelta;
          else if( step > (1-gamma)*stepdelta ) step = (1-gamma)*stepdelta;
          vw_out(DebugMessage, "math") << "After: " << step << std::endl;
          stepdelta = step;
        }
      }
    }
  };


  template <class FuncT, class StepT>
  typename FuncT::domain_type steepest_descent( FuncT const& func,
                                                typename FuncT::domain_type const& seed,
                                                StepT const& step,
                                                int numiters ) {
    typename FuncT::domain_type pos = seed;
    typename FuncT::result_type val = func(pos);
    vw_out(DebugMessage, "math") << "Initial: " << val << std::endl;
    for( int i=0; i<numiters; ++i ) {
      typename FuncT::gradient_type grad = func.gradient( pos );
      pos = step( func, pos, val, grad, -grad );
      val = func(pos);
      vw_out(DebugMessage, "math") << "Step " << i << ": " << val << std::endl;
    }
    return pos;
  }


  template <class FuncT, class StepT>
  typename FuncT::domain_type conjugate_gradient( FuncT const& func,
                                                  typename FuncT::domain_type const& seed,
                                                  StepT const& step,
                                                  int numiters ) {
    typename FuncT::domain_type pos = seed;
    typename FuncT::result_type val = func(pos);
    vw_out(DebugMessage, "math") << "Initial: " << val << std::endl;
    double last_grad_norm2 = 0;
    typename FuncT::gradient_type last_dir;
    for( int i=0; i<numiters; ++i ) {
      typename FuncT::gradient_type grad = func.gradient(pos);
      typename FuncT::gradient_type dir = -grad;
      double grad_norm2 = dot_prod(grad,grad);
      if( i % VW_CONJGRAD_MAX_ITERS_BETWEEN_SPACER_STEPS != 0 && func.dimension() != 0 )
        dir += (grad_norm2/last_grad_norm2)*last_dir;
      pos = step( func, pos, val, grad, dir );
      last_grad_norm2 = grad_norm2;
      last_dir = dir;
      val = func(pos);
      vw_out(DebugMessage, "math") << "Step " << i << ": " << val << std::endl;
    }
    return pos;
  }

} } // namespace vw::math

#endif // #ifndef __VW_MATH_CONJUGATEGRADIENT_H__
