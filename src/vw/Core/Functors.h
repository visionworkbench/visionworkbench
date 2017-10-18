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


/// \file vw/Core/Functors.h
///
/// General-purpose functors and functor support code.
///
/// This file provides a few mix-in base classes that you can use when
/// writing your own polymorphic functors.  They instruct the return
/// type deduction system about the return type of your functor in
/// some of the common cases.
///
/// It also provides a collection of functors implementing the
/// standard arithmetic operations and so forth, which also serve to
/// demonstrate how to use the mix-ins.
///
#ifndef __VW_CORE_FUNCTORS_H__
#define __VW_CORE_FUNCTORS_H__
#include <vw/config.h>

#include <vw/Core/TypeDeduction.h>

namespace vw {

  /// A mix-in specifying that a functor is an unary functor
  /// whose return type is the same as its argument type.
  struct VW_API UnaryReturnSameType {
    /// \cond INTERNAL
    template <class Args> struct result;
    template <class ChannelT, class FuncT>
    struct result<FuncT(ChannelT)> {
      typedef ChannelT type;
    };
    /// \endcond
  };

  /// A mix-in specifying that a functor is a unary functor
  /// whose return type is determined by the given traits
  /// template class.
  template <template<class> class ResultT>
  struct UnaryReturnTemplateType {
    /// \cond INTERNAL
    template <class Args> struct result {};
    template <class FuncT, class ArgT>
    struct result<FuncT(ArgT)> {
      typedef typename ResultT<ArgT>::type type;
    };
    /// \endcond
  };

  /// A mix-in specifying that a functor is a binary functor
  /// whose return type is determined by the given traits
  /// template class.
  template <template<class,class> class ResultT>
  struct BinaryReturnTemplateType {
    /// \cond INTERNAL
    template <class Args> struct result {};
    template <class FuncT, class Arg1T, class Arg2T>
    struct result<FuncT(Arg1T,Arg2T)> {
      typedef typename ResultT<Arg1T,Arg2T>::type type;
    };
    /// \endcond
  };

  /// A mix-in specifying that a functor always returns a
  /// particular type.
  template <class T>
  struct ReturnFixedType {
    /// \cond INTERNAL
    template <class Args>
    struct result {
      typedef T type;
    };
    /// \endcond
  };

  /// A mix-in specifying that a binary functor always returns
  /// the same type as its first argument.
  struct VW_API BinaryReturn1stType {
    /// \cond INTERNAL
    template <class Args> struct result;
    template <class FuncT, class Arg1T, class Arg2T>
    struct result<FuncT(Arg1T,Arg2T)> {
      typedef Arg1T type;
    };
    /// \endcond
  };

  /// A mix-in specifying that a binary functor always returns
  /// the same type as its second argument.
  struct VW_API BinaryReturn2ndType {
    /// \cond INTERNAL
    template <class Args> struct result;
    template<class FuncT, class Arg1T, class Arg2T>
    struct result<FuncT(Arg1T,Arg2T)> {
      typedef Arg2T type;
    };
    /// \endcond
  };

  /// A mix-in specifying that a functor is an unary functor
  /// whose return type is determined by the given binary
  /// type function with the given first argument.
  template <template<class,class> class ResultT, class ValT>
  struct UnaryReturnBinaryTemplateBind1st {
    /// \cond INTERNAL
    template <class Args> struct result { typedef void type; };
    template <class FuncT, class ArgT>
    struct result<FuncT(ArgT)> {
      typedef typename ResultT<ValT,ArgT>::type type;
    };
    /// \endcond
  };

  /// A mix-in specifying that a functor is an unary functor
  /// whose return type is determined by the given binary
  /// type function with the given second argument.
  template <template<class,class> class ResultT, class ValT>
  struct UnaryReturnBinaryTemplateBind2nd {
    /// \cond INTERNAL
    template <class Args> struct result { typedef void type; };
    template <class FuncT, class ArgT>
    struct result<FuncT(ArgT)> {
      typedef typename ResultT<ArgT,ValT>::type type;
    };
    /// \endcond
  };

  /// A mix-in specifying that a functor is a binary functor
  /// whose return type is determined by the given ternary
  /// type function with the given first argument.
  template <template<class,class,class> class ResultT, class ValT>
  struct BinaryReturnTernaryTemplateBind1st {
    /// \cond INTERNAL
    template <class Args> struct result { typedef void type; };
    template <class FuncT, class Arg1T, class Arg2T>
    struct result<FuncT(Arg1T,Arg2T)> {
      typedef typename ResultT<ValT,Arg1T,Arg2T>::type type;
    };
    /// \endcond
  };

  /// A mix-in specifying that a functor is an unary functor
  /// whose return type is determined by the given ternary
  /// type function with the given first and second arguments.
  template <template<class,class,class> class ResultT, class Val1T, class Val2T>
  struct UnaryReturnTernaryTemplateBind1st2nd {
    /// \cond INTERNAL
    template <class Args> struct result { typedef void type; };
    template <class FuncT, class ArgT>
    struct result<FuncT(ArgT)> {
      typedef typename ResultT<Val1T,Val2T,ArgT>::type type;
    };
    /// \endcond
  };

  /// A mix-in specifying that a functor is an unary functor
  /// whose return type is determined by the given ternary
  /// type function with the given first and third arguments.
  template <template<class,class,class> class ResultT, class Val1T, class Val3T>
  struct UnaryReturnTernaryTemplateBind1st3rd {
    /// \cond INTERNAL
    template <class Args> struct result { typedef void type; };
    template <class FuncT, class ArgT>
    struct result<FuncT(ArgT)> {
      typedef typename ResultT<Val1T,ArgT,Val3T>::type type;
    };
    /// \endcond
  };


  // ********************************************************************
  // Utility Functors
  // ********************************************************************

  // Unary no-operation functor
  struct VW_API NoOpFunctor : ReturnFixedType<void> {
    template <class ArgT>
    inline void operator()( ArgT const& /*arg*/ ) const {}
  };

  // Unary identity functor
  struct VW_API IdentityFunctor : UnaryReturnSameType {
    template <class ArgT>
    inline ArgT const& operator()( ArgT const& arg ) const {
      return arg;
    }
  };


  // ********************************************************************
  // Arithmetic Operator Functors
  // ********************************************************************

  // Unary negation of an argument
  struct VW_API ArgNegationFunctor : UnaryReturnSameType {
    template <class ArgT>
    typename result<ArgNegationFunctor(ArgT)>::type
    inline operator()( ArgT const& arg ) const { return -arg; }
  };

  // Binary sum of two arguments
  struct VW_API ArgArgSumFunctor : BinaryReturnTemplateType<SumType> {
    template <class Arg1T, class Arg2T>
    inline typename SumType<Arg1T,Arg2T>::type operator()( Arg1T const& arg1, Arg2T const& arg2 ) const { return arg1+arg2; }
  };

  // Unary sum of an argument and a value
  template <class ValT>
  struct ArgValSumFunctor : UnaryReturnBinaryTemplateBind2nd<SumType,ValT> {
  private:
    const ValT m_val;
  public:
    ArgValSumFunctor( ValT const& val ) : m_val(val) {}

    template <class ArgT>
    inline typename SumType<ArgT,ValT>::type operator()( ArgT const& arg ) const { return arg+m_val; }
  };

  // Unary sum of a value and an argument
  template <class ValT>
  struct ValArgSumFunctor : UnaryReturnBinaryTemplateBind1st<SumType,ValT> {
  private:
    const ValT m_val;
  public:
    ValArgSumFunctor( ValT const& val ) : m_val(val) {}

    template <class ArgT>
    inline typename SumType<ValT,ArgT>::type operator()( ArgT const& arg ) const { return m_val+arg; }
  };

  // Binary in-place sum of two arguments
  struct ArgArgInPlaceSumFunctor : BinaryReturn1stType {
    template <class Arg1T, class Arg2T>
    inline Arg1T& operator()( Arg1T& arg1, Arg2T const& arg2 ) const { return arg1=(Arg1T)(arg1+arg2); }
  };

  // Unary in-place sum of an argument and a value
  template <class ValT>
  struct ArgValInPlaceSumFunctor : UnaryReturnSameType {
  private:
    const ValT m_val;
  public:
    ArgValInPlaceSumFunctor( ValT const& val ) : m_val(val) {}

    template <class ArgT>
    inline ArgT& operator()( ArgT& arg ) const { return arg=(ArgT)(arg+m_val); }
  };

  // Binary difference of two arguments
  struct ArgArgDifferenceFunctor : BinaryReturnTemplateType<DifferenceType> {
    template <class Arg1T, class Arg2T>
    inline typename DifferenceType<Arg1T,Arg2T>::type operator()( Arg1T const& arg1, Arg2T const& arg2 ) const { return arg1-arg2; }
  };

  // Unary difference of an argument and a value
  template <class ValT>
  struct ArgValDifferenceFunctor : UnaryReturnBinaryTemplateBind2nd<DifferenceType,ValT> {
  private:
    const ValT m_val;
  public:
    ArgValDifferenceFunctor( ValT const& val ) : m_val(val) {}

    template <class ArgT>
    inline typename DifferenceType<ArgT,ValT>::type operator()( ArgT const& arg ) const { return arg-m_val; }
  };

  // Unary difference of a value and an argument
  template <class ValT>
  struct ValArgDifferenceFunctor : UnaryReturnBinaryTemplateBind1st<DifferenceType,ValT> {
  private:
    const ValT m_val;
  public:
    ValArgDifferenceFunctor( ValT const& val ) : m_val(val) {}

    template <class ArgT>
    inline typename DifferenceType<ValT,ArgT>::type operator()( ArgT const& arg ) const { return m_val-arg; }
  };

  // Binary in-place difference of two arguments
  struct VW_API ArgArgInPlaceDifferenceFunctor : BinaryReturn1stType {
    template <class Arg1T, class Arg2T>
    inline Arg1T& operator()( Arg1T& arg1, Arg2T const& arg2 ) const { return arg1=(Arg1T)(arg1-arg2); }
  };

  // Unary in-place difference of an argument and a value
  template <class ValT>
  struct ArgValInPlaceDifferenceFunctor : UnaryReturnSameType {
  private:
    const ValT m_val;
  public:
    ArgValInPlaceDifferenceFunctor( ValT const& val ) : m_val(val) {}

    template <class ArgT>
    inline ArgT& operator()( ArgT& arg ) const { return arg=(ArgT)(arg-m_val); }
  };

  // Binary product of two arguments
  struct VW_API ArgArgProductFunctor : BinaryReturnTemplateType<ProductType> {
    template <class Arg1T, class Arg2T>
    inline typename ProductType<Arg1T,Arg2T>::type operator()( Arg1T const& arg1, Arg2T const& arg2 ) const { return arg1*arg2; }
  };

  // Unary product of an argument and a value
  template <class ValT>
  struct ArgValProductFunctor : UnaryReturnBinaryTemplateBind2nd<ProductType,ValT> {
  private:
    const ValT m_val;
  public:
    ArgValProductFunctor( ValT const& val ) : m_val(val) {}

    template <class ArgT>
    inline typename ProductType<ArgT,ValT>::type operator()( ArgT const& arg ) const { return arg*m_val; }
  };

  // Unary product of a value and an argument
  template <class ValT>
  struct ValArgProductFunctor : UnaryReturnBinaryTemplateBind1st<ProductType,ValT> {
  private:
    const ValT m_val;
  public:
    ValArgProductFunctor( ValT const& val ) : m_val(val) {}

    template <class ArgT>
    inline typename ProductType<ValT,ArgT>::type operator()( ArgT const& arg ) const { return m_val*arg; }
  };

  // Binary in-place product of two arguments
  struct VW_API ArgArgInPlaceProductFunctor : BinaryReturn1stType {
    template <class Arg1T, class Arg2T>
    inline Arg1T& operator()( Arg1T& arg1, Arg2T const& arg2 ) const { return arg1=(Arg1T)(arg1*arg2); }
  };

  // Unary in-place product of an argument and a value
  template <class ValT>
  struct ArgValInPlaceProductFunctor : UnaryReturnSameType {
  private:
    const ValT m_val;
  public:
    ArgValInPlaceProductFunctor( ValT const& val ) : m_val(val) {}

    template <class ArgT>
    inline ArgT& operator()( ArgT& arg ) const { return arg=(ArgT)(arg*m_val); }
  };

  // Binary quotient of two arguments
  struct VW_API ArgArgQuotientFunctor : BinaryReturnTemplateType<QuotientType> {
    template <class Arg1T, class Arg2T>
    inline typename QuotientType<Arg1T,Arg2T>::type operator()( Arg1T const& arg1, Arg2T const& arg2 ) const { return arg1/arg2; }
  };

  // Unary quotient of an argument and a value
  template <class ValT>
  struct ArgValQuotientFunctor : UnaryReturnBinaryTemplateBind2nd<QuotientType,ValT> {
  private:
    const ValT m_val;
  public:
    ArgValQuotientFunctor( ValT const& val ) : m_val(val) {}

    template <class ArgT>
    inline typename QuotientType<ArgT,ValT>::type operator()( ArgT const& arg ) const { return arg/m_val; }
  };

  // Unary quotient of a value and an argument
  template <class ValT>
  struct ValArgQuotientFunctor : UnaryReturnBinaryTemplateBind1st<QuotientType,ValT> {
  private:
    const ValT m_val;
  public:
    ValArgQuotientFunctor( ValT const& val ) : m_val(val) {}

    template <class ArgT>
    inline typename QuotientType<ValT,ArgT>::type operator()( ArgT const& arg ) const { return m_val/arg; }
  };

  // Binary in-place quotient of two arguments
  struct VW_API ArgArgInPlaceQuotientFunctor : BinaryReturn1stType {
    template <class Arg1T, class Arg2T>
    inline Arg1T& operator()( Arg1T& arg1, Arg2T const& arg2 ) const { return arg1=(Arg1T)(arg1/arg2); }
  };

  // Unary in-place quotient of an argument and a value
  template <class ValT>
  struct ArgValInPlaceQuotientFunctor : UnaryReturnSameType {
  private:
    const ValT m_val;
  public:
    ArgValInPlaceQuotientFunctor( ValT const& val ) : m_val(val) {}

    template <class ArgT>
    inline ArgT& operator()( ArgT& arg ) const { return arg=(ArgT)(arg/m_val); }
  };

  // **** WARNING ****
  // The syntax in all three SafeQuotient functors is carefully
  // crafted to work around an issue with RedHat's gcc 3.2.3-56.
  // If you change this, be sure to test it on an RHEL3 box.
  // **** WARNING ****

  // Safe binary quotient of two arguments
  struct ArgArgSafeQuotientFunctor : BinaryReturnTemplateType<QuotientType> {
    template <class Arg1T, class Arg2T>
    inline typename QuotientType<Arg1T,Arg2T>::type operator()( Arg1T const& arg1, Arg2T const& arg2 ) const {
      if( arg2==Arg2T() ) return typename QuotientType<Arg1T,Arg2T>::type();
      else return (arg1/arg2);
    }
  };

  // Safe unary quotient of an argument and a value
  template <class ValT>
  struct ArgValSafeQuotientFunctor : UnaryReturnBinaryTemplateBind2nd<QuotientType,ValT> {
  private:
    const ValT m_val;
  public:
    ArgValSafeQuotientFunctor( ValT const& val ) : m_val(val) {}

    template <class ArgT>
    inline typename QuotientType<ArgT,ValT>::type operator()( ArgT const& arg ) const {
      if ( m_val==ValT() ) return typename QuotientType<ArgT,ValT>::type();
      else return (arg/m_val);
    }
  };

  // Safe unary quotient of a value and an argument
  template <class ValT>
  struct ValArgSafeQuotientFunctor : UnaryReturnBinaryTemplateBind1st<QuotientType,ValT> {
  private:
    const ValT m_val;
  public:
    ValArgSafeQuotientFunctor( ValT const& val ) : m_val(val) {}

    template <class ArgT>
    inline typename QuotientType<ValT,ArgT>::type operator()( ArgT const& arg ) const {
      if ( arg==ArgT() ) return typename QuotientType<ValT,ArgT>::type();
      else return (m_val/arg);
    }
  };

  // Safe binary in-place quotient of two arguments
  struct VW_API ArgArgInPlaceSafeQuotientFunctor : BinaryReturn1stType {
    template <class Arg1T, class Arg2T>
    inline Arg1T& operator()( Arg1T& arg1, Arg2T const& arg2 ) const {
      if ( arg2==Arg2T() ) return arg1=Arg1T();
      else return arg1=(Arg1T)(arg1/arg2);
    }
  };

  // Safe unary in-place quotient of an argument and a value
  template <class ValT>
  struct ArgValInPlaceSafeQuotientFunctor : UnaryReturnSameType {
  private:
    const ValT m_val;
  public:
    ArgValInPlaceSafeQuotientFunctor( ValT const& val ) : m_val(val) {}

    template <class ArgT>
    inline ArgT& operator()( ArgT& arg ) const {
      if ( m_val==ValT() ) return arg=ArgT();
      else return arg=(ArgT)(arg/m_val);
    }
  };

  // Binary equality operator of two arguments
  struct VW_API ArgArgEqualityFunctor : ReturnFixedType<bool> {
    template <class Arg1T, class Arg2T>
    inline bool operator()( Arg1T const& arg1, Arg2T const& arg2 ) const { return arg1==arg2; }
  };

  // Unary equality operator of an argument and a value
  template <class ValT>
  struct ArgValEqualityFunctor : ReturnFixedType<bool> {
  private:
    const ValT m_val;
  public:
    ArgValEqualityFunctor( ValT const& val ) : m_val(val) {}

    template <class ArgT>
    inline bool operator()( ArgT const& arg ) const { return arg==m_val; }
  };

  // Unary equality operator of a value and an argument
  template <class ValT>
  struct ValArgEqualityFunctor : ReturnFixedType<bool> {
  private:
    const ValT m_val;
  public:
    ValArgEqualityFunctor( ValT const& val ) : m_val(val) {}

    template <class ArgT>
    inline bool operator()( ArgT const& arg ) const { return m_val==arg; }
  };

  // Binary inequality operator of two arguments
  struct VW_API ArgArgInequalityFunctor : ReturnFixedType<bool> {
    template <class Arg1T, class Arg2T>
    inline bool operator()( Arg1T const& arg1, Arg2T const& arg2 ) const { return arg1!=arg2; }
  };

  // Unary inequality operator of an argument and a value
  template <class ValT>
  struct ArgValInequalityFunctor : ReturnFixedType<bool> {
  private:
    const ValT m_val;
  public:
    ArgValInequalityFunctor( ValT const& val ) : m_val(val) {}

    template <class ArgT>
    inline bool operator()( ArgT const& arg ) const { return arg!=m_val; }
  };

  // Unary inequality operator of a value and an argument
  template <class ValT>
  struct ValArgInequalityFunctor : ReturnFixedType<bool> {
  private:
    const ValT m_val;
  public:
    ValArgInequalityFunctor( ValT const& val ) : m_val(val) {}

    template <class ArgT>
    inline bool operator()( ArgT const& arg ) const { return m_val!=arg; }
  };

  // Binary less-than operator of two arguments
  struct VW_API ArgArgLessThanFunctor : ReturnFixedType<bool> {
    template <class Arg1T, class Arg2T>
    inline bool operator()( Arg1T const& arg1, Arg2T const& arg2 ) const { return arg1<arg2; }
  };

  // Unary less-than operator of an argument and a value
  template <class ValT>
  struct ArgValLessThanFunctor : ReturnFixedType<bool> {
  private:
    const ValT m_val;
  public:
    ArgValLessThanFunctor( ValT const& val ) : m_val(val) {}

    template <class ArgT>
    inline bool operator()( ArgT const& arg ) const { return arg<m_val; }
  };

  // Unary less-than operator of a value and an argument
  template <class ValT>
  struct ValArgLessThanFunctor : ReturnFixedType<bool> {
  private:
    const ValT m_val;
  public:
    ValArgLessThanFunctor( ValT const& val ) : m_val(val) {}

    template <class ArgT>
    inline bool operator()( ArgT const& arg ) const { return m_val<arg; }
  };

  // Binary less-than-or-equal operator of two arguments
  struct VW_API ArgArgLessThanOrEqualFunctor : ReturnFixedType<bool> {
    template <class Arg1T, class Arg2T>
    inline bool operator()( Arg1T const& arg1, Arg2T const& arg2 ) const { return arg1<=arg2; }
  };

  // Unary less-than-or-equal operator of an argument and a value
  template <class ValT>
  struct ArgValLessThanOrEqualFunctor : ReturnFixedType<bool> {
  private:
    const ValT m_val;
  public:
    ArgValLessThanOrEqualFunctor( ValT const& val ) : m_val(val) {}

    template <class ArgT>
    inline bool operator()( ArgT const& arg ) const { return arg<=m_val; }
  };

  // Unary less-than-or-equal operator of a value and an argument
  template <class ValT>
  struct ValArgLessThanOrEqualFunctor : ReturnFixedType<bool> {
  private:
    const ValT m_val;
  public:
    ValArgLessThanOrEqualFunctor( ValT const& val ) : m_val(val) {}

    template <class ArgT>
    inline bool operator()( ArgT const& arg ) const { return m_val<=arg; }
  };

  // Binary greater-than operator of two arguments
  struct VW_API ArgArgGreaterThanFunctor : ReturnFixedType<bool> {
    template <class Arg1T, class Arg2T>
    inline bool operator()( Arg1T const& arg1, Arg2T const& arg2 ) const { return arg1>arg2; }
  };

  // Unary greater-than operator of an argument and a value
  template <class ValT>
  struct ArgValGreaterThanFunctor : ReturnFixedType<bool> {
  private:
    const ValT m_val;
  public:
    ArgValGreaterThanFunctor( ValT const& val ) : m_val(val) {}

    template <class ArgT>
    inline bool operator()( ArgT const& arg ) const { return arg>m_val; }
  };

  // Unary greater-than operator of a value and an argument
  template <class ValT>
  struct ValArgGreaterThanFunctor : ReturnFixedType<bool> {
  private:
    const ValT m_val;
  public:
    ValArgGreaterThanFunctor( ValT const& val ) : m_val(val) {}

    template <class ArgT>
    inline bool operator()( ArgT const& arg ) const { return m_val>arg; }
  };

  // Binary greater-than-or-equal operator of two arguments
  struct VW_API ArgArgGreaterThanOrEqualFunctor : ReturnFixedType<bool> {
    template <class Arg1T, class Arg2T>
    inline bool operator()( Arg1T const& arg1, Arg2T const& arg2 ) const { return arg1>=arg2; }
  };

  // Unary greater-than-or-equal operator of an argument and a value
  template <class ValT>
  struct ArgValGreaterThanOrEqualFunctor : ReturnFixedType<bool> {
  private:
    const ValT m_val;
  public:
    ArgValGreaterThanOrEqualFunctor( ValT const& val ) : m_val(val) {}

    template <class ArgT>
    inline bool operator()( ArgT const& arg ) const { return arg>=m_val; }
  };

  // Unary greater-than-or-equal operator of a value and an argument
  template <class ValT>
  struct ValArgGreaterThanOrEqualFunctor : ReturnFixedType<bool> {
  private:
    const ValT m_val;
  public:
    ValArgGreaterThanOrEqualFunctor( ValT const& val ) : m_val(val) {}

    template <class ArgT>
    inline bool operator()( ArgT const& arg ) const { return m_val>=arg; }
  };

  struct NOP {
    template <typename T>
    void operator()(T) const {}
  };

} // namespace vw

#endif  // __VW_CORE_FUNCTORS_H__
