// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <gtest/gtest.h>
//#include <string>
#include <vw/Math.h>
//#include <vw/Math/ConjugateGradient.h>
#include <vw/Math/KDTree.h>
#include <vw/Math/PoseEstimation.h>
//#include <vw/Image/PixelTypes.h>
//#include <boost/type_traits.hpp>
//#include <boost/utility/enable_if.hpp>
#include <complex>

using namespace vw;
using namespace vw::math;

typedef int T;
typedef T SizeT;
typedef T ValT;
typedef T ElemT;
typedef T DataT;

const int DimN  = 2;
const int DimN1 = 2;
const int DimN2 = 2;

typedef BBox<double, 2> BBoxT;
typedef BBoxT BBoxT1;
typedef BBoxT BBoxT2;

typedef double RealT;
typedef double RealT1;
typedef double RealT2;
typedef double ScalarT;

typedef Vector3i VectorT;
typedef VectorT Vector1T;
typedef VectorT Vector2T;
typedef VectorT VectorT1;
typedef VectorT VectorT2;
typedef VectorT SingularValuesT;

typedef Vector3 BVectorT;

typedef Vector<std::complex<double>, 3> EigenvaluesT;

typedef Matrix3x3 MatrixT;
typedef MatrixT AMatrixT;
typedef MatrixT BMatrixT;
typedef MatrixT Matrix1T;
typedef MatrixT Matrix2T;
typedef MatrixT QMatrixT;
typedef MatrixT RMatrixT;
typedef MatrixT UMatrixT;
typedef MatrixT VTMatrixT;

typedef Matrix<std::complex<double>, 3, 3> VMatrixT;

typedef Quaternion<int> QuaternionT;

typedef Vector3 DomainT;
typedef DomainT ScaleT;

struct FuncT : math::LeastSquaresModelBase<FuncT> {
  typedef Vector<double> result_type;
  typedef Vector<double> domain_type;
  typedef Matrix<double> jacobian_type;
  typedef Vector<double> gradient_type;
  FuncT() {}
  FuncT(int) {}
  template <class T>
  T operator() (T a) const {return a;}
  template <class T>
  T operator() (T a, T b) const {return a+b;}
  gradient_type gradient(domain_type a) const {return a;}
  int dimension() const {return 3;}
  result_type difference(result_type a, result_type b) const {return a+b;}
  jacobian_type jacobian(domain_type a) const {return jacobian_type();}

  double operator() (DomainT a) const {return double();}
};
typedef FuncT ImplT;

typedef ConstantStepSize StepT;

struct ContainerT : public std::vector<int> {
  int rows() const {return 3;}
  int cols() const {return 3;}
};

typedef std::vector<ContainerT> FileT;
typedef FileT::iterator ForwardIterator;

#include "TestInstantiateFreeList.hh"

TEST(InstantiateMathFree, Inst) {}
