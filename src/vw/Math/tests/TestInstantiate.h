#include <vector>
#include <vw/Math.h>
#include <vw/Math/DisjointSet.h>
#include <vw/Math/KDTree.h>

using namespace vw;

const bool Arg1IsCompound = false;
const bool Arg2IsCompound = false;
const bool ArgIsCompound = false;

const int RowsN = 3;
const int SizeN = 5;
const int ColsN = 3;
const int DimN = 3;
const int Transpose1N = 3;
const int Transpose2N = 3;
const int TransposeN = 3;

typedef int AccumT;
typedef short Arg1T[];
typedef int Arg2T[];
typedef int Args;
typedef int ArgT[];
typedef int ElemT;
typedef int T;
typedef int T1;
typedef int T2;
typedef int ValT;
typedef double RealT;
typedef BBox<double, 2> BBoxT;
typedef int DstElemT;

typedef Matrix3x3i MatrixT;
typedef MatrixT Matrix1T;
typedef MatrixT Matrix2T;
typedef MatrixT MatT;
typedef MatrixT::iterator IterT;

typedef Vector3i VecT;
typedef VecT VectorT;
typedef VecT SrcVecT;
typedef VecT Vector1T;
typedef VecT Vector2T;
typedef Vector<int, 5> DstVecT;

typedef std::vector<int> ContainerT;
typedef std::vector<ContainerT> FileT;
typedef int DomainT;
typedef ContainerT RangeT;

struct FuncT : ReturnFixedType<int> {
  FuncT() {}
  FuncT(int) {}
  int operator() (int a) const {return a;}
  int operator() (int a, int b) const {return a+b;}
};
typedef FuncT F;


struct ImplT : math::LeastSquaresModelBase<ImplT> {

  typedef Vector<double> result_type;
  typedef Vector<double> domain_type;
  typedef Matrix<double> jacobian_type;

  ImplT() {}

  inline result_type operator()( domain_type const& x ) const {
    return x;
  }

};

typedef Quaternion<int> QuaternionT;

#include "TestInstantiateList.hh"

class TestInstantiateMath : public CxxTest::TestSuite
{
  public:
  void test_inst() {}
};
