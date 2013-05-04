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


// -*- c++ -*-
#ifndef vw_geometry_ATrans_h
#define vw_geometry_ATrans_h

#include <vw/Math/Vector.h>
#include <vw/Math/Matrix.h>

namespace vw
{
  namespace geometry {
    template<class ElemT, int SizeN>
    class ATrans;

    /** Helper function for inverting an affine transform. */
    template<class ElemT, int SizeN>
    ATrans<ElemT, SizeN> inverse(ATrans<ElemT, SizeN> const& atrans) throw();

    /** Output stream operator for affine transformation class. */
    template<class ElemT, int SizeN>
    std::ostream& operator<< (std::ostream& ostr, ATrans<ElemT, SizeN> const& rhs);


    /**
     * @brief Class representing an affine transform.
     * Affine transforms of dimension N consist of a NxN matrix and an
     * vector of size N.  An affine transform of dimension N can be
     * expressed as a matrix of dimension (N+1)x(N+1) with the last row
     * consisting of [0, 0, ..., 1]. While the later allows for a more
     * compact representation, the former is computationally more
     * efficient.
     *
     * This class represents affine transforms internally as matrix plus
     * vector. The user-visible part tries hard to follow the
     * (N+1)x(N+1) matrix notation.
     *
     * The primary customers of this class are the Frame &
     * FrameStore classes.  The affine transformation in this case holds
     * a position-vector and rotation matrix.
     *
     */
    template<class ElemT = double, int SizeN = 3>
    class ATrans
    {
    public:
      /** Vector type (translation part) */
      typedef vw::math::Vector<ElemT, SizeN> VectorT;
      /** Matrix type (rotation part) */
      typedef vw::math::Matrix<ElemT, SizeN, SizeN> MatrixT;

      /** Default constructor */
      ATrans() {}
      /** Initializing constructor */
      ATrans(VectorT const& p, MatrixT const& r) : m_p(p), m_r(r) {}
      /** Initialize from transformation matrix */
      ATrans(vw::math::Matrix<ElemT, SizeN + 1, SizeN + 1> const& m);
      /** Initialize from transformation matrix */
      ATrans(vw::math::Matrix<ElemT, 0, 0> const& m);


      /** Const accessor of the translation vector. */
      VectorT const& translation() const throw() {
        return m_p;
      }
      /** Accessor of the translation vector. */
      VectorT& translation() throw() {
        return m_p;
      }

      /** Const accessor of the (rotation) matrix. */
      MatrixT const& rotation() const throw() {
        return m_r;
      }
      /** Accessor of the (rotation) matrix. */
      MatrixT& rotation() throw() {
        return m_r;
      }

      /** Matrix N+1 style multiplication */
      ATrans operator*(ATrans const& rhs) const;
      /** In-place matrix N+1 style multiplication */
      void operator*=(ATrans const& rhs);
      /** Position vecotr multiplication */
      VectorT operator*(VectorT const& rhs) const;

      /** Conversion to matrix N+1 */
      operator vw::math::Matrix<ElemT, SizeN + 1, SizeN + 1>();

    protected:
      /** Position vector of the affine transform. */
      VectorT m_p;
      /** (Rotation) matrix of the affine transform. */
      MatrixT m_r;
    };

    template<class ElemT, int SizeN>
    inline
    ATrans<ElemT, SizeN> inverse(ATrans<ElemT, SizeN> const& atrans) throw()
    {
      typename ATrans<ElemT, SizeN>::MatrixT mI = math::inverse(atrans.rotation());
      return ATrans<ElemT, SizeN>(-(mI * atrans.translation()), mI);
    }

    template<class ElemT, int SizeN>
    std::ostream& operator<< (std::ostream& ostr, ATrans<ElemT, SizeN> const& rhs)
    {
      ostr << "xATrans(" << rhs.translation() << "), (" << rhs.rotation() << ")";
      return ostr;
    }


    template<class ElemT, int SizeN>
    inline
    ATrans<ElemT, SizeN>::ATrans(vw::math::Matrix<ElemT, SizeN + 1, SizeN + 1> const& m) :
        m_p(m(0, 3), m(1,3), m(2, 3))
    {
      for (int i = 0; i < SizeN; ++i) {
        for (int j = 0; j < SizeN; ++j) {
          m_r(i, j) = m(i, j);
        }
      }
    }

    template<class ElemT, int SizeN>
    inline
    ATrans<ElemT, SizeN>::ATrans(vw::math::Matrix<ElemT, 0,0> const& m) :
        m_p(m(0, 3), m(1,3), m(2, 3))
    {
      for (int i = 0; i < SizeN; ++i) {
        for (int j = 0; j < SizeN; ++j) {
          m_r(i, j) = m(i, j);
        }
      }
    }

    template<class ElemT, int SizeN>
    inline
    ATrans<ElemT, SizeN>
    ATrans<ElemT, SizeN>::operator*(ATrans const& rhs) const
    {
      return ATrans(m_r * rhs.m_p + m_p, m_r * rhs.m_r);
    }

    template<class ElemT, int SizeN>
    inline
    void
    ATrans<ElemT, SizeN>::operator*=(ATrans const& rhs)
    {
      *this = *this * rhs;
    }

    template<class ElemT, int SizeN>
    inline
    typename ATrans<ElemT, SizeN>::VectorT
    ATrans<ElemT, SizeN>::operator*(VectorT const& p) const
    {
      VectorT v = m_r * p;
      return v + m_p;
    }

    template<class ElemT, int SizeN>
    inline
    ATrans<ElemT, SizeN>::operator vw::math::Matrix<ElemT, SizeN+1, SizeN+1>()
    {
      vw::math::Matrix<ElemT, SizeN+1, SizeN+1> m;
      return m;
    }
  }

  typedef geometry::ATrans<double, 3> ATrans3;
  typedef geometry::ATrans<float, 3> ATrans3f;
}

#endif // vw_geometry_ATrans_h

