// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file BundleAdjust.h
/// 
/// Optimization classes for carrying out bundle adjustment of many
/// camera images.

#ifndef __VW_CAMERA_BUNDLE_ADJUST_H__
#define __VW_CAMERA_BUNDLE_ADJUST_H__

// Vision Workbench
#include <vw/Camera/ControlNetwork.h>
#include <vw/Math/Matrix.h>
#include <vw/Math/Vector.h>
#include <vw/Core/Log.h>

// Boost 
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION==103200
// Mapped matrix exist in 1.32, but Sparse Matrix does
#define boost_sparse_matrix boost::numeric::ublas::sparse_matrix
#define boost_sparse_vector boost::numeric::ublas::sparse_vector
#else
// Sparse Matrix was renamed Mapped Matrix in later editions
#define boost_sparse_matrix boost::numeric::ublas::mapped_matrix
#define boost_sparse_vector boost::numeric::ublas::mapped_vector
#endif

#include <string>

namespace vw {
namespace camera {
  
  // CRTP Base class for Bundle Adjustment functors.
  // 
  template <class ImplT, unsigned CameraParamsN, unsigned PointParamsN>
  struct BundleAdjustmentModelBase {

    static const unsigned camera_params_n = CameraParamsN;
    static const unsigned point_params_n = PointParamsN;

    /// \cond INTERNAL
    // Methods to access the derived type
    inline ImplT& impl() { return static_cast<ImplT&>(*this); }
    inline ImplT const& impl() const { return static_cast<ImplT const&>(*this); }
    /// \endcond

    virtual ~BundleAdjustmentModelBase() {}
    
    // Required access to camera
    virtual Vector2 operator() ( unsigned i, unsigned j,
				 Vector<double,camera_params_n> const& a_j,
				 Vector<double,point_params_n> const& b_j ) const = 0;

    // Approximate the jacobian for small variations in the a_j
    // parameters (camera parameters). 
    inline Matrix<double, 2, CameraParamsN> A_jacobian ( unsigned i, unsigned j,
                                                         Vector<double, CameraParamsN> const& a_j, 
                                                         Vector<double, PointParamsN> const& b_i ) {

      // Get nominal function value
      Vector2 h0 = impl()(i,j,a_j,b_i);

      // Jacobian is #outputs x #params
      Matrix<double, 2, CameraParamsN> J;

      // For each param dimension, add epsilon and re-evaluate h() to
      // get numerical derivative w.r.t. that parameter
      for ( unsigned n=0; n < CameraParamsN; ++n ){
        Vector<double, CameraParamsN> a_j_prime = a_j;

        // Variable step size, depending on parameter value
        double epsilon = 1e-7 + fabs(a_j(n)*1e-7);
        a_j_prime(n) += epsilon;

        // Evaluate function with this step and compute the derivative
        // w.r.t. parameter i
        Vector2 hi = impl()(i,j,a_j_prime, b_i);
        select_col(J,n) = (hi-h0)/epsilon;
      }

      return J;
    }

    // Approximate the jacobian for small variations in the b_i
    // parameters (3d point locations). 
    inline Matrix<double, 2, PointParamsN> B_jacobian ( unsigned i, unsigned j,
                                                        Vector<double, CameraParamsN> const& a_j, 
                                                        Vector<double, PointParamsN> const& b_i ) {
      // Get nominal function value
      Vector2 h0 = impl()(i,j,a_j, b_i);

      // Jacobian is #outputs x #params
      Matrix<double, 2, PointParamsN> J;

      // For each param dimension, add epsilon and re-evaluate h() to
      // get numerical derivative w.r.t. that parameter
      for ( unsigned n=0; n < PointParamsN; ++n ){
        Vector<double, PointParamsN> b_i_prime = b_i;

        // Variable step size, depending on parameter value
        double epsilon = 1e-7 + fabs(b_i(n)*1e-7);
        b_i_prime(n) += epsilon;

        // Evaluate function with this step and compute the derivative
        // w.r.t. parameter i
        Vector2 hi = impl()(i,j,a_j,b_i_prime);
        select_col(J,n) = (hi-h0)/epsilon;
      }
      return J;
    }

    // Report functions
    virtual std::string image_unit( void ) { return "pixels"; }
    virtual std::string camera_position_unit( void ) { return "meters"; }
    virtual std::string camera_pose_unit( void ) { return "radians"; }
    virtual std::string gcp_unit( void ) { return "meters"; }
    
    // Forced on the user to define
    virtual void image_errors(std::vector<double>&) = 0;
    virtual void camera_position_errors(std::vector<double>&) = 0;
    virtual void camera_pose_errors(std::vector<double>&) = 0;
    virtual void gcp_errors(std::vector<double>&) = 0;
    virtual boost::shared_ptr<ControlNetwork> control_network(void) = 0;
  };


  // Returns cholesky factor L, D in lower left hand corner and
  // diagonal modifies in place returns 1 if original matrix was
  // positive definite, 0 otherwise has verbose output if not positive
  // definite
  template<class T>
  inline unsigned mod_cholesky(Matrix<T>& M){
    
    unsigned n = M.rows();
    for(unsigned i = 0; i < n; i++){
      for(unsigned j = i; j < n; j++){
	T sum = M(i,j);
	for(int k = i-1; k >= 0; k--)
	  sum -= M(i,k) * M(j,k);
        
	if (i == j){
	  if (sum <= 0.0){
	    std::cout << " Not positive definite! " << "\n";
	    return 0; // not positive definite!
	  }
	  M(i,i) = sqrt(sum);
	}else
	  M(j,i) = sum/M(i,i);
      }
    }
    return 1;
  };
  


  // returns cholesky factor L, D in lower left hand corner and
  // diagonal modifies in place Assumes M is positive definite
  template<class T>
  inline void cholesky(Matrix<T>& M){
    unsigned n = M.rows();
    for (unsigned j = 0; j < n; j++){
      if(j > 0){
	submatrix(M, j, j, n-j, 1) -= submatrix(M, j, 0, n-j, j) * transpose(submatrix(M, j, 0, 1, j));
      }
      submatrix(M, j, j, n-j, 1)/= sqrt(M(j,j));
    }
  };

  
  // solves Lx = b
  // modifies in place 
  template<class T>
  inline void fsolve(Vector<T>& b, Matrix<T>& L){
    
    unsigned n = L.rows();
    
    b(0) /= L(0,0);
    for(unsigned i = 1; i < n; i++){
      Vector<double, 1> temp = submatrix(L, i, 0, 1, i) * subvector(b, 0, i);
      b(i) = (b(i) - temp(0))/L(i,i);
    }
    
  };
  
  // solves L^T y = x
  template<class T>
  inline void bsolve(Vector<T>& b, Matrix<T>& U){
    
    unsigned n = U.rows();
    
    b(n-1) /= U(n-1, n-1);
    for(unsigned i = b.size()-1; i>=1; i--){
      Vector<double, 1> temp = submatrix(U, i-1, i, 1, b.size()-i) * subvector(b, i, b.size()-i);
      b(i-1) = (b(i-1) - temp(0))/U(i-1,i-1);
    }
    
   };
  
  // solves Ax = b for positive definite A returns solution in b, and
  // cholesky decomp. of A in L returns 0 if A is not positive
  // definite
  template<class T>
  inline unsigned solve(Vector<T>&b, Matrix<T>& A){
     
     unsigned ret = mod_cholesky(A);
    
     if (ret == 0)
       return ret;

     fsolve(b, A);
     Matrix<T> U = transpose(A);
     bsolve(b, U);
     
     return ret;
   };
   

   // Replaces lower triangle of A with that of L^{-1} Replaces upper
   // triangle of A with zeros Returns 0 if A is not positive definite
   template<class T>
   inline unsigned chol_inverse(Matrix<T>& A){
     unsigned n = A.rows();
     unsigned ret = mod_cholesky(A);
     if (ret == 0)
       return ret;
     T sum;
     for(unsigned i=0; i<n; i++){
       A(i,i) = 1/A(i,i);
       for(unsigned j = i+1; j < n; j++){
	 sum = 0.0;
	 for(unsigned k = i; k < j; k++)
	   sum -= A(j,k) * A(k,i);
	 A(j,i) = sum/A(j,j);
	 A(i,j) = 0;	 
       }
     }
     return ret;
   };

  //------------------------------------------------------------------
  //                 Sparse Skyline Matrix 
  //
  // The hessian in the sparse LM algorithm takes the form of a Sparse
  // skyline matrix.  Under this sparseness condition, efficient
  // solutions to a linear system can be efficiently computed.  The
  // code below defines a class for a sparse skyline matrix as well as
  // methods for it's decomposition and for linear system solving.
  // -----------------------------------------------------------------

  /// An extremely simple sparse matrix class that wraps around a
  /// boost compessed_matrix<>, but keeps track of the first non-zero
  /// element in each row (i.e. the "skyline") for more efficient
  /// processing in some algorithms.
  template <class ElemT>
  class SparseSkylineMatrix {  
    typedef boost::numeric::ublas::compressed_matrix<ElemT> sparse_matrix_type;
    // Haven't decided yet which boost sparse matrix type is
    // fastest... this one might be good, too.
    //
    // typedef boost::numeric::ublas::generalized_vector_of_vector<ElemT, row_major, vector<coordinate_vector<ElemT> > > sparse_matrix_type;
    
    sparse_matrix_type m_matrix;
    std::vector<uint32> m_skyline;

  public:
    SparseSkylineMatrix(unsigned cols, unsigned rows) : 
      m_matrix(cols,rows), m_skyline(cols) {
      VW_ASSERT(cols == rows, ArgumentErr() << "SparseSkylineMatrix must be square and symmetric.\n");
      for (unsigned i = 0; i < cols; ++i)
        m_skyline[i] = i;
    }
  

    // Returns the "skyline" of the sparse, symmetric matrix.  That is,
    // each index i in the returned vector contains the index of the
    // first valid entry in row i (or equivelently, column i) of the
    // skyline matrix.
    const std::vector<uint32>& skyline() const { return m_skyline; }
  
    uint32 rows() const { return m_matrix.size1(); }
    uint32 cols() const { return m_matrix.size2(); }

    const ElemT* find_element (uint32 i, uint32 j) const { 
      if (j > i) 
        return m_matrix.find_element(j,i); 
      else
        return m_matrix.find_element(i,j); 
    }

    // Some boost sparse matrix types define this method...
    void push_back (uint32 i, uint32 j, ElemT const& val) {
      return m_matrix.push_back(i,j,val);
    }

    typename sparse_matrix_type::const_reference operator () (uint32 i, uint32 j) const {
      VW_DEBUG_ASSERT(i < 0 || i >= this->rows() || j < 0 || j >= this->cols(),
                      ArgumentErr() << "SparseSkylineMatrix: index " << i << " " << j << " out of bounds.");

      // Force symmetry by reflecting all points to the lower left
      // triangle.
      if (j > i) 
        return m_matrix(j,i);
      else
        return m_matrix(i,j);
    }
    
    typename sparse_matrix_type::reference operator () (uint32 i, uint32 j) {
      VW_DEBUG_ASSERT(i < 0 || i >= this->rows() || j < 0 || j >= this->cols(),
                      ArgumentErr() << "SparseSkylineMatrix: index " << i << " " << j << " out of bounds.");
      
      // Force symmetry by reflecting all points to the lower left
      // triangle.
      if (j > i) {
        uint32 temp = j; 
        j = i; 
        i = temp;
      }
      
      if (j < m_skyline[i])
        m_skyline[i] = j;
      return m_matrix(i,j);
    }
    
    // Handy for debugging...
    void print_sparse_structure() {
      std::cout << "SPARSE STRUCTURE: \n";
      for (unsigned i = 0; i < this->rows(); ++i) {
        for (unsigned j = 0; j < this->cols(); ++j) {
          const ElemT* e = this->find_element(i,j);
          if (e) 
            std::cout << "* ";
          else 
            std::cout << ". ";
        }
        std::cout << "\n";
      }
    }    

    // Handy for debugging...
    void print_matrix() {
      std::cout << "SPARSE SKYLINE MATRIX: \n";
      for (unsigned i = 0; i < this->rows(); ++i) {
        for (unsigned j = 0; j < this->cols(); ++j) {
          const ElemT* e = this->find_element(i,j);
          if (e) 
            std::cout << this->operator()(i,j) << "  ";
          else 
            std::cout << "0         ";
        }
        std::cout << "\n";
      }
    }    

  };

  //--------------------------------------------------------------------
  //        L*D*L^T Decompostion for Symmetric, Spares, Skyline Matrices
  //--------------------------------------------------------------------
  
  // Perform L*D*L^T decomposition on a sparse, skyline symmetric
  // semi-definite matrix.  WARNING: The results are stored in place,
  // so this operation destroys the previous contents of A.
  //
  // Once this operation is complete, the diagonal entries of A
  // contain the values from D, and the lower left block diagonal of A
  // contains L.  (The diagonal entries of L are always 1, so those
  // are assumed here...)
  template <class ElemT>
  void ldl_decomposition(SparseSkylineMatrix<ElemT>& A) {
    VW_ASSERT(A.cols() == A.rows(), ArgumentErr() << "ldl_decomposition: argument must be square and symmetric.\n");
    
    const std::vector<uint32>& skyline = A.skyline();
    
    for (unsigned j = 0; j < A.cols(); ++j) {

      // Compute v(1:j)
      std::vector<double> v(j+1);
      v[j] = A(j,j);
      for (unsigned i = skyline[j]; i < j; ++i) {
        v[i] = A(j,i)*A(i,i);
        v[j] -= A(j,i)*v[i];
      }
      
      // Store d(j) and compute L(j+1:n,j)
      A(j,j) = v[j];
      for (unsigned i = j+1; i < A.cols(); ++i) {
        double row_sum = 0;
        for (unsigned jj = skyline[i]; jj < j; ++jj) 
          row_sum += A(i,jj)*v[jj];
        if (j >= skyline[i])
          A(i,j) = ( A(i,j)-row_sum ) / v[j];
      }
    }
  }

  //--------------------------------------------------------------
  //            Solve Spare Skyline Linear System: Ax=b
  //--------------------------------------------------------------

  /// Perform L*D*L^T decomposition on a sparse skyline symmetric
  /// semi-definite matrix A (in place) and then solves an equation of the
  /// form Ax=b using forward and backward substitution.
  /// 
  /// WARNING: Modifies the contents of the matrix A.
  template <class ElemT, class VectorT>
  vw::Vector<double> sparse_solve(SparseSkylineMatrix<ElemT>& A, VectorT const& b) {
    VW_ASSERT(A.cols() == A.rows(), ArgumentErr() << "sparse_solve: matrix must be square and symmetric.\n");


    // Compute the L*D*L^T decomposition of A
    ldl_decomposition(A);
    
    
    const std::vector<uint32>& skyline = A.skyline();
    std::vector<uint32> inverse_skyline(skyline.size());
    
    // Construct the inverse skyline, which is used to optimize the final
    // back substitution step below.
    for (unsigned j = 0; j < inverse_skyline.size(); ++j) {
      inverse_skyline[j] = 0;
      for (int i = skyline.size()-1; i>=0; --i) {
        if (j < skyline[i]) 
          ++(inverse_skyline[j]);
        else
          break; // Break out of the inner loop
      }
    }
   
    // Forward Substitution Step ( L*x'=b )
    vw::Vector<double> x_prime(A.cols());
    for (unsigned i = 0; i < x_prime.size(); ++i) {
      double sum = 0;
      for (unsigned j = skyline[i]; j < i; ++j) 
        sum += A(i,j)*x_prime(j);
      x_prime(i) = b(i)-sum;
    }
    
    // Divide by D ( D*x''=x' )
    vw::Vector<double> x_doubleprime(A.cols());
    for (unsigned i = 0; i < x_doubleprime.size(); ++i) {
      x_doubleprime(i) = x_prime(i)/A(i,i);
    }

    // Back Substitution step ( L^T*x=x'' )
    vw::Vector<double> x(A.cols());
    for (int32 i = x.size()-1; i >= 0; --i) {
      double sum = 0;
      for (unsigned j = i+1; j < A.cols()-inverse_skyline[i]; ++j) 
        sum += A(j,i)*x(j);
      x(i) = x_doubleprime(i) - sum;
    }
    return x;
  }

  
  //----------------------------------------------------------------
  // Robust cost functions.  These cost function can help to reduce
  // the impact of outliers in the bundle adjustment.
  //----------------------------------------------------------------
  
  struct PseudoHuberError { 
    double m_b;
    double m_b_2;
    PseudoHuberError(double b) : m_b(b) {
      m_b_2 = b*b;
    }

    double operator() (double delta_norm) {
      return 2.0f * m_b_2* (sqrt(1.0f + delta_norm*delta_norm/m_b_2) - 1.0f);
    }

    std::string name_tag (void) const { return "PsudeoHuberError"; }
    double threshold(void) const { return m_b; }
  };

  struct HuberError { 
    double m_b;
    HuberError(double b) : m_b(b) {}

    double operator() (double delta_norm) {
      if (delta_norm < m_b)
        return delta_norm*delta_norm;
      else
        return 2*m_b*delta_norm - m_b*m_b;
    }

    std::string name_tag (void) const { return "HuberError"; }
    double threshold(void) const { return m_b; }
  };

  struct L1Error {
    double operator() (double delta_norm) { return fabs(delta_norm); }
    
    std::string name_tag (void) const { return "L1Error"; }
    double threshold(void) const { return 0.0; }
  };


  struct L2Error { 
    double operator() (double delta_norm) { return delta_norm*delta_norm; }

    std::string name_tag (void) const { return "L2Error"; }
    double threshold(void) const { return 0.0; }
  };

  // Our Cauchy Error is missing a sigma^2 at the front of the log,
  // this is okay. Now when increasing sigma the weighting on higher error
  // is relaxed, increasing sigma increases the tail thickness of our PDF.
  struct CauchyError { 
    double m_sigma;
    double m_sigma_2;
    CauchyError(double sigma) : m_sigma(sigma) {
      m_sigma_2 = m_sigma*m_sigma;
    }

    double operator() (double delta_norm) { 
      return log(1+delta_norm*delta_norm/m_sigma_2);
    }

    std::string name_tag (void) const { return "CauchyError"; }
    double threshold(void) const { return m_sigma; }
  };

  //--------------------------------------------------------------
  //                     BundleAdjustment
  //--------------------------------------------------------------

  template <class BundleAdjustModelT, class RobustCostT>
  class BundleAdjustment {   
    
  private:
    boost::shared_ptr<ControlNetwork> m_control_net;
    BundleAdjustModelT &m_model;
    unsigned m_control; //0 for new control, 1 for old control
    double m_lambda;
    double m_nu;
    double g_tol;
    double d_tol;

    int m_iterations;
    RobustCostT m_robust_cost_func;
    
  public:
    
    BundleAdjustment(BundleAdjustModelT &model, 
                     RobustCostT const& robust_cost_func) : 
      m_model(model), m_robust_cost_func(robust_cost_func) {

      m_iterations = 0;

      // Extracting a control network;
      m_control_net = m_model.control_network();

      // Set an initial value for lambda, lambda_low, and is.  These values vary as the
      // algorithm proceeds.
      m_lambda = 1e-3;
      m_control = 0; //use fast control (Fletcher), 1 = Traditional lambda / or * by 10
      m_nu = 2;
      g_tol = 1e-10;
      d_tol = 1e-10;

      // Compute the initial error

      // Here are some useful variable declarations that make the code
      // below more readable.
      unsigned num_cam_params = BundleAdjustModelT::camera_params_n;
      unsigned num_pt_params = BundleAdjustModelT::point_params_n;

      unsigned num_cameras = m_model.num_cameras();
      unsigned num_ground_control_points = m_control_net->num_ground_control_points();
      unsigned num_observations = 2*m_model.num_pixel_observations() +
                                  num_cameras*num_cam_params + 
                                  num_ground_control_points*num_pt_params;

      Vector<double> epsilon(num_observations);                   // Error vector

      int idx = 0;
      for (unsigned i = 0; i < m_control_net->size(); ++i) {       // Iterate over control points
        for (unsigned m = 0; m < (*m_control_net)[i].size(); ++m) {  // Iterate over control measures
          int camera_idx = (*m_control_net)[i][m].image_id();

          // Apply robust cost function weighting and populate the error vector
          Vector2 unweighted_error = (*m_control_net)[i][m].dominant() - m_model(i, camera_idx, 
                                                                              m_model.A_parameters(camera_idx), 
                                                                              m_model.B_parameters(i));
	  double mag = norm_2(unweighted_error);
          double weight = sqrt(m_robust_cost_func(mag)) / mag;
          subvector(epsilon,2*idx,2) = unweighted_error * weight;

          ++idx;
        }
      }

      // Add rows to J and epsilon for a priori position/pose constraints...
      for (unsigned j=0; j < num_cameras; ++j) {
        subvector(epsilon,
                  2*m_model.num_pixel_observations() + j*num_cam_params,
                  num_cam_params) = m_model.A_initial(j)-m_model.A_parameters(j);
      }

      // ... and the position of the 3D points to J and epsilon ...
      idx = 0;
      for (unsigned i=0; i < m_model.num_points(); ++i) {
        if ((*m_control_net)[i].type() == ControlPoint::GroundControlPoint) {
          subvector(epsilon,
                    2*m_model.num_pixel_observations() + num_cameras*num_cam_params + idx*num_pt_params,
                    num_pt_params) = m_model.B_initial(i)-m_model.B_parameters(i);
          ++idx;
        }
      }

      // Summarize the stats from this step in the iteration
      double overall_norm = transpose(epsilon) *  epsilon;
      
    }
  
    /// Set/Read Controls

    double lambda() const { return m_lambda; }
    void set_lambda(double lambda) { m_lambda = lambda; }
    unsigned control() const { return m_control; }
    void set_control(unsigned control){m_control = control;}
    
    /// Information Reporting

    int iterations() const { return m_iterations; }
    RobustCostT costfunction() const { return m_robust_cost_func; }

    /// Compare matrices methods (FOR DEBUG PURPOSES)

    void compare_matrices(SparseSkylineMatrix<double> &a,
                          Matrix<double> &b, 
                          std::string debug_text,
                          double tol) {
      std::cout << "Checking " << debug_text << "...\n";
      
      std::cout << "A : " << a.rows() << " x " << a.cols() << "    B : " << b.rows() << " x " << b.cols() << "\n";
      for (unsigned r = 0; r < a.rows(); ++r) {    // rows
        for (unsigned c = 0; c < a.cols(); ++c) {  // cols
          double s = a(r,c);
          double n = b(r,c);
          if ( fabs(s - n) > tol )
            std::cout << "Mismatch SPSKLN MATRIX at " << r << " " << c << " " << "  :  " << s << " vs. " << n << "   diff: " << (s-n) << "\n";
        }
      }
    }

    template <class ElemT>
    void compare_matrices(boost_sparse_matrix<ElemT> &sparse,
                          Matrix<double> normal, 
                          std::string debug_text,
                          double tol) {
      std::cout << "Checking " << debug_text << "...\n";
      unsigned s1 = sparse.size1();
      unsigned s2 = sparse.size2();
      
      std::cout << "Sparse size : " << s1 << " x " << s2 << "    Non-sparse size: " << normal.rows() << " x " << normal.cols() << "\n";
      for (unsigned r = 0; r < sparse.size1(); ++r) {    // rows
        for (unsigned c = 0; c < sparse.size2(); ++c) {  // cols
          int n_cols = sparse(r,c).ref().cols();
          int n_rows = sparse(r,c).ref().rows();
          
          for (int rr = 0; rr < n_rows; ++rr) {
            for (int cc = 0; cc < n_cols; ++cc) {
              double s = sparse(r,c).ref()(rr,cc);
              double n = normal(r*n_rows+rr, c*n_cols+cc);
              if ( fabs(s - n) > tol )
                std::cout << "Mismatch MATRIX at " << r << " " << c << " " << rr << " " << cc << " : " << s << " vs. " << n << "   diff: " << (s-n) << "\n";
            }
          }
        }
      }
    }

    // For comparing diagonal sparse matrices
    template <class ElemT>
    void compare_matrices(boost_sparse_vector<ElemT> &sparse,
                          Matrix<double> normal, 
                          std::string debug_text,
                          double tol) {
      std::cout << "Checking " << debug_text << "...\n";
      unsigned sz = sparse.size();
      std::cout << "Sparse size : " << sz << " x " << sz << "     Non-sparse size: " << normal.rows() << " x " << normal.cols() << "\n";

      for (unsigned c = 0; c < sparse.size(); ++c) { 
        int n_cols = sparse(c).ref().cols();
        int n_rows = sparse(c).ref().rows();
        
        for (int rr = 0; rr < n_rows; ++rr) {
          for (int cc = 0; cc < n_cols; ++cc) {
            double s = sparse(c).ref()(rr,cc);
            double n = normal(c*n_rows+rr, c*n_cols+cc);
            if ( fabs(s - n) > tol )
              std::cout << "Mismatch VECTOR at " << c << " " << rr << " " << cc << " : " << s << " vs. " << n << "    diff: " << (s-n) << "\n";
          }
        }
      }
    }

    // For comparing diagonal sparse matrices
    void compare_vectors(Vector<double> &a, Vector<double> &b, 
                         std::string debug_text, double tol) {
      std::cout << "Checking " << debug_text << "...\n";
      std::cout << "A size : " << a.size() << "     B size: " << b.size() << "\n";

      for (unsigned c = 0; c < a.size(); ++c) { 
        double s = a[c];
        double n = b[c];
        if ( fabs(s - n) > tol )
          std::cout << "Mismatch VW VECTOR at " << c << " : " << s << " vs. " << n << "    diff: " << (s-n) << "\n";
      }
    }


    void debug_impl(Matrix<double> &J, Matrix<double> &sigma, Vector<double> &epsilon) {

      // Here are some useful variable declarations that make the code
      // below more readable.
      unsigned num_cam_params = BundleAdjustModelT::camera_params_n;
      unsigned num_pt_params = BundleAdjustModelT::point_params_n;
      unsigned num_points = m_model.num_points();
      unsigned num_model_parameters = m_model.num_cameras()*num_cam_params + m_model.num_points()*num_pt_params;

      unsigned num_cameras = m_model.num_cameras();
      unsigned num_ground_control_points = m_control_net->num_ground_control_points();
      unsigned num_observations = 2*m_model.num_pixel_observations() +
                                  num_cameras*num_cam_params + 
                                  num_ground_control_points*num_pt_params;


      // --- SETUP STEP ----
      // Add rows to J and epsilon for the imaged pixel observations
      int idx = 0;
      for (unsigned i = 0; i < m_control_net->size(); ++i) {       // Iterate over control points
        for (unsigned m = 0; m < (*m_control_net)[i].size(); ++m) {  // Iterate over control measures
          int camera_idx = (*m_control_net)[i][m].image_id();

          Matrix<double> J_a = m_model.A_jacobian(i,camera_idx, 
                                                  m_model.A_parameters(camera_idx), 
                                                  m_model.B_parameters(i));
          Matrix<double> J_b = m_model.B_jacobian(i,camera_idx, 
                                                  m_model.A_parameters(camera_idx), 
                                                  m_model.B_parameters(i));

          // Populate the Jacobian Matrix
          submatrix(J, 2*idx, num_cam_params*camera_idx, 2, num_cam_params) = J_a;
          submatrix(J, 2*idx, num_cam_params*num_cameras + i*num_pt_params, 2, num_pt_params) = J_b;

          // Apply robust cost function weighting and populate the error vector
          Vector2 unweighted_error = (*m_control_net)[i][m].dominant() - m_model(i, camera_idx, 
                                                                              m_model.A_parameters(camera_idx), 
                                                                              m_model.B_parameters(i));
	  double mag = norm_2(unweighted_error);
          double weight = sqrt(m_robust_cost_func(mag)) / mag;
          subvector(epsilon,2*idx,2) = unweighted_error * weight;

          // Fill in the entries of the sigma matrix with the uncertainty of the observations.
          Matrix2x2 inverse_cov;
          Vector2 pixel_sigma = (*m_control_net)[i][m].sigma();
          inverse_cov(0,0) = 1/(pixel_sigma(0)*pixel_sigma(0));
          inverse_cov(1,1) = 1/(pixel_sigma(1)*pixel_sigma(1));
          submatrix(sigma, 2*idx, 2*idx, 2, 2) = inverse_cov;

          ++idx;
        }
      }
      
      // Add rows to J and epsilon for a priori camera parameters...
      for (unsigned j=0; j < num_cameras; ++j) {
        Matrix<double> id(num_cam_params, num_cam_params);
        id.set_identity();
        submatrix(J,
                  2*m_model.num_pixel_observations() + j*num_cam_params,
                  j*num_cam_params,
                  num_cam_params,
                  num_cam_params) = id;
        Vector<double> unweighted_error = m_model.A_initial(j)-m_model.A_parameters(j);
        subvector(epsilon,
                  2*m_model.num_pixel_observations() + j*num_cam_params,
                  num_cam_params) = unweighted_error;
        submatrix(sigma,
                  2*m_model.num_pixel_observations() + j*num_cam_params,
                  2*m_model.num_pixel_observations() + j*num_cam_params,
                  num_cam_params, num_cam_params) = m_model.A_inverse_covariance(j);
      }

      // ... and the position of the 3D points to J and epsilon ...
      idx = 0;
      for (unsigned i=0; i < m_model.num_points(); ++i) {
        if ((*m_control_net)[i].type() == ControlPoint::GroundControlPoint) {
          Matrix<double> id(num_pt_params,num_pt_params);
          id.set_identity();
          submatrix(J,2*m_model.num_pixel_observations() + num_cameras*num_cam_params + idx*num_pt_params,
                    num_cameras*num_cam_params + idx*num_pt_params,
                    num_pt_params,
                    num_pt_params) = id;
          Vector<double> unweighted_error = m_model.B_initial(i)-m_model.B_parameters(i);
          subvector(epsilon,
                    2*m_model.num_pixel_observations() + num_cameras*num_cam_params + idx*num_pt_params,
                    num_pt_params) = unweighted_error;
          submatrix(sigma,
                    2*m_model.num_pixel_observations() + num_cameras*num_cam_params + idx*num_pt_params,
                    2*m_model.num_pixel_observations() + num_cameras*num_cam_params + idx*num_pt_params,
                    num_pt_params, num_pt_params) = m_model.B_inverse_covariance(i);
          ++idx;
        }
      }
    }


    //----------------------------------------------------------------
    //              Levenberg Marquardt Update Methods
    //----------------------------------------------------------------


    // This is a simple, non-sparse, unoptimized implementation of LM
    // bundle adjustment.  It is primarily used for validation and
    // debugging.
    //
    // Each entry in the outer vector corresponds to a distinct 3D
    // point.  The inner vector contains a list of image IDs and
    // pixel coordinates where that point was imaged.
    double update_reference_impl(double &abs_tol, double &rel_tol) { 
      ++m_iterations;

      // Here are some useful variable declarations that make the code
      // below more readable.
      unsigned num_cam_params = BundleAdjustModelT::camera_params_n;
      unsigned num_pt_params = BundleAdjustModelT::point_params_n;
      unsigned num_points = m_model.num_points();
      unsigned num_model_parameters = m_model.num_cameras()*num_cam_params + m_model.num_points()*num_pt_params;

      unsigned num_cameras = m_model.num_cameras();
      unsigned num_ground_control_points = m_control_net->num_ground_control_points();
      unsigned num_observations = 2*m_model.num_pixel_observations() +
                                  num_cameras*num_cam_params + 
                                  num_ground_control_points*num_pt_params;

      // The core LM matrices and vectors
      Matrix<double> J(num_observations, num_model_parameters);   // Jacobian Matrix
      Vector<double> epsilon(num_observations);                   // Error vector
      Matrix<double> sigma(num_observations, num_observations);   // Sigma (uncertainty) matrix

      // --- SETUP STEP ----
      // Add rows to J and epsilon for the imaged pixel observations
      int idx = 0;
      for (unsigned i = 0; i < m_control_net->size(); ++i) {       // Iterate over control points
        for (unsigned m = 0; m < (*m_control_net)[i].size(); ++m) {  // Iterate over control measures
          int camera_idx = (*m_control_net)[i][m].image_id();

          Matrix<double> J_a = m_model.A_jacobian(i,camera_idx, 
                                                  m_model.A_parameters(camera_idx), 
                                                  m_model.B_parameters(i));
          Matrix<double> J_b = m_model.B_jacobian(i,camera_idx, 
                                                  m_model.A_parameters(camera_idx), 
                                                  m_model.B_parameters(i));

          // Populate the Jacobian Matrix
          submatrix(J, 2*idx, num_cam_params*camera_idx, 2, num_cam_params) = J_a;
          submatrix(J, 2*idx, num_cam_params*num_cameras + i*num_pt_params, 2, num_pt_params) = J_b;

          // Apply robust cost function weighting and populate the error vector
          Vector2 unweighted_error = (*m_control_net)[i][m].dominant() - m_model(i, camera_idx, 
                                                                              m_model.A_parameters(camera_idx), 
                                                                              m_model.B_parameters(i));
	  double mag = norm_2(unweighted_error);
          double weight = sqrt(m_robust_cost_func(mag)) / mag;
          subvector(epsilon,2*idx,2) = unweighted_error * weight;

          // Fill in the entries of the sigma matrix with the uncertainty of the observations.
          Matrix2x2 inverse_cov;
          Vector2 pixel_sigma = (*m_control_net)[i][m].sigma();
          inverse_cov(0,0) = 1/(pixel_sigma(0)*pixel_sigma(0));
          inverse_cov(1,1) = 1/(pixel_sigma(1)*pixel_sigma(1));
          submatrix(sigma, 2*idx, 2*idx, 2, 2) = inverse_cov;

          ++idx;
        }
      }
      
      // Add rows to J and epsilon for a priori camera parameters...
      for (unsigned j=0; j < num_cameras; ++j) {
        Matrix<double> id(num_cam_params, num_cam_params);
        id.set_identity();
        submatrix(J,
                  2*m_model.num_pixel_observations() + j*num_cam_params,
                  j*num_cam_params,
                  num_cam_params,
                  num_cam_params) = id;
        Vector<double> unweighted_error = m_model.A_initial(j)-m_model.A_parameters(j);
        subvector(epsilon,
                  2*m_model.num_pixel_observations() + j*num_cam_params,
                  num_cam_params) = unweighted_error;
        submatrix(sigma,
                  2*m_model.num_pixel_observations() + j*num_cam_params,
                  2*m_model.num_pixel_observations() + j*num_cam_params,
                  num_cam_params, num_cam_params) = m_model.A_inverse_covariance(j);
      }

      // ... and the position of the 3D points to J and epsilon ...
      idx = 0;
      for (unsigned i=0; i < m_model.num_points(); ++i) {
        if ((*m_control_net)[i].type() == ControlPoint::GroundControlPoint) {
          Matrix<double> id(num_pt_params,num_pt_params);
          id.set_identity();
          submatrix(J,2*m_model.num_pixel_observations() + num_cameras*num_cam_params + idx*num_pt_params,
                    num_cameras*num_cam_params + idx*num_pt_params,
                    num_pt_params,
                    num_pt_params) = id;
          Vector<double> unweighted_error = m_model.B_initial(i)-m_model.B_parameters(i);
          subvector(epsilon,
                    2*m_model.num_pixel_observations() + num_cameras*num_cam_params + idx*num_pt_params,
                    num_pt_params) = unweighted_error;
          submatrix(sigma,
                    2*m_model.num_pixel_observations() + num_cameras*num_cam_params + idx*num_pt_params,
                    2*m_model.num_pixel_observations() + num_cameras*num_cam_params + idx*num_pt_params,
                    num_pt_params, num_pt_params) = m_model.B_inverse_covariance(i);
          ++idx;
        }
      }

      // --- UPDATE STEP ----

      // Build up the right side of the normal equation...
      Vector<double> del_J = -1.0 * (transpose(J) * sigma * epsilon);

      // ... and the left side.  (Remembering to rescale the diagonal
      // entries of the approximated hessian by lambda)
      Matrix<double> hessian = transpose(J) * sigma * J;

      // initialize m_lambda on first iteration, ignore if user has
      // changed it.
      double max = 0.0;
      if (m_iterations == 1 && m_lambda == 1e-3){
	for (unsigned i = 0; i < hessian.rows(); ++i){
	  if (fabs(hessian(i,i)) > max)
	    max = fabs(hessian(i,i));
	}
	m_lambda = 1e-10 * max;
      }

      for ( unsigned i=0; i < hessian.rows(); ++i )
        hessian(i,i) +=  m_lambda;

      //Cholesky decomposition. Returns Cholesky matrix in lower left hand corner.
      Vector<double> delta = del_J;

      // Here we want to make sure that if we apply Schur methods as on p. 604, we can get the same answer as in the general delta.  
      unsigned num_cam_entries = num_cam_params * num_cameras;
      unsigned num_pt_entries = num_pt_params * num_points;

      Matrix<double> U = submatrix(hessian, 0, 0, num_cam_entries, num_cam_entries);
      Matrix<double> W = submatrix(hessian, 0, num_cam_entries,  num_cam_entries, num_pt_entries);
      Matrix<double> Vinv = submatrix(hessian, num_cam_entries, num_cam_entries, num_pt_entries, num_pt_entries); 
      chol_inverse(Vinv);
      Matrix<double> Y = W * transpose(Vinv) * Vinv;
      Vector<double> e = subvector(delta, 0, num_cam_entries) - W * transpose(Vinv) * Vinv * subvector(delta, num_cam_entries, num_pt_entries); 
      Matrix<double> S = U - Y * transpose(W);
      solve(e, S); // using cholesky
      solve(delta, hessian);
      
      double nsq_x = 0;
      for (unsigned j=0; j<m_model.num_cameras(); ++j){
	Vector<double> vec = m_model.A_parameters(j);
	nsq_x += norm_2(vec);
      }
      for (unsigned i=0; i<m_model.num_points(); ++i){
	Vector<double> vec =  m_model.B_parameters(i);
	nsq_x += norm_2(vec);
      }


      // Solve for update

      // --- EVALUATE UPDATE STEP ---
      Vector<double> new_epsilon(num_observations);                  // Error vector

      idx = 0;
      for (unsigned i = 0; i < m_control_net->size(); ++i) {         // Iterate over control points
        for (unsigned m = 0; m < (*m_control_net)[i].size(); ++m) {  // Iterate over control measures
          int camera_idx = (*m_control_net)[i][m].image_id();

          Vector<double> cam_delta = subvector(delta, num_cam_params*camera_idx, num_cam_params);
          Vector<double> pt_delta = subvector(delta, num_cam_params*num_cameras + num_pt_params*i, num_pt_params);

          // Apply robust cost function weighting and populate the error vector
          Vector2 unweighted_error = (*m_control_net)[i][m].dominant() - m_model(i, camera_idx,
                                                                              m_model.A_parameters(camera_idx)-cam_delta, 
                                                                              m_model.B_parameters(i)-pt_delta);
	  double mag = norm_2(unweighted_error);
          double weight = sqrt(m_robust_cost_func(mag)) / mag;
          subvector(new_epsilon,2*idx,2) = unweighted_error * weight;

          ++idx;
        }
      }

      // Add rows to J and epsilon for a priori position/pose constraints...
      for (unsigned j=0; j < num_cameras; ++j) {
        Vector<double> cam_delta = subvector(delta, num_cam_params*j, num_cam_params);
        subvector(new_epsilon,
                  2*m_model.num_pixel_observations() + j*num_cam_params,
                  num_cam_params) = m_model.A_initial(j)-(m_model.A_parameters(j)-cam_delta);
      }

      // ... and the position of the 3D points to J and epsilon ...
      idx = 0;
      for (unsigned i=0; i < m_model.num_points(); ++i) {
        if ((*m_control_net)[i].type() == ControlPoint::GroundControlPoint) {
          Vector<double> pt_delta = subvector(delta, num_cam_params*num_cameras + num_pt_params*i, num_pt_params);
          subvector(new_epsilon,
                    2*m_model.num_pixel_observations() + num_cameras*num_cam_params + idx*num_pt_params,
                    num_pt_params) = m_model.B_initial(i)-(m_model.B_parameters(i)-pt_delta);
          ++idx;
        }
      }

      //Fletcher modification
      double Splus = .5*transpose(new_epsilon) * sigma * new_epsilon; //Compute new objective
      double SS = .5*transpose(epsilon) * sigma * epsilon;            //Compute old objective
      double dS = .5 * transpose(delta) *(m_lambda * delta + del_J);

      // WARNING: will want to replace dS later

      double R = (SS - Splus)/dS; // Compute ratio
      unsigned ret = 0;
      
      if (R > 0){
	ret = 1;	
	for (unsigned j=0; j<m_model.num_cameras(); ++j) 
          m_model.set_A_parameters(j, m_model.A_parameters(j) - subvector(delta, num_cam_params*j, num_cam_params));
        for (unsigned i=0; i<m_model.num_points(); ++i)
          m_model.set_B_parameters(i, m_model.B_parameters(i) - subvector(delta, num_cam_params*num_cameras + num_pt_params*i, num_pt_params));
	
	// Summarize the stats from this step in the iteration
        double overall_norm = sqrt(.5 * transpose(new_epsilon) * sigma * new_epsilon);
	
        double overall_delta = sqrt(.5 * transpose(epsilon) * sigma * epsilon) - sqrt(.5 * transpose(new_epsilon) * sigma * new_epsilon) ; 
	
	abs_tol = overall_norm;
        rel_tol = overall_delta;	

	if (m_control==0){
	  double temp = 1 - pow((2*R - 1),3);		
	  if (temp < 1.0/3.0)
	    temp = 1.0/3.0;
	  
	  m_lambda *= temp;
	  m_nu = 2;
	} else if (m_control == 1)
	  m_lambda /= 10;
	
	return overall_delta;
      
      } else { // R <= 0
	
	if (m_control == 0){
	  m_lambda *= m_nu; 
	  m_nu*=2; 
	} else if (m_control == 1)
	  m_lambda *= 10;

	double overall_delta = sqrt(.5 * transpose(epsilon) * sigma * epsilon) - sqrt(.5 * transpose(new_epsilon) * sigma * new_epsilon);
 
	return ScalarTypeLimits<double>::highest();
      }
    }


      /////////////////////////////////////////////////////////////
     /////                 SPARSE LM IMPLEMENTATION    ///////////
    /////////////////////////////////////////////////////////////
  
    // This is the sparse levenberg marquardt update step.  Returns
    // the average improvement in the cost function.
    double update(double &abs_tol, double &rel_tol) { 
      ++m_iterations;
      
      VW_DEBUG_ASSERT(m_control_net->size() == m_model.num_points(), LogicErr() << "BundleAdjustment::update() : Number of bundles does not match the number of points in the bundle adjustment model.");

      // Jacobian Matrices and error values
      boost_sparse_matrix<Matrix<double, 2, BundleAdjustModelT::camera_params_n> > A(m_model.num_points(), m_model.num_cameras());
      boost_sparse_matrix<Matrix<double, 2, BundleAdjustModelT::point_params_n> > B(m_model.num_points(), m_model.num_cameras());
      boost_sparse_matrix<Vector2> epsilon(m_model.num_points(), m_model.num_cameras());
      boost_sparse_matrix<Vector2> new_epsilon(m_model.num_points(), m_model.num_cameras());
      
      // Data structures necessary for Fletcher modification
      // boost::numeric::ublas::mapped_matrix<Vector2> Jp(m_model.num_points(), m_model.num_cameras());
           
      // Intermediate Matrices and vectors
      boost_sparse_vector< Matrix<double,BundleAdjustModelT::camera_params_n,BundleAdjustModelT::camera_params_n> > U(m_model.num_cameras());
      boost_sparse_vector< Matrix<double,BundleAdjustModelT::point_params_n,BundleAdjustModelT::point_params_n> > V(m_model.num_points());  
      boost_sparse_matrix<Matrix<double,BundleAdjustModelT::camera_params_n,BundleAdjustModelT::point_params_n> > W(m_model.num_cameras(), m_model.num_points());     
     
      // Copies of Intermediate Marices
      boost_sparse_vector< Vector<double,BundleAdjustModelT::camera_params_n> > epsilon_a(m_model.num_cameras());
      boost_sparse_vector< Vector<double,BundleAdjustModelT::point_params_n > > epsilon_b(m_model.num_points());
      boost_sparse_matrix<Matrix<double,BundleAdjustModelT::camera_params_n,BundleAdjustModelT::point_params_n> > Y(m_model.num_cameras(), m_model.num_points());


      unsigned num_cam_params = BundleAdjustModelT::camera_params_n;
      unsigned num_pt_params = BundleAdjustModelT::point_params_n;

      unsigned delta_length = U.size() * num_cam_params + V.size() * num_pt_params;
      
      Vector<double> g(delta_length);  
      unsigned current_g_length = 0;
      
      Vector<double> delta(delta_length);
      unsigned current_delta_length = 0;

      // Fletcher LM parameteres
      double dS = 0; //Predicted improvement for Fletcher modification
      
      // Populate the Jacobian, which is broken into two sparse
      // matrices A & B, as well as the error matrix and the W 
      // matrix.
      vw_out(DebugMessage, "bundle_adjustment") << "Image Error: " << std::endl;
      unsigned i = 0;
      double error_total = 0; // assume this is r^T\Sigma^{-1}r
      for (typename ControlNetwork::const_iterator iter = m_control_net->begin(); iter != m_control_net->end(); ++iter) {
        for (typename ControlPoint::const_iterator measure_iter = (*iter).begin(); measure_iter != (*iter).end(); ++measure_iter) {

          unsigned j = (*measure_iter).image_id();
          VW_DEBUG_ASSERT(j >=0 && j < m_model.num_cameras(), ArgumentErr() << "BundleAdjustment::update() : image index out of bounds.");
	  
          // Store jacobian values
          A(i,j) = m_model.A_jacobian(i,j,m_model.A_parameters(j),m_model.B_parameters(i));
          B(i,j) = m_model.B_jacobian(i,j,m_model.A_parameters(j),m_model.B_parameters(i));

          // Apply robust cost function weighting
          Vector2 unweighted_error = measure_iter->dominant() - m_model(i,j,m_model.A_parameters(j),m_model.B_parameters(i));
	  double mag = norm_2(unweighted_error);
          double weight = sqrt(m_robust_cost_func(mag)) / mag;
          epsilon(i,j) = unweighted_error * weight;
	 	            
	  Matrix2x2 inverse_cov;
          Vector2 pixel_sigma = measure_iter->sigma();
          inverse_cov(0,0) = 1/(pixel_sigma(0)*pixel_sigma(0));
          inverse_cov(1,1) = 1/(pixel_sigma(1)*pixel_sigma(1));
	 
	  error_total += .5 * transpose(static_cast<Vector2>(epsilon(i,j))) * 
	    inverse_cov * static_cast<Vector2>(epsilon(i,j));
          
          // Store intermediate values
          U(j) += transpose(A(i,j).ref()) * inverse_cov * A(i,j).ref();
          V(i) += transpose(B(i,j).ref()) * inverse_cov * B(i,j).ref();
          W(j,i) = transpose(A(i,j).ref()) * inverse_cov * B(i,j).ref();

	  epsilon_a(j) += transpose(A(i,j)) * inverse_cov * epsilon(i,j);
          epsilon_b(i) += transpose(B(i,j)) * inverse_cov * epsilon(i,j);

	  // If GCP debug?
	  if ((*iter).type() == ControlPoint::GroundControlPoint) {
	    vw_out(DebugMessage, "bundle_adjustment") << "\t>" << i << " " << j << " error: " << unweighted_error << std::endl;
	  }
	}
        ++i;
      }

      // Add in the camera position and pose constraint terms and covariances.
      vw_out(DebugMessage, "bundle_adjustment") << "Camera Constraint Error:" << std::endl;
      for (unsigned j = 0; j < U.size(); ++j) {
        
	Matrix<double,BundleAdjustModelT::camera_params_n,BundleAdjustModelT::camera_params_n> inverse_cov;
        inverse_cov = m_model.A_inverse_covariance(j);

        Matrix<double, BundleAdjustModelT::camera_params_n, BundleAdjustModelT::camera_params_n> C;
        C.set_identity();
        U(j) += transpose(C) * inverse_cov * C;

        Vector<double, BundleAdjustModelT::camera_params_n> eps_a = m_model.A_initial(j)-m_model.A_parameters(j);
 	
	error_total += .5  * transpose(eps_a) * inverse_cov * eps_a;
	
	epsilon_a(j) += transpose(C) * inverse_cov * eps_a;
        
        // For debugging
	vw_out(DebugMessage, "bundle_adjustment") << "\t>" << j << " error: " << eps_a << std::endl;
      }

      // Add in the 3D point position constraint terms and
      // covariances. We only add constraints for Ground Control
      // Points (GCPs), not for 3D tie points.
      vw_out(DebugMessage, "bundle_adjustment") << "GCP Error: " << std::endl;
      for (unsigned i = 0; i < V.size(); ++i) {
        if ((*m_control_net)[i].type() == ControlPoint::GroundControlPoint) {

          Matrix<double,BundleAdjustModelT::point_params_n,BundleAdjustModelT::point_params_n> inverse_cov;
          inverse_cov = m_model.B_inverse_covariance(i);
          
          Matrix<double, BundleAdjustModelT::point_params_n, BundleAdjustModelT::point_params_n> D;
          D.set_identity();
          V(i) +=  transpose(D) * inverse_cov * D;
          
          Vector<double, BundleAdjustModelT::point_params_n> eps_b = m_model.B_initial(i)-m_model.B_parameters(i);
	  
	  error_total += .5 * transpose(eps_b) * inverse_cov * eps_b;

	  epsilon_b(i) += transpose(D) * inverse_cov * eps_b;

          // For debugging
	  vw_out(DebugMessage, "bundle_adjustment") << "\t>" << i << " error: " << eps_b << std::endl;
        }
      }

      
      // flatten both epsilon_b and epsilon_a into a vector
      for (unsigned j = 0; j < U.size(); j++){
	subvector(g, current_g_length, num_cam_params) = epsilon_a(j).ref();
	current_g_length += num_cam_params;
      }

      for (unsigned i = 0; i < V.size(); i++){
	subvector(g, current_g_length, num_pt_params) = epsilon_b(i).ref();
	current_g_length += num_pt_params;
      }


      //e at this point should be -g_a
      
      // set initial lambda, and ignore if the user has touched it
      if (m_iterations == 1 && m_lambda == 1e-3){

	double max = 0.0;
	for (unsigned i = 0; i < U.size(); ++i)  
	  for (unsigned j = 0; j < BundleAdjustModelT::camera_params_n; ++j){
	    if (fabs(U(i).ref()(j,j)) > max)
	      max = fabs(U(i).ref()(j,j));
	  }
	for (unsigned i = 0; i < V.size(); ++i) 
	  for (unsigned j = 0; j < BundleAdjustModelT::point_params_n; ++j) {
	    if ( fabs(V(i).ref()(j,j)) > max)
	      max = fabs(V(i).ref()(j,j));
	  }
	m_lambda = max * 1e-10;
      }
      
     
      // "Augment" the diagonal entries of the U and V matrices with
      // the parameter lambda.
      for (i = 0; i < U.size(); ++i) {
        for (unsigned j = 0; j < BundleAdjustModelT::camera_params_n; ++j) {
	  U(i).ref()(j,j) += m_lambda;
	   //U(i).ref()(j,j) *= (1 + m_lambda); 
	}
      }
      for (i = 0; i < V.size(); ++i) {
        for (unsigned j = 0; j < BundleAdjustModelT::point_params_n; ++j) {
	  V(i).ref()(j,j) += m_lambda;
	  //V(i).ref()(j,j) *= (1 + m_lambda);
	}
      }
      
      // Create the 'e' vector in S * delta_a = e.  The first step is
      // to "flatten" our block structure to a vector that contains
      // scalar entries.
      Vector<double> e(m_model.num_cameras() * BundleAdjustModelT::camera_params_n);
      for (unsigned j = 0; j < epsilon_a.size(); ++j) {
        subvector(e, j*BundleAdjustModelT::camera_params_n, BundleAdjustModelT::camera_params_n) = epsilon_a(j).ref();
      }
               
      //Second Pass.  Compute Y and finish constructing e.
      i = 0; 
      for (typename ControlNetwork::const_iterator iter = m_control_net->begin(); iter != m_control_net->end(); ++iter) { 
        for (typename ControlPoint::const_iterator measure_iter = (*iter).begin(); measure_iter != (*iter).end(); ++measure_iter) { 
          unsigned j = measure_iter->image_id(); 
	 
          // Compute the blocks of Y 
 	  Matrix<double> V_temp = V(i).ref(); 
 	  chol_inverse(V_temp); 
 	  Y(j,i) = W(j,i).ref() * transpose(V_temp) * V_temp; 

          // "Flatten the block structure to compute 'e'. 
          Vector<double, BundleAdjustModelT::camera_params_n> temp = Y(j,i).ref()*epsilon_b(i).ref(); 
 	  subvector(e, j*BundleAdjustModelT::camera_params_n, BundleAdjustModelT::camera_params_n) -= temp;
        } 
        ++i; 
      } 
      
      // The S matrix is a m x m block matrix with blocks that are
      // camera_params_n x camera_params_n in size.  It has a sparse
      // skyline structure, which makes it more efficient to solve
      // through L*D*L^T decomposition and forward/back substitution
      // below.
      SparseSkylineMatrix<double> S(m_model.num_cameras()*BundleAdjustModelT::camera_params_n, 
                                    m_model.num_cameras()*BundleAdjustModelT::camera_params_n);
      
      i = 0;
      for (typename ControlNetwork::const_iterator iter = m_control_net->begin(); iter != m_control_net->end(); ++iter) {
        for (typename ControlPoint::const_iterator j_measure_iter = (*iter).begin(); j_measure_iter != (*iter).end(); ++j_measure_iter) {
          unsigned j = j_measure_iter->image_id();
          
          for (typename ControlPoint::const_iterator k_measure_iter = (*iter).begin(); k_measure_iter != (*iter).end(); ++k_measure_iter) {
            unsigned k = k_measure_iter->image_id();
            
            // Compute the block entry...
            Matrix<double, BundleAdjustModelT::camera_params_n, BundleAdjustModelT::camera_params_n> temp = -Y(j,i).ref() * transpose( W(k,i).ref() );
            // ... and "flatten" this matrix into the scalar entries of S
            for (unsigned aa = 0; aa < BundleAdjustModelT::camera_params_n; ++aa) {
              for (unsigned bb = 0; bb < BundleAdjustModelT::camera_params_n; ++bb) {
                // FIXME: This if clause is required at the moment to
                // ensure that we do not use the += on the symmetric
                // entries of the SparseSkylineMatrix.  These
                // symmetric entries are shallow, hence this code
                // would add the value twice if we're not careful
                // here.
                if (k*BundleAdjustModelT::camera_params_n + bb <=
                    j*BundleAdjustModelT::camera_params_n + aa) {
                  S(j*BundleAdjustModelT::camera_params_n + aa,
                    k*BundleAdjustModelT::camera_params_n + bb) += temp(aa,bb);

		  //  S_old(j*BundleAdjustModelT::camera_params_n + aa,
                  //  k*BundleAdjustModelT::camera_params_n + bb) += temp_old(aa,bb);
		  
                }
              }
            }
          }
        }
        ++i;
      }
      
      // Augment the diagonal entries S(i,i) with U(i)
      for (unsigned i = 0; i < m_model.num_cameras(); ++i) {
        // ... and "flatten" this matrix into the scalar entries of S
        for (unsigned aa = 0; aa < BundleAdjustModelT::camera_params_n; ++aa) {
          for (unsigned bb = 0; bb < BundleAdjustModelT::camera_params_n; ++bb) {
            // FIXME: This if clause is required at the moment to
            // ensure that we do not use the += on the symmetric
            // entries of the SparseSkylineMatrix.  These
            // symmetric entries are shallow, hence this code
            // would add the value twice if we're not careful
            // here.
            if (i*BundleAdjustModelT::camera_params_n + bb <= 
                i*BundleAdjustModelT::camera_params_n + aa) {
              S(i*BundleAdjustModelT::camera_params_n + aa,
                i*BundleAdjustModelT::camera_params_n + bb) += U(i).ref()(aa,bb);
	    }
          }
        }
      } 

      // Compute the LDL^T decomposition and solve using sparse methods.
      Vector<double> delta_a = sparse_solve(S, e);

      subvector(delta, current_delta_length, e.size()) = delta_a;
      current_delta_length += e.size();

      boost_sparse_vector<Vector<double, BundleAdjustModelT::point_params_n> > delta_b(m_model.num_points());
      boost_sparse_vector<Vector<double, BundleAdjustModelT::camera_params_n> > delta_a_aux(m_model.num_cameras());

      i = 0;
      for (typename ControlNetwork::const_iterator iter = m_control_net->begin(); iter != m_control_net->end(); ++iter) {
        Vector<double, BundleAdjustModelT::point_params_n> temp;
        for (typename ControlPoint::const_iterator j_measure_iter = (*iter).begin(); j_measure_iter != (*iter).end(); ++j_measure_iter) {
          unsigned j = j_measure_iter->image_id();
          delta_a_aux(j) =  subvector(delta_a, j*BundleAdjustModelT::camera_params_n, BundleAdjustModelT::camera_params_n);

	  temp += transpose( W(j,i).ref() ) * delta_a_aux(j).ref();
	}

	Vector<double> delta_temp = epsilon_b(i).ref() - temp;
	
	Matrix<double> hessian = V(i).ref();
		
	solve(delta_temp, hessian);
	delta_b(i) = delta_temp;

	subvector(delta, current_delta_length, num_pt_params) = delta_temp;
	current_delta_length += num_pt_params;

        ++i;
      }
      
      i = 0;
      double nsq_x = 0.0;
      for (typename ControlNetwork::const_iterator iter = m_control_net->begin(); iter != m_control_net->end(); ++iter) {
        for (typename ControlPoint::const_iterator measure_iter = (*iter).begin(); measure_iter != (*iter).end(); ++measure_iter) {
          unsigned j = measure_iter->image_id();
	  nsq_x += norm_2( m_model.A_parameters(j));
	}nsq_x += norm_2(m_model.B_parameters(i));
      }
      
      dS = .5 * transpose(delta) *(m_lambda * delta + g);

      // -------------------------------
      // Compute the update error vector and predicted change
      // -------------------------------
      i = 0;
      double new_error_total = 0;
      for (typename ControlNetwork::const_iterator iter = m_control_net->begin(); iter != m_control_net->end(); ++iter) {
        for (typename ControlPoint::const_iterator measure_iter = (*iter).begin(); measure_iter != (*iter).end(); ++measure_iter) {

          unsigned j = measure_iter->image_id();
          
          // Compute error vector
          Vector<double, BundleAdjustModelT::camera_params_n> new_a = m_model.A_parameters(j) + subvector(delta_a, BundleAdjustModelT::camera_params_n*j, BundleAdjustModelT::camera_params_n);
	  Vector<double> del_a = subvector(delta_a, BundleAdjustModelT::camera_params_n*j, BundleAdjustModelT::camera_params_n);

          Vector<double, BundleAdjustModelT::point_params_n> new_b = m_model.B_parameters(i) + delta_b(i).ref();

	 
          // Apply robust cost function weighting
          Vector2 unweighted_error = measure_iter->dominant() - m_model(i,j,new_a,new_b);
	  double mag = norm_2(unweighted_error);
          double weight = sqrt(m_robust_cost_func(mag)) / mag;
          new_epsilon(i,j) = weight * unweighted_error;
          
	  Matrix2x2 inverse_cov;
          Vector2 pixel_sigma = measure_iter->sigma();
          inverse_cov(0,0) = 1/(pixel_sigma(0)*pixel_sigma(0));
          inverse_cov(1,1) = 1/(pixel_sigma(1)*pixel_sigma(1));
	 
	  new_error_total += .5 * transpose(static_cast<Vector2>(new_epsilon(i,j))) * 
	    inverse_cov * static_cast<Vector2>(new_epsilon(i,j));
        }
        ++i;
      }

      for (unsigned j = 0; j < U.size(); ++j) {

        Vector<double, BundleAdjustModelT::camera_params_n> new_a = m_model.A_parameters(j) + subvector(delta_a, BundleAdjustModelT::camera_params_n*j, BundleAdjustModelT::camera_params_n);
        Vector<double, BundleAdjustModelT::camera_params_n> eps_a = m_model.A_initial(j)-new_a;
        
	Matrix<double,BundleAdjustModelT::camera_params_n,BundleAdjustModelT::camera_params_n> inverse_cov;
	inverse_cov = m_model.A_inverse_covariance(j);

	new_error_total += .5 * transpose(eps_a) * inverse_cov * eps_a;
      }
      
      // We only add constraints for Ground Control Points (GCPs), not for 3D tie points.
      for (unsigned i = 0; i < V.size(); ++i) {
        if ((*m_control_net)[i].type() == ControlPoint::GroundControlPoint) {

          Vector<double, BundleAdjustModelT::point_params_n> new_b = m_model.B_parameters(i) + delta_b(i).ref();
          Vector<double, BundleAdjustModelT::point_params_n> eps_b = m_model.B_initial(i)-new_b;
	  
	  Matrix<double,BundleAdjustModelT::point_params_n,BundleAdjustModelT::point_params_n> inverse_cov;
          inverse_cov = m_model.B_inverse_covariance(i);
          
	  new_error_total += .5 * transpose(eps_b) * inverse_cov * eps_b;
        }
      }
      
      //Fletcher modification
      double Splus = new_error_total;     //Compute new objective
      double SS = error_total;            //Compute old objective
      double R = (SS - Splus)/dS;         // Compute ratio
      
      if (R>0){

	for (unsigned j=0; j<m_model.num_cameras(); ++j)
          m_model.set_A_parameters(j, m_model.A_parameters(j) + subvector(delta_a,
                                                                          BundleAdjustModelT::camera_params_n*j,
                                                                          BundleAdjustModelT::camera_params_n));
        for (unsigned i=0; i<m_model.num_points(); ++i)
          m_model.set_B_parameters(i, m_model.B_parameters(i) + delta_b(i).ref());
	
        // Summarize the stats from this step in the iteration
        double overall_norm = sqrt(new_error_total);
        double overall_delta = sqrt(error_total) - sqrt(new_error_total);

	/* Taken over mostly by Bundle Adjust Report
        std::cout << "\n" << "Sparse LM Iteration " << m_iterations << ":     "
                  << "  Overall: " << overall_norm << "  delta: " << overall_delta 
		  << "  lambda: " << m_lambda <<"   Ratio:  " << R << "\n";
	*/

	abs_tol = overall_norm;
        rel_tol = fabs(overall_delta);
	
	if(m_control == 0){
	  double temp = 1 - pow((2*R - 1),3);
	  if (temp < 1.0/3.0)
	    temp = 1.0/3.0;
	  
	  m_lambda *= temp;
	  m_nu = 2;

	} else if (m_control == 1){
	  m_lambda /= 10;
	}

	return overall_delta;

      } else { // here we didn't make progress
	
	if (m_control == 0){
	  m_lambda *= m_nu;
	  m_nu*=2;
	} else if (m_control == 1){
	  m_lambda *= 10;
	}
	double overall_delta = sqrt(error_total) - sqrt(new_error_total);

	/* // Taken over mostly by Bundle Adjust Report
	std::cout <<"\n" << "Sparse LM Iteration " << m_iterations << "  delta  "<< overall_delta << "  lambda: " << m_lambda << "   Ratio:  " << R <<"\n";
	*/

	return ScalarTypeLimits<double>::highest();
      }

      return 0;
    }
    
    BundleAdjustModelT& bundle_adjust_model() { return m_model; }
  };
  
}} // namespace vw::camera

#endif // __VW_CAMERA_BUNDLE_ADJUST_H__
