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

#include <vw/Core/FundamentalTypes.h>
#include <vw/Math/Geometry.h>
#include <vw/Math/LinearAlgebra.h>
#include <vw/Math/Matrix.h>
#include <vw/Math/Vector.h>
#include <vw/Math/RANSAC.h>

#include <math.h>

namespace vw{
namespace math{

double normalize_longitude(double lon, bool center_on_zero) {
  const double MULTIPLE = 360.0;
  double min, max;
  if (center_on_zero) {
    min=-180.0;
    max= 180.0;
  } else {
    min=  0.0;
    max=360.0;
  }
 
  if (lon < min){
    double factor = (min - lon) / MULTIPLE;
    lon += MULTIPLE*ceil(factor);
  }
  if (lon > max){
    double factor = (lon - max) / MULTIPLE;
    lon -= MULTIPLE*ceil(factor);
  }
   
  return lon;
}

  double degree_diff(double d1, double d2) {
  double diff  = fabs(d1 - d2);
  double diff2 = fabs(d1 - (d2+360));
  double diff3 = fabs(d1 - (d2-360));
  if (diff2 < diff)
    diff = diff2;
  if (diff3 < diff)
    diff = diff3;
  return diff;
}


}} // End namespace vw::math

using namespace vw;
using namespace vw::math;


/// Solve for Normalization Similarity Matrix used for noise rej.
Matrix3x3
HomographyFittingFunctor::NormSimilarity( std::vector<Vector3> const& pts ) const {
  size_t num_points = pts.size();
  size_t dimension = 3;

  Vector2 translation;
  for ( size_t i = 0; i < num_points; i++ )
    translation+=subvector(pts[i],0,dimension-1);
  translation /= num_points;

  double scale = 0;
  for ( size_t i = 0; i < num_points; i++ )
    scale += norm_2( subvector(pts[i],0,dimension-1) - translation );
  scale = num_points*sqrt(2.)/scale;

  Matrix3x3 t;
  t(2,2) = 1;
  t(0,0) = scale;
  t(1,1) = scale;
  t(0,2) = -scale*translation[0];
  t(1,2) = -scale*translation[1];
  return t;
}

vw::Matrix3x3
HomographyFittingFunctor::BasicDLT( std::vector<Vector3 > const& input,
                                    std::vector<Vector3 > const& output )  const {
  VW_ASSERT( input.size() == 4 && output.size() == 4,
             vw::ArgumentErr() << "DLT in this implementation expects to have only 4 inputs." );
  VW_ASSERT( input[0][input[0].size()-1] == 1,
             vw::ArgumentErr() << "Input data doesn't seem to be normalized.");
  VW_ASSERT( output[0][output[0].size()-1] == 1,
             vw::ArgumentErr() << "Secondary input data doesn't seem to be normalized.");
  VW_ASSERT( input[0].size() == 3,
             vw::ArgumentErr() << "BasicDLT only supports homogeneous 2D vectors.");

  vw::Matrix<double,8,9> A;
  for ( uint8 i = 0; i < 4; i++ )
    for ( uint8 j = 0; j < 3; j++ ) {
      // Filling in -wi'*xi^T
      A(i,j+3) = -output[i][2]*input[i][j];
      // Filling in yi'*xi^T
      A(i,j+6) = output[i][1]*input[i][j];
      // Filling in wi'*xi^T
      A(i+4,j) = output[i][2]*input[i][j];
      // Filling in -xi'*xi^T
      A(i+4,j+6) = -output[i][0]*input[i][j];
    }

  Matrix<double> nullsp = nullspace(A);
  nullsp /= nullsp(8,0);
  Matrix3x3 H;
  for ( uint8 i = 0; i < 3; i++ )
    for ( uint8 j = 0; j < 3; j++ )
      H(i,j) = nullsp(i*3+j,0);
  return H;
}

// Given a set of points in 3D, heuristically estimate what it means
// for two points to be "not far" from each other. The logic is to
// find a bounding box of an inner cluster and multiply that by 0.2.
// This will produce an inlier threshold. The input matrices must
// have 3 rows and N columns where N is equal to the number of points.
double estimate3DTransInlierThresh(vw::Matrix<double> const& points) {
  
  VW_ASSERT(points.rows() == 3 && points.cols() >= 3, 
            vw::ArgumentErr() << "Too few points in estimate3DTransInlierThresh().\n");
  
  vw::Vector3 range(0.0, 0.0, 0.0);
  int num_pts = points.cols();
  
  std::vector<double> vals(num_pts);
  for (int it = 0; it < range.size(); it++) {  // iterate in each coordinate
  
    // Sort all values in given coordinate
    for (int p = 0; p < num_pts; p++)
      vals[p] = points(it, p);
    std::sort(vals.begin(), vals.end());
    
    // Find some percentiles
    int min_p = round(num_pts*0.25);
    int max_p = round(num_pts*0.75);
    if (min_p >= num_pts) min_p = num_pts - 1;
    if (max_p >= num_pts) max_p = num_pts - 1;
    double min_val = vals[min_p], max_val = vals[max_p];
    range[it] = 0.2*(max_val - min_val);
  }

  // Find the average of all ranges
  double range_val = 0.0;
  for (int it = 0; it < range.size(); it++)
    range_val += range[it];
  range_val /= range.size();

  return range_val;
}

// This fitting functor attempts to find a rotation + translation +
// scale transformation between two vectors of points or
// just a rotation + translation or just a translation.
struct Similarity3DFittingFunctor {

  typedef vw::Matrix<double, 4, 4> result_type;
  std::string m_transform_type;

  Similarity3DFittingFunctor(std::string const & transform_type):
    m_transform_type(transform_type){}
  
  /// A transformation requires 3 inputs and 3 outputs to make a fit.
  size_t min_elements_needed_for_fit(vw::Vector3 const& /* point */) const { return 3; }

  result_type operator()(std::vector<vw::Vector3> const& in_vec,
                         std::vector<vw::Vector3> const& out_vec,
                         vw::Matrix<double> const& /* initial_guess */
                         = vw::Matrix<double>() ) const {
    // check consistency
    if (in_vec.size() != out_vec.size())
      vw_throw( vw::ArgumentErr() << "There must be as many inputs as outputs to be "
                << "able to compute a transform between them.\n");
    if (in_vec.size() < min_elements_needed_for_fit(vw::Vector3()))
      vw_throw( vw::ArgumentErr() << "Cannot compute a transformation. Insufficient data.\n");

    // Convert to expected format
    int num_pts = in_vec.size();
    vw::Matrix<double> in, out;
    in.set_size(3, num_pts);
    out.set_size(3, num_pts);
    for (int col = 0; col < in.cols(); col++) {
      for (int row = 0; row < in.rows(); row++) {
        in(row, col)  = in_vec[col][row];
        out(row, col) = out_vec[col][row];
      }
    }

    // Find the 3D transform for these inputs
    vw::Matrix<double, 3, 3> rotation;
    vw::Vector<double, 3>    translation;
    double                   scale;
    vw::math::find_3D_transform_aux(in, out, rotation, translation, scale, m_transform_type);

    // Convert to expected format
    result_type out_trans;
    out_trans.set_identity();
    submatrix(out_trans, 0, 0, 3, 3) = scale * rotation;
    for (int it = 0; it < 3; it++) 
      out_trans(it, 3) = translation[it];
    
    return out_trans;
  }
  
};

/// This metric can be used to measure the error between a 3D 
/// point p and a 3D point p1 that is transformed by a
/// 4x4 matrix H corresponding to a 3D affine transform (or its
/// particular cases).
typedef HomogeneousL2NormErrorMetric<3> Affine3DErrorMetric;

void vw::math::find_3D_transform(vw::Matrix<double> const & in, 
                                 vw::Matrix<double> const & out,
                                 vw::Matrix<double,3,3>   & rotation,
                                 vw::Vector<double,3>     & translation,
                                 double                   & scale,
                                 std::string        const & transform_type,
                                 bool                       filter_outliers,
                                 vw::Vector2        const & ransac_params) {

  // Initialize the outputs
  rotation.set_identity();
  translation = vw::Vector3();
  scale = 1.0;
  
  // Sanity checks
  VW_ASSERT((in.rows() == 3) && (in.rows() == out.rows()) &&
            (in.cols() == out.cols()), 
            vw::ArgumentErr() << "find_3D_transform(): input data "
            << "has incorrect size.\n");
  VW_ASSERT(in.cols() >= 3, 
            vw::ArgumentErr() << "find_3D_transform(): Must have at least "
            << "three data points.\n");

  if (!filter_outliers) {
    find_3D_transform_aux(in, out, rotation, translation, scale, transform_type);
    return;
  }

  // Copy the data to vectors.
  // TODO(oalexan1): Convert all uses of this function to an interface
  // using vectors and not matrices for the points.
  int num_pts = in.cols();
  std::vector<vw::Vector3> in_pts(num_pts), out_pts(num_pts);
  for (int col = 0; col < in.cols(); col++) {
    for (int row = 0; row < in.rows(); row++) {
      in_pts[col][row]  = in(row, col);
      out_pts[col][row] = out(row, col);
    }
  }

  int num_ransac_iterations = ransac_params[0];
  double outlier_factor     = ransac_params[1];
  if (num_ransac_iterations < 1 || outlier_factor <= 0.0)
    vw_throw(ArgumentErr() << "Invalid parameters were provided for outlier filtering.\n");

  // Find the inlier threshold based on the distribution of points in the output
  double inlier_threshold = outlier_factor * estimate3DTransInlierThresh(out);
  int min_num_output_inliers = std::max(num_pts/2, 3);
  bool reduce_min_num_output_inliers_if_no_fit = true;

  vw_out() << "Starting RANSAC.\n";
  vw::Matrix<double, 4, 4> transform;
  std::vector<size_t> inlier_indices;
  try {
    // Must first create the functor and metric, then pass these to ransac. If
    // created as inline arguments to ransac, these may go go out
    // of scope prematurely, which will result in incorrect behavior.
    Similarity3DFittingFunctor fitting_functor(transform_type);
    Affine3DErrorMetric error_metric;
    RandomSampleConsensus<Similarity3DFittingFunctor, Affine3DErrorMetric>
      ransac(fitting_functor, error_metric, num_ransac_iterations, inlier_threshold,
                min_num_output_inliers, reduce_min_num_output_inliers_if_no_fit);
    transform = ransac(in_pts, out_pts);
    inlier_indices = ransac.inlier_indices(transform, in_pts, out_pts);
  } catch (const vw::math::RANSACErr& e ) {
    vw_out() << "RANSAC failed: " << e.what() << "\n";
    return;
  }
  vw_out() << "Found " << inlier_indices.size() << " / " << num_pts << " inliers.\n";

  // Create the outputs
  rotation = vw::math::submatrix(transform, 0, 0, 3, 3);
  scale = pow(vw::math::det(rotation), 1.0/3.0);
  rotation /= scale;
  for (int it = 0; it < 3; it++) 
    translation[it] = transform(it, 3);
  
  return;
}
  
void vw::math::find_3D_transform_aux(vw::Matrix<double> const & in_vec, 
                                     vw::Matrix<double> const & out_vec,
                                     vw::Matrix<double,3,3>   & rotation,
                                     vw::Vector<double,3>     & translation,
                                     double                   & scale,
                                     std::string        const & transform_type) {
  
  if (transform_type != "similarity" && transform_type != "rigid" &&
      transform_type != "translation") {
    vw_throw(vw::ArgumentErr() << "find_3D_transform_aux: Expecting to compute a "
             << "transform which is either similarity, or rigid, or translation." );
  }
    
  // Make copies that we can modify inline
  vw::Matrix<double> in  = in_vec; 
  vw::Matrix<double> out = out_vec;
    
  // Default output
  rotation.set_identity();
  translation.set_all(0.0);
  scale = 1.0;
    
  VW_ASSERT((in.rows() == 3) && (in.rows() == out.rows()) && (in.cols() == out.cols()), 
            vw::ArgumentErr() << "find_3D_transform(): input data "
            << "has incorrect size.\n");
  VW_ASSERT(in.cols() >= 3, 
            vw::ArgumentErr() << "find_3D_transform(): Must have at least "
            << "three data points.\n");
    
  typedef vw::math::MatrixCol<vw::Matrix<double>> ColView;

  // First find the scale, by finding the ratio of sums of some distances,
  // then bring the datasets to the same scale.
  if (transform_type == "similarity") {
    int num_segments = 0;
    double dist_in = 0, dist_out = 0;
    for (int col = 0; col < in.cols() - 1; col++) {
      
      int next_col = col + 1;
      
      num_segments++;
      
      ColView inCol1 (in,  col), inCol2 (in,  next_col);
      ColView outCol1(out, col), outCol2(out, next_col);
      dist_in  += vw::math::norm_2(inCol2  - inCol1 );
      dist_out += vw::math::norm_2(outCol2 - outCol1);
    }

    if (num_segments < 1 || dist_in <= 0 || dist_out <= 0) 
      vw_throw(vw::ArgumentErr() << "find_3D_transform(): not enough distinct points "
               << "to find the scale.\n");
    
    scale = dist_out/dist_in;
    out /= scale;
  }

  // Find the centroids then shift to the origin
  vw::Vector3 in_ctr;
  vw::Vector3 out_ctr;
  int col_count = 0;
  for (size_t col = 0; col < in.cols(); col++) {

    ColView inCol (in,  col);
    ColView outCol(out, col);
    in_ctr  += inCol;
    out_ctr += outCol;
    col_count++;
  }

  // Get the mean
  in_ctr  /= col_count;
  out_ctr /= col_count;

  // Subtract mean from in and out
  for (size_t col = 0; col < in.cols(); col++) { 
    ColView inCol (in,  col);
    ColView outCol(out, col);
    inCol  -= in_ctr;
    outCol -= out_ctr;
  }

  if (transform_type != "translation") {
    // SVD
    vw::Matrix<double> cov = in * vw::math::transpose(out);
    vw::Matrix<float> U, VT;
    vw::Vector<float> s;
    vw::math::svd(cov, U, s, VT);
    
    // Find the rotation
    double d = vw::math::det(vw::math::transpose(VT) * vw::math::transpose(U));
    if (d > 0)
      d = 1.0;
    else
      d = -1.0;
    vw::Matrix3x3 I;
    I.set_identity();
    I(2, 2) = d;
    rotation = vw::math::transpose(VT) * I * vw::math::transpose(U);
  }
  
  translation = scale*(out_ctr - rotation*in_ctr);

  return;
}
