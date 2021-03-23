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

void vw::math::find_3D_affine_transform(vw::Matrix<double> const & in_vec, 
                                        vw::Matrix<double> const & out_vec,
                                        vw::Matrix<double,3,3>   & rotation,
                                        vw::Vector<double,3>     & translation,
                                        double                   & scale,
                                        std::string        const & transform_type,
                                        bool                       filter_outliers,
                                        vw::Vector2        const & outlier_removal_params) {
  
  std::vector<bool> is_outlier;
  if (!filter_outliers) {
    find_3D_affine_transform_aux(in_vec, out_vec, rotation, translation, scale,
                                 transform_type, filter_outliers, outlier_removal_params,
                                 is_outlier);
    return;
  }

  // Filter outliers using 3*(75-th percentile).

  // TODO: Need to do an honest RANSAC, but then the question
  // becomes how to estimate a good outlier factor for it which will
  // likely require RANSAC to be done twice with the factor being
  // figured at the first pass, like the error for half of the
  // measurements multiplied by 2 or so. The current approach could
  // be good enough if outliers are not many or not large.
    
  int num_attempts      = 5;
  double percentile     = outlier_removal_params[0];
  double outlier_factor = outlier_removal_params[1];
  int num_pts           = in_vec.cols();

  // Start with no outliers
  is_outlier.resize(num_pts);
  for (int ipt = 0; ipt < num_pts; ipt++) 
    is_outlier[ipt] = false;

  std::vector<double> errors(num_pts);
  for (int attempt = 0; attempt < num_attempts; attempt++) {
      
    find_3D_affine_transform_aux(in_vec, out_vec,  rotation, translation,
                                 scale, transform_type, filter_outliers, outlier_removal_params,
                                 is_outlier);

    for (int col = 0; col < num_pts; col++) {
      Vector3 src, ref;
      for (int row = 0; row < 3; row++) {
        src[row] = in_vec(row, col);
        ref[row] = out_vec(row, col);
      }
      Vector3 trans_src = scale*rotation*src + translation;
      errors[col] = norm_2(ref - trans_src);
    }

    std::sort(errors.begin(), errors.end());

    int cutoff = round(num_pts * (percentile/100.0)) - 1;
    if (cutoff < 0) cutoff = 0;

    double thresh = outlier_factor * errors[cutoff];
    
    // Flag the outliers. Must recompute the errors since they were sorted
    for (int col = 0; col < num_pts; col++) {
      Vector3 src, ref;
      for (int row = 0; row < 3; row++) {
        src[row] = in_vec(row, col);
        ref[row] = out_vec(row, col);
      }
      Vector3 trans_src = scale*rotation*src + translation;
      errors[col] = norm_2(ref - trans_src);
      is_outlier[col] = (errors[col] > thresh);
    }
  }
    
}
  
void vw::math::find_3D_affine_transform_aux(vw::Matrix<double> const & in_vec, 
                                            vw::Matrix<double> const & out_vec,
                                            vw::Matrix<double,3,3>   & rotation,
                                            vw::Vector<double,3>     & translation,
                                            double                   & scale,
                                            std::string        const & transform_type,
                                            bool                       filter_outliers,
                                            vw::Vector2        const & outlier_removal_params,
                                            std::vector<bool>        & is_outlier) {

  
  if (transform_type != "similarity" && transform_type != "rigid" &&
      transform_type != "translation") {
    vw_throw( vw::ArgumentErr() << "find_3D_affine_transform_aux: Expecting to compute a "
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
            vw::ArgumentErr() << "find_3D_affine_transform(): input data has incorrect size.\n");
  VW_ASSERT((in.cols() >= 3), 
            vw::ArgumentErr() << "find_3D_affine_transform(): Must have at least "
            << "three data points.\n");
    
  typedef vw::math::MatrixCol<vw::Matrix<double> > ColView;

  // First find the scale, by finding the ratio of sums of some distances,
  // then bring the datasets to the same scale.
  if (transform_type == "similarity") {
    double dist_in = 0, dist_out = 0;
    for (size_t col = 0; col < in.cols() - 1; col++) {
      
      if (filter_outliers && is_outlier[col]) continue;
      
      size_t next_col = col + 1;
      if (filter_outliers) {
        // Find the next column that is not an outlier
        bool success = false;
        while (next_col < in.cols()){
          if (!is_outlier[next_col]) {
            success = true;
            break;
          }
          next_col++;
        }
        if (!success) continue;
      }
      
      ColView inCol1 (in,  col), inCol2 (in,  next_col);
      ColView outCol1(out, col), outCol2(out, next_col);
      dist_in  += vw::math::norm_2(inCol2  - inCol1 );
      dist_out += vw::math::norm_2(outCol2 - outCol1);
    }
      
    if (dist_in <= 0 || dist_out <= 0)
      return;

    scale = dist_out/dist_in;
    out /= scale;
  }

  // Find the centroids then shift to the origin
  vw::Vector3 in_ctr;
  vw::Vector3 out_ctr;
  int good_col = 0;
  for (size_t col = 0; col < in.cols(); col++) {

    if (filter_outliers && is_outlier[col]) continue;

    ColView inCol (in,  col);
    ColView outCol(out, col);
    in_ctr  += inCol;
    out_ctr += outCol;
    good_col++;
  }

  // Get the mean
  in_ctr  /= good_col;
  out_ctr /= good_col;

  // Subtract mean from in and out
  for (size_t col = 0; col < in.cols(); col++) { 

    if (filter_outliers && is_outlier[col]) continue;

    ColView inCol (in,  col);
    ColView outCol(out, col);
    inCol  -= in_ctr;
    outCol -= out_ctr;
  }

  // Wipe columns with outliers
  good_col = 0;
  if (filter_outliers) {
    for (size_t col = 0; col < in.cols(); col++) {

      if (filter_outliers && is_outlier[col]) continue;
        
      for (size_t row = 0; row < in.rows(); row++) {
        in(row, good_col)  = in(row, col);
        out(row, good_col) = out(row, col);
      }
      good_col++;
    }
   
    VW_ASSERT((good_col >= 1), 
              vw::ArgumentErr() << "find_3D_affine_transform(): no data "
              << "left after outlier filtering.\n");
      
    in.set_size(in.rows(), good_col);
    out.set_size(out.rows(), good_col);
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
