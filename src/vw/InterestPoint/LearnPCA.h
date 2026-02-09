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


/// \file LearnPCA.h
///
/// Learn an interest point descriptor set via PCA.
///

#ifndef __VW_INTERESTPOINT_LEARNPCA_H__
#define __VW_INTERESTPOINT_LEARNPCA_H__

#include <vw/Math/LinearAlgebra.h>
#include <vw/Image/ImageMath.h>
#include <vw/Image/ImageView.h>
#include <vw/FileIO/DiskImageView.h>

#include <vw/InterestPoint/InterestPoint.h>
#include <vw/InterestPoint/DetectorBase.h>
#include <vw/InterestPoint/Descriptor.h>
#include <vw/FileIO/MatrixIO.h>

#define PCA_BASIS_SIZE  20
#define MAX_POINTS_TO_DRAW 1000

// Use DescriptorGeneratorBase::operator() and compute_descriptor methods to
// help us populate the training data matrix
// DescriptorGeneratorBase takes care of finding the support region for
// each interest point
class LearnPCADataFiller : public vw::ip::DescriptorGeneratorBase<LearnPCADataFiller> {

private:
  vw::Matrix<double> *data;
  int column;

public:
  // *** TODO: Pass by ref without copy?
  LearnPCADataFiller(vw::Matrix<double>* training_data, unsigned int start_column) {
    data = training_data;
    column = start_column;
  }

  template <class ViewT, class IterT>
  void compute_descriptor (vw::ImageViewBase<ViewT> const& support,
                           IterT /*first*/, IterT /*last*/ ) {
    // This secretly does not create a descriptor

    int support_squared = support.impl().cols() * support.impl().rows();

    double norm_const = 0.0;
    // Copy support region into training data matrix
    int row = 0;
    for (int j = 0; j < support.impl().rows(); j++) {
      for (int i = 0; i < support.impl().cols(); i++) {
        (*data)(row++, column) = support.impl()(i, j);
        norm_const += support.impl()(i, j) * support.impl()(i, j);
      }
    }

    // Normalize values in column
    norm_const = sqrt(norm_const);
    for (int i = 0; i < support_squared; i++) {
      (*data)(i, column) /= norm_const;
    }

    // Increment column number for next call to compute_descriptor
    column++;
  }

  int descriptor_size() { return 0; }
};

class LearnPCA {

private:
  std::string basis_filename;
  std::string avg_filename;

  vw::Matrix<double> training_data;
  unsigned int training_data_size;

  vw::Matrix<float> pca_basis;
  vw::Vector<float> pca_avg;

  int support_squared;

public:

  static const int DEFAULT_SUPPORT_SIZE = 41;

  LearnPCA(const std::string& pcabasis_filename,
           const std::string& pcaavg_filename)
    : basis_filename(pcabasis_filename), avg_filename(pcaavg_filename) {

    support_squared = DEFAULT_SUPPORT_SIZE * DEFAULT_SUPPORT_SIZE;
    training_data_size = 0;
  }

  template <class T>
  vw::ImageView<T> bin_subsample(const vw::ImageView<T> &image) {
    vw::ImageView<T> ret(image.cols()/2, image.rows()/2);
    for (int y= 0; y< ret.rows(); y++) {
      for (int x= 0; x< ret.cols(); x++) {
        ret(x,y) = T(
                     (image(x*2, y*2  )/2 + image(x*2+1, y*2  )/2)/2 +
                     (image(x*2, y*2+1)/2 + image(x*2+1, y*2+1)/2)/2);
      }
    }
    return ret;
  }

  void processImage(vw::DiskImageView<vw::PixelRGB<vw::uint8> > &dimage) {

    float log_threshold = 0.01;
    int tile_size = 2048;

    // *** TODO: let user choose interest point detector
    // Find interest points in image
    vw::ip::InterestPointList ipl;

                // *** TODO: set this somewhere else
    int max_x_dim = 1000;
    vw::ImageView<vw::PixelRGB<vw::uint8> > image = dimage;
    while(image.cols() > max_x_dim) {
      image = bin_subsample(image);
      std::cout << "Reduced image to " << image.cols() << "x" << image.rows() << std::endl;
    }
    image = vw::gaussian_filter(image, 1);

    vw::ip::LogInterestOperator interest_operator(log_threshold);
    vw::ip::ScaledInterestPointDetector<vw::ip::LogInterestOperator> detector(interest_operator);
    //ScaledInterestPointDetector<LogInterestOperator> detector;
    std::cout << "Running interest point detector on " << dimage.filename() << std::endl;
    ipl = detector(image, tile_size);
    write_point_image("ip_" + dimage.filename(), image, ipl);

    // Resize the training data matrix to accommodate new interest points
    unsigned int interest_point_index = training_data_size;
    training_data_size += ipl.size();
    training_data.set_size(support_squared, training_data_size, true);

    // Populate training data matrix
    LearnPCADataFiller fill_matrix(&training_data, interest_point_index);

    std::cout << "Populating training matrix" << std::endl;
    std::cout << "  Starting fill at column " << interest_point_index << std::endl;
    std::cout << "  " << ipl.size() << " interest points" << std::endl;
    std::cout << "  " << training_data_size << " total interest points\n" << std::endl;
    fill_matrix(image, ipl);
  }

  void runPCA() {
  std::cout << "Running PCA on training data" << std::endl;

    // Compute average
    std::cout << "  Computing average" << std::endl;
    pca_avg.set_size(support_squared);
    for (int i = 0; i < support_squared; i++) {
      for (unsigned j = 0; j < training_data.cols(); j++) {
        pca_avg(i) += training_data(i, j);
      }
      pca_avg(i) /= training_data.cols();
    }

    // Subtract average
    std::cout << "  Subtracting average" << std::endl;
    for (int i = 0; i < support_squared; i++) {
      for (unsigned j = 0; j < training_data.cols(); j++) {
        training_data(i, j) -= pca_avg(i);
      }
    }

    // Compute SVD of training data matrix
    std::cout << "  Computing SVD" << std::endl;
    vw::Matrix<float> U;
    vw::Vector<float> E;
    vw::Matrix<float> VT;
    svd(training_data, U, E, VT);

    assert(E.size() >= PCA_BASIS_SIZE);
    std::cout << "  Top " << PCA_BASIS_SIZE << " singular values from " << E.size()
         << " total singular values" << std::endl;
    for (int i = 0; i < PCA_BASIS_SIZE; i++) {
      std::cout << "  " << i << ": " << E[i] << std::endl;
    }

    // Take top n eignvectors of U as PCA basis
    vw::Matrix<float> pca_basis = submatrix(U, 0, 0, U.rows(),
                                        PCA_BASIS_SIZE);

    std::cout << "pca_basis: " << pca_basis.rows() << " x " << pca_basis.cols() << std::endl;
    std::cout << "pca_avg: " << pca_avg.size() << std::endl;

    vw::write_matrix(basis_filename, pca_basis);
    vw::write_vector(avg_filename, pca_avg);
  }

  // Draw the interest points and write as an image.
  template <class ViewT>
  void write_point_image(std::string out_file_name,
                         vw::ImageViewBase<ViewT> const& src,
                         vw::ip::InterestPointList const& points) {

    vw::ImageView<vw::PixelRGB<vw::uint8> > viz = vw::pixel_cast<vw::PixelRGB<vw::uint8> >(vw::channel_cast_rescale<vw::uint8>(src));

    // Draw points into color planes
    int n = 0;
    for (vw::ip::InterestPointList::const_iterator pt = points.begin();
         pt != points.end() && n < MAX_POINTS_TO_DRAW; ++pt, ++n) {
      // Draw a red line from the point outward along the orientation
      for (int r=0; r<(int)(8*(*pt).scale); ++r){
        int i = (int)(0.5 + (*pt).x + r*cos((*pt).orientation));
        int j = (int)(0.5 + (*pt).y + r*sin((*pt).orientation));
        // red, green, blue
        viz(i,j) = vw::PixelRGB<vw::uint8>(255, 0, 0);
      }
      // Draw a green 3x3 filled square at the point to indicate center
      int i0 = (int)(0.5 + (*pt).x);
      int j0 = (int)(0.5 + (*pt).y);
      for (int j=j0-1; j<=j0+1; ++j){
        for (int i=i0-1; i<=i0+1; ++i){
          // red, green, blue
          viz(i,j) = vw::PixelRGB<vw::uint8>(0, 255, 0);
        }
      }
    }

    vw::write_image(out_file_name, viz);
  }

};

#endif // __VW_INTERESTPOINT_LEARNPCA_H__
