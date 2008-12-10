// __BEGIN_LICENSE__
// 
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
// 
// Copyright 2006 Carnegie Mellon University. All rights reserved.
// 
// This software is distributed under the NASA Open Source Agreement
// (NOSA), version 1.3.  The NOSA has been approved by the Open Source
// Initiative.  See the file COPYING at the top of the distribution
// directory tree for the complete NOSA document.
// 
// THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
// KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
// LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
// SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
// A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
// THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
// DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
// 
// __END_LICENSE__

/// \file LearnPCA.h
///
/// Learn an interest point descriptor set via PCA.
///

#ifndef __VW_INTERESTPOINT_LEARNPCA_H__
#define __VW_INTERESTPOINT_LEARNPCA_H__

#include <vw/InterestPoint/InterestData.h>
#include <vw/InterestPoint/Detector.h>
#include <vw/InterestPoint/Descriptor.h>
#include <vw/InterestPoint/MatrixIO.h>
#include <vw/InterestPoint/VectorIO.h>
#include <vw/Math/Matrix.h>
#include <vw/Math/LinearAlgebra.h>
#include <vw/FileIO/DiskImageView.h>

#define PCA_BASIS_SIZE	20
#define MAX_POINTS_TO_DRAW 1000

using namespace std;
using namespace vw;
using namespace vw::ip;

// Use DescriptorGeneratorBase::operator() and compute_descriptor methods to 
// help us populate the training data matrix
// DescriptorGeneratorBase takes care of finding the support region for
// each interest point
class LearnPCADataFiller : public DescriptorGeneratorBase<LearnPCADataFiller> {

private:
  Matrix<double> *data;
  int column;

public:
  // *** TODO: Pass by ref without copy?
  LearnPCADataFiller(Matrix<double>* training_data, unsigned int start_column) {
    data = training_data;
    column = start_column;
  }
  
  template <class ViewT>
  Vector<float> compute_descriptor (ImageViewBase<ViewT> const& support) {
    // Create dummy return vector
    Vector<float> result(0.0f);
    
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

    return result;
  }

  template <class ViewT>
  inline ImageView<typename ViewT::pixel_type> get_support( float x, float y, float scale, float ori,
                                                              ImageViewBase<ViewT> const& source, int size=DEFAULT_SUPPORT_SIZE ) {
    float half_size = ((float)(size - 1)) / 2.0f;
    float scaling = 1.0f / scale;

    // This is mystifying - why won't the four-arg compose work?
    return transform(source.impl(),
		     compose(TranslateTransform(half_size, half_size),
			     compose(ResampleTransform(scaling, scaling),
				     RotateTransform(-ori),
				     TranslateTransform(-x, -y))),
		     size, size);
  }
};

class LearnPCA {

private:
  std::string basis_filename;
  std::string avg_filename;

  Matrix<double> training_data;
  unsigned int training_data_size;

  Matrix<float> pca_basis;
  Vector<float> pca_avg;

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

  void processImage(DiskImageView<PixelRGB<uint8> > &dimage) {
    
    float log_threshold = 0.01;
    int tile_size = 2048;
    
    // *** TODO: let user choose interest point detector
    // Find interest points in image
    InterestPointList ipl;

		// *** TODO: set this somewhere else
    int max_x_dim = 1000;
    ImageView<PixelRGB<uint8> > image = dimage;
    while(image.cols() > max_x_dim) {
      image = bin_subsample(image);
      cout << "Reduced image to " << image.cols() << "x" << image.rows() << endl;
    }
    image = gaussian_filter(image, 1);

    LogInterestOperator interest_operator(log_threshold);
    ScaledInterestPointDetector<LogInterestOperator> detector(interest_operator);
    //ScaledInterestPointDetector<LogInterestOperator> detector;
    cout << "Running interest point detector on " << dimage.filename() << endl;
    ipl = detector(image, tile_size);
    write_point_image("ip_" + dimage.filename(), image, ipl);

    // Resize the training data matrix to accommodate new interest points
    unsigned int interest_point_index = training_data_size;
    training_data_size += ipl.size();
    training_data.set_size(support_squared, training_data_size, true);
    
    // Populate training data matrix
    LearnPCADataFiller fill_matrix(&training_data, interest_point_index);
    
    cout << "Populating training matrix" << endl;
    cout << "  Starting fill at column " << interest_point_index << endl;
    cout << "  " << ipl.size() << " interest points" << endl;
    cout << "  " << training_data_size << " total interest points\n" << endl;
    fill_matrix(image, ipl);
  }

  void runPCA() {
    cout << "Running PCA on training data" << endl;

    // Compute average
    cout << "  Computing average" << endl;
    pca_avg.set_size(support_squared);
    for (int i = 0; i < support_squared; i++) {
      for (unsigned j = 0; j < training_data.cols(); j++) {
	pca_avg(i) += training_data(i, j);
      }
      pca_avg(i) /= training_data.cols();
    }

    // Subtract average
    cout << "  Subtracting average" << endl;
    for (int i = 0; i < support_squared; i++) {
      for (unsigned j = 0; j < training_data.cols(); j++) {
	training_data(i, j) -= pca_avg(i);
      }
    }

    // Compute SVD of training data matrix
    cout << "  Computing SVD" << endl;
    Matrix<float> U;
    Vector<float> E;
    Matrix<float> VT;
    svd(training_data, U, E, VT);

    assert(E.size() >= PCA_BASIS_SIZE);
    cout << "  Top " << PCA_BASIS_SIZE << " singular values from " << E.size()
	 << " total singular values" << endl;
    for (int i = 0; i < PCA_BASIS_SIZE; i++) {
      cout << "  " << i << ": " << E[i] << endl;
    }

    // Take top n eignvectors of U as PCA basis
    Matrix<float> pca_basis = submatrix(U, 0, 0, U.rows(), 
					PCA_BASIS_SIZE);

    cout << "pca_basis: " << pca_basis.rows() << " x " << pca_basis.cols() << endl;
    cout << "pca_avg: " << pca_avg.size() << endl;

    write_matrix(basis_filename, pca_basis);
    write_vector(avg_filename, pca_avg);
  }

  // Draw the interest points and write as an image.
  template <class ViewT>
  void write_point_image(std::string out_file_name, 
			 ImageViewBase<ViewT> const& src,
			 InterestPointList const& points) {

    ImageView<PixelRGB<uint8> > viz = pixel_cast<PixelRGB<uint8> >(channel_cast_rescale<uint8>(src));

    // Draw points into color planes
    int n = 0;
    for (InterestPointList::const_iterator pt = points.begin();
	 pt != points.end() && n < MAX_POINTS_TO_DRAW; ++pt, ++n) {
      // Draw a red line from the point outward along the orientation
      for (int r=0; r<(int)(8*(*pt).scale); ++r){
	int i = (int)(0.5 + (*pt).x + r*cos((*pt).orientation));
	int j = (int)(0.5 + (*pt).y + r*sin((*pt).orientation));
	// red, green, blue
	viz(i,j) = PixelRGB<uint8>(255, 0, 0);
      }
      // Draw a green 3x3 filled square at the point to indicate center
      int i0 = (int)(0.5 + (*pt).x);
      int j0 = (int)(0.5 + (*pt).y);
      for (int j=j0-1; j<=j0+1; ++j){
	for (int i=i0-1; i<=i0+1; ++i){
	  // red, green, blue
	  viz(i,j) = PixelRGB<uint8>(0, 255, 0);
	}
      }
    }

    write_image(out_file_name, viz);
  }

};

#endif // __VW_INTERESTPOINT_LEARNPCA_H__
