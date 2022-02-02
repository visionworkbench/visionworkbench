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


/// \file Statistics.tcc

template <class ViewT>
void find_image_min_max(const ImageViewBase<ViewT> &view, double &min_val, double &max_val){

  max_val = -std::numeric_limits<double>::max();
  min_val = -max_val;
  for (int row = 0; row < view.impl().rows(); row++){
    for (int col = 0; col < view.impl().cols(); col++){
      if ( !is_valid(view.impl()(col, row)) ) 
        continue;
      double val = view.impl()(col, row);
      if (val < min_val) min_val = val;
      if (val > max_val) max_val = val;
    }
  }
}

// The image bounds min_val and max_val must be found before invoking this function
template <class ViewT>
void histogram(const ImageViewBase<ViewT> &view, int num_bins, double min_val, double max_val,
               math::Histogram &hist){
  
  VW_ASSERT(num_bins > 0, ArgumentErr() << "histogram: number of input bins must be positive");
  
  if (max_val == min_val) 
    max_val = min_val + 1.0;
  
  hist.initialize(num_bins, min_val, max_val);
  for (int row = 0; row < view.impl().rows(); row++){
    for (int col = 0; col < view.impl().cols(); col++){
      if ( !is_valid(view.impl()(col, row)) )
        continue;
      double val = view.impl()(col, row);
      hist.add_value(val);
    }
  }
  return;
}

// See Statistics.h for documentation
template <class ViewT>
double otsu_threshold(const ImageViewBase<ViewT> &view){

  math::Histogram hist;
  int num_bins = 256;
  double max_val = -std::numeric_limits<double>::max();
  double min_val = -max_val;

  // First get the min and max values
  find_image_min_max(view, min_val, max_val);

  if (max_val == min_val) 
    max_val++;

  // Compute the input image histogram
  histogram(view, num_bins, min_val, max_val, hist);
  
  double sum = static_cast<double>(hist.get_total_num_values());
  if (sum == 0.0)
    return 0.0;

  // Get a weighted total of normalized histogram values 
  double totalAccum = 0.0;
  for (int i = 0; i < num_bins; i++)
    totalAccum += i*(hist.get_bin_value(i)/sum);

  // Find the variance between classes
  std::vector<double> V;
  V.assign(num_bins, 0.0);

  double leftProb = 0.0, leftAccum = 0.0, rightProb = 0.0, rightAccum = 0.0;
  for (int i = 0; i < num_bins; i++){

    leftProb  += hist.get_bin_value(i)/sum;
    leftAccum += i*(hist.get_bin_value(i)/sum);

    rightProb  = 1.0 - leftProb;
    rightAccum = totalAccum - leftAccum;

    if (leftProb == 0 || rightProb == 0.0) continue;
    
    double leftMean = leftAccum/leftProb;
    double rightMean = rightAccum/rightProb;
    V[i] = leftProb*rightProb*(leftMean-rightMean)*(leftMean-rightMean);
  }

  double maxV = *std::max_element(V.begin(), V.end());

  // If the maximum is reached at more than one index, find the average index
  double indexSum = 0, numIndices = 0;
  for (size_t i = 0; i < V.size(); i++){
    if (V[i] == maxV){
      indexSum += i;
      numIndices++;
    }
  }
  double meanIndex = indexSum/numIndices;

  // Normalize the value
  double normalized_index =  meanIndex/(num_bins - 1.0);

  // Scale 
  double thresh = min_val + normalized_index * (max_val - min_val);

  return thresh;
}

// See Statistics.h for documentation
template <class ViewT>
double otsu_threshold(const ImageViewBase<ViewT> &view,
                      int num_sample_rows, int num_sample_cols,
                      int num_bins){

  double row_ratio = double(view.impl().rows() - 1.0)/double(num_sample_rows - 1.0);
  double col_ratio = double(view.impl().cols() - 1.0)/double(num_sample_cols - 1.0);

  // Find min and max vals
  double max_val = -std::numeric_limits<double>::max();
  double min_val = -max_val;
  for (int sample_row = 0; sample_row < num_sample_rows; sample_row++){
    int row = round(sample_row * row_ratio);
    
    for (int sample_col = 0; sample_col < num_sample_cols; sample_col++){
      int col = round(sample_col * col_ratio);
      
      if (!is_valid(view.impl()(col, row))) 
        continue;
      double val = view.impl()(col, row);
      if (val < min_val) min_val = val;
      if (val > max_val) max_val = val;
    }
  }
  
  // Build the histogram
  std::vector<double> hBuf(num_bins, 0.0);
  for (int sample_row = 0; sample_row < num_sample_rows; sample_row++){
    int row = round(sample_row * row_ratio);
    
    for (int sample_col = 0; sample_col < num_sample_cols; sample_col++){
      int col = round(sample_col * col_ratio);

      double scaled = (view.impl()(col, row) - min_val) / (max_val - min_val);
      int bin = static_cast<int>(round((num_bins - 1) * scaled));
      if (bin < 0) 
        bin = 0;
      if (bin > num_bins - 1) 
        bin = num_bins - 1;
      
      hBuf[bin]++;
    }
  }
  
  // Note how we cast to double to avoid integer overflow when doing the product
  double product = double(view.impl().rows()) * double(view.impl().cols());
  double mu = 0, scale = 1.0/product;
  for(int i = 0; i < num_bins; i++ ) {
    mu += i*(double)hBuf[i];
  }
  
  mu *= scale;
  double mu1 = 0, q1 = 0;
  double max_sigma = 0, max_pos = 0;
  
  for (int i = 0; i < num_bins; i++ ) {
    
    double p_i, q2, mu2, sigma;
    
    p_i = hBuf[i]*scale;
    mu1 *= q1;
    q1 += p_i;
    q2 = 1.0 - q1;
    
    if (std::min(q1, q2) < FLT_EPSILON || std::max(q1,q2) > 1.0 - FLT_EPSILON )
      continue;
    
    mu1 = (mu1 + i*p_i)/q1;
    mu2 = (mu - q1*mu1)/q2;
    sigma = q1*q2*(mu1 - mu2)*(mu1 - mu2);
    if (sigma > max_sigma) {
      max_sigma = sigma;
      max_pos = i;
    }
  }
  
  // Undo the normalization
  return (max_val - min_val) * (max_pos / (num_bins - 1.0)) + min_val;
}

template <class ViewT>
void percentile_scale_convert(ImageViewBase<ViewT> const& input_image,
                              ImageView<vw::uint8> &output_image,
                              double low_percentile, double high_percentile, int num_bins) {

  // First get the min and max values
  double min_val = -1.0, max_val = -1.0;
  find_image_min_max(input_image, min_val, max_val);

  // Compute the input image histogram
  math::Histogram hist;
  histogram(input_image, num_bins, min_val, max_val, hist);

  // Find the bins at the input percentiles
  size_t low_bin  = hist.get_percentile(low_percentile );
  size_t high_bin = hist.get_percentile(high_percentile);

  // Find the input values that correspond to the bin indices
  double bin_width  = hist.get_bin_width();
  double low_value  = (low_bin +1) * bin_width + min_val;
  double high_value = (high_bin+1) * bin_width + min_val;

  // Scale the image using the computed values and convert to uint8
  output_image = pixel_cast<vw::uint8>(normalize( clamp(input_image, low_value, high_value),
                                                        low_value, high_value, 0.0, 255.0 ));
}



