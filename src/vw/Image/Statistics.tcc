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
void histogram( const ImageViewBase<ViewT> &view, int num_bins, math::Histogram &hist){
  
  VW_ASSERT(num_bins > 0, ArgumentErr() << "histogram: number of input bins must be positive");
  
  // Find the maximum and minimum
  double max_val = -std::numeric_limits<double>::max(), min_val = -max_val;
  for (int row = 0; row < view.impl().rows(); row++){
    for (int col = 0; col < view.impl().cols(); col++){
      if ( !is_valid(view.impl()(col, row)) )
        continue;
      double val = view.impl()(col, row);
      if (val < min_val) min_val = val;
      if (val > max_val) max_val = val;
    }
  }
  if (max_val == min_val)
    max_val = min_val + 1.0;
    
  hist.initialize(num_bins, min_val, max_val);
  for (int row = 0; row < view.impl().rows(); row++){
    for (int col = 0; col < view.impl().cols(); col++){
      if ( !is_valid(view.impl()(col, row)) )
        continue;
      double val = view.impl()(col, row);
      hist.add_value_no_check(val);
    }
  }

  return;
}


template <class ViewT>
void find_image_min_max( const ImageViewBase<ViewT> &view, double &min_val, double &max_val){

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

template <class ViewT>
void histogram( const ImageViewBase<ViewT> &view, int num_bins, double min_val, double max_val,
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


template <class ViewT>
double optimal_threshold( const ImageViewBase<ViewT> &view){

  math::Histogram hist;
  int num_bins = 256;
  histogram(view, num_bins, hist);

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
  return meanIndex/(num_bins - 1.0);
}

template <class ViewT>
void percentile_scale_convert(ImageViewBase<ViewT> const& input_image,
                              ImageView<vw::uint8> &output_image,
                              double low_percentile, double high_percentile, int num_bins) {

  // First get the min and max values
  double min_val, max_val;
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



