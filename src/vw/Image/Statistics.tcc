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
void histogram( const ImageViewBase<ViewT> &view, int num_bins, std::vector<double> &hist){
  
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
  if (max_val == min_val) max_val = min_val + 1.0;
    
  hist.assign(num_bins, 0.0);
  for (int row = 0; row < view.impl().rows(); row++){
    for (int col = 0; col < view.impl().cols(); col++){
      if ( !is_valid(view.impl()(col, row)) )
        continue;
      double val = view.impl()(col, row);
      int bin = (int)round( (num_bins - 1) * ( (val - min_val)/(max_val - min_val) ) );
      hist[bin]++;
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
                std::vector<double> &hist){
  
  VW_ASSERT(num_bins > 0, ArgumentErr() << "histogram: number of input bins must be positive");
  
  if (max_val == min_val) 
    max_val = min_val + 1.0;
  
  hist.assign(num_bins, 0.0);
  for (int row = 0; row < view.impl().rows(); row++){
    for (int col = 0; col < view.impl().cols(); col++){
      if ( !is_valid(view.impl()(col, row)) )
        continue;
      double val = view.impl()(col, row);
      int bin = (int)round( (num_bins - 1) * ( (val - min_val)/(max_val - min_val) ) );
      
      // Saturate the bin assignment to prevent a memory exception.
      if (bin < 0)
        bin = 0;
      if (bin > num_bins-1)
        bin = num_bins-1;
        
      hist[bin]++;
    }
  }
  return;
}

inline size_t get_histogram_percentile(std::vector<double> const& hist, double percentile) {

  // Verify the input percentile is in the legal range
  if ((percentile < 0) || (percentile > 1.0)) {
    vw_throw(ArgumentErr() << "get_histogram_percentile: illegal percentile request: " << percentile << "\n");
  }

  // Get the total pixel count
  const size_t num_bins = hist.size();
  double num_pixels = 0;
  for (size_t i=0; i<num_bins; ++i)
    num_pixels += hist[i];

  // Now go through and find the requested percentile
  double running_percentile = 0;
  for (size_t i=0; i<num_bins; ++i) {
    double this_percent = hist[i] / num_pixels;
    running_percentile += this_percent;
    //std::cout << "p - " << running_percentile << std::endl;
    if (running_percentile >= percentile)
      return i;
  }
  vw_throw(LogicErr() << "get_histogram_percentile: Illegal histogram encountered!");
}


template <class ViewT>
double optimal_threshold( const ImageViewBase<ViewT> &view){

  std::vector<double> hist;
  int num_bins = 256;
  histogram(view, num_bins, hist);

  double sum = 0.0;
  for (size_t i = 0; i < hist.size(); i++) sum += hist[i];
  if (sum == 0.0) return 0.0;

  // Normalize the histogram
  for (size_t i = 0; i < hist.size(); i++) hist[i] /= sum;
  sum = 1.0;
  
  double totalAccum = 0.0;
  for (size_t i = 0; i < hist.size(); i++) totalAccum += i*hist[i];

  // Find the variance between classes
  std::vector<double> V;
  V.assign(num_bins, 0.0);

  double leftProb = 0.0, leftAccum = 0.0, rightProb = 0.0, rightAccum = 0.0;
  for (size_t i = 0; i < hist.size(); i++){

    leftProb  += hist[i];
    leftAccum += i*hist[i];

    rightProb  = sum - leftProb;
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
			ImageView<PixelGray<vw::uint8> > &output_image,
			double low_percentile, double high_percentile, int num_bins) {

  // First get the min and max values
  double min_val, max_val;
  find_image_min_max(input_image, min_val, max_val);

  // Compute the input image histogram
  std::vector<double> hist;
  histogram(input_image, num_bins, min_val, max_val, hist);

  // Find the bins at the input percentiles
  size_t low_bin  = get_histogram_percentile(hist, low_percentile );
  size_t high_bin = get_histogram_percentile(hist, high_percentile);

  // Find the input values that correspond to the bin indices
  double bin_width = (max_val - min_val) / (static_cast<double>(num_bins - 1));
  double low_value  = (low_bin +1) * bin_width + min_val;
  double high_value = (high_bin+1) * bin_width + min_val;

  // Scale the image using the computed values and convert to uint8
  output_image = pixel_cast<vw::uint8>(normalize( clamp(input_image, low_value, high_value),
					    low_value, high_value, 0.0, 255.0 ));
}



