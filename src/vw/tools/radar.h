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


#ifndef __VW_RADAR_H__
#define __VW_RADAR_H__

#include <stdlib.h>
#include <boost/filesystem.hpp>
#include <vw/Core/Functors.h>
#include <vw/Math/Functors.h>
#include <vw/Math/Statistics.h>
#include <vw/Image/Algorithms.h>
#include <vw/Image/ImageSurface.h>
#include <vw/Image/ImageChannels.h>
#include <vw/Image/Filter.h>
#include <vw/Image/BlobIndex.h>
#include <vw/Image/Statistics.h>
#include <vw/Image/PixelMask.h>
#include <vw/Image/Manipulation.h>
#include <vw/Image/Transform.h>
#include <vw/Image/ImageIO.h>
#include <vw/Image/PixelMask.h>
#include <vw/Image/MaskViews.h>
#include <vw/Image/UtilityViews.h>
#include <vw/Image/PerPixelViews.h>
#include <vw/FileIO/DiskImageView.h>
#include <vw/Cartography/GeoReferenceUtils.h>
#include <vw/Cartography/GeoTransform.h>
#include <vw/tools/flood_common.h>

/**
  Tools for processing radar data.

  Currently only Sentinel-1 data is supported.
*/

namespace vw {
namespace radar {


//=========================================================================================

// TODO: Move these into vw/Math?

/// Standard Z shaped fuzzy math membership function between a and b
template <typename T>
class FuzzyMembershipZFunctor : public ReturnFixedType<T> {
  float m_a, m_b, m_c, m_dba; ///< Constants
public:
  /// Constructor
  FuzzyMembershipZFunctor(float a, float b) : m_a(a), m_b(b), m_c((a+b)/2.0), m_dba(b-a) {}

  /// Apply Z function
  T operator()(T const& v) const {
    if (v < m_a)
      return T(1.0);
    if (v < m_c)
      return T(1.0 - 2.0*pow(((v-m_a)/m_dba), 2.0));
    if (v < m_b)
      return T(2.0*pow(((v-m_b)/m_dba), 2.0));
    return T(0.0);
  }
}; // End class FuzzyMembershipZFunctor


/// Standard S shaped fuzzy math membership function between a and b
template <typename T>
class FuzzyMembershipSFunctor : public ReturnFixedType<T> {
  float m_a, m_b, m_c, m_dba; ///< Constants
public:
  /// Constructor
  FuzzyMembershipSFunctor(float a, float b) : m_a(a), m_b(b), m_c((a+b)/2.0), m_dba(b-a) {}

  /// Apply S function
  T operator()(T const& v) const {
    if (v < m_a)
        return T(0.0);
    if (v < m_c)
        return T(2*pow(((v-m_a)/m_dba), 2.0));
    if (v < m_b)
        return T(1.0 - 2*pow(((v-m_b)/m_dba), 2.0));
    return T(1.0);
  }
}; // End class FuzzyMembershipSFunctor



//=========================================================================================

// TODO: Consolidate these functions in vw/math/statistics

/// As part of the Kittler/Illingworth method, compute J(T) for a given T
/// - Input should be a percentage histogram
double compute_kittler_illingworth_jt(std::vector<double> const& histogram,
                   int num_bins, double min_val, double max_val, int bin) {

    double FAIL_VAL = 999999.0; // Just returning a high error should work fine

    // For convenience, generate a list of bin values.
    std::vector<double> bin_values(num_bins);
    double bin_width = (max_val - min_val) / num_bins;
    for (int i=0; i<num_bins; ++i)
      bin_values[i] = min_val + i*bin_width;

    // Compute the total value in the two classes and the mean value.
    double P1 = 0, P2 = 0;
    double weighted_sum1=0, weighted_sum2=0;
    for (int i=0; i<=bin; ++i) {
      P1 += histogram[i];
      weighted_sum1 += histogram[i]*bin_values[i];
    }
    for (int i=bin+1; i<num_bins; ++i) {
      P2 += histogram[i];
      weighted_sum2 += histogram[i]*bin_values[i];
    }

    // Only continue if both classes contain at least one pixel.
    if ((P1 <= 0) or (P2 <= 0))
        return FAIL_VAL;
    double mean1 = weighted_sum1 / P1;
    double mean2 = weighted_sum2 / P2;

    // Compute the standard deviations of the classes.
    double sigma1=0, sigma2=0;
    for (int i=0; i<=bin; ++i)
      sigma1 += pow(bin_values[i] - mean1, 2.0) * histogram[i];
    for (int i=bin+1; i<num_bins; ++i)
      sigma2 += pow(bin_values[i] - mean2, 2.0) * histogram[i];
    sigma1 /= P1;
    sigma2 /= P2;

    // Make sure both classes contain at least two intensity values.
    if ((sigma1 <= 0) or (sigma2 <= 0))
        return FAIL_VAL;

    // Compute J(T).
    double J = 1.0 + 2.0*(P1*log(sigma1) + P2*log(sigma2)) 
                   - 2.0*(P1*log(P1    ) + P2*log(P2    ));
    return J;
} // End function compute_kittler_illingworth_jt

/// Tries to compute an optimal histogram threshold using the Kittler/Illingworth method
double split_histogram_kittler_illingworth(std::vector<double> const& histogram,
                                     int num_bins, double min_val, double max_val) {
  double bin_width = (max_val - min_val) / num_bins;

  // Normalize the histogram (each bin is now a percentage)
  std::vector<double> histogram_percentages(num_bins);
  double sum = 0;
  for (int i=0; i<num_bins; ++i) // Compute sum
    sum += histogram[i];
  //std::cout << "== Normalized histogram ==\n";
  for (int i=0; i<num_bins; ++i) {// Divide
    histogram_percentages[i] = histogram[i] / sum;
    //std::cout << histogram_percentages[i] << "\n";  
  }
  //std::cout << "\n";

  // Try out every bin value in the histogram and pick the best one
  //  - Skip the first bin due to computation below.
  //  - For more resolution, use more bins!
  
  // Compute score for each possible answer.
  std::vector<double> scores(num_bins-1);
  for (int i=0; i<num_bins-1; ++i)
    scores[i] = compute_kittler_illingworth_jt(histogram_percentages, num_bins, min_val, max_val, i+1);

  // Find the lowest score      
  double min_score = scores[0];
  int    min_index = 0;
  for (int i=0; i<num_bins-1; ++i) {
    if (scores[i] < min_score) {
      min_score = scores[i];
      min_index = i+1;
    }
  }
  //std::cout << "Min index = " << min_index << std::endl;
  //std::cout << "bin_width = " << bin_width << std::endl;
  // Compute the final threshold which is below the current bin value
  double threshold = min_val + bin_width*(static_cast<double>(min_index) - 0.5);
  
  return threshold;
} // End function split_histogram_kittler_illingworth


//=========================================================================================

typedef uint16 Sentinel1Type;
typedef float  RadarType;
typedef PixelMask<RadarType> RadarTypeM;


/// Functor to convert a sentinel1 pixel from digital numbers (DN) to decibels (DB)
struct Sentinel1DnToDb : public ReturnFixedType<PixelMask<float> > {
  Sentinel1DnToDb(){}
  PixelMask<float> operator()( PixelMask<Sentinel1Type> p ) const {
    PixelMask<float> output = PixelMask<float>(10*log10(p.child()));
    if (!is_valid(p))
      invalidate(output);
    return output;
  }
};

/// Function to convert a sentinel1 image from digital numbers (DN) to decibels (DB)
template <class ImageT>
UnaryPerPixelView<ImageT,Sentinel1DnToDb>
inline sentinel1_dn_to_db( ImageViewBase<ImageT> const& image) {
  return UnaryPerPixelView<ImageT,Sentinel1DnToDb>( image.impl(), Sentinel1DnToDb() );
}


/// Crop and preprocess the input image in preparation for the sar_martinis algorithm.
template <class ImageT>
void preprocess_sentinel1_image(ImageT const& input_image, BBox2i const& roi,
                                RadarType &global_min, RadarType &global_max,
                                vw::GdalWriteOptions const& write_options,
                                std::string const& temporary_path,
                                double nodata_value,
                                ImageViewRef<RadarTypeM>      & processed_image) {

  // We write the preprocessed image to disk to save on processing time later.

  // The images all appear to fall in the 0-1000 value range.  If this changes more
  //  preprocessing will be needed.


  global_min = 0.0;  // This can be kept constant
  global_max = 35.0; // Computing this exact value has proven to be more expensize than it is worth.

  const double PROC_MIN = 0;   // The paper reccomends these scale values.
  const double PROC_MAX = 400;

  // Preprocess image and write it to disk.
  // - Perform median filter to correct speckles (see section 2.1.4)
  int kernel_size = 3;
  cartography::block_write_gdal_image(temporary_path,
                                      apply_mask(normalize(median_filter_view(sentinel1_dn_to_db(crop(input_image, roi)), 
                                                                              Vector2i(kernel_size, kernel_size)
                                                                             ),
                                                 global_min, global_max, PROC_MIN, PROC_MAX),
                                                 nodata_value
                                                ),
                                      nodata_value,
                                      write_options,
                                      TerminalProgressCallback("vw", "\t--> Preprocessing:"));

  // Return a view of the image on disk so it never gets recomputed.
  processed_image = create_mask(DiskImageView<RadarType>(temporary_path), nodata_value);

  // Update these to reflect the scaled values
  global_min = PROC_MIN;
  global_max = PROC_MAX; 

} // End preprocess_sentinel1_image

/// Class to compute the mean and standard deviation of the four sub-tiles from each input tile.
/// - See the source paper from Martinis et al for details.
/// - Each input tile produces a single pixel in the output image.
template <class ImageT>
class ImageTileMeansView : public ImageViewBase<ImageTileMeansView<ImageT> > {

public: // Definitions

  typedef PixelMask<Vector2f> pixel_type;
  typedef pixel_type          result_type;

private: // Variables

  ImageT const& m_input_image;
  int m_width, m_height, m_tile_size;

public: // Functions

  // Constructor
  ImageTileMeansView( ImageT  const& input_image, int tile_size)
                  : m_input_image(input_image), 
                    m_width (input_image.cols()/tile_size), 
                    m_height(input_image.rows()/tile_size), 
                    m_tile_size(tile_size){}

  inline int32 cols  () const { return m_width;  }
  inline int32 rows  () const { return m_height; }
  inline int32 planes() const { return 1; }

  inline result_type operator()( int32 i, int32 j, int32 p=0 ) const
  { return 0; } // NOT IMPLEMENTED!
 
  typedef ProceduralPixelAccessor<ImageTileMeansView<ImageT> > pixel_accessor;
  inline pixel_accessor origin() const { return pixel_accessor( *this, 0, 0 ); }

  /// Compute the mean and percentage of valid pixels in an entire image
  double mean_and_validity(CropView<ImageView<RadarTypeM> > const& image, double &mean) const {
    mean = 0;
    double count = 0, sum = 0;
    for (int r=0; r<image.rows(); ++r) {
      for (int c=0; c<image.cols(); ++c) {
        if (!is_valid(image(c,r)))
          continue;
        sum   += image(c,r);
        count += 1.0;
      }
    }
    if (count > 0)
      mean = sum / count;
    return count / static_cast<double>(image.rows() * image.cols());
  }
                           

  typedef CropView<ImageView<result_type> > prerasterize_type;
  inline prerasterize_type prerasterize( BBox2i const& bbox ) const {
    
    // This function is only set up for 1x1 output tiles.
    if ((bbox.width() != 1) || (bbox.height() != 1))
      vw_throw(LogicErr() << "ImageTileMeansView only works with 1x1 output tiles!");
    
    // Compute the bbox in the input image
    BBox2i input_bbox = bbox;
    input_bbox *= m_tile_size;

    // Compute the four sub-ROIs
    const int NUM_SUB_ROIS = 4;
    int hw = input_bbox.width() /2;
    int hh = input_bbox.height()/2;
    std::vector<BBox2i> sub_rois(NUM_SUB_ROIS); // ROIs relative to the whole tile
    std::vector<double> means;
    means.reserve(NUM_SUB_ROIS);
    sub_rois[0] = BBox2i(0,  0,  hw, hh); // Top left
    sub_rois[1] = BBox2i(hw, 0,  hw, hh); // Top right
    sub_rois[2] = BBox2i(hw, hh, hw, hh); // Bottom right
    sub_rois[3] = BBox2i(0,  hh, hw, hh); // Bottom left
    
    ImageView<RadarTypeM> section = crop(m_input_image, input_bbox);

    // Don't compute statistics from regions with a lot of bad pixels
    const double MIN_PERCENT_VALID = 0.95;
    
    // Compute the mean in each of the four sub-rois
    double percent_valid;
    for (int i=0; i<NUM_SUB_ROIS; ++i) {
      double mean;
      percent_valid = mean_and_validity(crop(section, sub_rois[i]), mean);
      if (percent_valid >= MIN_PERCENT_VALID)
        means.push_back(mean);
    }
    bool is_valid = false;
    double mean_of_means   = 0;
    double stddev_of_means = 0;
    if (means.size() > 0) {
      // Compute the standard deviation of the means
      // - Currently both are set to zero if all pixels are invalid.
      mean_of_means = math::mean(means);
      is_valid = (mean_of_means > 0);
      if (is_valid)
        stddev_of_means = math::standard_deviation(means, mean_of_means);
    }

    // Set up the output image tile - only 1 by 1 pixel!
    ImageView<result_type> tile(1, 1);
    tile(0, 0)[0] = mean_of_means;
    tile(0, 0)[1] = stddev_of_means;
    if (!is_valid)
      invalidate(tile(0, 0));
    else
      validate(tile(0, 0));

    // Return the tile we created with fake borders to make it look the size of the entire output image
    return prerasterize_type(tile, -bbox.min().x(), -bbox.min().y(), cols(), rows() );

  } // End prerasterize function

 template <class DestT>
 inline void rasterize( DestT const& dest, BBox2i const& bbox ) const {
   vw::rasterize( prerasterize(bbox), dest, bbox );
 }

}; // End class ImageTileMeansView


/// Function to apply ImageTileMeansView on an image.
template <class ImageT>
void generate_tile_means(ImageT input_image, int tile_size,
                         ImageView<PixelMask<Vector2f> > &tile_means_stddevs) {

  // Set up computation object 
  ImageTileMeansView<ImageT> tile_mean_generator(input_image, tile_size);
 
  // Use multiple threads to process the image.
  // - Each block generates one output pixel frome tile_sizeXtile_size input pixels.
  tile_means_stddevs = block_rasterize(tile_mean_generator, Vector2i(1,1));

} // End function generate_tile_means




/// Combine the four fuzzy scores into one final score
/// - Inputs are expected to be in the range 0-1
template <class PixelT1, class PixelT2, class PixelT3, class PixelT4>
struct DefuzzFourFunctor : public ReturnFixedType<float> {

  typedef PixelMask<float> result_type;

  result_type operator()(PixelT1 const& p1, PixelT2 const& p2, PixelT3 const& p3, PixelT4 const& p4) const {
    // If any input score is zero, the output score is zero.
    if ((p1 == 0) || (p2 == 0) || (p3 == 0) || (p4 == 0))
      return result_type(0.0);

    float mean = (p1 + p2 + p3 + p4) / 4.0;
    return result_type(mean);
  }
}; // End class DefuzzFourFunctor

/// Combine only two fuzzy scores into one final score
/// - Inputs are expected to be in the range 0-1
template <class PixelT1, class PixelT2>
struct DefuzzTwoFunctor : public ReturnFixedType<float> {

  typedef PixelMask<float> result_type;

  result_type operator()(PixelT1 const& p1, PixelT2 const& p2) const {
    // If any input score is zero, the output score is zero.
    if ((p1 == 0) || (p2 == 0))
      return result_type(0.0);

    float mean = (p1 + p2) / 2.0;
    return result_type(mean);
  }
}; // End class DefuzzTwoFunctor

//============================================================================

// TODO: Move these up into VW.


// Define comparison function (just compare the first elements)
template <typename T>
bool less_than_function_pair_first( const T& a, const T& b) { 
  return a.first < b.first; 
}


/// Generate a sorted list of indices into an input vector.
/// - Intended for simple vectors of ints or floats.
template <typename T>
void sort_vector_indices(std::vector<T> const& v, std::vector<size_t> &indices) {

  // Copy data into a new vector with indices
  typedef std::pair<T, size_t> PairV;
  std::vector<PairV> v_new(v.size());
  for (size_t i=0; i<v.size(); ++i) {
    v_new[i].first  = v[i];
    v_new[i].second = i;
  }

  // Run sort and store the sorted indices
  std::sort(v_new.begin(), v_new.end(), less_than_function_pair_first<PairV>);
  indices.resize(v.size());
  for (size_t i=0; i<v.size(); ++i)
    indices[i] = v_new[i].second;
}


/// Converts a normal vector into a slope angle in degrees.
struct GetAngleFunc : public ReturnFixedType<PixelMask<float> > {
  PixelMask<float> operator() (PixelMask<Vector3f> const& pix) const {
    if (is_transparent(pix))
      return PixelMask<float>();
    return PixelMask<float>(RAD_TO_DEG*acos(fabs(dot_prod(pix.child(),Vector3(0,0,1)))));
  }
};
template <class ViewT>
UnaryPerPixelView<ViewT, GetAngleFunc> get_angle(ImageViewBase<ViewT> const& view) {
  return UnaryPerPixelView<ViewT, GetAngleFunc>(view.impl(), GetAngleFunc());
}



/// Select the best N tiles to use for computing the global water threshold
/// - This follows the process in the source paper.
size_t select_best_tiles(ImageView<PixelMask<Vector2f> >    & tile_means_stddevs,
                         std::vector<int>                   & kept_tile_indices,
                         vw::GdalWriteOptions const& write_options,
                         bool debug) {

  std::cout << "Computing global tile statistics...\n";

  // For convenient reference of the two input channels
  ImageViewRef<RadarTypeM> tile_means   = copy_mask(select_channel(tile_means_stddevs, 0), tile_means_stddevs);
  ImageViewRef<RadarTypeM> tile_stddevs = copy_mask(select_channel(tile_means_stddevs, 1), tile_means_stddevs);

  // Compute the global mean
  double global_mean = mean_channel_value(tile_means);
  if (debug)
    std::cout << "Computed global mean: " << global_mean << "\n";

  // Compute 95% quantile standard deviation
  double stddev_min, stddev_max;
  find_image_min_max(tile_stddevs, stddev_min, stddev_max);
  
  int    num_bins = 255;
  std::vector<double> hist;
  histogram(tile_stddevs, num_bins, stddev_min, stddev_max, hist);
  const double TILE_STDDEV_PERCENTILE_CUTOFF = 0.95;
  int    bin = get_histogram_percentile(hist, TILE_STDDEV_PERCENTILE_CUTOFF);
  double bin_width      = (stddev_max - stddev_min)/static_cast<double>(num_bins);
  double std_dev_cutoff = stddev_min + bin_width*bin;
  if (debug)
    std::cout << "Computed standard deviation cutoff: " << std_dev_cutoff << "\n";

  // Select the tiles with the highest STD values (N')
  ImageView<uint8> kept_tile_display;
  kept_tile_display.set_size(tile_means.cols(), tile_means.rows());
  std::vector<Vector2i> n_prime_tiles;
  std::vector<double  > n_prime_std_dev, n_prime_mean;
  double mean_of_selected = 0;
  for (int r=0; r<tile_stddevs.rows(); ++r) {
    for (int c=0; c<tile_stddevs.cols(); ++c) {
      // The tile must have a high stddev and also be below the global mean
      //  since water tends to be darker than land.
      if ( (tile_stddevs(c,r) > std_dev_cutoff) &&
           (tile_means  (c,r) < global_mean   )   ) {
        n_prime_tiles.push_back(Vector2i(c,r));
        n_prime_std_dev.push_back(tile_stddevs(c,r));
        n_prime_mean.push_back   (tile_means  (c,r));
        mean_of_selected += tile_means(c,r);

        kept_tile_display(c,r) = 255;
      }
    }
  }

  size_t num_tiles_kept = n_prime_tiles.size();
  mean_of_selected /= static_cast<double>(num_tiles_kept);
  std::cout << "Selected " << num_tiles_kept << " initial tiles with mean value "
            << mean_of_selected << std::endl;
  
  if (debug)
    block_write_gdal_image("initial_kept_tiles.tif", kept_tile_display, write_options); // DEBUG

  if (n_prime_tiles.empty())
    vw_throw(LogicErr() << "No tiles left after std_dev filtering!");

  // Cap the number of selected tiles
  const size_t MAX_NUM_TILES = 5; // From the paper  
  
  // If already at/below the cap we are finished.
  if (num_tiles_kept <= MAX_NUM_TILES) {
    // Copy linear indices of kept tiles
    kept_tile_indices.resize(n_prime_tiles.size());
    for (size_t i=0; i<n_prime_tiles.size(); ++i)
      kept_tile_indices[i] = n_prime_tiles[i][1]*tile_stddevs.cols() + n_prime_tiles[i][0];
    return n_prime_tiles.size();
  }
  // Reset this image
  fill(kept_tile_display, 0);
  
  // Keep the N tiles with the highest standard deviation
  std::vector<size_t> indices;
  sort_vector_indices(n_prime_std_dev, indices); // Sorts by STD low to high
  
  kept_tile_indices.resize(MAX_NUM_TILES);
  size_t max_index = num_tiles_kept - 1;
  size_t index_out=0;
  for (size_t i=max_index; i>0; --i) {
    size_t index_to_keep = indices[i];
    
    // Don't keep tiles above the mean of the initial set of kept tiles.
    // - This was in the paper but in our case it made things worse
    //double this_mean = n_prime_mean[index_to_keep];
    //if (this_mean > mean_of_selected)
    //  continue;
      
    Vector2i index_2d = n_prime_tiles[index_to_keep];
    kept_tile_indices[index_out] = index_2d[1]*tile_stddevs.cols() + index_2d[0]; // Get linear index
    kept_tile_display(index_2d[0], index_2d[1]) = 255;
    
    // Quit when we have kept the desired number of tiles.
    ++index_out;
    if (index_out >= MAX_NUM_TILES)
      break;
  }

  if (debug)
    block_write_gdal_image("final_kept_tiles.tif", kept_tile_display, write_options); // DEBUG
  
  std::cout << "Reduced to " << kept_tile_indices.size() << " kept tiles.\n";
  
  return kept_tile_indices.size();
} // End function select_best_tiles


/// Use the kept tiles to compute a global water threshold for the entire image.
bool compute_global_threshold(ImageViewRef<RadarTypeM> const& preprocessed_image, 
                              std::vector<int>         const& kept_tile_indices,
                              std::vector<BBox2i>      const& large_tile_boxes,
                              float global_min, float global_max,
                              double &threshold_mean) {

  // For each selected tile, find optimal threshold using Kittler-Illingworth method.
  // - Could compute in parallel, but this is not a big bottleneck.
  int num_bins = 255;
  const size_t num_tiles = kept_tile_indices.size();
  std::vector<double> optimal_tile_thresholds(num_tiles);
  double dmin = static_cast<double>(global_min);
  double dmax = static_cast<double>(global_max);
  //std::cout << "Histogram min = " << dmin << std::endl;
  //std::cout << "Histogram max = " << dmax << std::endl;
  threshold_mean = 0;
  for (size_t i=0; i<num_tiles; ++i) {
  
    // Compute tile histogram     
    BBox2i   roi = large_tile_boxes[kept_tile_indices[i]]; // The ROIs are stored row-first.
    std::vector<double> hist;
    histogram(crop(preprocessed_image, roi), num_bins, dmin, dmax, hist);
    
    // Compute optimal split
    optimal_tile_thresholds[i] = split_histogram_kittler_illingworth(hist, num_bins, global_min, global_max);
    //std::cout << "For ROI " << roi << ", computed threshold " << optimal_tile_thresholds[i] << std::endl;
    threshold_mean += optimal_tile_thresholds[i];
  }
  threshold_mean /= static_cast<double>(num_tiles);

  // This number is in the rescaled units of the preprocessed image.
  const double MAX_THRESHOLD_STDDEV = 10.0;
  
  double threshold_stddev = math::standard_deviation(optimal_tile_thresholds, threshold_mean);
  
  std::cout << "Mean of tile thresholds: " << threshold_mean   << std::endl;
  std::cout << "STD  of tile thresholds: " << threshold_stddev << std::endl;

  // If our selected tiles are too dissimilar our final answer is probably not trustworthy.
  if (threshold_stddev > MAX_THRESHOLD_STDDEV) {
    std::cout << "WARNING: Standard deviation of computed thresholds exceeds maximum value of "
              << MAX_THRESHOLD_STDDEV << std::endl;
    return false;
  }
  // The paper has some backup methods of computing the threshold but it is probably safer just
  //  to declare failure in those instances.

  return true;  

} // End compute_global_threshold

//============================================================================



/** Main function of algorithm from:
      Martinis, Sandro, Jens Kersten, and Andre Twele. 
      "A fully automated TerraSAR-X based flood service." 
      ISPRS Journal of Photogrammetry and Remote Sensing 104 (2015): 203-212.
*/
void sar_martinis(std::string const& input_image_path, std::string const& output_path,
                  vw::GdalWriteOptions const& write_options,
                  std::string dem_path="", bool debug=false, int tile_size = 512,
                  double sensitivity = 1.0) {

  // Set up paths to temporary files in the output folder.
  boost::filesystem::path fs_path(output_path); 
  fs_path = fs_path.parent_path().parent_path();
  fs_path /= "preprocessed_image.tif";
  std::string preprocessed_image_path = fs_path.string();
  
  BBox2i roi = bounding_box(DiskImageView<Sentinel1Type>(input_image_path));
 
  // Load the georeference from the input image
  // - The input image won't be georeferenced unless it goes through gdalwarp.
  cartography::GeoReference georef;
  bool have_georef = cartography::read_georeference(georef, input_image_path);
  if (!have_georef)
    vw_throw(ArgumentErr() << "Failed to read image georeference!");
  georef = crop(georef, roi); // Account for the input ROI
  if (debug)
    std::cout << "Read georeference: " << georef << std::endl;
  
  double input_meters_per_pixel = cartography::get_image_meters_per_pixel(roi.width(), roi.height(), georef);
  std::cout << "Computed image pixel resolution in meters: " << input_meters_per_pixel << std::endl;
  
  // Read nodata value
  double input_nodata_value = 0;
  read_nodata_val(input_image_path, input_nodata_value);
   
  // Compute the min and max values of the image
  std::cout << "Preprocessing...\n";

  // Apply any needed preprocessing to the image
  double internal_nodata_value = -32767.0; // For use with our temporary files.
  ImageViewRef<RadarTypeM> preprocessed_image; // This is a view of the image on disk.
  RadarType global_min, global_max;
  preprocess_sentinel1_image(create_mask(DiskImageView<Sentinel1Type>(input_image_path), input_nodata_value), 
                             roi, global_min, global_max, 
                             write_options, preprocessed_image_path, internal_nodata_value, preprocessed_image);
  
  // Perform the automatic threshold computation step.
  // - If it fails, try a second time with half-size tiles.
  const int MAX_TILE_ATTEMPTS = 2;
  bool   tile_thresh_success = false;
  double threshold_mean;
  for (int num_tile_attempts = 0; num_tile_attempts < MAX_TILE_ATTEMPTS; ++num_tile_attempts) {
    
    // Generate vector of BBoxes for each tile in the input image (S+)
    std::vector<BBox2i> large_tile_boxes = subdivide_bbox(bounding_box(preprocessed_image), 
                                                          tile_size, tile_size, false);
    
    std::cout << "Computing tile means...\n";

    // For each tile compute the mean value and the standard deviation of the four sub-tiles.
    ImageView<PixelMask<Vector2f> > tile_means_stddevs; // These are much smaller than the input image
    generate_tile_means(preprocessed_image, tile_size, tile_means_stddevs);

    if (debug) {
      std::cout << "Writing DEBUG images...\n";
      float debug_nodata = -32768.0;
      block_write_gdal_image("tile_means.tif", 
        apply_mask(select_channel(tile_means_stddevs, 0), debug_nodata), debug_nodata, write_options);
      block_write_gdal_image("tile_stddevs.tif", 
        apply_mask(select_channel(tile_means_stddevs, 1), debug_nodata), debug_nodata, write_options);
    }


    // Select the tiles that we will use to compute the optimal global threshold.
    std::vector<int> kept_tile_indices;
    size_t num_tiles_kept = select_best_tiles(tile_means_stddevs, kept_tile_indices, write_options, debug);

    if (num_tiles_kept > 0) {
      // Use the selected tiles to compute the optimal image threshold.
      tile_thresh_success = compute_global_threshold(preprocessed_image, kept_tile_indices, large_tile_boxes,
                                                     global_min, global_max, threshold_mean);
    }
    
    // Exit the loop if we succeeded, otherwise try one more time with half the tile size.
    if (tile_thresh_success)
      break;
    tile_size /= 2;
    if (num_tile_attempts == 0)
      std::cout << "Making one more attempt to auto-compute the threshold using smaller image tiles...\n";
  } // End of loop to attemp to automaticall compute threshold.
  
  // Quit if we did not compute a good threshold.
  if (!tile_thresh_success) {
    vw_throw(ArgumentErr() << "Unable to compute a good water threshold for this image!");
  }
  
  // This will mask the water pixels, setting water pixels to 255, land pixels to 1, and invalid pixels to 0.
  ImageViewRef<RadarTypeM> raw_water = threshold(preprocessed_image, threshold_mean, 
                                                 FLOOD_DETECT_WATER, FLOOD_DETECT_LAND);

  // Apply the initial threshold to the image and save it to disk!
  std::string initial_water_detect_path = "initial_water_detect.tif";
  block_write_gdal_image(initial_water_detect_path,
                         pixel_cast<uint8>(apply_mask(raw_water, FLOOD_DETECT_NODATA)),
                         have_georef, georef,
                         true, FLOOD_DETECT_NODATA, // Choose the nodata value
                         write_options,
                         TerminalProgressCallback("vw", "\t--> Applying initial threshold:"));


  // Get information needed for fuzzy logic results filtering

  // Write out an image containing the water blob size at each pixel, then read it back in
  // as needed to avoid recomputing the expensive blob computations.
  // - In order to parallelize this step, blob computations are approximated.
  //   Setting TILE_EXPAND (in pixels) to a larger number improves the approximation.
  
  const double MIN_BLOB_SIZE_METERS = 1000; // These numbers are larger than in the paper but 
  const double MAX_BLOB_SIZE_METERS = 5000; //  seem to work better.
  const int    TILE_EXPAND   = 256;
  
  uint32 min_blob_size = MIN_BLOB_SIZE_METERS / input_meters_per_pixel;
  uint32 max_blob_size = MAX_BLOB_SIZE_METERS / input_meters_per_pixel;

  std::cout << "Image meters per pixel: " << input_meters_per_pixel << std::endl;
  if (debug) {
    std::cout << "Min blob size pixels = " << min_blob_size << std::endl;
    std::cout << "Max blob size pixels = " << max_blob_size << std::endl;
  }

  // Compute the blob sizes and record directly to disk so they are not recomputed.
  std::string blobs_path = "blob_sizes.tif";
  const uint32 BLOBS_NODATA = 0;
  block_write_gdal_image(blobs_path,
                         get_blob_sizes(create_mask_less_or_equal(DiskImageView<uint8>(initial_water_detect_path), FLOOD_DETECT_LAND),
                                        TILE_EXPAND, max_blob_size),
                         have_georef, georef,
                         true, BLOBS_NODATA,
                         write_options,
                         TerminalProgressCallback("vw", "\t--> Counting blob sizes:"));
  // Get handle of image on disk.
  DiskImageView<uint32> blob_sizes(blobs_path);

  // Compute these statistics at low resolution for speed.
  int DEM_STATS_SUBSAMPLE_FACTOR = 10;

  std::cout << "Computing mean of flooded regions...\n";

  // Load a low-res version of our initial water results.
  ImageViewRef<RadarTypeM> low_res_raw_water =  
          copy_mask(subsample(preprocessed_image, DEM_STATS_SUBSAMPLE_FACTOR),
                    subsample(create_mask_less_or_equal(DiskImageView<uint8>(initial_water_detect_path), 
                                                                                     FLOOD_DETECT_LAND), 
                                                                 DEM_STATS_SUBSAMPLE_FACTOR)
                   );
  cartography::GeoReference low_res_georef = resample(georef, 1.0/DEM_STATS_SUBSAMPLE_FACTOR);

  // Compute mean radar value of pixels under initial water threshold
  // - This is also computed at a lower resolution to increase speed.
  // - Could do full res with a multi-threaded implementation.
  double mean_raw_water_value = mean_channel_value(low_res_raw_water);

  std::cout << "Mean value of flooded regions = " << mean_raw_water_value << std::endl;

  //block_write_gdal_image("low_res_water.tif", low_res_raw_water, 
  //                       have_georef, low_res_georef, true, -9999,
  //                       write_options, TerminalProgressCallback("vw", "\t--> DEBUG:"));

  // Go ahead and set up the fuzzy logic results that can be done without a DEM

  // Compute fuzzy classifications on four categories
  typedef PixelMask<float> FuzzyPixelType;
  typedef FuzzyMembershipSFunctor<RadarTypeM> FuzzyFunctorS;
  typedef FuzzyMembershipZFunctor<RadarTypeM> FuzzyFunctorZ;
  ImageViewRef<FuzzyPixelType> defuzzed;
 
  // SAR
  FuzzyFunctorZ radar_fuzz_functor(mean_raw_water_value, threshold_mean);
  ImageViewRef<FuzzyPixelType> radar_fuzz = per_pixel_view(preprocessed_image, radar_fuzz_functor);

  // Body size 
  FuzzyFunctorS blob_fuzz_functor(min_blob_size, max_blob_size);
  ImageViewRef<FuzzyPixelType> blob_fuzz = per_pixel_view(blob_sizes, blob_fuzz_functor);

  if (debug) {
    const double temp_nodata = -1;
    block_write_gdal_image("radar_fuzz.tif", apply_mask(radar_fuzz, temp_nodata),
                           have_georef, georef, true, temp_nodata,
                           write_options, TerminalProgressCallback("vw", "\t--> radar_fuzz:"));
    block_write_gdal_image("blob_fuzz.tif", apply_mask(blob_fuzz, temp_nodata),
                           have_georef, georef, true, temp_nodata,
                           write_options, TerminalProgressCallback("vw", "\t--> blob_fuzz:"));
  }


  if (dem_path.empty()) {
    // Handle the case where a DEM was not provided
  
    std::cout << "No input DEM file provided, finishing calculations without one.\n";
  
    // Defuzz only two fuzzy classifiers and compare to a fixed threshold in the 0-1 range.
    typedef DefuzzTwoFunctor<FuzzyPixelType, FuzzyPixelType> DefuzzFunctorType;
    defuzzed = per_pixel_view(radar_fuzz, blob_fuzz, DefuzzFunctorType());
  
  } else {
    // Handle the case where a DEM was provided

    typedef PixelMask<float> DemPixelType;

    // Should be safe to use this as a DEM nodata value!
    double dem_nodata_value = -3.4028234663852886e+38;
    read_nodata_val(dem_path, dem_nodata_value); 

    DiskImageView<float> dem(dem_path);

    cartography::GeoReference dem_georef;
    if (!cartography::read_georeference(dem_georef, dem_path))
      vw_throw(ArgumentErr() << "Failed to read DEM georeference!");

    std::cout << "Input DEM file found, using it to complete the calculations.\n";

    // Generate a low-resolution DEM masked by the initial flood detection
    // - This is used to compute image-wide statistics in a more reasonable amount of time
    ImageView<DemPixelType> low_res_dem = subsample(create_mask(dem, dem_nodata_value), DEM_STATS_SUBSAMPLE_FACTOR);
    cartography::GeoReference low_res_dem_georef = resample(dem_georef, 1.0/DEM_STATS_SUBSAMPLE_FACTOR);

    ImageViewRef<PixelMask<float> > low_res_dem_in_image_coords = 
      cartography::geo_transform(low_res_dem, low_res_dem_georef, low_res_georef,
                                 low_res_raw_water.cols(), low_res_raw_water.rows(), ConstantEdgeExtension());
    ImageViewRef<PixelMask<float> > dem_in_image_coords = 
      cartography::geo_transform(create_mask(dem, dem_nodata_value), dem_georef, georef,
                                 preprocessed_image.cols(), preprocessed_image.rows(), ConstantEdgeExtension());

    if (debug) {
      block_write_gdal_image("geotrans_dem.tif",
                              apply_mask(dem_in_image_coords, dem_nodata_value),
                             have_georef, low_res_georef, true, dem_nodata_value,
                             write_options, TerminalProgressCallback("vw", "\t--> geotrans:"));
    }

    // Now go through and compute statistics across the water covered locations of the DEM
    StdDevAccumulator<float> stddev_functor;
    FunctorMaskWrapper<StdDevAccumulator<float>, PixelMask<float> > dem_stats_functor(stddev_functor);
    for_each_pixel(copy_mask(low_res_dem_in_image_coords, low_res_raw_water), dem_stats_functor);
    
    float mean_water_height   = dem_stats_functor.child().value();
    float stddev_water_height = dem_stats_functor.child().mean();

    std::cout << "Mean height of flooded regions = " << mean_water_height 
              << ", and stddev = " << stddev_water_height << std::endl;

    // Elevation
    // - The max value looks a little weird but it comes straight from the paper.
    const double high_height = mean_water_height + stddev_water_height*(stddev_water_height + 3.5);
    FuzzyFunctorZ height_fuzz_functor(mean_water_height, high_height);
    ImageViewRef<FuzzyPixelType> height_fuzz = per_pixel_view(dem_in_image_coords, height_fuzz_functor);
    
    // Slope
    const double degrees_low  = 0;
    const double degrees_high = 15;
    FuzzyFunctorZ slope_fuzz_functor(degrees_low, degrees_high);
    ImageViewRef<FuzzyPixelType> slope_fuzz = per_pixel_view(get_angle(compute_normals(dem_in_image_coords, 1.0, 1.0)), slope_fuzz_functor);

    if (debug) {
      block_write_gdal_image("height_fuzz.tif", apply_mask(height_fuzz, dem_nodata_value),
                             have_georef, georef, true, dem_nodata_value,
                             write_options, TerminalProgressCallback("vw", "\t--> height_fuzz:"));
      block_write_gdal_image("slope_fuzz.tif", apply_mask(slope_fuzz, dem_nodata_value),
                             have_georef, georef, true, dem_nodata_value,
                             write_options, TerminalProgressCallback("vw", "\t--> slope_fuzz:"));
    }
                           
    // Defuzz all four fuzzy classifiers and compare to a fixed threshold in the 0-1 range.
    typedef DefuzzFourFunctor<FuzzyPixelType, FuzzyPixelType, FuzzyPixelType, FuzzyPixelType> DefuzzFunctorType;
    defuzzed = per_pixel_view(radar_fuzz, height_fuzz, slope_fuzz, blob_fuzz, DefuzzFunctorType());
    
  } // End of case where a DEM was provided

  if (debug) {                         
    block_write_gdal_image("defuzzed.tif", apply_mask(defuzzed, -999),
                           have_georef, georef, false, 999,
                           write_options, TerminalProgressCallback("vw", "\t--> Defuzz:"));
  }

  // Perform two-level flood fill of the defuzzed image and write it to disk.
  // - The mask is added back in at this point.
  const double FINAL_FLOOD_THRESHOLD = sensitivity * 0.6;
  const double WATER_GROW_THRESHOLD  = sensitivity * 0.45;
  block_write_gdal_image(output_path,
                         apply_mask(
                           copy_mask(
                             two_threshold_fill(defuzzed, TILE_EXPAND, FINAL_FLOOD_THRESHOLD, 
                                                WATER_GROW_THRESHOLD, FLOOD_DETECT_LAND, FLOOD_DETECT_WATER),
                             create_mask(DiskImageView<uint8>(initial_water_detect_path, FLOOD_DETECT_NODATA))
                           ),
                           FLOOD_DETECT_NODATA
                         ),      
                         true, georef,
                         true, FLOOD_DETECT_NODATA,
                         write_options,
                         TerminalProgressCallback("vw", "\t--> Generating final output:"));

  if (debug == false) {
    // Clean up temporary image files
    std::remove(initial_water_detect_path.c_str());
    std::remove(preprocessed_image_path.c_str());
    std::remove(blobs_path.c_str());
  }

} // End function sar_martinis


}} // end namespace vw/radar
#endif

