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
//#include <vw/Math/Statistics.h>
#include <vw/Image/Statistics.h>
#include <vw/Image/PixelMask.h>
#include <vw/Image/Manipulation.h>
#include <vw/Image/Transform.h>
#include <vw/Image/ImageIO.h>
#include <vw/Image/PixelMask.h>
#include <vw/Image/MaskViews.h>
#include <vw/Image/UtilityViews.h>
#include <vw/FileIO/DiskImageView.h>
#include <vw/Cartography/GeoReferenceUtils.h>

/**
  Tools for processing radar data.
  Add more files as needed.
  
  TODO: Break this file up?
*/

namespace vw {
namespace radar {







// TODO: Move all of this!


/// Generic class for implementing a function that generates a pixel value based
///  on a window around a pixel.
/// - TODO: Could ConvolutionView use this class?
template <class ImageT, class FuncT, class EdgeT>
class WindowFunctionView : public ImageViewBase<WindowFunctionView<ImageT,FuncT,EdgeT> >
{
private:
  ImageT   m_image;
  EdgeT    m_edge;     ///< Edge extension type
  FuncT    m_functor;  ///< Functor that operates on each window.
  Vector2i m_window_size;
  int      m_half_width;
  int      m_half_height;

public:
  typedef typename ImageT::pixel_type pixel_type;  ///< The pixel type of the image view.
  typedef pixel_type                  result_type; ///< We compute the result, so we return by value.
  typedef ProceduralPixelAccessor<WindowFunctionView<ImageT, FuncT, EdgeT> > 
                                      pixel_accessor; ///< The view's pixel_accessor type.

  /// Constructs a ConvolutionView with the given image and kernel and 
  /// with the origin of the kernel located at the point (ci,cj).
  WindowFunctionView( ImageT const& image, Vector2i window_size,
                      FuncT  const& functor, 
                      EdgeT  const& edge = EdgeT() )
    : m_image(image), m_edge(edge), m_functor(functor), m_window_size(window_size) {
    m_half_width  = m_window_size[0]/2;
    m_half_height = m_window_size[1]/2;
  }

  inline int32 cols  () const { return m_image.cols  (); }
  inline int32 rows  () const { return m_image.rows  (); }
  inline int32 planes() const { return m_image.planes(); }

  /// Returns a pixel_accessor pointing to the origin.
  inline pixel_accessor origin() const { return pixel_accessor( *this ); }

  /// Returns the pixel at the given position in the given plane.
  inline result_type operator()( int32 x, int32 y, int32 p=0 ) const {
    BBox2i roi(x-m_half_width, y-m_half_height, m_window_size[0], m_window_size[1]);
    return m_functor(edge_extend(m_image, roi, m_edge));
  }

  // Edge extension is done in the prerasterize function so the returned type does not need edge extension
  typedef WindowFunctionView<CropView<ImageView<typename ImageT::pixel_type> >, 
                             FuncT, NoEdgeExtension> prerasterize_type;
  inline prerasterize_type prerasterize( BBox2i const& bbox ) const {
    // Compute the required base of support for the input bounding box
    BBox2i src_bbox( bbox.min().x() - m_half_width, 
                     bbox.min().y() - m_half_height,
                     bbox.width () + m_window_size[0]-1, 
                     bbox.height() + m_window_size[1]-1 );
    // Take an edge extended image view of the input support region
    ImageView<typename ImageT::pixel_type> src = edge_extend(m_image, src_bbox, m_edge);
    // Use the crop trick to fake that the support region is the same size as the entire image.
    return prerasterize_type( crop(src, -src_bbox.min().x(), -src_bbox.min().y(), m_image.cols(), m_image.rows()),
                              m_window_size, m_functor, NoEdgeExtension() );
  }

  template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const {
    vw::rasterize( prerasterize(bbox), dest, bbox );
  }
}; // End class WindowFunctionView

/// Functor to find the median of an input image.
/// - Usually used with WindowFunctionView.
template <typename ImageT>
struct WindowMedianFunctor {

  mutable std::vector<double> m_values; ///< Persistent storage location

  /// Constructor
  WindowMedianFunctor(Vector2i window_size) {
    m_values.resize(window_size[0]*window_size[1]);
  }

  /// Returns the median pixel of the provided input image.
  /// - Generally the input image will be a cropped view of a whole image.
  template <class T>
  typename ImageT::pixel_type operator()(ImageViewBase<T> const& image) const {

    // Loop through the kernel and collect the values
    int index = 0;
    for (int r=0; r<image.impl().rows(); ++r) {
      for (int c=0; c<image.impl().cols(); ++c) {
        if (is_valid(image.impl()(c,r))) {
          m_values[index] = image.impl()(c,r);
          ++index;
        }
      }
    }
    if (index == 0) { // All pixels invalid!
      typename ImageT::pixel_type result(0);
      invalidate(result);
      return result;
    }

    // Now that we have all the values, compute the median.
    double median;
    if (index == static_cast<int>(m_values.size())-1) // No invalid pixels
      median = math::destructive_median(m_values);
    else { // Invalid pixels
      // Resize the vector twice so we can call the median function
      size_t full_size = m_values.size();
      m_values.resize(index + 1);
      median = math::destructive_median(m_values);
      m_values.resize(full_size);
    }
    return typename ImageT::pixel_type(median);
  }
};

/// Apply a median filter to an input image
template <class ImageT, class EdgeT>
WindowFunctionView<ImageT, WindowMedianFunctor<ImageT>, EdgeT> 
median_view(ImageT const& image, Vector2i window_size, EdgeT edge) {
  typedef WindowFunctionView<ImageT, WindowMedianFunctor<ImageT>, EdgeT> return_type;
  WindowMedianFunctor<ImageT> functor(window_size);
  return return_type(image, window_size, functor, edge);
}


/*

/// Standard Z shaped fuzzy math membership function between a and b
double fuzzMemZ(double x, double a, double b) {
    
    double c = (a+b)/2.0;
    
    if (x < a)
        return 1.0;
    if (x < c)
        return 1.0 - 2*pow(((x-a)/(b-a)), 2.0);
    if (x < b)
        return 2*pow(((x-b)/(b-a)), 2.0);
    return 0.0;
}

/// Standard S shaped fuzzy math membership function between a and b
double fuzzMemS(double x, double a, double b) {

    double c = (a+b)/2.0;

    if (x < a)
        return 0.0;
    if (x < c)
        return 2*pow(((x-a)/(b-a)), 2.0);
    if (x < b)
        return 1.0 - 2*pow(((x-b)/(b-a)), 2.0);
    return 1.0;
}


/// Changes the scaling of a number from one range to a new one.
double rescaleNumber(double num, double currMin, double currMax, double newMin, double newMax) {
    
    double currRange = currMax - currMin;
    double newRange  = newMax - newMin;
    double scaled    = (num - currMin) / currRange;
    double output    = scaled*newRange + newMin;
    return output;
}

*/
// TODO: Consolidate these functions in vw/math/statistics
// TODO: Really need a histogram class!

/// As part of the Kittler/Illingworth method, compute J(T) for a given T
/// - Input should be a percentage histogram
double computeKIJT(std::vector<double> const& histogram,
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
} // End function computeKIJT

/// Tries to compute an optimal histogram threshold using the Kittler/Illingworth method
double splitHistogramKittlerIllingworth(std::vector<double> const& histogram,
                                     int num_bins, double min_val, double max_val) {

  // Normalize the histogram (each bin is now a percentage)
  std::vector<double> histogram_percentages(num_bins);
  double sum = 0;
  for (int i=0; i<num_bins; ++i) // Compute sum
    sum += histogram[i];
  for (int i=0; i<num_bins; ++i) // Divide
    histogram_percentages[i] = histogram[i] / sum;


  // Try out every bin value in the histogram and pick the best one
  //  - Skip the first bin due to computation below.
  //  - For more resolution, use more bins!
  //  - TODO: write up a smarter solver for this.
  
  // Compute score for each possible answer.
  std::vector<double> scores(num_bins-1);
  for (int i=0; i<num_bins-1; ++i)
    scores[i] = computeKIJT(histogram_percentages, num_bins, min_val, max_val, i+1);

  // Find the lowest score      
  double min_score = scores[0];
  int    min_index = 0;
  for (int i=0; i<num_bins-1; ++i) {
    if (scores[i] < min_score) {
      min_score = scores[i];
      min_index = i+1;
    }
  }
  // Compute the final threshold which is below the current bin value
  double bin_width = (max_val - min_val) / num_bins;
  double threshold = min_val + bin_width*(static_cast<double>(min_index) - 0.5);
  
  return threshold;
} // End function splitHistogramKittlerIllingworth





//=========================================================================================

// TODO: What data type to use?
typedef uint16 SarType;
typedef PixelMask<SarType> SarTypeM;

/// Crop and preprocess the input image in preparation for the sar_martinis algorithm.
void preprocess_sentinel1_image(ImageView<SarType> const& input_image, 
                                double nodata_value, BBox2i const& roi,
                                SarType global_min, SarType global_max,
                                cartography::GdalWriteOptions const& write_options,
                                ImageViewRef<SarTypeM>      & processed_image) {

  // Currently we write the preprocessed image to disk, but maybe in the future
  // we should not.

  // TODO: Apply processing!
  // TODO: If we use a disk image, need to clean it up later!

  // --> For now, leave the image as-is knowing that most pixels fall in the 0-1000 range.
  //  ---> Later, apply a log scale or something.

  // TODO: Detect if image already exists, don't reprocess!
  std::cout << "Skipping preprocess already on disk!!!!\n";
/*

  TODO: Handle the mask in the filter!
  // Perform median filter to correct speckles (see section 2.1.4)
  int kernel_size = 3;
  std::string temp_image_path = "preprocessed_image.tif";
  cartography::block_write_gdal_image(temp_image_path,
                                      median_view(crop(input_image, roi), 
                                                  Vector2i(kernel_size, kernel_size),
                                                  ConstantEdgeExtension()),
                                      write_options,
                                      TerminalProgressCallback("vw", "\t--> Preprocessing:"));
*/
  // Return a view of the image on disk for easy access
  processed_image = create_mask(DiskImageView<SarType>("preprocessed_image.tif"), nodata_value);

/*  
  // EE does most of the same preprocessing as the paper but we still need to
  //  duplicate the 0 to 400 scale they used.
  PROC_global_min =   0.0
  PROC_global_max = 400.0
  minmax = rawImage.reduceRegion(ee.Reducer.minMax(), domain.bounds, scale=BASE_RES).getInfo();
  minVal = [value for key, value in minmax.items() if 'min' in key.lower()][0]
  maxVal = [value for key, value in minmax.items() if 'max' in key.lower()][0]
  radarImage = rawImage.unitScale(minVal, maxVal).multiply(PROC_global_max)
*/
} // End preprocess_sentinel1_image

/// Splits up one large BBox into a grid of smaller BBoxes.
/// - Returns the number of BBoxes created.
/// - If include_partials is set to false, incomplete BBoxes will be discarded.
int divide_roi(BBox2i const& full_roi, int size,
               std::vector<std::vector<BBox2i> > &new_rois, bool include_partials = true) {
               
  // Compute the number of boxes
  double dsize = static_cast<double>(size);
  int num_boxes_x, num_boxes_y;
  if (include_partials) {
    num_boxes_x = ceil(static_cast<double>(full_roi.width())  / dsize);
    num_boxes_y = ceil(static_cast<double>(full_roi.height()) / dsize);
  } else { // Discard partial boxes
    num_boxes_x = floor(static_cast<double>(full_roi.width())  / dsize);
    num_boxes_y = floor(static_cast<double>(full_roi.height()) / dsize);
  }

  std::cout << "Dividing up ROI " << full_roi << std::endl;
  printf("Num boxes: %d, %d\n", num_boxes_x, num_boxes_y);

  // Generate all of the boxes
  int max_x = full_roi.max()[0];
  int max_y = full_roi.max()[1];
  int num_boxes = num_boxes_y*num_boxes_x;
  new_rois.resize(num_boxes_y);
  for (int r=0; r<num_boxes_y; ++r) {
    new_rois[r].resize(num_boxes_x);
    for (int c=0; c<num_boxes_x; ++c) {
      // Set up each box being mindful of partial boxes at the edges
      int x      = c*size;
      int y      = r*size;
      int height = size;
      int width  = size;
      if (y+height > max_y) height = max_y - y;
      if (x+width  > max_x) width  = max_x - x;
      new_rois[r][c] = BBox2i(x, y, width, height);
      //std::cout << "--> New ROI: " << new_rois[r][c] << std::endl;
    }
  }
  return num_boxes;
}



// TODO: MOVE THESE

/// Compute the mean of a vector.
template <typename T>
double mean(std::vector<T> const& values) {

  double count = 0;
  double result = 0;
  for (size_t i = 0; i < values.size(); ++i) {
    if (is_valid(values[i])) {
      result += values[i];
      count  += 1.0;
    }
  }
  if (count < 0) // TODO: Exception?
    return 0;

  return result / count;
}

/// Compute the standard deviation of a vector.
template <typename T>
double standard_deviation(std::vector<T> const& values, double mean_value) {

  double result = 0, temp = 0, count = 0;
  for (size_t i = 0; i < values.size(); ++i) {
    if (is_valid(values[i])) {
      temp = values[i] - mean_value;
      result += temp * temp;
      count += 1.0;
    }
  }

  return sqrt(result / count);
}


// TODO: Replace with a simpler multi-threaded processing method to get the means!
template <class ImageT>
class ImageTileMeansView : public ImageViewBase<ImageTileMeansView<ImageT> > {

public: // Definitions

  // This only controls junk written to disk so use a small type.
  typedef uint8      pixel_type;
  typedef pixel_type result_type;

private: // Variables

  ImageT const& m_input_image;
  mutable ImageView<SarTypeM> m_tile_means, m_tile_stddevs;
  int m_tile_size;

public: // Functions

  // Constructor
  ImageTileMeansView( ImageT  const& input_image, int num_boxes_x, int num_boxes_y, int tile_size)
                  : m_input_image(input_image), m_tile_size(tile_size){
    // Init the true output objects
    m_tile_means.set_size(num_boxes_x,   num_boxes_y);
    m_tile_stddevs.set_size(num_boxes_x, num_boxes_y);
  }

  inline int32 cols  () const { return m_input_image.cols(); }
  inline int32 rows  () const { return m_input_image.rows(); }
  inline int32 planes() const { return 1; }

  inline result_type operator()( int32 i, int32 j, int32 p=0 ) const
  { return 0; } // NOT IMPLEMENTED!

  // Accessors to grab the results
  ImageView<SarTypeM> const& tile_means  () const {return m_tile_means;}
  ImageView<SarTypeM> const& tile_stddevs() const {return m_tile_stddevs;}
 
  typedef ProceduralPixelAccessor<ImageTileMeansView<ImageT> > pixel_accessor;
  inline pixel_accessor origin() const { return pixel_accessor( *this, 0, 0 ); }


  typedef CropView<ImageView<result_type> > prerasterize_type;
  inline prerasterize_type prerasterize( BBox2i const& bbox ) const {

    // Figure out which tile this is
    int this_col = bbox.min()[0] / m_tile_size;
    int this_row = bbox.min()[1] / m_tile_size;

    // Skip processing for this tile if it falls out of bounds (can happen on borders)
    if ( (this_col < m_tile_means.cols()) && (this_row < m_tile_means.rows())) {
      // Compute the four sub-ROIs
      const int NUM_SUB_ROIS = 4;
      int hw = bbox.width() /2;
      int hh = bbox.height()/2;
      std::vector<BBox2i> sub_rois(NUM_SUB_ROIS);
      std::vector<PixelMask<double> > means   (NUM_SUB_ROIS);
      /*
      sub_rois[0] = BBox2i(bbox.min()[0],    bbox.min()[1],    hw, hh); // Top left
      sub_rois[1] = BBox2i(bbox.min()[0]+hw, bbox.min()[1],    hw, hh); // Top right
      sub_rois[2] = BBox2i(bbox.min()[0]+hw, bbox.min()[1]+hh, hw, hh); // Bottom right
      sub_rois[3] = BBox2i(bbox.min()[0],    bbox.min()[1]+hh, hw, hh); // Bottom left
      */
      sub_rois[0] = BBox2i(0,  0,  hw, hh); // Top left
      sub_rois[1] = BBox2i(hw, 0,  hw, hh); // Top right
      sub_rois[2] = BBox2i(hw, hh, hw, hh); // Bottom right
      sub_rois[3] = BBox2i(0,  hh, hw, hh); // Bottom left
      
      //printf("Means for tile: %d, %d\n", this_col, this_row);
      
      ImageView<SarTypeM> section = crop(m_input_image, bbox);
      
      // Compute the mean in each of the four sub-rois
      for (int i=0; i<NUM_SUB_ROIS; ++i) {
        try {
          means[i] = mean_channel_value(crop(section, sub_rois[i]));
        } catch(ArgumentErr) {
          invalidate(means[i]); // No valid pixels in this region.
        }
        //std::cout << "Mean for sub-roi " << sub_rois[i] << " = " << means[i] << std::endl;
      }

      // Compute the standard deviation of the means
      // - Currently both are set to zero if all pixels are invalid.
      double mean_of_means   = mean(means);
      bool is_valid = (mean_of_means > 0);
      double stddev_of_means = 0;
      if (is_valid) {
        stddev_of_means = standard_deviation(means, mean_of_means);
        //printf("Mean, stddev for tile (%d, %d) = %lf, %lf\n", this_col, this_row, mean_of_means, stddev_of_means);
      }

      // Assign the REAL outputs
      if (is_valid) {
        m_tile_means  (this_col, this_row) = SarTypeM(mean_of_means);
        m_tile_stddevs(this_col, this_row) = SarTypeM(stddev_of_means);
      } else {
        invalidate(m_tile_means  (this_col, this_row));
        invalidate(m_tile_stddevs(this_col, this_row));
      }
    } // End mean/stddev calculation

    // Set up the output image tile - This is junk we don't care about!
    ImageView<result_type> tile(bbox.width(), bbox.height());

    // Return the tile we created with fake borders to make it look the size of the entire output image
    return prerasterize_type(tile,
                             -bbox.min().x(), -bbox.min().y(),
                             cols(), rows() );

  } // End prerasterize function

 template <class DestT>
 inline void rasterize( DestT const& dest, BBox2i const& bbox ) const {
   vw::rasterize( prerasterize(bbox), dest, bbox );
 }

}; // End class ImageTileMeansView



template <class ImageT>
void generate_tile_means(ImageT input_image, int tile_size, int num_boxes_x, int num_boxes_y,
                         ImageView<SarTypeM> &tile_means, ImageView<SarTypeM> &tile_stddevs) {

  // These tiles must be written at this exact size to get the correct results!
  //Vector2i block_size(tile_size, tile_size);
  
  ImageTileMeansView<ImageT> tile_mean_generator(input_image, num_boxes_x, num_boxes_y, tile_size);
  
  std::string dummy_path = "dummy.tif";
  
  // TODO: Using this View object is a hack to get multi-threaded processing
  //       and it should be replaced with a simpler implementation!!!
  
  //DiskImageResourceGDAL mean_resource(dummy_path, input_image.format(), block_size);
  
  // Write dummy image to force multi-threaded processing
  //block_write_image(mean_resource, tile_mean_generator, 
  //                  TerminalProgressCallback("vw", "\t--> Computing tile means:"));

  // TODO: Load input options!
  cartography::GdalWriteOptions write_options;
  write_options.raster_tile_size = Vector2i(tile_size, tile_size);
  block_write_gdal_image(dummy_path, tile_mean_generator, write_options,
                         TerminalProgressCallback("vw", "\t--> Computing tile means:"));


  // Grab results from the view object
  tile_means   = tile_mean_generator.tile_means();
  tile_stddevs = tile_mean_generator.tile_stddevs();

} // End function generate_tile_means





/* TODO: Convert this for local use!

// Erode blobs from given image by iterating through tiles, biasing
// each tile by a factor of blob size, removing blobs in the tile,
// then shrinking the tile back. The bias is necessary to help avoid
// fragmenting (and then unnecessarily removing) blobs.
template <class ImageT>
class PerTileErode: public ImageViewBase<PerTileErode<ImageT> >{
  ImageT m_img;
public:
  PerTileErode( ImageViewBase<ImageT>   const& img):
    m_img(img.impl()){}

  // Image View interface
  typedef typename ImageT::pixel_type pixel_type;
  typedef pixel_type                  result_type;
  typedef ProceduralPixelAccessor<PerTileErode> pixel_accessor;

  inline int32 cols  () const { return m_img.cols(); }
  inline int32 rows  () const { return m_img.rows(); }
  inline int32 planes() const { return 1; }

  inline pixel_accessor origin() const { return pixel_accessor( *this, 0, 0 ); }

  inline pixel_type operator()( double i, double j, int32 p = 0 ) const {
    vw_throw(NoImplErr() << "PerTileErode::operator()(...) is not implemented");
    return pixel_type();
  }

  typedef CropView<ImageView<pixel_type> > prerasterize_type;
  inline prerasterize_type prerasterize(BBox2i const& bbox) const {

    int area = stereo_settings().erode_max_size;

    // We look a beyond the current tile, to avoid cutting blobs
    // if possible. Skinny blobs will be cut though.
    int bias = 2*int(ceil(sqrt(double(area))));

    BBox2i bbox2 = bbox;
    bbox2.expand(bias);
    bbox2.crop(bounding_box(m_img));
    ImageView<pixel_type> tile_img = crop(m_img, bbox2);

    int tile_size = max(bbox2.width(), bbox2.height()); // don't subsplit
    BlobIndexThreaded smallBlobIndex(tile_img, area, tile_size);
    ImageView<pixel_type> clean_tile_img = applyErodeView(tile_img,
                                                          smallBlobIndex);
    return prerasterize_type(clean_tile_img,
                             -bbox2.min().x(), -bbox2.min().y(),
                             cols(), rows() );
  }

  template <class DestT>
  inline void rasterize(DestT const& dest, BBox2i bbox) const {
    vw::rasterize(prerasterize(bbox), dest, bbox);
  }
};

template <class ImageT>
PerTileErode<ImageT>
per_tile_erode( ImageViewBase<ImageT> const& img) {
  typedef PerTileErode<ImageT> return_type;
  return return_type( img.impl() );
}
*/




/** Main function of algorithm from:
      Martinis, Sandro, Jens Kersten, and Andre Twele. 
      "A fully automated TerraSAR-X based flood service." 
      ISPRS Journal of Photogrammetry and Remote Sensing 104 (2015): 203-212.
*/
void sar_martinis(std::string const& input_image_path, 
                  cartography::GdalWriteOptions const& write_options) {

  // TODO: How to specify the ROI?
  
  BBox2i roi = bounding_box(DiskImageView<SarType>(input_image_path));
 
  // Load the georeference from the input image
  cartography::GeoReference georef;
  bool have_georef = cartography::read_georeference(georef, input_image_path);
  if (have_georef) {
    georef = crop(georef, roi); // Account for the input ROI
    std::cout << "Read georeference: " << georef << std::endl;
  }
  else
    std::cout << "FAILED to read georeference!\n";
  
  // Read nodata value
  const double DEFAULT_NODATA = 0; // TODO!!!!
  boost::scoped_ptr<SrcImageResource> rsrc(DiskImageResource::open(input_image_path));
  double nodata_value = DEFAULT_NODATA;
  if (rsrc->has_nodata_read()) {
    nodata_value = rsrc->nodata_read();
    std::cout << "Read nodata: " << nodata_value << std::endl;
  }
  else
    std::cout << "Failed to read nodata value, using default value of 0.\n";
  const bool have_nodata = true;
  
  //std::cout << "MinMax...\n";
  
  // Compute the min and max values of the image
  // TODO: How to handle the preprocessing?  Inputs are 16 bit but almost everything is < 1000.
  //       actually computing the real max seems worthless.
  SarType global_min=0, global_max=1000;
  //min_max_channel_values(crop(DiskImageView<SarType>(input_image_path), roi), global_min, global_max);

  std::cout << "Preprocessing...\n";

  // Apply any needed preprocessing to the image TODO
  ImageViewRef <SarTypeM> preprocessed_image;
  preprocess_sentinel1_image(DiskImageView<SarType>(input_image_path), nodata_value, roi, 
                             global_min, global_max, write_options, preprocessed_image);

  int tile_size = 512; // TODO
  
  // Generate vector of BBoxes for each tile in the input image (S+)
  std::vector<std::vector<BBox2i> > large_tile_boxes;
  divide_roi(bounding_box(preprocessed_image), tile_size, large_tile_boxes, false);
  int num_boxes_y = large_tile_boxes.size();
  int num_boxes_x = large_tile_boxes[0].size();
  //printf("Num boxes: %d, %d\n", num_boxes_x, num_boxes_y);
  std::cout << "Computing tile means...\n";

  // For each tile compute the mean value and the standard deviation of the four sub-tiles.
  ImageView<SarTypeM> tile_means, tile_stddevs; // These are much smaller than the input image
  generate_tile_means(preprocessed_image, tile_size, num_boxes_x, num_boxes_y,
                      tile_means, tile_stddevs);

  std::cout << "Writing DEBUG images...\n";
  block_write_gdal_image("tile_means.tif",   tile_means,   write_options);
  block_write_gdal_image("tile_stddevs.tif", tile_stddevs, write_options);


  std::cout << "Computing global tile statistics...\n";

  // Compute the global mean
  double global_mean = mean_channel_value(tile_means);

  std::cout << "Computing global mean = " << global_mean << "\n";

  // Compute 95% quantile standard deviation
  double stddev_min, stddev_max;
  find_image_min_max(tile_stddevs, stddev_min, stddev_max);
  int num_bins = 255;
  std::vector<double> hist;
  histogram(tile_stddevs, num_bins, stddev_min, stddev_max, hist);
  const double TILE_STDDEV_PERCENTILE_CUTOFF = 0.95;
  int bin = get_histogram_percentile(hist, TILE_STDDEV_PERCENTILE_CUTOFF);
  double bin_width      = (stddev_max - stddev_min)/static_cast<double>(num_bins);
  double std_dev_cutoff = stddev_min + bin_width*bin;
  std::cout << "std_dev_cutoff " << std_dev_cutoff << "\n";

  // Select the tiles with the highest STD values (N')
  ImageView<uint8> kept_tile_display;
  kept_tile_display.set_size(tile_means.cols(), tile_means.rows());
  std::vector<Vector2i> n_prime_tiles;
  std::vector<double  > n_prime_std_dev;
  for (int r=0; r<tile_stddevs.rows(); ++r) {
    for (int c=0; c<tile_stddevs.cols(); ++c) {
      // The tile must have a high stddev and also be below the global mean
      //  since water tends to be darker than land.
      if ( (tile_stddevs(c,r) > std_dev_cutoff) &&
           (tile_means  (c,r) < global_mean   )   ) {
        n_prime_tiles.push_back(Vector2i(c,r));
        n_prime_std_dev.push_back(tile_stddevs(c,r));
        std::cout << "Keeping tile " << Vector2i(c,r) << std::endl;
        kept_tile_display(c,r) = 255;
      }
    }
  }
  block_write_gdal_image("kept_tiles.tif", kept_tile_display, write_options); // DEBUG
  
  int num_tiles_kept = static_cast<int>(n_prime_tiles.size());
  std::cout << "Selected " << num_tiles_kept << " tiles.\n";
  
  // TODO: When this happens, repeat the earlier steps with the tile size cut in half.
  if (n_prime_tiles.empty())
    vw_throw(LogicErr() << "No tiles left after std_dev filtering!");

  // TODO: Cap the number of selected tiles!!!!!!!!!


  // For each selected tile, find optimal threshold using Kittler-Illingworth method.
  // - TODO: Does any rescaling need to be applied to these thresholds values?
  //         splitValDb = rescaleNumber(splitVal, PROC_global_min, PROC_global_max, minVal, maxVal)        
  // - TODO: Compute this in parallel!
  std::vector<double> optimal_tile_thresholds(num_tiles_kept);
  double dmin = static_cast<double>(global_min);
  double dmax = static_cast<double>(global_max);
  double threshold_mean = 0;
  for (int i=0; i<num_tiles_kept; ++i) {
  
    // Compute tile histogram     
    Vector2i tile_index = n_prime_tiles[i];
    BBox2i roi = large_tile_boxes[tile_index[1]][tile_index[0]]; // The ROIs are stored row-first.
    std::vector<double> hist;
    histogram(crop(preprocessed_image, roi), num_bins, dmin, dmax, hist);
    
    // Compute optimal split
    optimal_tile_thresholds[i] = splitHistogramKittlerIllingworth(hist, num_bins, global_min, global_max);
    std::cout << "For ROI " << roi << ", computed threshold " << optimal_tile_thresholds[i] << std::endl;
    threshold_mean += optimal_tile_thresholds[i];
  }
  threshold_mean /= static_cast<double>(num_tiles_kept);

  // TODO: Check statistics of the computed tile thresholds
  // TODO: Some method to discard outliers

  // If the standard deviation of the local thresholds in DB are greater than this,
  //  the result is probably bad (number from the paper)
  double MAX_STD_DB = 5.0;

  // The maximum allowed value, from the paper.
  double MAX_THRESHOLD_DB = 10.0;

  double threshold_stddev = standard_deviation(optimal_tile_thresholds, threshold_mean);
  
  // TODO // Recompute these values in the original DB units.
  //tileThresholdsDb = [rescaleNumber(x, PROC_global_min, PROC_global_max, minVal, maxVal) for x in tileThresholds]
  //threshMeanDb = numpy.mean(tileThresholdsDb)
  //threshStdDb  = numpy.std(tileThresholdsDb, ddof=1)

  std::cout << "Mean of tile thresholds: " << threshold_mean   << std::endl;
  std::cout << "STD  of tile thresholds: " << threshold_stddev << std::endl;
  //print 'Mean of tile thresholds (DB): ' + str(threshMeanDb)
  //print 'STD  of tile thresholds (DB): ' + str(threshStdDb)

  // TODO: Verify that the computed threshold looks reasonable!

  //// TODO: Use an alternate method of computing the threshold like they do in the paper!
  //if threshStdDb > MAX_STD_DB:
  //    raise Exception('Failed to compute a good threshold! STD = ' + str(threshStdDb))
  // TODO
  //if threshMeanDb > MAX_THRESHOLD_DB:
  //    threshMean = rescaleNumber(MAX_THRESHOLD_DB, minVal, maxVal, PROC_global_min, PROC_global_max)
  //    print 'Saturating the computed threshold at 10 DB!'
  
  
  // This will mask the water pixels
  ImageViewRef<SarTypeM> raw_water = create_mask_less_or_equal(preprocessed_image, threshold_mean);

  // DEBUG: Apply the initial threshold to the image and save it to disk!
  block_write_gdal_image("initial_water_detect.tif",
                         apply_mask(copy_mask(constant_view(uint8(0), preprocessed_image), 
                                              raw_water),
                                    uint8(255)),
                         have_georef, georef,
                         have_nodata, nodata_value,
                         write_options,
                         TerminalProgressCallback("vw", "\t--> Applying initial threshold:"));

  // Get information needed for fuzzy logic results filtering
return;
  std::cout << "Computing mean of flooded regions...\n";

  // Compute mean value of pixels under initial water threshold
  // - TODO: Is this in parallel?
  double mean_raw_water_value = mean_channel_value(raw_water);

  std::cout << "Mean value of flooded regions = " << mean_raw_water_value << std::endl;

  std::cout << "TODO: Implement DEM, blobs, fuzzy logic, expansion!\n";

  // No point going through the fuzzy logic stuff when most of the inputs are missing!

/*
  // TODO: Form image blobs. --> This is complicated for large images so come back later!
LOCAL HOPEFULLY - Local blobs plus a blob consolidation step???
  // TODO: Compute the number of pixels in each blob, up to the maximum we care about (1000m*m)
  maxBlobSize = 1000/BASE_RES
  blobSizes   = rawWater.mask(rawWater).connectedPixelCount(maxBlobSize)
*/

/*
LOCAL
  // TODO: Get a DEM to go with the image!!!  --> Come back to this later!!!!!
  // TODO: Fill in DEM holes with zero elevation
  // TODO: Compute slope at each DEM pixel!
  slopeImage = ee.Terrain.slope(dem)

GLOBAL
  // TODO: Compute mean and std of the elevations of the pixels we marked as water
  waterHeights = dem.mask(rawWater)
  meanWaterHeight = waterHeights.reduceRegion(ee.Reducer.mean(),   domain.bounds, scale=BASE_RES).getInfo()['elevation']
  stdWaterHeight  = waterHeights.reduceRegion(ee.Reducer.stdDev(), domain.bounds, scale=BASE_RES).getInfo()['elevation']
*/


  /*
  TODO: Come back to these later!

  // For each pixel: 
    
  // TODO: Compute fuzzy classifications on four categories:
  
  // SAR
  sarFuzz = fuzzMemZ(radarImage, mean_raw_water_value,  threshold_mean)
     
  // Elevation
  // TODO: The max value seems a little strange.
  heightFuzz = fuzzMemZ(dem, meanWaterHeight, meanWaterHeight + stdWaterHeight*(stdWaterHeight + 3.5))
  
  // Slope
  slopeFuzz = fuzzMemZ(slopeImage, 0, 15) // Values in degrees
  
  
  // Body size
  minBlobSize = 250/BASE_RES
  blobFuzz = fuzzMemS(blobSizes, minBlobSize, maxBlobSize).mask(blobSizes)


  // Combine fuzzy classification into a single fuzzy classification
  zeroes    = sarFuzz.Not().Or(heightFuzz.Not().Or(slopeFuzz.Not().Or(blobFuzz.Not())))
  meanFuzz  = sarFuzz.add(heightFuzz).add(slopeFuzz).add(blobFuzz).divide(ee.Image(4.0))
  finalFuzz = meanFuzz.And(zeroes.Not())


  double final_flood_threshold = 0.6;
  double water_grow_threshold  = 0.45;
  
  // Apply fixed threshold to get updated flooded pixels
  defuzz = finalFuzz.gt(ee.Image(final_flood_threshold))
  
  // TODO: Region growing on water-classified pixels using a lower threshold.
  
  // Expand water classification outwards using a lower fuzzy threshold.
  // - This step is a little tough for EE so we approximate using a dilation step.
  dilatedWater = defuzz.focal_max(radius=1000, units='meters')
  finalWater   = dilatedWater.And(finalFuzz.gt(ee.Image(water_grow_threshold)))
  
*/
} // End function sar_martinis













}} // end namespace vw/radar
#endif

