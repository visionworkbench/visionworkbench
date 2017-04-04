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

/**
  Program to exercise VW's water detection tools.
*/

#include <queue>
#include <vw/Image/ImageIO.h>
#include <vw/Math/Statistics.h>
#include <vw/tools/modis_utilities.h>
#include <vw/tools/modis_water_detection.h>
//#include <vw/tools/radar.h>

#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

using namespace vw;


// TODO: MOVE
/// Applies a median filter to a vector with a window size
/// - TODO: Better implementation!
template <typename T>
void window_median_filter(std::vector<T> const& v_in, std::vector<T> &v_out, const int width) {

  // Setup
  const int length = static_cast<int>(v_in.size());
  const int half_width = width / 2;
  v_out = v_in;

  // Handle small input sizes
  if ((width < 3) || (length < width))
    return;
  
//  std::cout << "length = " << length << std::endl; 
//  std::cout << "half_width = " << half_width << std::endl; 
  
  // Just copy the first few elements
  for (int i=0; i<half_width; ++i)
    v_out[i] = v_in[i];

  // Init a queue of neigboring values
  std::list<T> locals;    
  for (int i=0; i<width-1; ++i)    
    locals.push_back(v_in[i]);

  for (int i=half_width; i<length-half_width; ++i) {
    // Load next value into queue
    locals.push_back(v_in[i+half_width]);
  
    // Compute the median of the neigbors
    v_out[i] = math::median(locals);
  
    // Drop oldest value
    locals.pop_front();
  }
  // Just copy the last few elements
  for (int i=length-half_width; i<length; ++i) {
    v_out[i] = v_in[i];
  }
} // End median_filter


/// Returns the percentage of values equal to a given value
template <typename T>
double percent_equal(T const& values, double equal_to) {

  double count = 0;
  typename T::const_iterator iter;
  for (iter = values.begin(); iter != values.end(); ++iter) {
    if (static_cast<double>(*iter) == equal_to)
      count += 1.0;
  }

  return count / static_cast<double>(values.size());
}

/// Class to locate the best position (if any) of a jump from low pixel values to higher pixel values.
class JumpFinder {
public:

  JumpFinder(int lead_size) 
    : m_lead_size(static_cast<size_t>(lead_size)), 
      m_best_position(0), m_best_score(15), m_zero_count(0), m_count(0), m_trail_sum(0) {
  }

  /// Add a new value and return true if a jump was detected
  bool add_value(double value) {
  
    const int    BUFFER_WIDTH = 1; // Keep this many numbers out of the statistics
    const double CLOSE_ZERO   = 8; // Treat values this small as basically zero
  
    // Add new value to the leading list
    m_lead_list.push_front(value);
    
    if (m_lead_list.size() > m_lead_size) {
      // Move back of lead list to the trailing list
      m_buffer_list.push_front(m_lead_list.back());
      m_lead_list.pop_back();
      
      if (m_buffer_list.size() > BUFFER_WIDTH) {
        double trail_val = m_buffer_list.back();
        if (trail_val <= CLOSE_ZERO)
          ++m_zero_count;
        ++m_count;
        m_trail_list.push_front(trail_val);
        m_buffer_list.pop_back();
        m_trail_sum += trail_val;
      }
    }
    
    if (m_trail_list.empty())
      return false;

    // Compute statistics
    double mean_trail=0, mean_lead=0, dev_trail=0;
    double percent_zero = m_zero_count / m_count;
    mean_lead = math::mean(m_lead_list);
    if (percent_zero < 1.0) { // Save time on big blank regions
      mean_trail = m_trail_sum / m_count;
      dev_trail  = math::standard_deviation(m_trail_list, mean_trail);
    }
    
    const double MAX_DEV_TRAIL    = 5;    // After STD_DEV of trail portion is higher, not a good jump
    const double MIN_PERCENT_ZERO = 0.30; // Must be this percent or higher zeroes in trail for a jump
    const double DEV_TRAIL_BREAK  = 20;   // Quit if the trail deviation exceeds this threshold.

    // Evaluate the value of an edge being at this location
    double score = (mean_lead - mean_trail);
    if ((dev_trail > MAX_DEV_TRAIL) ||
        (percent_zero < MIN_PERCENT_ZERO))
      score = 0;
    //std::cout << "Lead  stats: " << mean_lead  << ", " << dev_lead  << std::endl;
    //std::cout << "Trail stats: " << mean_trail << ", " << dev_trail 
    //          << ", size = " << m_trail_list.size() << ", percent = " << percent_zero<< std::endl;
    if (score > m_best_score) { // Record the best edge score to date
      m_best_score    = score;
      m_best_position = get_jump_index();
    }
    
    // If one of these things is true, very unlikely to see the edge beyond this point.
    return ((percent_zero < MIN_PERCENT_ZERO) || (dev_trail > DEV_TRAIL_BREAK)); 
  }
  
  /// Return the index at which the jump occurs.
  int get_jump_index() const {
    return m_trail_list.size() + 1; // Advance into the buffer pixel
  }
  
  int get_best_jump() const {
    //printf("Best score = %lf and position %d\n", m_best_score, m_best_position);
    return m_best_position;
  }

private:
  size_t m_lead_size;
  std::list<double> m_lead_list;
  std::list<double> m_trail_list;
  std::list<double> m_buffer_list;
  int    m_best_position;
  double m_best_score;
  double m_zero_count;
  double m_count;
  double m_trail_sum;

}; // End class JumpFinder



template <class ImageT>
class Sentinel1CleanBordersView : public ImageViewBase<Sentinel1CleanBordersView<ImageT> > {

public: // Definitions

  typedef uint16     pixel_type;
  typedef pixel_type result_type;

private: // Variables

  ImageT const& m_input_image;

public: // Functions

  // Constructor
  Sentinel1CleanBordersView( ImageT  const& input_image)
                  : m_input_image(input_image){}

  inline int32 cols  () const { return m_input_image.cols();  }
  inline int32 rows  () const { return m_input_image.rows(); }
  inline int32 planes() const { return 1; }

  inline result_type operator()( int32 i, int32 j, int32 p=0 ) const
  { return 0; } // NOT IMPLEMENTED!
 
  typedef ProceduralPixelAccessor<Sentinel1CleanBordersView<ImageT> > pixel_accessor;
  inline pixel_accessor origin() const { return pixel_accessor( *this, 0, 0 ); }

  void check_bbox_edges(BBox2i const& bbox, bool &left, bool &right, bool &top, bool &bottom) const {
    left   = (bbox.min().x() == 0);
    top    = (bbox.min().y() == 0);
    right  = (bbox.max().x() == cols()); // TODO: CHECK
    bottom = (bbox.max().y() == rows());
  }

  typedef CropView<ImageView<result_type> > prerasterize_type;
  inline prerasterize_type prerasterize( BBox2i const& bbox ) const {

    ImageView<result_type> output_tile = crop(m_input_image, bbox);

    // Check if this is a border tile.  If not, don't do any processing.
    
    // Determine if this tile sits on the border of the tiles.
    bool left_edge, right_edge, top_edge, bottom_edge;
    check_bbox_edges(bbox, left_edge, right_edge, top_edge, bottom_edge);

    const unsigned short NODATA_VALUE        = 0; // Output nodata flag
    const int            MEDIAN_FILTER_WIDTH = 9; // Width of median filter applied to edge detect results
    const int            JUMP_WIDTH          = 8; // Width of lead region in the JumpFinder class
    
    // For each line coming in from the border, choose the best break point.
    int num_rows = bbox.height();
    int num_cols = bbox.width();
    
    if (left_edge) {
      //std::cout << "LEFT bbox = " << bbox << std::endl;
      std::vector<int> breaks(num_rows), breaks_filtered;
      for (int r=0; r<num_rows; ++r) { // For each row in this tile
        breaks[r] = 0; // Default if no jump found
        JumpFinder edge_detector(JUMP_WIDTH);
        // Keep iterating to the right until we detect the start of valid pixels
        for (int c=0; c<num_cols; ++c) { // Move along the row     
          int value = output_tile(c,r);
          if (edge_detector.add_value(value))
            break;
        } // End col loop
        breaks[r] = edge_detector.get_best_jump();
      } // End row loop

      // Smooth out the borders
      window_median_filter(breaks, breaks_filtered, MEDIAN_FILTER_WIDTH);

      // Update the output image
      for (int r=0; r<num_rows; ++r) { // For each row in this tile
        for (int c=0; c<=breaks_filtered[r]; ++c) // Move along the row
          output_tile(c,r) = NODATA_VALUE;
      }
    } // End left edge case

    if (right_edge) {
      //std::cout << "RIGHT bbox = " << bbox << std::endl;
      std::vector<int> breaks(num_rows), breaks_filtered;
      for (int r=0; r<num_rows; ++r) { // For each row in this tile
        //std::cout << "Row = " << r+bbox.min().y() << std::endl;
        breaks[r] = num_cols-1; // Default if no jump found
        JumpFinder edge_detector(JUMP_WIDTH);
        // Keep iterating to the right until we detect the start of valid pixels
        for (int c=num_cols-1; c>0; --c) { // Move along the row     
          int value = output_tile(c,r);
          if (edge_detector.add_value(value))
            break;
        } // End col loop
        breaks[r] = num_cols-1 - edge_detector.get_best_jump();
      } // End row loop

      // Smooth out the borders
      window_median_filter(breaks, breaks_filtered, MEDIAN_FILTER_WIDTH);

      // Update the output image
      for (int r=0; r<num_rows; ++r) { // For each row in this tile
        for (int c=num_cols-1; c>=breaks_filtered[r]; --c) // Move along the row
          output_tile(c,r) = NODATA_VALUE;
      }
    } // End right edge case

    if (top_edge) {
      //std::cout << "TOP bbox = " << bbox << std::endl;
      std::vector<int> breaks(num_cols), breaks_filtered;
      for (int c=0; c<num_cols; ++c) { // For each column in this tile
        breaks[c] = 0; // Default if no jump found
        JumpFinder edge_detector(JUMP_WIDTH);
        // Keep iterating to the right until we detect the start of valid pixels
        for (int r=0; r<num_rows; ++r) { // Move along the column
          int value = output_tile(c,r);
          if (edge_detector.add_value(value))
            break;
        } // End row loop
        breaks[c] = edge_detector.get_best_jump();
        //std::cout << "breaks[c] = " << breaks[c] << std::endl;
      } // End col loop

      // Smooth out the borders
      window_median_filter(breaks, breaks_filtered, MEDIAN_FILTER_WIDTH);

      // Update the output image
      for (int c=0; c<num_cols; ++c) { // For each column in this tile
        //std::cout << "breaks_filtered[c] = " << breaks_filtered[c] << std::endl;
        for (int r=0; r<breaks_filtered[c]; ++r) // Move along the column
          output_tile(c,r) = NODATA_VALUE;
      }
    } // End top edge case

    if (bottom_edge) {
      //std::cout << "BOTTOM bbox = " << bbox << std::endl;
      std::vector<int> breaks(num_cols), breaks_filtered;
      for (int c=0; c<num_cols; ++c) { // For each column in this tile
        breaks[c] = num_rows-1; // Default if no jump found
        JumpFinder edge_detector(JUMP_WIDTH);
        // Keep iterating to the right until we detect the start of valid pixels
        for (int r=num_rows-1; r>0; --r) { // Move along the column
          int value = output_tile(c,r);
          if (edge_detector.add_value(value))
            break;
        } // End row loop
        breaks[c] = num_rows-1 - edge_detector.get_best_jump();
      } // End col loop

      // Smooth out the borders
      window_median_filter(breaks, breaks_filtered, MEDIAN_FILTER_WIDTH);

      // Update the output image
      for (int c=0; c<num_cols; ++c) { // For each column in this tile
        for (int r=num_rows-1; r>breaks_filtered[c]; --r) // Move along the column
          output_tile(c,r) = NODATA_VALUE;
      }
    } // End bottom edge case

    // Return the tile we created with fake borders to make it look the size of the entire output image
    return prerasterize_type(output_tile, -bbox.min().x(), -bbox.min().y(), cols(), rows() );

  } // End prerasterize function

 template <class DestT>
 inline void rasterize( DestT const& dest, BBox2i const& bbox ) const {
   vw::rasterize( prerasterize(bbox), dest, bbox );
 }
}; // End class Sentinel1CleanBordersView

/// Helper function
template <class T>
Sentinel1CleanBordersView<T> clean_sentinel1_borders(T const& input) {
  return Sentinel1CleanBordersView<T>(input);
}

// TODO: MOVE!
/// Compute an optimal tile size close to the input tile size
int compute_new_tile_size(int input_size, int image_size) {

  // Define the safe range of percent size of the last tile used
  double MIN_TILE_PERCENTAGE = 0.75;
  double MAX_TILE_PERCENTAGE = 0.99;

  // The current tile usage
  double num_tiles_used       = (double)image_size / (double)input_size;
  double last_tile_percentage = num_tiles_used - floor(num_tiles_used);
  
  // If the input is in the acceptable range, keep it!
  if ((last_tile_percentage >= MIN_TILE_PERCENTAGE) && (last_tile_percentage <= MAX_TILE_PERCENTAGE))
    return input_size;

  printf("Adjusting input tile size:\n");
  printf("num_tiles = %lf, input tile size: %d\n", num_tiles_used, input_size);

  // Otherwise the tile size should be adjusted.
  
  // Decide how many tiles to shoot for
  double target_num_tiles = floor(num_tiles_used);
  if (last_tile_percentage > 0.75)
    target_num_tiles = ceil(num_tiles_used);
  
  // Compute new tile size, rounding up to ensure we don't end up with a tile sliver.
  int tile_size = ceil((double)image_size / target_num_tiles);

  double new_num_tiles = (double)image_size / (double)tile_size;
  printf("new_num_tiles = %lf, Computed new tile size: %d\n", new_num_tiles, tile_size);

  // GDAL block write sizes must be a multiple to 16 so if the input value is
  //  not a multiple of 16 increase it until it is.
  const int TILE_MULTIPLE = 16;
  if (tile_size % TILE_MULTIPLE != 0)
    tile_size = ((tile_size / TILE_MULTIPLE) + 1) * TILE_MULTIPLE;
    
  new_num_tiles = (double)image_size / (double)tile_size;
  printf("new_num_tiles = %lf, Computed new tile size: %d\n", new_num_tiles, tile_size);
    
  return tile_size;
}


int main(int argc, char **argv) {

  std::string input_file;
  std::string output_path;
  int num_threads = 0;
  int tile_size   = 2048;

  cartography::GdalWriteOptions write_options;

  po::options_description general_options("Clean up the borders of a Sentinel1 radar image.\n\nGeneral Options");
  general_options.add_options()
    ("output-path,o",    po::value<std::string>(&output_path), "The output file path")
    ("num-threads",      po::value<int>(&num_threads)->default_value(0), 
                         "Number of threads to use for writing")
    ("tile-size",        po::value<int>(&tile_size)->default_value(2048), 
                         "This is the farthest in that bad edges can be corrected.")
    ("help,h", "Display this help message");

  po::options_description hidden_options("");
  hidden_options.add_options()
    ("input-file", po::value<std::string>(&input_file));

  po::options_description options("Allowed Options");
  options.add(general_options).add(hidden_options);

  po::positional_options_description p;
  p.add("input-file", -1);

  std::ostringstream usage;
  usage << "Usage: clean_sentinel1_borders [options] <image file>" <<std::endl << std::endl;
  usage << general_options << std::endl;

  po::variables_map vm;
  try {
    po::store( po::command_line_parser( argc, argv ).options(options).positional(p).run(), vm );
    po::notify( vm );
  } catch (const po::error& e) {
    std::cout << "An error occured while parsing command line arguments.\n";
    std::cout << "\t" << e.what() << "\n\n";
    std::cout << usage.str();
    return 1;
  }

  if( vm.count("help") ) {
    std::cout << usage.str();
    return 0;
  }

  // Modify tile size to make sure good coverage of the borders.
  ImageFormat input_format = image_format(input_file);
  int tile_size_h = compute_new_tile_size(tile_size, input_format.cols);
  int tile_size_v = compute_new_tile_size(tile_size, input_format.rows); 


  // TODO: Clean up settings usage!
  write_options.raster_tile_size = Vector2i(tile_size_h, tile_size_v);
  vw_settings().set_default_tile_size(tile_size);
  if (num_threads > 0) {
    write_options.num_threads = num_threads;
    vw_settings().set_default_num_threads(write_options.num_threads);
  }

  if( vm.count("output-path") < 1 ) {
    std::cerr << "Error: must specify the output file!" << std::endl << std::endl;
    std::cout << usage.str();
    return 1;
  }

  // Generate the output image
  std::remove(output_path.c_str());
  block_write_gdal_image(output_path,
                         clean_sentinel1_borders(DiskImageView<uint16>(input_file)), 
                         0, write_options);
  
}
