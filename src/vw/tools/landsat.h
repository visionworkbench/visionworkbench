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


#ifndef __VW_LANDSAT_H__
#define __VW_LANDSAT_H__

#include <stdlib.h>
#include <algorithm>
#include <vw/Image/PixelMask.h>
#include <vw/Image/Manipulation.h>
#include <vw/Image/Transform.h>
#include <vw/FileIO/DiskImageView.h>
#include <vw/Cartography/GeoReferenceUtils.h>

/**
  Tools for processing LANDSAT data.
  Landsat images are split into multiple geotif files and are medium sized.
*/

namespace vw {

  // Used to save the landsat image as stored in memory
  template<> struct PixelFormatID<PixelMask<Vector<float, 7> > > { static const PixelFormatEnum value = VW_PIXEL_GENERIC_8_CHANNEL; };

namespace landsat {


//----------------------------------------------------
// Landsat types

typedef PixelMask<Vector<uint16, 7> > LandsatPixelType;
typedef ImageViewRef<LandsatPixelType> LandsatImage;

// These are the Landsat channels we use, available on both 7 and 8.
enum LANDSAT_CHANNEL_INDICES { BLUE  = 0, 
                               GREEN = 1, 
                               RED   = 2, 
                               NIR   = 3, 
                               SWIR1 = 4, 
                               TEMP  = 5, // = Temperature
                               SWIR2 = 6};
                               
//----------------------------------------------------

/// View for treating multiple single channel images on disk as one multi-channel image.
/// - T is the data type, N is the number of channels.
template <typename T, int N>
class SplitChannelFileView : public ImageViewBase<SplitChannelFileView<T, N> > {

public: // Definitions

  typedef Vector<T, N> pixel_type;
  typedef pixel_type   result_type;

private: // Variables

  std::vector<std::string> const m_input_files;
  int m_cols, m_rows;

public: // Functions

  // Constructor
  SplitChannelFileView( std::vector<std::string> const& input_files)
                  : m_input_files(input_files){
    // Verify that the correct number of input files were passed in.
    if (input_files.size() != N)
      vw_throw( ArgumentErr() << "SplitChannelFileView created with incorrect number of input files!\n");
      
    // Record the image size to save time later
    boost::scoped_ptr<SrcImageResource> resource1(DiskImageResource::open(input_files[0]));
    m_cols = resource1->format().cols;
    m_rows = resource1->format().rows;
    
    // Make sure all images are the same size
    for (size_t i=0; i<input_files.size(); ++i) {
      boost::scoped_ptr<SrcImageResource> resource2(DiskImageResource::open(input_files[i]));
      if ((m_cols != static_cast<int>(resource2->format().cols)) || 
          (m_rows != static_cast<int>(resource2->format().rows))   )
        vw_throw( ArgumentErr() << "SplitChannelFileView: Input files must all be the same size!\n");
    }
  }

  inline int32 cols  () const { return m_cols; }
  inline int32 rows  () const { return m_rows; }
  inline int32 planes() const { return 1; }

  inline result_type operator()( int32 i, int32 j, int32 p=0 ) const { 
    vw_throw( NoImplErr() << "SplitChannelFileView: () function not implemented!\n");
    return result_type();
  }

 
  typedef ProceduralPixelAccessor<SplitChannelFileView<T, N> > pixel_accessor;
  inline pixel_accessor origin() const { return pixel_accessor( *this, 0, 0 ); }

  typedef CropView<ImageView<result_type> > prerasterize_type;
  inline prerasterize_type prerasterize( BBox2i const& bbox ) const {
    // Init the tile
    ImageView<result_type> tile(bbox.width(), bbox.height());
    
    // Load info into the file one file at a time
    for (size_t channel=0; channel<m_input_files.size(); ++channel) {
      ImageView<T> temp = crop(DiskImageView<T>(m_input_files[channel]), bbox);
      select_channel(tile, channel) = temp;
    }

    // Return the tile we created with fake borders to make it look the size of the entire output image
    return prerasterize_type(tile,
                             -bbox.min().x(), -bbox.min().y(),
                             cols(), rows() );

  } // End prerasterize function

 template <class DestT>
 inline void rasterize( DestT const& dest, BBox2i const& bbox ) const {
   vw::rasterize( prerasterize(bbox), dest, bbox );
 }
}; // End class SplitChannelFileView


/*  This function is for the wrong file format??
/// Process MODIS input files into a single image object
void load_landsat_image(LandsatImage &image, std::vector<std::string> const& image_files,
                        int landsat_type, float &sun_elevation) {

  // TODO: Are all bands available?
  const int NUM_BANDS_OF_INTEREST = 7;
  const std::string LS5_BAND_LOCATIONS[NUM_BANDS_OF_INTEREST] = {"0", "1", "2", "3", "4", "5", "6"};
  const std::string LS7_BAND_LOCATIONS[NUM_BANDS_OF_INTEREST] = {"0", "1", "2", "3", "4", "5", "7"};
  const std::string LS8_BAND_LOCATIONS[NUM_BANDS_OF_INTEREST] = {"1", "2", "3", "4", "5", "cloud", "6"};

  const std::string BAND_PREFIX = "sr_band";
  
  sun_elevation = 45; // The default sun elevation angle

  // TODO: Verify always have the same nodata value!          
  const double DEFAULT_NODATA = -9999;
  double nodata_value = DEFAULT_NODATA;
  //const bool have_nodata = true;
  
  // Check that all the required bands are present
  for (int chan=0; chan<NUM_BANDS_OF_INTEREST; ++chan) {
  
    // For the given output channel, determine the input channel containing it
    std::string text;
    switch(landsat_type) {
      case 8:  text = LS8_BAND_LOCATIONS[chan]; break;
      case 7:  text = LS7_BAND_LOCATIONS[chan]; break;
      default: text = LS5_BAND_LOCATIONS[chan]; break;
    };
    if (text.size() == 1)
      text = "band" + text;
    text = "_sr_" + text;

    //printf("Looking for data %d in channel %d\n", chan, input_channel);
    std::cout << "Looking for text: " << text << std::endl;
    
    // Look for the string in the input file names
    bool found = false;
    for (size_t f=0; f<image_files.size(); ++f) {
      if (image_files[f].find(text) != std::string::npos) {
        found = true;
        DiskImageView<int16> input_image(image_files[f]);
        // Init the output image
        if (image.rows() == 0) {
          image.set_size(input_image.cols(), input_image.rows());


          // Read nodata value
          boost::scoped_ptr<SrcImageResource> rsrc(DiskImageResource::open(image_files[f]));
          if (rsrc->has_nodata_read()) {
            nodata_value = rsrc->nodata_read();
            std::cout << "Read nodata: " << nodata_value << std::endl;
          }
          else
            std::cout << "Failed to read nodata value, using default value of -9999.\n";
          
          // TODO: Verify all masks are the same
          // Copy the mask
          image = copy_mask(image, create_mask(input_image, nodata_value));
        }
        // Load this file
        select_channel(image, chan) = pixel_cast<float>(input_image);
        break;
      }
    }
    if (!found)
      vw_throw( ArgumentErr() << "Error: No input file contained landsat channel " << chan+1<< "\n");   
  } // End loop for loading input channels

  // At this point we should have loaded all the required channels.

  // Try to load any useful metadata
  // - Find the metadata file
  for (size_t f=0; f<image_files.size(); ++f) {
    if (image_files[f].find("_MTL.txt") != std::string::npos) {
      // Search the file for the metadata
      std::ifstream handle(image_files[f].c_str());
      std::string line;
      while (std::getline(handle, line)) {
        if (line.find("SUN_ELEVATION") == std::string::npos)
          continue;
        // Parse the line containing the info
        size_t eqpos = line.find("=");
        std::string num = line.substr(eqpos+1);
        sun_elevation = atof(num.c_str());
        break; // Done searching for info
      }
      handle.close();
      break;
    }
  }

  // TODO: Set the image mask!

} // End function load_landsat_image
*/


/// Process MODIS input files into a single image object
void load_landsat_image(LandsatImage &image, std::vector<std::string> const& image_files,
                        int landsat_type, float &sun_elevation) {

  // The bands we want are: BLUE, GREEN, RED, NIR, SWIR1, TEMP, SWIR2
  const int NUM_BANDS_OF_INTEREST = 7;
  const std::string LS5_BAND_LOCATIONS[NUM_BANDS_OF_INTEREST] = {"1", "2", "3", "4", "5", "6", "7"};
  const std::string LS7_BAND_LOCATIONS[NUM_BANDS_OF_INTEREST] = {"1", "2", "3", "4", "5", "6", "8"};
  const std::string LS8_BAND_LOCATIONS[NUM_BANDS_OF_INTEREST] = {"2", "3", "4", "5", "6", "10", "7"};

  const std::string BAND_PREFIX = "_B";
  
  sun_elevation = 45; // The default sun elevation angle
  
  std::vector<std::string> sorted_input_files(NUM_BANDS_OF_INTEREST);
  
  // Check that all the required bands are present
  for (int chan=0; chan<NUM_BANDS_OF_INTEREST; ++chan) {
  
    // For the given output channel, determine the input channel containing it
    std::string text;
    switch(landsat_type) {
      case 8:  text = LS8_BAND_LOCATIONS[chan]; break;
      case 7:  text = LS7_BAND_LOCATIONS[chan]; break;
      default: text = LS5_BAND_LOCATIONS[chan]; break;
    };
    if (text.size() == 1)
      text = "0" + text;
    text = BAND_PREFIX + text + ".TIF";

    //printf("Looking for data %d in channel %d\n", chan, input_channel);
    std::cout << "Looking for text: " << text << std::endl;
    
    // Look for the string in the input file names
    bool found = false;
    for (size_t f=0; f<image_files.size(); ++f) {
      if (image_files[f].find(text) != std::string::npos) {
        found = true;
        sorted_input_files[chan] = image_files[f];
        break;
      }
    }
    if (!found)
      vw_throw( ArgumentErr() << "Error: No input file contained landsat channel " << chan+1<< "\n");
  } // End loop for finding input channels

  // Set up image view object which reads from the input files
  // - Currently we assume that the pixel is invalid when all values are zero!
  image.reset(create_mask(SplitChannelFileView<uint16, NUM_BANDS_OF_INTEREST>(sorted_input_files),
                          Vector<uint16, NUM_BANDS_OF_INTEREST>()));

  // At this point we should have loaded all the required channels.
/*
  // Try to load any useful metadata
  // - Find the metadata file
  for (size_t f=0; f<image_files.size(); ++f) {
    if (image_files[f].find(".txt") != std::string::npos) {
      // Search the file for the metadata
      std::ifstream handle(image_files[f].c_str());
      std::string line;
      while (std::getline(handle, line)) {
        if (line.find("SUN_ELEVATION") == std::string::npos)
          continue;
        // Parse the line containing the info
        size_t eqpos = line.find("=");
        std::string num = line.substr(eqpos+1);
        sun_elevation = atof(num.c_str());
        break; // Done searching for info
      }
      handle.close();
      break;
    }
  }
*/
  // TODO: Set the image mask!

} // End function load_landsat_image




/// Loads the georeference for a Landsat file
void load_landsat_georef(std::vector<std::string> const& image_files,
                         cartography::GeoReference &georef) {

  /// Each channel should contain the same georeference so
  ///  we probably won't need to run more iterations of this loop.
  size_t num_images = image_files.size();
  for (size_t i=0; i<num_images; ++i) {
    std::string path = image_files[i];
    DiskImageResourceGDAL disk_resource(path);
    bool success = read_georeference(georef, disk_resource);
    if (success)
      return;
  }
  
  vw_throw( ArgumentErr() << "Failed to load georef from any input Landsat channel!\n");
}


/// Compute the Floating Algae Index
struct detect_water_fai_functor {
  bool operator()( LandsatPixelType const& pixel) const {
    return pixel[NIR] + (pixel[RED] + (pixel[SWIR1] - pixel[RED])*(170/990));
  }
};

/// Compute NDTI (Turbidity) index
struct detect_water_ndti_functor {
  bool operator()( LandsatPixelType const& pixel) const {
    if (pixel[RED] + pixel[GREEN] == 0)
      return 0; // Avoid divide-by-zero
    return (pixel[RED] - pixel[GREEN]) / (pixel[RED] + pixel[GREEN]);
  }
};

/// Compute NDSI index
struct detect_water_ndsi_functor {
  bool operator()( LandsatPixelType const& pixel) const {
    if (pixel[GREEN] + pixel[SWIR1] == 0)
      return 0; // Avoid divide-by-zero
    return (pixel[GREEN] - pixel[SWIR1]) / (pixel[GREEN] + pixel[SWIR1]);
  }
};

/// Scale the given input range to the 0 to 1 range.
float rescale_to_01(float value, float min, float max) {
  float rng = max - min;
  return ((value - min) / rng);
}

/// Guess if a pixel shows a cloud or not.
bool detect_clouds(LandsatPixelType const& pixel) {
 
  // Compute several indicators of cloudiness and take the minimum of them.
  
  float score = 1.0;
  
  // Clouds are reasonably bright in the blue band.
  score = std::min(score, rescale_to_01(pixel[BLUE], 0.1, 0.3));

  // Clouds are reasonably bright in all visible bands.
  score = std::min(score, rescale_to_01(pixel[RED]+pixel[GREEN]+pixel[BLUE], 0.2, 0.8));

  // Clouds are reasonably bright in all infrared bands.
  score = std::min(score, rescale_to_01(pixel[NIR]+pixel[SWIR1]+pixel[SWIR2], 0.3, 0.8));

  // Clouds are reasonably cool in temperature.
  score = std::min(score, rescale_to_01(pixel[TEMP], 300, 290));

  // However, clouds are not snow.
  detect_water_ndsi_functor f;
  score = std::min(score, rescale_to_01(f(pixel), 0.8, 0.6)); // TODO: Is this correct?

  const float CLOUD_THRESHOLD = 0.35;
  return (score > CLOUD_THRESHOLD);
}

/// Function to scale water detection sensitivity based on sun angle.
float compute_water_threshold(float sun_angle_degrees) {
  return ((0.6 / 54.0) * (62.0 - sun_angle_degrees)) + .05;
}


// TODO: MOVE THIS
template <typename T>
T clamp(T value, T min, T max) {
  if (value > max) return max;
  if (value < min) return min;
  return value;
}

/// Shortcut to clamp between 0 and 1
template <typename T>
T clamp01(T value) {
  return clamp<T>(value, 0, 1);
}


/// Classifies one pixel as water or not.
bool detect_water(LandsatPixelType const& pixel_in, float sun_elevation_degrees=45) {

    return false; // DEBUG

    // Cast the input pixel to float
    const int NUM_USED_BANDS = 7;
    Vector<float, NUM_USED_BANDS> pixel = remove_mask(pixel_in);

    // Check this first!
    if (detect_clouds(pixel))
      return false;
    
    // Compute several indicators of water and take the minimum of them.
    float score = 1.0;

    // Water tends to be dark
    
    // Bands for shadow masking
    float shadow_sum = pixel[NIR] + pixel[SWIR1] + pixel[SWIR2];
    shadow_sum = clamp01(rescale_to_01(shadow_sum, 0.35, 0.2));
    score      = std::min(score, shadow_sum);

    // It also tends to be relatively bright in the blue band
    std::vector<float> dark_values(5);
    dark_values[0] = pixel[GREEN];
    dark_values[1] = pixel[RED];
    dark_values[2] = pixel[NIR];
    dark_values[3] = pixel[SWIR2];
    dark_values[4] = pixel[SWIR1];
    float mean  = math::mean(dark_values); // TODO: Check calls and divide-by-zero
    float std   = math::standard_deviation(dark_values, mean);
    float z;
    if (mean == 0)
      z = 0; // TODO: 0 or 1?
    else {
      z = (pixel[BLUE] - std) / mean;
      z = rescale_to_01(z, 0.0, 1.0);
      z = clamp01(z);
    }
    score = std::min(score, z);

    // Water is at or above freezing
    score = std::min(score, rescale_to_01(pixel[TEMP], 273, 275));

    // Water is nigh in ndsi (aka mndwi)
    float mndwi = 0;
    if (pixel[GREEN] + pixel[SWIR1] != 0)
      mndwi = (pixel[GREEN] - pixel[SWIR1]) / (pixel[GREEN] + pixel[SWIR1]);
    score = clamp01(std::min(score, rescale_to_01(mndwi, 0.3, 0.8)));

    // Select water pixels from the raw score
    float water_thresh = compute_water_threshold(sun_elevation_degrees);
    
    // TODO: Detect snow also

    return (score > water_thresh);
}

/// Use this to call detect_water on each pixel like this:
/// --> = UnaryPerPixelView(landsat_image, landsat::DetectWaterLandsatFunctor());
class DetectWaterLandsatFunctor  : public ReturnFixedType<uint8> {
  float m_sun_elevation_degrees;
public:
  DetectWaterLandsatFunctor(float sun_elevation_degrees=45)
   : m_sun_elevation_degrees(sun_elevation_degrees) {}
  
  uint8 operator()( LandsatPixelType const& pixel) const {
    if (is_valid(pixel)){
      //if (detect_water(pixel, m_sun_elevation_degrees))
        return 255;
      return 1;
    }
    else
      return 0; // Nodata value
  }
};

void detect_water(std::vector<std::string> const& image_files, 
                  cartography::GdalWriteOptions const& write_options) {

  LandsatImage ls_image;
  int landsat_type= 8; // TODO: Where to get this!
  float sun_elevation_degrees;
  load_landsat_image(ls_image, image_files, landsat_type, sun_elevation_degrees);

  // Throws on failure
  cartography::GeoReference georef;
  load_landsat_georef(image_files, georef);

  const uint8 nodata_out = 0; // TODO: Put all of these constants somewhere

/*
  typedef UnaryPerPixelView<LandsatImage, DetectWaterLandsatFunctor> LandsatWaterDetectView;
  block_write_gdal_image("landsat_input.tif",
                         ls_image,
                         true, georef,
                         true, -9999.0,
                         write_options,
                         TerminalProgressCallback("vw", "\t--> DEBUG write input image"));
*/

  Stopwatch timer;
  timer.start();

  // TODO: Setting!
  std::string output_path = "landsat_output.tif";
  typedef UnaryPerPixelView<LandsatImage, DetectWaterLandsatFunctor> LandsatWaterDetectView;
  block_write_gdal_image(output_path,
                         crop(LandsatWaterDetectView(ls_image, DetectWaterLandsatFunctor(sun_elevation_degrees)),
                         BBox2(3496, 6010, 335, 464)),
                         true, georef,
                         true, nodata_out,
                         write_options,
                         TerminalProgressCallback("vw", "\t--> Classifying Landsat:"));

  timer.stop();
  std::cout << "TT time = " << timer.elapsed_seconds() << std::endl;

}


/*

/// Estimates the cloud cover percentage in a Landsat image
def getCloudPercentage(image, region):
    
    
    # The function will attempt the calculation in these ranges
    # - Native Landsat resolution is 30
    MIN_RESOLUTION = 60
    MAX_RESOLUTION = 1000
    resolution = MIN_RESOLUTION

    while True:
        try:
            cloudScore = detect_clouds(image)
            reducedimage = image.reduce(ee.Reducer.allNonZero())
            areaCount = reducedimage.reduceRegion(ee.Reducer.sum(), region, resolution)
            cloudCount = cloudScore.reduceRegion(ee.Reducer.sum(), region, resolution)
            percentage = safe_get_info(cloudCount)['constant'] / safe_get_info(areaCount)['all']
            return percentage
        except Exception as e:
            # Keep trying with lower resolution until we succeed
            resolution = 2*resolution
            #print resolution, MAX_RESOLUTION
            if resolution > MAX_RESOLUTION:
                raise e


*/



}} // end namespace vw::landsat

#endif

