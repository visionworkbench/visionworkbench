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
#include <vw/Image/ImageChannels.h>
#include <vw/Image/Transform.h>
#include <vw/FileIO/DiskImageView.h>
#include <vw/Cartography/GeoReferenceUtils.h>
#include <vw/tools/flood_common.h>

/**
  Tools for processing Landsat data.
  Landsat images are split into multiple geotif files and are medium sized.
  
  Some support for Landsat 5/7 is provided but only Landsat 8 is fully supported.
*/

namespace vw {

  // These definitions are needed to save the landsat image as stored in memory
  template<> struct PixelFormatID<PixelMask<Vector<float,  7> > > { static const PixelFormatEnum value = VW_PIXEL_GENERIC_8_CHANNEL; };
  template<> struct PixelFormatID<PixelMask<Vector<uint16, 7> > > { static const PixelFormatEnum value = VW_PIXEL_GENERIC_8_CHANNEL; };
  template<> struct PixelFormatID<Vector<uint16, 7> > { static const PixelFormatEnum value = VW_PIXEL_GENERIC_7_CHANNEL; };
  template<> struct PixelFormatID<Vector<float, 7> > { static const PixelFormatEnum value = VW_PIXEL_GENERIC_7_CHANNEL; };



namespace landsat {


//----------------------------------------------------
// Landsat type definitions

const int NUM_BANDS_OF_INTEREST = 7;
typedef PixelMask<Vector<uint16, NUM_BANDS_OF_INTEREST> > LandsatPixelType;
typedef PixelMask<Vector<float,  NUM_BANDS_OF_INTEREST> > LandsatToaPixelType;
typedef ImageViewRef<LandsatPixelType> LandsatImage;

// These are the Landsat channels we use, available on both Landsat 7 and Landsat 8.
enum LANDSAT_CHANNEL_INDICES { BLUE  = 0, 
                               GREEN = 1, 
                               RED   = 2, 
                               NIR   = 3, 
                               SWIR1 = 4, 
                               TEMP  = 5, // = Temperature
                               SWIR2 = 6};

// These are the band indices where this data is located in each Landsat sensor.
// - Note that these are zero-based indices.
const int LS5_BAND_LOCATIONS[NUM_BANDS_OF_INTEREST] = {0, 1, 2, 3, 4, 5, 6};
const int LS7_BAND_LOCATIONS[NUM_BANDS_OF_INTEREST] = {0, 1, 2, 3, 4, 5, 7};
const int LS8_BAND_LOCATIONS[NUM_BANDS_OF_INTEREST] = {1, 2, 3, 4, 5, 9, 6};

/// Given an input channel, see what the matching output channel is.
/// - Returns -1 if it does not have an output channel.
int get_output_channel(int input_channel, int landsat_type) {

  // Load the corrrect list of input channels
  const int* list = 0;
  switch(landsat_type) {
    case 8:  list = LS8_BAND_LOCATIONS; break;
    case 7:  list = LS7_BAND_LOCATIONS; break;
    default: list = LS5_BAND_LOCATIONS; break;
  };
  
  // Iterate through the input channels
  for (int i=0; i<NUM_BANDS_OF_INTEREST; ++i) {
    if (list[i] == input_channel)
      return i; // Found channel in the output list
  }
  return -1; // Not in the output list!
}
   
//----------------------------------------------------

/// View for treating multiple single channel images on disk as one multi-channel image.
/// - T is the data type, N is the number of channels.
/// - TODO: Move this into VW
template <typename T, int N>
class SplitChannelFileView : public ImageViewBase<SplitChannelFileView<T, N> > {

public: // Definitions

  typedef Vector<T, N> pixel_type;
  typedef pixel_type   result_type;

public: // Functions

  /// Constructor
  SplitChannelFileView( std::vector<std::string> const& input_files)
                  : m_input_files(input_files){
    // Verify that the correct number of input files were passed in.
    if (input_files.size() != N)
      vw_throw( ArgumentErr() << "SplitChannelFileView created with incorrect number of input files!\n");
      
    // Record the image size to save time later
    boost::shared_ptr<vw::DiskImageResource> resource1(vw::DiskImageResourcePtr(input_files[0]));
    m_cols = resource1->format().cols;
    m_rows = resource1->format().rows;
    
    // Make sure all images are the same size
    for (size_t i=0; i<input_files.size(); ++i) {
      boost::shared_ptr<vw::DiskImageResource> resource2(vw::DiskImageResourcePtr(input_files[i]));
            
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
    return prerasterize_type(tile, -bbox.min().x(), -bbox.min().y(), cols(), rows());

  } // End prerasterize function

 template <class DestT>
 inline void rasterize( DestT const& dest, BBox2i const& bbox ) const {
   vw::rasterize( prerasterize(bbox), dest, bbox );
 }
 
private: // Variables

  std::vector<std::string> const m_input_files;
  int m_cols, m_rows;
 
}; // End class SplitChannelFileView


/// Process Landsat input files into a single image object
void load_landsat_image(LandsatImage &image, std::vector<std::string> const& image_files,
                        int landsat_type) {

  const std::string BAND_PREFIX = "_B";
  
  std::vector<std::string> sorted_input_files(NUM_BANDS_OF_INTEREST);
  
  // Check that all the required bands are present
  for (int chan=0; chan<NUM_BANDS_OF_INTEREST; ++chan) {
  
    // For the given output channel, determine the input channel containing it
    // - Convert index to 1-based string.
    std::string text1, text2;
    switch(landsat_type) {
      case 8:  text1 = num2str(LS8_BAND_LOCATIONS[chan]+1); break;
      case 7:  text1 = num2str(LS7_BAND_LOCATIONS[chan]+1); break;
      default: text1 = num2str(LS5_BAND_LOCATIONS[chan]+1); break;
    };
    // Complete the search pattern, and if needed add an alternate search pattern with a zero padding.
    text2 = "";
    if (text1.size() == 1) {
      text2 = "0" + text1;
      text2 = BAND_PREFIX + text2 + ".TIF";
    }
    text1 = BAND_PREFIX + text1 + ".TIF";
    
    //printf("Looking for data %d in channel %d\n", chan, input_channel);
    
    // Look for the string in the input file names
    bool found = false;
    for (size_t f=0; f<image_files.size(); ++f) {
      // Look for all applicable search strings
      found = (image_files[f].find(text1) != std::string::npos);
      if (!found && !text2.empty())
        found = (image_files[f].find(text2) != std::string::npos);

      // If we got a match, this is the channel we need!
      if (found) {
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

} // End function load_landsat_image


/// Extract numeric the value from one line of the metadata file.
int get_band_number_from_line(std::string const& line) {
  // Parse the line containing the info
  size_t us_pos   = line.rfind("_")+1;
  size_t stop_pos = line.rfind("=")-1;
  std::string num = line.substr(us_pos, stop_pos-us_pos);
  return atoi(num.c_str()) - 1;
}

/// Parses a landsat metadata line to find the value and the location to store it in.
template<size_t N>
bool update_vector_from_line(std::string const& line, 
                             std::string const& prefix, int landsat_type,
                             Vector<float, N> &vec) {
  // See if line contains the desired prefix
  if (line.find(prefix) == std::string::npos)
    return false;
  
  // Check if we are using this band and if so in which output channel.
  int band = get_band_number_from_line(line);
  int output_band = get_output_channel(band, landsat_type);
  if (output_band >= 0) {
    float value = parse_metadata_line(line);
    vec[output_band] = value; // Store result
  }
  return true;
}

/// Convenience structure for storing Landsat metadata information.
struct LandsatMetadataContainer {

  typedef Vector<float, NUM_BANDS_OF_INTEREST> CoefficientVector;
  
  // These are all values from the config used to convert to TOA values.
  CoefficientVector rad_mult;
  CoefficientVector rad_add;
  CoefficientVector toa_mult;
  CoefficientVector toa_add;
  Vector<float, 4>  k_constants;
  float sun_elevation_degrees;
};

/// Loads miscellaneous information from the metadata file
/// - Tuned for USGS format Landsat 8 metadata file
void load_landsat_metadata(std::vector<std::string> const& image_files,
                           int landsat_type,
                           LandsatMetadataContainer & metadata) {
  // Find the metadata file
  std::string metadata_path = find_string_in_list(image_files,".txt");
  if (metadata_path.empty())
    vw_throw( ArgumentErr() << "Error: Landsat metadata file not found!\n");

  metadata.sun_elevation_degrees = 0; // Init to flag value

  // Search the file for the metadata
  std::ifstream handle(metadata_path.c_str());
  std::string line;
  while (std::getline(handle, line)) {
    // Check line for elevation angle
    if (line.find("SUN_ELEVATION") != std::string::npos) {
      metadata.sun_elevation_degrees = parse_metadata_line(line);
      continue;
    }

    // Check line for the radiance conversion constants
    if (update_vector_from_line(line, "RADIANCE_MULT_BAND_", landsat_type, metadata.rad_mult))
      continue;
    if (update_vector_from_line(line, "RADIANCE_ADD_BAND_", landsat_type, metadata.rad_add))
      continue;

    // Chek line for TOA conversion constants
    if (update_vector_from_line(line, "REFLECTANCE_MULT_BAND_", landsat_type, metadata.toa_mult))
      continue;
    if (update_vector_from_line(line, "REFLECTANCE_ADD_BAND_", landsat_type, metadata.toa_add))
      continue;

    // Read in the four TIR temperature constants
    if (line.find("K1_CONSTANT_BAND_10") != std::string::npos)
      metadata.k_constants[0] = parse_metadata_line(line);
    if (line.find("K1_CONSTANT_BAND_11") != std::string::npos)
      metadata.k_constants[1] = parse_metadata_line(line);
    if (line.find("K2_CONSTANT_BAND_10") != std::string::npos)
      metadata.k_constants[2] = parse_metadata_line(line);
    if (line.find("K2_CONSTANT_BAND_11") != std::string::npos)
      metadata.k_constants[3] = parse_metadata_line(line);

  } // End of search through metadata lines.
  handle.close();

  if ((metadata.sun_elevation_degrees == 0) || (metadata.toa_mult[0] == 0))
    vw_throw( ArgumentErr() << "Error: Failed to read in required metadata!\n");
    
  // Correct the TOA conversions to account for sun elevation angle
  double sun_elevation_radians = DEG_TO_RAD*metadata.sun_elevation_degrees;
  metadata.toa_mult /= sin(sun_elevation_radians);
  metadata.toa_add  /= sin(sun_elevation_radians);
} // End load_landsat_metadata


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

/// Convert an input landsat pixel to top-of-atmosphere.
/// - Currently only works for Landsat 8
LandsatToaPixelType convert_to_toa(LandsatPixelType const& pixel_in,
                                   LandsatMetadataContainer const& metadata)
{             
  // This will convert the non-temperature bands to TOA
  LandsatToaPixelType pixel_f = pixel_cast<LandsatToaPixelType>(pixel_in);
  LandsatToaPixelType pixel   = pixel_f;
  for (int i=0; i<NUM_BANDS_OF_INTEREST; ++i)
    pixel[i] = pixel_f[i]*metadata.toa_mult[i] + metadata.toa_add[i];
  
  // Now convert the temperature band
  // - This is hard coded to use the K constants from LS8 Band 10!
  float temp_rad = pixel_f[TEMP]*metadata.rad_mult[TEMP] + metadata.rad_add[TEMP];
  pixel[TEMP] = metadata.k_constants[2] / log(metadata.k_constants[0]/temp_rad + 1.0);

  return pixel;
}

/// Functor wrapper for TOA conversion function
class LandsatToaFunctor  : public ReturnFixedType<LandsatToaPixelType > {
  LandsatMetadataContainer m_metadata;
public:
  LandsatToaFunctor(LandsatMetadataContainer const& metadata)
   : m_metadata(metadata) {}
  
  LandsatToaPixelType operator()( LandsatPixelType const& pixel) const {
    return convert_to_toa(pixel, m_metadata);
  }
};


/// Compute the Floating Algae Index
struct detect_water_fai_functor {
  bool operator()( LandsatToaPixelType const& pixel) const {
    return pixel[NIR] + (pixel[RED] + (pixel[SWIR1] - pixel[RED])*(170/990));
  }
};

/// Compute NDTI (Turbidity) index
struct detect_water_ndti_functor {
  bool operator()( LandsatToaPixelType const& pixel) const {
    return compute_index(pixel[RED], pixel[GREEN]);
  }
};

/// Compute NDSI index
struct detect_water_ndsi_functor {
  bool operator()( LandsatToaPixelType const& pixel) const {
    return compute_index(pixel[GREEN], pixel[SWIR1]);
  }
};

// The following water detection functions are copied from our
//  Earth Engine files and have a bunch of hard coded constants.

/// Guess if a pixel shows a cloud or not.
bool detect_clouds(LandsatToaPixelType const& pixel) {
 
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
  score = std::min(score, rescale_to_01(f(pixel), 0.8, 0.6));

  const float CLOUD_THRESHOLD = 0.35;
  return (score > CLOUD_THRESHOLD);
}

/// Function to scale water detection sensitivity based on sun angle.
float compute_water_threshold(float sun_angle_degrees) {
  float thresh = ((0.6 / 54.0) * (62.0 - sun_angle_degrees));
  
  // Unfortunately the computation above can return a negative value!
  // - For now just set a lower bound on the threshold to avoid this problem.
  const float MIN_THRESH = 0.05;
  if (thresh < MIN_THRESH)
    return MIN_THRESH;
  return thresh;
}


/// Classifies one pixel as water or not.
float detect_water(LandsatToaPixelType const& pixel) {

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
    float mean  = math::mean(dark_values);
    float std   = math::standard_deviation(dark_values, mean);
    float z;
    if (mean == 0)
      z = 1.0;
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

    return score;
}

/// Use this to call detect_water on each pixel like this:
/// --> = per_pixel_view(landsat_image, landsat::DetectWaterLandsatFunctor());
class DetectWaterLandsatFunctor  : public ReturnFixedType<uint8> {
  float m_water_thresh;
public:
  DetectWaterLandsatFunctor(float sun_elevation_degrees=45, double sensitivity=1.0)
   : m_water_thresh(compute_water_threshold(sun_elevation_degrees)*sensitivity) {}
  
  uint8 operator()( LandsatToaPixelType const& pixel) const {
    if (is_valid(pixel)){
      if (detect_water(pixel) > m_water_thresh)
        return FLOOD_DETECT_WATER;
      return FLOOD_DETECT_LAND;
    }
    else
      return FLOOD_DETECT_NODATA;
  }
};


/// DEBUG class to output the raw water detection score.
class DetectWaterLandsatFunctorScore  : public ReturnFixedType<float> {
public: 
  float operator()( LandsatToaPixelType const& pixel) const {
    if (is_valid(pixel))
      return detect_water(pixel);
    else
      return -1.0;
  }
};


/// The sensitivity option multiplies the water detection threshold and is
/// an easy way to adjust the gross sensitivity of the algorithm.
void detect_water(std::vector<std::string> const& image_files, std::string const& output_path,
                  vw::GdalWriteOptions const& write_options, 
                  double sensitivity=1.0, bool debug=false) {

  const int landsat_type= 8; // This could be a parameter but for now it is the only option.
  LandsatImage ls_image;
  load_landsat_image(ls_image, image_files, landsat_type);

  cartography::GeoReference georef;
  load_landsat_georef(image_files, georef);

  
  LandsatMetadataContainer metadata;
  load_landsat_metadata(image_files, landsat_type, metadata);

  if (debug) {
    std::cout << "Read metadata:\n";
    std::cout << "rad_mult "    << metadata.rad_mult    << std::endl;
    std::cout << "rad_add "     << metadata.rad_add     << std::endl;
    std::cout << "toa_mult "    << metadata.toa_mult    << std::endl;
    std::cout << "toa_add "     << metadata.toa_add     << std::endl;
    std::cout << "k_constants " << metadata.k_constants << std::endl;
    std::cout << "sun_elevation " <<  metadata.sun_elevation_degrees << std::endl;
    std::cout << "water threshold = " << 
              compute_water_threshold(metadata.sun_elevation_degrees)*sensitivity << std::endl;

    block_write_gdal_image("landsat_toa_input.tif",
                           apply_mask(per_pixel_view(ls_image, LandsatToaFunctor(metadata)), 
                                      Vector<float,  NUM_BANDS_OF_INTEREST>()),
                           true, georef,
                           false, 0,
                           write_options,
                           TerminalProgressCallback("vw", "\t--> Writing TOA input:"));
                             
    block_write_gdal_image("landsat_raw_output.tif",
                           apply_mask(
                             per_pixel_view(
                                per_pixel_view(
                                  ls_image,
                                  LandsatToaFunctor(metadata)
                                ),
                                DetectWaterLandsatFunctorScore()
                             ),
                             -1.0
                           ),
                           true, georef,
                           true, -1.0,
                           write_options,
                           TerminalProgressCallback("vw", "\t--> Writing raw output:"));
  }


  block_write_gdal_image(output_path,
                         apply_mask(
                           per_pixel_view(
                              per_pixel_view(
                                ls_image,
                                LandsatToaFunctor(metadata)
                              ),
                              DetectWaterLandsatFunctor(metadata.sun_elevation_degrees, sensitivity)
                           ),
                           FLOOD_DETECT_NODATA
                         ),
                         true, georef,
                         true, FLOOD_DETECT_NODATA,
                         write_options,
                         TerminalProgressCallback("vw", "\t--> Classifying Landsat:"));
} // End detect_water


}} // end namespace vw::landsat

#endif

