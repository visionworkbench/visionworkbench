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


#ifndef __VW_MULTISPECTRAL_H__
#define __VW_MULTISPECTRAL_H__

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
  Tools for processing multispectral image data.  Only data from the 
  WorldView and SPOT6/7 satellites are supported.
  
  SPOT was not fully developed because the header files were missing key metadata.
*/

namespace vw {

namespace multispectral {

const int NUM_SPOT67_BANDS    = 5;
const int NUM_WORLDVIEW_BANDS = 8;

/// Band-averaged solar spectral irradiance
/// - Copied from "Radiometric Use of WorldView-2 Imagery"
const float WORLDVIEW_ESUN[NUM_WORLDVIEW_BANDS] = {
  //1580.8140, // PAN
  1758.2229, // Coastal
  1974.2416, // Blue
  1856.4104, // Green
  1738.4791, // Yellow
  1559.4555, // Red
  1342.0695, // Red Edge
  1069.7302, // NIR 1
  861.2866   // NIR 2
};

// Channel indices specify which channels in the file contain which spectra.

/*  This is nominally the correct order, but the HDDS files seem to have a
     different configuration!
enum SPOT67_CHANNEL_INDICES { SPOT_PAN   = 0, 
                              SPOT_BLUE  = 1, 
                              SPOT_GREEN = 2, 
                              SPOT_RED   = 3, 
                              SPOT_NIR   = 4};
*/

enum SPOT67_CHANNEL_INDICES { SPOT_BLUE  = 0, 
                              SPOT_GREEN = 1, 
                              SPOT_RED   = 2, 
                              SPOT_NIR   = 3};

enum WORLDVIEW23_CHANNEL_INDICES { COASTAL  = 0, 
                                   BLUE     = 1, 
                                   GREEN    = 2, 
                                   YELLOW   = 3, 
                                   RED      = 4, 
                                   RED_EDGE = 5,
                                   NIR1     = 6,
                                   NIR2     = 7};

// Define the input and TOA (Top Of Atmosphere) pixel types.
typedef PixelMask<Vector<uint8,  NUM_SPOT67_BANDS   > > Spot67PixelType;
typedef PixelMask<Vector<uint16, NUM_WORLDVIEW_BANDS> > WorldView23PixelType;
typedef PixelMask<Vector<float,  NUM_SPOT67_BANDS   > > Spot67ToaPixelType;
typedef PixelMask<Vector<float,  NUM_WORLDVIEW_BANDS> > WorldView23ToaPixelType;

// We use Ref (delayed action) types for these images in case they are very large.
typedef ImageViewRef<Spot67PixelType     > Spot67Image;
typedef ImageViewRef<WorldView23PixelType> WorldView23Image;

/// Loads an image from either Spot6 or Spot7 (they are the same format)
void load_spot67_image(std::vector<std::string> const& input_paths,
                       Spot67Image & image,
                       cartography::GeoReference &georef) {
  
  std::string image_path = find_string_in_list(input_paths, ".tif");
  if (image_path.empty())
    vw_throw( ArgumentErr() << "Error: SPOT image file not found!\n");
    
  // TODO: Is zero always the nodata value?
  image = create_mask(
            planes_to_channels<Vector<uint8, NUM_SPOT67_BANDS> >(DiskImageView<uint8>(image_path)));
  
  DiskImageResourceGDAL disk_resource(image_path);
  if (!read_georeference(georef, disk_resource)) 
    vw_throw( ArgumentErr() << "Failed to read georeference from image " << image_path <<"\n");
}


/// Load a Worldview 2 or 3 multispectral image
void load_worldview23_image(std::vector<std::string> const& input_paths,
                            WorldView23Image & image,
                            cartography::GeoReference &georef) {

  // Find the image file
  std::string image_path = find_string_in_list(input_paths, ".tif");
  if (image_path.empty())
    vw_throw( ArgumentErr() << "Error: WorldView image file not found!\n");

  // Load 8 bands from one image
  // - The band order is Coastal, Blue, Green, Yellow, Red, Red-Edge, Near-IR1, Near-IR2

  // WorldView 3 can capture some additional bands but that does not affect the
  //  multispectral file format.

  // TODO: Is zero the standard nodata value?
  // - The image is uint16 but only 11 bits are used (max value 2047)
  // - Need to use the planes_to_channels converter to get the image to read properly.
  typedef DiskImageView<uint16> WvDiskView;
  image = create_mask(planes_to_channels<Vector<uint16, NUM_WORLDVIEW_BANDS>,  WvDiskView>(WvDiskView(image_path)));

  DiskImageResourceGDAL disk_resource(image_path);
  if (!read_georeference(georef, disk_resource)) 
    vw_throw( ArgumentErr() << "Failed to read georeference from image " << image_path <<"\n");
}
/*
void read_spot67_metadata() {

TODO: Where are the gain/offset values we need for each band to convert to TOA??

}
*/
/// Convenience structure for storing WorldView metadata information
struct WorldViewMetadataContainer {

  typedef Vector<float, NUM_WORLDVIEW_BANDS> CoefficientVector;
  
  CoefficientVector abs_cal_factor;
  CoefficientVector effective_bandwidth;
  float             mean_sun_elevation; // Units = degrees
  float             earth_sun_distance; // Units = AU
  std::string       datetime;
  
  /// Populate derived values from input values
  void populate_derived_values() {
    // Extract numbers
    // - Input format: "2016-10-23T17:46:54.796950Z;"
    int   year   = atoi(datetime.substr( 0,4).c_str());
    int   month  = atoi(datetime.substr( 5,2).c_str());
    int   day    = atoi(datetime.substr( 8,2).c_str());
    int   hour   = atoi(datetime.substr(11,2).c_str());
    int   minute = atoi(datetime.substr(14,2).c_str());
    float second = atof(datetime.substr(17,8).c_str());
    
    earth_sun_distance = compute_earth_sun_distance(year, month, day, hour, minute, second);
  }
}; // End struct WorldViewMetadataContainer

/// Populate a WorldViewMetadataContainer object from a WorldView metadata file.
void load_worldview23_metadata(std::vector<std::string> const& input_paths,
                               WorldViewMetadataContainer &metadata) {
  // Find the metadata file
  std::string metadata_path = find_string_in_list(input_paths, ".IMD");
  if (metadata_path.empty())
    vw_throw( ArgumentErr() << "Error: WorldView metadata file not found!\n");
    
  // Search the file for the metadata
  std::ifstream handle(metadata_path.c_str());
  std::string line;
  int channel_index = -1, found_count = 0;
  while (std::getline(handle, line)) {
  
    // Check for new group
    if (line.find("BEGIN_GROUP") != std::string::npos) {
      size_t eqpos = line.find("=");
      std::string name = line.substr(eqpos+2);
      channel_index = -1;
      if (name == "BAND_C" ) channel_index = 0;
      if (name == "BAND_B" ) channel_index = 1;
      if (name == "BAND_G" ) channel_index = 2;
      if (name == "BAND_Y" ) channel_index = 3;
      if (name == "BAND_R" ) channel_index = 4;
      if (name == "BAND_RE") channel_index = 5;
      if (name == "BAND_N" ) channel_index = 6;
      if (name == "BAND_N2") channel_index = 7;
      continue;
    }
    
    // Check for values
    if (line.find("absCalFactor") != std::string::npos) {
      if (channel_index < 0)
        vw_throw( ArgumentErr() << "Error reading absCalFactor in metadata file!\n");
      metadata.abs_cal_factor[channel_index] = parse_metadata_line(line);
      ++found_count;
      continue;
    }
    if (line.find("effectiveBandwidth") != std::string::npos) {
      if (channel_index < 0)
        vw_throw( ArgumentErr() << "Error reading effectiveBandwidth in metadata file!\n");
      metadata.effective_bandwidth[channel_index] = parse_metadata_line(line);
      ++found_count;
      continue;
    } 
    if (line.find("meanSunEl") != std::string::npos) {
      metadata.mean_sun_elevation = parse_metadata_line(line);
      ++found_count;
      continue;
    }
    if (line.find("firstLineTime") != std::string::npos) {
      size_t eqpos = line.find("=");
      metadata.datetime = line.substr(eqpos+1);
      ++found_count;
      continue;
    }
  }
  handle.close();

  // Check that we got what we need
  if (found_count != 2*NUM_WORLDVIEW_BANDS+2)
    vw_throw( ArgumentErr() << "Failed to find all required metadata!\n");

  // Compute derived metadata values
  metadata.populate_derived_values();
  
} // End function read_worldview3_metadata



/// Convert an input WorldView pixel to top-of-atmosphere value.
WorldView23ToaPixelType convert_to_toa(WorldView23PixelType       const& pixel_in,
                                       WorldViewMetadataContainer const& metadata)
{             
  // First convert to radiance values
  WorldView23ToaPixelType rad_pixel = pixel_cast<WorldView23ToaPixelType>(pixel_in);
  for (int i=0; i<NUM_WORLDVIEW_BANDS; ++i)
    rad_pixel[i] = rad_pixel[i]*(metadata.abs_cal_factor[i]/metadata.effective_bandwidth[i]);

  // Now convert to reflectance values
  float scale_factor = metadata.earth_sun_distance*metadata.earth_sun_distance*M_PI /
                       cos(DEG_TO_RAD*(90.0 - metadata.mean_sun_elevation));
  WorldView23ToaPixelType pixel = rad_pixel;
  for (int i=0; i<NUM_WORLDVIEW_BANDS; ++i)
    pixel[i] = rad_pixel[i] * scale_factor / WORLDVIEW_ESUN[i];

  return pixel;
}

/// Functor wrapper for TOA conversion function
class WorldView23ToaFunctor : public ReturnFixedType<WorldView23ToaPixelType > {
  WorldViewMetadataContainer m_metadata;
public:
  WorldView23ToaFunctor(WorldViewMetadataContainer const& metadata)
   : m_metadata(metadata) {}
  
  WorldView23ToaPixelType operator()( WorldView23PixelType const& pixel) const {
    return convert_to_toa(pixel, m_metadata);
  }
};

/// Compute NDVI index
float compute_ndvi( WorldView23ToaPixelType const& pixel) {
  return compute_index(pixel[RED], pixel[NIR2]);
}

/// Compute NDWI index
float compute_ndwi( WorldView23ToaPixelType const& pixel) {
  return compute_index(pixel[BLUE], pixel[NIR1]);
}

/// Compute NDWI2 index
/// - Both of these calculations are sometimes listed as "NDWI"
float compute_ndwi2( WorldView23ToaPixelType const& pixel) {
  return compute_index(pixel[COASTAL], pixel[NIR2]);
}

/// Compute SDI index
float compute_sdi( WorldView23ToaPixelType const& pixel) {
  float denom = pixel[NIR2] + pixel[BLUE];
  if (denom == 0)
    return 10; // Avoid divide-by-zero
  return ((pixel[NIR2] - pixel[BLUE]) / denom) - pixel[NIR1];
}

// SPOT functions - note that TOA is not available!

/// Compute NDVI index - SPOT
float compute_ndvi_spot( Spot67PixelType const& pixel) {
  return compute_index(pixel[SPOT_RED], pixel[SPOT_NIR]);
}

/// Compute NDWI index - SPOT
float compute_ndwi_spot( Spot67PixelType const& pixel) {
  return compute_index(pixel[SPOT_BLUE], pixel[SPOT_NIR]);
}

// Functor wrappers for the compute functions above

class NdviFunctor  : public ReturnFixedType<PixelMask<float> > {
public:
  PixelMask<float> operator()( WorldView23ToaPixelType const& pixel) const {
    PixelMask<float> output(compute_ndvi(pixel));
    if (!is_valid(pixel))
      invalidate(output);
    return output;
  }
};
class NdwiFunctor : public ReturnFixedType<PixelMask<float> > {
public:
  PixelMask<float> operator()( WorldView23ToaPixelType const& pixel) const {
    PixelMask<float> output(compute_ndwi(pixel));
    if (!is_valid(pixel))
      invalidate(output);
    return output;
  }
};
class Ndwi2Functor : public ReturnFixedType<PixelMask<float> > {
public:
  PixelMask<float> operator()( WorldView23ToaPixelType const& pixel) const {
    PixelMask<float> output(compute_ndwi2(pixel));
    if (!is_valid(pixel))
      invalidate(output);
    return output;
  }
};
class SdiFunctor : public ReturnFixedType<PixelMask<float> > {
public:
  PixelMask<float> operator()( WorldView23ToaPixelType const& pixel) const {
    PixelMask<float> output(compute_sdi(pixel));
    if (!is_valid(pixel))
      invalidate(output);
    return output;
  }
};

class NdviFunctorSpot67 : public ReturnFixedType<PixelMask<float> > {
public:
  PixelMask<float> operator()( Spot67PixelType const& pixel) const {
    PixelMask<float> output(compute_ndvi_spot(pixel));
    if (!is_valid(pixel))
      invalidate(output);
    return output;
  }
};
class NdwiFunctorSpot67 : public ReturnFixedType<PixelMask<float> > {
public:
  PixelMask<float> operator()( Spot67PixelType const& pixel) const {
    PixelMask<float> output(compute_ndwi_spot(pixel));
    if (!is_valid(pixel))
      invalidate(output);
    return output;
  }
};


/// Use this to call detect_water on each pixel like this:
/// --> = per_pixel_view(landsat_image, landsat::DetectWaterWorldView23Functor());
/// - This seems to work fairly well except that it can confuse cloud shadows with water.
///   The NDVI index seems like the best way to discriminate between them but it is still
///   not very good.  To actually do a good job, would probably need to do something with
///   cloud detection (easier) and sun geometry.
class DetectWaterWorldView23Functor  : public ReturnFixedType<uint8> {
  double m_sensitivity;
public:
  DetectWaterWorldView23Functor(double sensitivity=1.0) : m_sensitivity(sensitivity) {}
  
  uint8 operator()( WorldView23ToaPixelType const& pixel) const {
    if (is_valid(pixel)){

      float ndvi  = compute_ndvi(pixel);
      float ndwi2 = compute_ndwi2(pixel);

      if ((ndvi < 0.1*m_sensitivity) || (ndwi2 < 0.3*m_sensitivity))
        return FLOOD_DETECT_LAND;
      if ((ndvi > 0.5*m_sensitivity) || (ndwi2 > 0.5*m_sensitivity))
        return FLOOD_DETECT_WATER;
      return FLOOD_DETECT_LAND;
    }
    else
      return FLOOD_DETECT_NODATA;
  }
};

/// Very ad-hoc detector for SPOT67.  TOA pixel conversion is not available,
///  and only BGR + NIR bands are available.
class DetectWaterSpot67Functor  : public ReturnFixedType<uint8> {
  double m_sensitivity;
public:
  DetectWaterSpot67Functor(double sensitivity=1.0) : m_sensitivity(sensitivity) {}
  
  uint8 operator()( Spot67PixelType const& pixel) const {
    if (is_valid(pixel)){
      // Very simple way to look for water!
      float ndvi  = compute_ndvi_spot(pixel);
      float ndwi = compute_ndwi_spot(pixel);
      //std::cout << "ndvi " << ndvi << ", ndwi " << ndwi <<std::endl;
      if ((ndwi > 0.3*m_sensitivity) || ((ndvi+ndwi) > 0.6*m_sensitivity))
        return FLOOD_DETECT_WATER;
      return FLOOD_DETECT_LAND;
    }
    else
      return FLOOD_DETECT_NODATA;
  }
};

// High level water detection function for WorldView.
void detect_water_worldview23(std::vector<std::string> const& image_files, 
                              std::string const& output_path,
                              vw::GdalWriteOptions const& write_options,
                              double sensitivity = 1.0, bool debug = false) {

  // Load the image and metadata
  WorldView23Image wv_image;
  cartography::GeoReference georef;
  load_worldview23_image(image_files, wv_image, georef);
 
  WorldViewMetadataContainer metadata;
  load_worldview23_metadata(image_files, metadata);

  if (debug) {
    std::cout << "Loaded metadata: \n";
    std::cout << "abs_cal_factor     "  << metadata.abs_cal_factor      << std::endl;
    std::cout << "effective_bandwidth " << metadata.effective_bandwidth << std::endl;
    std::cout << "mean_sun_elevation "  << metadata.mean_sun_elevation  << std::endl;
    std::cout << "earth_sun_distance "  << metadata.earth_sun_distance  << std::endl;
    std::cout << "datetime "            << metadata.datetime            << std::endl;
    
    block_write_gdal_image("ndvi.tif",
                           apply_mask(
                             per_pixel_view(
                                per_pixel_view(wv_image, WorldView23ToaFunctor(metadata)),
                                NdviFunctor()
                             ),
                             -999
                           ),
                           true, georef, true, -999, write_options,
                           TerminalProgressCallback("vw", "\t--> NDVI"));
/*
    block_write_gdal_image("ndwi.tif",
                           apply_mask(
                             per_pixel_view(
                                per_pixel_view(wv_image, WorldView3ToaFunctor(metadata)),
                                NdwiFunctor()
                             ),
                             -999
                           ),
                           true, georef, true, -999, write_options,
                           TerminalProgressCallback("vw", "\t--> NDWI"));
*/
    block_write_gdal_image("ndwi2.tif",
                           apply_mask(
                             per_pixel_view(
                                per_pixel_view(wv_image, WorldView23ToaFunctor(metadata)),
                                Ndwi2Functor()
                             ),
                             -999
                           ),
                           true, georef, true, -999, write_options,
                           TerminalProgressCallback("vw", "\t--> NDWI2"));    

    block_write_gdal_image("sdi.tif",
                           apply_mask(
                             per_pixel_view(
                                per_pixel_view(wv_image, WorldView23ToaFunctor(metadata)),
                                SdiFunctor()
                             ),
                             -999
                           ),
                           true, georef, true, -999, write_options,
                           TerminalProgressCallback("vw", "\t--> SDI"));
    
  }

  // All the work happens when this call executes
  block_write_gdal_image(output_path,
                         apply_mask(
                           per_pixel_view(
                              per_pixel_view(
                                wv_image,
                                WorldView23ToaFunctor(metadata)
                              ),
                              DetectWaterWorldView23Functor(sensitivity)
                           ),
                           FLOOD_DETECT_NODATA
                         ),
                         true, georef,
                         true, FLOOD_DETECT_NODATA,
                         write_options,
                         TerminalProgressCallback("vw", "\t--> Classifying WorldView:"));
}



// High level water detection function for SPOT.
void detect_water_spot67(std::vector<std::string> const& image_files, 
                         std::string const& output_path,
                         vw::GdalWriteOptions const& write_options,
                         double sensitivity = 1.0, bool debug = false) {

  // Load image, no useful metadata available.
  Spot67Image spot_image;
  cartography::GeoReference georef;
  load_spot67_image(image_files, spot_image, georef);

  if (debug) {
    block_write_gdal_image("ndvi.tif",
                           apply_mask(
                             per_pixel_view(spot_image, NdviFunctorSpot67()),
                             -999
                           ),
                           true, georef, true, -999, write_options,
                           TerminalProgressCallback("vw", "\t--> NDVI"));

    block_write_gdal_image("ndwi.tif",
                           apply_mask(
                             per_pixel_view(spot_image, NdwiFunctorSpot67()),
                             -999
                           ),
                           true, georef, true, -999, write_options,
                           TerminalProgressCallback("vw", "\t--> NDWI"));
  }

  // All the work happens when this call executes
  block_write_gdal_image(output_path,
                         apply_mask(
                           per_pixel_view(
                             spot_image,
                             DetectWaterSpot67Functor(sensitivity)
                           ),
                           FLOOD_DETECT_NODATA
                         ),
                         true, georef,
                         true, FLOOD_DETECT_NODATA,
                         write_options,
                         TerminalProgressCallback("vw", "\t--> Classifying Spot:"));
}




}} // end namespace vw::multispectral

#endif

