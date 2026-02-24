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


#ifndef __VW_MODIS_UTILITIES_H__
#define __VW_MODIS_UTILITIES_H__

#include <stdlib.h>
#include <vw/Image/PixelMask.h>
#include <vw/Image/Manipulation.h>
#include <vw/Image/ImageChannels.h>
#include <vw/Image/Transform.h>
#include <vw/FileIO/DiskImageView.h>
#include <vw/Cartography/GeoReferenceUtils.h>

/**
  Tools for processing MODIS data.
  MODIS images are in HDF format and are fairly small.
  
  Development of MODIS code was stopped early and this file is incomplete.
  All MODIS algorithms were originally developed using Google Earth Engine as part
  of the Crisis Mapping project.
*/

namespace vw {
namespace modis {


  //template<> struct PixelFormatID<float> { static const PixelFormatEnum value = VW_PIXEL_GENERIC_1_CHANNEL; };

//----------------------------------------------------
// MODIS types

typedef PixelMask<Vector<float, 9> > ModisPixelType;
typedef ImageView<ModisPixelType> ModisImage;

enum MODIS_CHANNEL_INDICES { B1 = 0, 
                             B2 = 1, 
                             B3 = 2, 
                             B4 = 3, 
                             B5 = 4, 
                             B6 = 5, 
                             B7 = 6,
                             QC250 = 7,
                             QC500 = 8};

// Types to conveniently package MODIS derived products
typedef PixelMask<Vector<float, 4> > ModisProductPixelType;
typedef ImageView<ModisProductPixelType> ModisProductImage;
enum MODIS_PRODUCT_NAME { NDVI=0, NDWI=1, EVI=2, LSWI=3 };

//----------------------------------------------------


/// Determine what scale factor needs to be applied to a file 
///  in order to get it to the maximum MODIS resolution.
/// - Currently this requires fixed file naming conventions and is not robust!
double get_modis_file_scale_factor(std::string const& path) {

  if ((path.find("sur_refl_b01") != std::string::npos) || 
      (path.find("sur_refl_b02") != std::string::npos) || 
      (path.find("QC_250m"     ) != std::string::npos) ) {
    return 1.0; // 250m
  }
  return 2.0; // 500m
  // TODO: Handle 1000m
}

/// Process MODIS input files into a single image object
void load_modis_image(ModisImage &image, std::vector<std::string> const& image_files) {

  // This is the maximum number of channels we will load
  const size_t NUM_MODIS_CHANNELS = 9;

  // These are the hard coded indices where we will store the names.
  std::vector<std::string> channel_names(NUM_MODIS_CHANNELS);
  channel_names[B1] = "sur_refl_b01_1";
  channel_names[B2] = "sur_refl_b02_1";
  channel_names[B3] = "sur_refl_b03_1";
  channel_names[B4] = "sur_refl_b04_1";
  channel_names[B5] = "sur_refl_b05_1";
  channel_names[B6] = "sur_refl_b06_1";
  channel_names[B7] = "sur_refl_b07_1";
  channel_names[QC250] = "QC_250m_1";
  channel_names[QC500] = "QC_500m_1";
      

  // Load each input image and store in the correct channel
  Vector2i image_size;
  ImageView<float> input_image;

  for (size_t i=0; i<image_files.size(); ++i) {
    std::string path = image_files[i];
    double scale = get_modis_file_scale_factor(path);

   // TODO: Manual type conversion?
   
   try {
     // Load and upsample the input image
     input_image = resample(DiskImageView<int>(path), scale);
     
     //printf("image size: %d %d\n", input_image.rows(), input_image.cols());
     
     // Init the output image
     if (image.rows() == 0) {
       image_size[0] = input_image.cols();
       image_size[1] = input_image.rows();
       image.set_size(image_size[0], image_size[1]);
     }
     
     // Find the channel this image goes in
     size_t channel = 999;
     for (size_t c=0; c<NUM_MODIS_CHANNELS; ++c) {
       if (path.find(channel_names[c]) != std::string::npos) {
         channel = c;
         break;
       }
     }
     if (channel >= NUM_MODIS_CHANNELS) {
       std::cout << "Failed to find channel for input " << path << std::endl;
       continue;
     }
     // Copy input image data to output image
     std::cout << "Storing image " << path << " in channel " << channel
               << " with scale " << scale << std::endl;
     select_channel(image, channel) = input_image;
   } catch(...) {
   }
  }

  // TODO: Set the image mask!

} // End function load_modis_image


/// Loads the georeference for a MODIS file
void load_modis_georef(std::vector<std::string> const& image_files,
                       cartography::GeoReference &georef) {

  size_t num_images = image_files.size();
  for (size_t i=0; i<num_images; ++i) {
    std::string path = image_files[i];
    DiskImageResourceGDAL disk_resource(path);
    bool success = read_georeference(georef, disk_resource);
    if (!success)
      continue;

    // Scale the georef of lower resolution channels    
    double scale = get_modis_file_scale_factor(path);
    georef = cartography::resample(georef, scale);
      
    return;
  }
  
  vw_throw( ArgumentErr() << "Failed to load georef from any input MODIS channel!\n");
}


//===========================================================================
// The next set of functions compute basic MODIS products
/*
void compute_ndvi(ModisImage const& modis, ImageView<float> &ndvi) {
  ndvi.set_size(modis.cols(), modis.rows());
  for (int r=0; r<modis.rows(); ++r) {
    for (int c=0; c<modis.cols(); ++c) {
      ModisPixelType pixel = modis(r,c);
      ndvi(c,r) = (pixel[1] - pixel[0]) / (pixel[1] + pixel[0]);
    }
  }
}

/// Normalized difference water index
void compute_ndwi(ModisImage const& modis, ImageView<float> &ndwi) {
  ndwi.set_size(modis.cols(), modis.rows());
  for (int r=0; r<modis.rows(); ++r) {
    for (int c=0; c<modis.cols(); ++c) {
      ModisPixelType pixel = modis(r,c);
      ndwi(c,r) = (pixel[0] - pixel[5]) / (pixel[0] + pixel[5]);
    }
  }
}

/// Enhanced vegetation index
void compute_evi(ModisImage const& modis, ImageView<float> &evi) {
  evi.set_size(modis.cols(), modis.rows());
  for (int r=0; r<modis.rows(); ++r) {
    for (int c=0; c<modis.cols(); ++c) {
      ModisPixelType pixel = modis(r,c);
      evi(c,r) = (pixel[1] - pixel[0])*2.5 / (pixel[1] + 6.0*pixel[0] - 7.5*pixel[2] + 1.0);
    }
  }
}

/// Land surface water index
void compute_lswi(ModisImage const& modis, ImageView<float> &lswi) {
  lswi.set_size(modis.cols(), modis.rows());
  for (int r=0; r<modis.rows(); ++r) {
    for (int c=0; c<modis.cols(); ++c) {
      ModisPixelType pixel = modis(r,c);
      lswi(c,r) = (pixel[1] - pixel[5]) / (pixel[1] + pixel[5]);
    }
  }
}
*/

/// Forms a convenience image consisting of all the modis image products
void form_modis_product_image(ModisImage const& modis, ModisProductImage &products) {

  // Loop through the output image and compute each product at each pixel
  products.set_size(modis.cols(), modis.rows());
  for (int r=0; r<modis.rows(); ++r) {
    for (int c=0; c<modis.cols(); ++c) {
      ModisPixelType pixel = modis(c,r);
      ModisProductPixelType pixel_out;
      pixel_out[NDVI] = 0;
      pixel_out[NDWI] = 0;
      pixel_out[EVI ] = 0;
      pixel_out[LSWI] = 0;
      
      if (pixel[B2] + pixel[B1] != 0)
        pixel_out[NDVI] = (pixel[B2] - pixel[B1]) / (pixel[B2] + pixel[B1]);
      if (pixel[B1] + pixel[B6] != 0)
        pixel_out[NDWI] = (pixel[B1] - pixel[B6]) / (pixel[B1] + pixel[B6]);
      if ((6.0*pixel[B1] + pixel[B2] - 7.5*pixel[B3] + 1.0) != 0)
        pixel_out[EVI ] = 2.5*(pixel[B2] - pixel[B1]) / (6.0*pixel[B1] + pixel[B2] - 7.5*pixel[B3] + 1.0);
      if (pixel[B2] + pixel[B6] != 0)
        pixel_out[LSWI] = (pixel[B2] - pixel[B6]) / (pixel[B2] + pixel[B6]);
      products(c,r) = pixel_out;
    }
  }

}} // end namespace vw::modis


//===========================================================================

} // end namespace modis

/*
def compute_modis_indices(domain):
    '''Compute several common interpretations of the MODIS bands'''
    
    band1 = domain.modis.sur_refl_b01 # pRED
    band2 = domain.modis.sur_refl_b02 # pNIR

    # Other bands must be used at lower resolution
    band3 = domain.modis.sur_refl_b03 # pBLUE
    band4 = domain.modis.sur_refl_b04
    band5 = domain.modis.sur_refl_b05.mask(ee.Image(1)) # band5 masks zero pixels...
    band6 = domain.modis.sur_refl_b06 # pSWIR
    
    # Convenience measure
    DVEL = EVI.subtract(LSWI)



def getQABits(image, start, end, newName):
    '''Extract bits from positions "start" to "end" in the image'''
    # Create a bit mask of the bits we need
    pattern = 0
    for i in range(start,end):
       pattern += 2**i
    # Extract the bits, shift them over, and rename the channel.
    temp = ee.Image(pattern)
    return image.select([0], [newName]).bitwise_and(temp).rightShift(start)

def getModisBadPixelMask(lowResModis):
    '''Retrieves the 1km MODIS bad pixel mask (identifies clouds)'''

    # Select the QA band
    qaBand = lowResModis.select('state_1km').uint16()
   
    # Get the cloud_state bits and find cloudy areas.
    cloudBits = getQABits(qaBand, 0, 1, 'cloud_state')
    cloud = cloudBits.eq(1).Or(cloudBits.eq(2))

    return cloud # The second part of this, the land water flag, does not work well at all.
    
    ## Get the land_water_flag bits.
    #landWaterFlag = getQABits(qaBand, 3, 5, 'land_water_flag')
    #
    ## Create a mask that filters out deep ocean and cloudy areas.
    #mask = landWaterFlag.neq(7).And(cloud.Not())
    #return mask

def getCloudPercentage(lowResModis, region):
    '''Returns the percentage of a region flagged as clouds by the MODIS metadata'''

    MODIS_CLOUD_RESOLUTION = 1000 # Clouds are flagged at this resolution

    # Divide the number of cloud pixels by the total number of pixels
    oneMask    = ee.Image(1.0) 
    cloudMask  = getModisBadPixelMask(lowResModis)
    areaCount  = oneMask.reduceRegion(  ee.Reducer.sum(), region, MODIS_CLOUD_RESOLUTION)
    cloudCount = cloudMask.reduceRegion(ee.Reducer.sum(), region, MODIS_CLOUD_RESOLUTION)
    percentage = safe_get_info(cloudCount)['cloud_state'] / safe_get_info(areaCount)['constant']

    return percentage

def get_permanent_water_mask():
    return ee.Image("MODIS/MOD44W/MOD44W_005_2000_02_24").select(['water_mask'], ['b1'])



# If mixed_thresholds is true, we find the thresholds that contain 0.05 land and 0.95 water,
#  otherwise we find the threshold that most accurately splits the training data.
def compute_binary_threshold(valueImage, classification, bounds, mixed_thresholds=False):
    '''Computes a threshold for a value given examples in a classified binary image'''
    
    # Build histograms of the true and false labeled values
    valueInFalse   = valueImage.mask(classification.Not())
    valueInTrue    = valueImage.mask(classification)
    NUM_BINS       = 128
    SCALE          = 250 # In meters
    histogramFalse = safe_get_info(valueInFalse.reduceRegion(ee.Reducer.histogram(NUM_BINS, None, None), bounds, SCALE))['b1']
    histogramTrue  = safe_get_info(valueInTrue.reduceRegion( ee.Reducer.histogram(NUM_BINS, None, None), bounds, SCALE))['b1']
    
    # Get total number of pixels in each histogram
    false_total = sum(histogramFalse['histogram'])
    true_total  = sum(histogramTrue[ 'histogram'])
    
    # WARNING: This method assumes that the false histogram is composed of greater numbers than the true histogram!!
    #        : This happens to be the case for the three algorithms we are currently using this for.
    
    false_index       = 0
    false_sum         = false_total
    true_sum          = 0.0
    threshold_index   = None
    lower_mixed_index = None
    upper_mixed_index = None
    for i in range(len(histogramTrue['histogram'])): # Iterate through the bins of the true histogram
        # Add the number of pixels in the current true bin
        true_sum += histogramTrue['histogram'][i]
        
        # Set x equal to the max end of the current bin
        x = histogramTrue['bucketMin'] + (i+1)*histogramTrue['bucketWidth']
        
        # Determine the bin of the false histogram that x falls in
        # - Also update the number of 
        while ( (false_index < len(histogramFalse['histogram'])) and
                (histogramFalse['bucketMin'] + false_index*histogramFalse['bucketWidth'] < x) ):
            false_sum   -= histogramFalse['histogram'][false_index] # Remove the pixels from the current false bin
            false_index += 1 # Move to the next bin of the false histogram
    
        percent_true_under_thresh = true_sum/true_total
        percent_false_over_thresh = false_sum/false_total
            
        if mixed_thresholds:
            if (false_total - false_sum) / float(true_sum) <= 0.05:
                lower_mixed_index = i
            if upper_mixed_index == None and (true_total - true_sum) / float(false_sum) <= 0.05:
                upper_mixed_index = i
        else:
            if threshold_index == None and (percent_false_over_thresh < percent_true_under_thresh) and (percent_true_under_thresh > 0.5):
                break

    
    if mixed_thresholds:
        if (not lower_mixed_index) or (not upper_mixed_index):
            raise Exception('Failed to compute mixed threshold values!')
        lower = histogramTrue['bucketMin'] + lower_mixed_index * histogramTrue['bucketWidth'] + histogramTrue['bucketWidth']/2
        upper = histogramTrue['bucketMin'] + upper_mixed_index * histogramTrue['bucketWidth'] + histogramTrue['bucketWidth']/2
        if lower > upper:
            temp  = lower
            lower = upper
            upper = temp
        print 'Thresholds (%g, %g) found.' % (lower, upper)
        return (lower, upper)
    else:
        # Put threshold in the center of the current true histogram bin/bucket
        threshold = histogramTrue['bucketMin'] + i*histogramTrue['bucketWidth'] + histogramTrue['bucketWidth']/2
        print 'Threshold %g Found. %g%% of water pixels and %g%% of land pixels separated.' % \
            (threshold, true_sum / true_total * 100.0, false_sum / false_total * 100.0)
        return threshold


def compute_dem_slope_degrees(dem, resolution):
    '''Computes a slope in degrees for each pixel of the DEM'''
    
    deriv = dem.derivative()
    dZdX    = deriv.select(['elevation_x']).divide(resolution)
    dZdY    = deriv.select(['elevation_y']).divide(resolution)
    slope = dZdX.multiply(dZdX).add(dZdY.multiply(dZdY)).sqrt().reproject("EPSG:4269", None, resolution); 
    RAD2DEG = 180 / 3.14159
    slopeAngle = slope.atan().multiply(RAD2DEG);
    return slopeAngle


def apply_dem(domain, water_fraction, speckleFilter=False):
    '''Convert a low resolution water fraction map into a higher resolution DEM-based flooded pixel map'''
    
    MODIS_PIXEL_SIZE_METERS = 250
    ## Treating the DEM values contained in the MODIS pixel as a histogram, find the N'th percentile
    ##  where N is the water fraction computed by DNNS.  That should be the height of the flood water.
    #modisPixelKernel = ee.Kernel.square(MODIS_PIXEL_SIZE_METERS, 'meters', False)
    #dem.mask(water_fraction).reduceNeighborhood(ee.Reducer.percentile(), modisPixelKernel)
    # --> We would like to compute a percentile here, but this would require a different reducer input for each pixel!
    
    # Get whichever DEM is loaded in the domain
    dem = domain.get_dem().image
    
    water_present = water_fraction.gt(0.0)#.mask(water_high)
    
    if speckleFilter: # Filter out isolated pixels
        water_present_despeckle = water_present.focal_min(500, 'circle', 'meters').focal_max(500, 'circle', 'meters')
        water_fraction = water_fraction.multiply(water_present_despeckle)
    
    # Get min and max DEM height within each water containing pixel
    # - If a DEM pixel contains any water then the water level must be at least that high.    
    dem_min = dem.mask(water_fraction).focal_min(MODIS_PIXEL_SIZE_METERS, 'square', 'meters')
    dem_max = dem.mask(water_fraction).focal_max(MODIS_PIXEL_SIZE_METERS, 'square', 'meters')
    
    # Approximation, linearize each tile's fraction point
    # - The water percentage is used as a percentage between the two elevations
    # - Don't include full or empty pixels, they don't give us clues to their height.
    water_high = dem_min.add(dem_max.subtract(dem_min).multiply(water_fraction))
    water_high = water_high.multiply(water_fraction.lt(1.0)).multiply(water_fraction.gt(0.0)) 
    
    
    # Problem: Averaging process spreads water way out to pixels where it was not detected!!
    #          Reducing the averaging is a simple way to deal with this and probably does not hurt results at all
    
    #dilate_kernel = ee.Kernel.circle(250, 'meters', False)
    #allowed_water_mask = water_fraction.gt(0.0)#.convolve(dilate_kernel)
    
    ## Smooth out the water elevations with a broad kernel; nearby pixels probably have the same elevation!
    water_dem_kernel        = ee.Kernel.circle(5000, 'meters', False)
    num_nearby_water_pixels = water_high.gt(0.0).convolve(water_dem_kernel)
    average_high            = water_high.convolve(water_dem_kernel).divide(num_nearby_water_pixels)
    
    # Alternate smoothing method using available tool EE
    #connected_water = ee.Algorithms.ConnectedComponentLabeler(water_present, water_dem_kernel, 256) # Perform blob labeling
    #addToMap(water_present,   {'min': 0, 'max':   1}, 'water present', False);
    #addToMap(connected_water, {'min': 0, 'max': 256}, 'labeled blobs', False);

    #addToMap(water_high, {'min': 0, 'max': 100}, 'water_high', False);    
    #addToMap(water_fraction, {'min': 0, 'max':   1}, 'Water Fraction', False);
    #addToMap(dem.subtract(average_high), {min : -0, max : 10}, 'Water Difference', false);
    #addToMap(dem.lte(average_high).and(domain.groundTruth.not()));
    
    #addToMap(allowed_water_mask, {'min': 0, 'max': 1}, 'allowed_water', False);
    #addToMap(average_high, {'min': 0, 'max': 100}, 'average_high', False);
    #addToMap(dem, {'min': 0, 'max': 100}, 'DEM', False);
    
    # Classify DEM pixels as flooded based on being under the local water elevation or being completely flooded.
    #return dem.lte(average_high).Or(water_fraction.eq(1.0)).select(['elevation'], ['b1'])
    dem_water = dem.lte(average_high).mask(ee.Image(1))#multiply(water_fraction.gt(0.0)).#.mask(water_fraction) # Mask prevents pixels with 0% water from being labeled as water
    return dem_water.Or(water_fraction.eq(1.0)).select(['elevation'], ['b1'])
*/

#endif

