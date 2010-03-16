#include <iostream>
#include <fstream>
#include <vw/Image.h>
#include <vw/Image/PixelMath.h>
#include <vw/Image/PixelMask.h>
#include <vw/Image/MaskViews.h>
#include <vw/FileIO.h>
#include <math.h>
#include <time.h>

//using namespace std;
//using namespace vw;

#include <vw/Core.h>
#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Cartography.h>
#include <vw/Math.h>
using namespace vw;
using namespace vw::math;
using namespace vw::cartography;
#include <vw/Photometry/Reflectance.h>
#include <vw/Photometry/Reconstruct.h>
#include <vw/Photometry/Misc.h>
#include <vw/Photometry/Weights.h>

//determines the best guess for the exposure time from the reflectance model
//forces unit exposure time for the first frame
//Ara Nefian
void Init_Exposure(float *reflectanceAvgArray, int numImages, float *exposureTimeArray)
{
  int i;
  exposureTimeArray[0] = 1.0;
  for (i = 0; i < numImages; i++){
       exposureTimeArray[i] = (exposureTimeArray[0]*reflectanceAvgArray[0])/reflectanceAvgArray[i];
       printf("init exposure time[%d]=%f\n", i, exposureTimeArray[i]);
  }  
}

//Ara Nefian
float ComputeGradient_Exposure(float T, float albedo)
{
  float grad;
  grad = T*albedo; 

  return grad;
}

//Ara Nefian
float ComputeError_Exposure(float intensity, float T, float albedo, float reflectance, Vector3 xyz, Vector3 xyz_prior)
{
  float error;
  error = (intensity-T*albedo*reflectance);
  return error;
}

float ComputeError_Exposure(float intensity, float T, float albedo, float reflectance)
{
  float error;
  error = (intensity-T*albedo*reflectance);
  return error;
}

void AppendExposureInfoToFile(string exposureFilename, string currInputFile, modelParams currModelParams)
{
  FILE *fp;
  
  fp = fopen(exposureFilename.c_str(), "a");
  
  fprintf(fp, "%s %f\n", currInputFile.c_str(), currModelParams.exposureTime);

  fclose(fp);
}

std ::vector<float> ReadExposureInfoFile(string exposureFilename, int numEntries)
{
  FILE *fp;
  std::vector<float> exposureTimeVector(numEntries);

  fp = fopen(exposureFilename.c_str(), "r");
  float x, y, z;

  for (unsigned int i = 0; i < numEntries; i++){
       char *filename = new char[500];
       float exposureTime;
       fscanf(fp, "%s %f\n", filename, &exposureTime);
       printf("%f\n", exposureTime);
       delete filename;
       exposureTimeVector[i] = exposureTime;
  }
  fclose(fp);

  return exposureTimeVector;
}


//computes the exposure time for image mosaicing (no reflectance model)
void ComputeExposure(std::string curr_input_file,  
                     std::string curr_albedo_file,
                     modelParams *currModelParams,
                     GlobalParams globalParams)
{

  
    DiskImageView<PixelMask<PixelGray<uint8> > > curr_image(curr_input_file);
    DiskImageView<PixelMask<PixelGray<uint8> > > curr_albedo(curr_albedo_file);

    GeoReference curr_geo;
    read_georeference(curr_geo, curr_input_file);

    float currReflectance;

    printf("init exposure time = %f, file = %s\n", currModelParams->exposureTime, curr_input_file.c_str());

    float delta_nominator = 0.0;
    float delta_denominator = 0.0;
    
    for (unsigned k=0; k < curr_image.rows(); ++k) {
      for (unsigned l=0; l < curr_image.cols(); ++l) {
	 
         Vector2 curr_sample_pix(l,k);
         
         if ( is_valid(curr_image(l,k)) ) {
          
	     currReflectance = 1;
             float error = ComputeError_Exposure((float)curr_image(l,k), currModelParams->exposureTime, 
                                                 (float)curr_albedo(l,k), currReflectance);
             float gradient = ComputeGradient_Exposure( (float)curr_albedo(l,k), currReflectance); 
	       
             delta_nominator = delta_nominator + error*gradient;
             delta_denominator = delta_denominator + gradient*gradient;  
             
	 }
       }
    }
     
    float delta = delta_nominator/delta_denominator;
    currModelParams->exposureTime = currModelParams->exposureTime+delta;
    printf("updated exposure time = %f\n", currModelParams->exposureTime);
 
}


//computes the exposure time for albedo mosaicing (uses a reflectance model)
void ComputeExposure(std::string curr_input_file,  
                     std::string curr_albedo_file,
                     std::string DEM_file,
                     modelParams *currModelParams,
                     GlobalParams globalParams)
{

  
    DiskImageView<PixelMask<PixelGray<uint8> > > curr_image(curr_input_file);
    DiskImageView<PixelMask<PixelGray<uint8> > > curr_albedo(curr_albedo_file);

    GeoReference curr_geo;
    read_georeference(curr_geo, curr_input_file);

    float currReflectance;

    //read the DEM file   
    DiskImageView<PixelGray<float> >  dem_image(DEM_file); 
    GeoReference curr_dem_geo;
    read_georeference(curr_dem_geo, DEM_file);

    printf("init exposure time = %f, file = %s\n", currModelParams->exposureTime, curr_input_file.c_str());

    float delta_nominator = 0.0;
    float delta_denominator = 0.0;
    
    for (unsigned k=0; k < curr_image.rows(); ++k) {
      for (unsigned l=0; l < curr_image.cols(); ++l) {
	 
         Vector2 curr_sample_pix(l,k);
         
         if ( is_valid(curr_image(l,k)) ) {
          
	   Vector2 lon_lat = curr_geo.pixel_to_lonlat(curr_sample_pix);
           Vector2 sample_pix_dem = curr_dem_geo.lonlat_to_pixel(curr_geo.pixel_to_lonlat(curr_sample_pix));

           int x = (int)sample_pix_dem[0];
	   int y = (int)sample_pix_dem[1];
         
           //check for valid DEM coordinates 
           if ((x>=0) && (x < dem_image.cols()) && (y>=0) && (y< dem_image.rows())){
                       
	     Vector3 longlat3(lon_lat(0),lon_lat(1),(dem_image)(x, y));
	     Vector3 xyz = curr_geo.datum().geodetic_to_cartesian(longlat3);
       
	     Vector2 sample_pix_dem_left;
	     sample_pix_dem_left(0) = x-1;
	     sample_pix_dem_left(1) = y;
	  
	     Vector2 sample_pix_dem_top;
	     sample_pix_dem_top(0) = x;
	     sample_pix_dem_top(1) = y-1;
	   

             //check for valid DEM pixel value and valid left and top coordinates
	     if ((sample_pix_dem_left(0) >= 0) && (sample_pix_dem_top(1) >= 0) && (dem_image(x,y) != -10000)){
 
	       Vector2 lon_lat_left = curr_dem_geo.pixel_to_lonlat(sample_pix_dem_left);
	       Vector3 longlat3_left(lon_lat_left(0),lon_lat_left(1),(dem_image)(sample_pix_dem_left(0), sample_pix_dem_left(1)));
	       Vector3 xyz_left = curr_geo.datum().geodetic_to_cartesian(longlat3_left);
           
	       Vector2 lon_lat_top= curr_dem_geo.pixel_to_lonlat(sample_pix_dem_top);
	       Vector3 longlat3_top(lon_lat_top(0),lon_lat_top(1),(dem_image)(sample_pix_dem_top(0), sample_pix_dem_top(1)));
	       Vector3 xyz_top = curr_geo.datum().geodetic_to_cartesian(longlat3_top);
                
	       //Vector3 normal = computeNormalFrom3DPoints(xyz, xyz_left, xyz_top);
               Vector3 normal = computeNormalFrom3DPointsGeneral(xyz, xyz_left, xyz_top);

               currReflectance = ComputeReflectance(normal, xyz, *currModelParams, globalParams);
               float error = ComputeError_Exposure((float)curr_image(l,k), currModelParams->exposureTime, 
                                                   (float)curr_albedo(l,k), currReflectance, xyz, xyz);

               float gradient = ComputeGradient_Exposure(currModelParams->exposureTime, (float)curr_albedo(l,k)); 
            
               if (globalParams.useWeights == 0){
                   delta_nominator = delta_nominator + error*gradient;
                   delta_denominator = delta_denominator + gradient*gradient;  
	
               }
               else{
		   //float weight = ComputeWeights(curr_sample_pix, currModelParams->center2D, currModelParams->maxDistance);
                   float weight = ComputeLineWeights(curr_sample_pix, currModelParams->centerLine, currModelParams->maxDistArray);
		   delta_nominator = delta_nominator + error*gradient*weight;
                   delta_denominator = delta_denominator + gradient*gradient*weight;  
               }
	       
               //delta_nominator = delta_nominator + error*gradient;
               //delta_denominator = delta_denominator + gradient*gradient;  
          
	     }
	   }
	 }
       }
    }
     
    float delta = delta_nominator/delta_denominator;
    printf("delta = %f\n", delta);
    currModelParams->exposureTime = currModelParams->exposureTime+delta;
    printf("updated exposure time = %f\n", currModelParams->exposureTime);
 
}



//--------------------------------------------------------------------------------------------------------------------------------------------------------
// Functions by Taemin Kim
// Written by Taemin Kim
void save_exposure_times(char * output_file, Vector<float> exposure_times) {
	std::cout << "Writing exposure times to " << output_file << std::endl;
	
	float buffer[SIZE_OF_BUFFER];
	for (unsigned i = 0; i < exposure_times.size(); ++i) buffer[i] = exposure_times(i);
	save_binary_file(buffer,exposure_times.size(),output_file);	
}

Vector<float> save_exposure_times(char * output_file, std::vector<std::string> radiance_files, std::vector<std::string> response_files, 
								  std::vector<std::string> index_files, char * weight_file) {
	Vector<float> exposure_times(response_files.size());
	Vector<uint8> inverse_weight = load_inverse_weight(weight_file, 256);
	
	std::cout << std::endl << "Estimating exposure times." << std::endl;
	for (unsigned i = 0; i < response_files.size(); ++i) {
		DiskImageView<PixelMask<PixelGray<uint8> > > index(index_files[i]);
		DiskImageView<PixelMask<PixelGray<float> > > tm_response(response_files[i]);
		DiskImageView<PixelMask<PixelGray<float> > > tm_radiance(radiance_files[i]);
		float sum_radiance = 0, sum_response = 0;
		for (unsigned x=0; x<tm_response.cols(); ++x)
			for (unsigned y=0; y<tm_response.rows(); ++y)
				if ( is_valid(tm_response(x,y)) ) {
					sum_response += tm_response(x,y)/inverse_weight(index(x,y));
					sum_radiance += tm_radiance(x,y)/inverse_weight(index(x,y));
				}
		exposure_times[i] = sum_response/sum_radiance;
		std::cout << i << ": " << exposure_times[i] << ", s " << sum_radiance << " ";
	}
	std::cout << std::endl;
	
	float buffer[SIZE_OF_BUFFER];
	for (unsigned i = 0; i < response_files.size(); ++i) buffer[i] = exposure_times(i);
	save_binary_file(buffer,response_files.size(),output_file);	
	
	return exposure_times;
}

// Written by Taemin Kim
// Preserve output_files if output_files is modified later than input_files and exposure_times.
Vector<float> save_exposure_times(char * output_file, std::vector<std::string> radiance_files, std::vector<std::string> response_files, 
				  std::vector<std::string> index_files, std::vector<std::string> weight_files, char * weight_file) {
	
	Vector<float> exposure_times(response_files.size());
	Vector<uint8> inverse_weight = load_inverse_weight(weight_file, 256);
	
	std::cout << std::endl << "Estimating exposure times." << std::endl;
	for (unsigned i = 0; i < response_files.size(); ++i) {
		DiskImageView<PixelMask<PixelGray<uint8> > > index(index_files[i]);
		DiskImageView<PixelMask<PixelGray<float> > > response(response_files[i]);
		DiskImageView<PixelMask<PixelGray<float> > > radiance(radiance_files[i]);
		DiskImageView<PixelMask<PixelGray<float> > > weight(weight_files[i]);
		ImageView<PixelMask<PixelGray<float> > > tm_response = weight*response;
		ImageView<PixelMask<PixelGray<float> > > tm_radiance = weight*radiance;
		float sum_radiance = 0, sum_response = 0;
		for (unsigned x=0; x<tm_response.cols(); ++x)
			for (unsigned y=0; y<tm_response.rows(); ++y)
				if ( is_valid(tm_response(x,y)) ) {
					sum_response += tm_response(x,y)/inverse_weight(index(x,y));
					sum_radiance += tm_radiance(x,y)/inverse_weight(index(x,y));
				}
		exposure_times[i] = sum_response/sum_radiance;
		std::cout << i << ": " << exposure_times[i] << ", s " << sum_radiance << " ";
	}
	std::cout << std::endl;
	
	float buffer[SIZE_OF_BUFFER];
	for (unsigned i = 0; i < response_files.size(); ++i) buffer[i] = exposure_times(i);
	save_binary_file(buffer,response_files.size(),output_file);	
	
	return exposure_times;
}

Vector<float> load_exposure_times(char * input_file, int num) {
	Vector<float> exposure_times(num);
	float buffer[SIZE_OF_BUFFER];
	FILE *fp;
	
	struct stat file_stat;
	if ( stat(input_file, &file_stat) ) { 
		std::cout << "Writing uniform exposure times to  " << input_file << std::endl;
		for (unsigned i = 0; i < num; ++i) {
			exposure_times(i) = 1;
			buffer[i] = 1;
		}
		save_binary_file(buffer,num,input_file);	
	} else {
		std::cout << "Reading exposure times from " << input_file << std::endl;
		load_binary_file(buffer, num, input_file);	
		for (unsigned i = 0; i < num; ++i) exposure_times(i) = buffer[i];
		std::cout << exposure_times << std::endl;
	}
	
	return exposure_times;
}

// Preserve output_files if output_files is modified later than input_files and exposure_times.
float normalize_exposures(std::vector<std::string> input_files, char * exp_time_file) {
	float max_value = 0, scaling_factor;
	float lo, hi;
	struct stat input_stat, output_stat, exp_time_stat;
	for (int i = 0; i < input_files.size(); ++i) {
		DiskImageView<PixelMask<PixelGray<float> > > image(input_files[i]);
		std::cout << " Checking " << i << "th radiance image: " << input_files[i] << "." << std::endl;
		min_max_channel_values(image, lo, hi);
		std::cout << " Max: " << hi << ", " << "Min: " << lo << "." << std::endl;
		if (max_value < hi) max_value = hi;
	}
	
	scaling_factor = (DYNAMIC_RANGE-1)/max_value;
	for (int i = 0; i < input_files.size(); ++i) {
		GeoReference geo;
		read_georeference(geo, input_files[i]);
		std::cout << " Scaling " << i << "th radiance image: " << input_files[i] << " by " << scaling_factor << "." << std::endl;
		DiskImageView<PixelMask<PixelGray<float> > > image(input_files[i]);
		ImageView<PixelMask<PixelGray<float> > > tm_image = scaling_factor*image;
		write_georeferenced_image(input_files[i], tm_image, geo, TerminalProgressCallback("{Core}","Processing:"));
 	}
	Vector<float> exposure_times = load_exposure_times(exp_time_file, input_files.size());
	exposure_times /= scaling_factor;
	save_exposure_times(exp_time_file, exposure_times);	
	
	return scaling_factor;
}


