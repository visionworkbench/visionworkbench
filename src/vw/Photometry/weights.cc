// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <string>
#include <fstream>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>

#include <boost/operators.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/fstream.hpp>
namespace fs = boost::filesystem;

#include <vw/Core.h>
#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Cartography.h>
#include <vw/Math.h>
using namespace vw;
using namespace vw::math;
using namespace vw::cartography;

#include "shape_from_shading.h"
#include "shape.h"
#include "albedo.h"
#include "exposure.h"
#include "reflectance.h"
#include "reconstruct.h"
#include "shadow.h"
#include "index.h"
#include <math.h>
//#include <cv.h>
//#include <highgui.h>

#include "reconstruct.h"


/*
Vector2 ComputeImageCenter2D(std::string input_img_file, float *maxDist)
{
 
    //compute the center of the image
    DiskImageView<PixelMask<PixelGray<uint8> > >  input_img(input_img_file); 
    GeoReference input_geo;
    read_georeference(input_geo, input_img_file);
 
    int k, l;
    int numSamples = 0;

    Vector2 C;
    C(0) = 0;
    C(1) = 0;
  
    //initialize  output_img, and numSamples
    for (k = 0 ; k < input_img.rows(); ++k) {
        for (l = 0; l < input_img.cols(); ++l) {
           
           Vector2 input_image_pix(l,k);
           
           if ( is_valid(input_img(l,k)) ) {
           
              C(0) = C(0) + l; 
              C(1) = C(1) + k; 
            
              numSamples++;
	  }
      }
   }
   C(0) = C(0)/numSamples;
   C(1) = C(1)/numSamples;
 

   float l_maxDist = 0;
    for (k = 0 ; k < input_img.rows(); ++k) {
        for (l = 0; l < input_img.cols(); ++l) {
           
           Vector2 input_image_pix(l,k);
           
           if ( is_valid(input_img(l,k)) ) {
              
              Vector2 dist_t;
              dist_t[0] = l-C(0); 
              dist_t[1] = k-C(1); 
          
              
              float dist = (dist_t[0]*dist_t[0]) + (dist_t[0]*dist_t[0]) ;
              if (dist > l_maxDist){
		 l_maxDist = dist;
	      }
	  }
      }
   }

    l_maxDist = sqrt(l_maxDist);
   
    printf("file=%s, i = %f, j = %f, maxDist = %f\n",input_img_file.c_str(), C(0), C(1), l_maxDist);
   
    *maxDist = l_maxDist;
    return C;
}
float ComputeWeights(Vector2 pix, Vector2 C, float maxDistance)
{

  float a = -1.0/maxDistance;
  float b = 1;

  float weight;
  Vector2 dist;

  dist[0]= pix[0]-C[0];
  dist[1]= pix[1]-C[1];

  float t_dist = sqrt(dist[0]*dist[0] + dist[1]*dist[1]);
  weight = a*t_dist + b; 

  return weight;
}
*/

int* ComputeImageCenterLine(std::string input_img_file, int **r_maxDistArray)
{
 
    //compute the center of the image
    DiskImageView<PixelMask<PixelGray<uint8> > >  input_img(input_img_file); 
 
    int k, l;
  
    int *centerLine = new int[input_img.rows()];  
    int *maxDistArray = new int[input_img.rows()];  

    int minVal, maxVal;

    printf("file=%s\n",input_img_file.c_str());
    //initialize  output_img, and numSamples
    for (k = 0 ; k < input_img.rows(); ++k) {
        
        minVal = input_img.cols();
        maxVal = 0;
        for (l = 0; l < input_img.cols(); ++l) {
           
           Vector2 input_image_pix(l,k);
          
     
           if ( is_valid(input_img(l,k)) ) {
           
	     if (l < minVal){
	         minVal = l;
             }
             if (l > maxVal){
	         maxVal = l; 
             }  
	  }
      }
      centerLine[k] = (minVal + maxVal)/2;   
      maxDistArray[k] = maxVal - minVal;
      if (maxDistArray[k] < 0){
	  maxDistArray[k]=0;
      }
      //printf("cl[%d] = %d, maxDist[%d] = %d\n", k, centerLine[k], k, maxDistArray[k]); 
   }
  
  
 
   *r_maxDistArray = maxDistArray; 
   return centerLine;
}


int* ComputeImageHorCenterLine(std::string input_img_file, int **r_maxDistArray)
{
 
    //compute the center of the image
    DiskImageView<PixelMask<PixelGray<uint8> > >  input_img(input_img_file); 
 
    int k, l;
  
    int *centerLine = new int[input_img.cols()];  
    int *maxDistArray = new int[input_img.cols()];  

    int minVal, maxVal;

    printf("file=%s\n",input_img_file.c_str());
    //initialize  output_img, and numSamples
    for (k = 0 ; k < input_img.cols(); ++k) {
        
        minVal = input_img.rows();
        maxVal = 0;
        for (l = 0; l < input_img.rows(); ++l) {
           
           Vector2 input_image_pix(l,k);
          
     
           if ( is_valid(input_img(l,k)) ) {
           
	     if (l < minVal){
	         minVal = l;
             }
             if (l > maxVal){
	         maxVal = l; 
             }  
	  }
      }
      centerLine[k] = (minVal + maxVal)/2;   
      maxDistArray[k] = maxVal - minVal;
      if (maxDistArray[k] < 0){
	  maxDistArray[k]=0;
      }
      //printf("cl[%d] = %d, maxDist[%d] = %d\n", k, centerLine[k], k, maxDistArray[k]); 
   }
 
   *r_maxDistArray = maxDistArray; 
    return centerLine;
}


int* ComputeDEMCenterLine(std::string input_DEM_file, int **r_maxDistArray)
{
 
    //compute the center of the image
    DiskImageView<PixelGray<float> >   input_DEM(input_DEM_file); 
 
    int k, l;
  
    int *centerLine = new int[input_DEM.rows()];  
    int *maxDistArray = new int[input_DEM.rows()];  

    int minVal, maxVal;

    printf("file=%s\n",input_DEM_file.c_str());
    //initialize  output_img, and numSamples
    for (k = 0 ; k < input_DEM.rows(); ++k) {
        
        minVal = input_DEM.cols();
        maxVal = 0;
        for (l = 0; l < input_DEM.cols(); ++l) {
           
           if ( input_DEM(l,k) != -10000 ) {
	     if (l < minVal){
	         minVal = l;
             }
             if (l > maxVal){
	         maxVal = l; 
             }  
	  }
      }
      centerLine[k] = (minVal + maxVal)/2;   
      maxDistArray[k] = maxVal - minVal;
      if (maxDistArray[k] < 0){
	  maxDistArray[k]=0;
      }
      //printf("cl[%d] = %d, maxDist[%d] = %d\n", k, centerLine[k], k, maxDistArray[k]); 
   }
  
  
 
   *r_maxDistArray = maxDistArray; 
   return centerLine;
}

int* ComputeDEMHorCenterLine(std::string input_DEM_file, int **r_maxDistArray)
{
 
    //compute the center of the image
    DiskImageView<PixelGray<float> >   input_DEM(input_DEM_file); 
 
    int k, l;
  
    int *centerLine = new int[input_DEM.cols()];  
    int *maxDistArray = new int[input_DEM.cols()];  

    int minVal, maxVal;

    printf("file=%s\n",input_DEM_file.c_str());
    //initialize  output_img, and numSamples
    for (k = 0 ; k < input_DEM.cols(); ++k) {
        
        minVal = input_DEM.rows();
        maxVal = 0;
        for (l = 0; l < input_DEM.rows(); ++l) {
           
           if ( input_DEM(l,k) != -10000 ) {
	     if (l < minVal){
	         minVal = l;
             }
             if (l > maxVal){
	         maxVal = l; 
             }  
	  }
      }
      centerLine[k] = (minVal + maxVal)/2;   
      maxDistArray[k] = maxVal - minVal;
      if (maxDistArray[k] < 0){
	  maxDistArray[k]=0;
      }
      printf("cl[%d] = %d, maxDist[%d] = %d\n", k, centerLine[k], k, maxDistArray[k]); 
   }
  
  
 
   *r_maxDistArray = maxDistArray; 
   return centerLine;
}


float ComputeLineWeights(Vector2 pix, int *centerLine, int *maxDistArray)
{
  int maxDist = maxDistArray[(int)pix[1]]/2.0;
  int center = centerLine[(int)pix[1]];
  float dist = fabs((int)pix[0]-center);
  float a;
  float b; 
  
  if (maxDist == 0){
      maxDist = 1;
  }

  a = -1.0/(float)maxDist;
  b = 1;
  float weight = a*dist + b;
  //printf("i = %d, j = %d, weight = %f\n", (int)pix[1], (int)pix[0], weight); 
  return weight;
}

float ComputeLineWeightsV(Vector2 pix, int *centerLine, int *maxDistArray)
{
  int maxDist = maxDistArray[(int)pix[0]]/2.0;
  int center = centerLine[(int)pix[0]];
  float dist = fabs((int)pix[1]-center);
  float a;
  float b; 
  
  if (maxDist == 0){
      maxDist = 1;
  }

  a = -1.0/(float)maxDist;
  b = 1;
  float weight = a*dist + b;
  //printf("i = %d, j = %d, weight = %f\n", (int)pix[1], (int)pix[0], weight); 
  return weight;
}
float ComputeWeightsVH(Vector2 pix, modelParams imgParams)
{
  float weightV = ComputeLineWeightsV(pix, imgParams.horCenterLine, imgParams.maxVerDistArray);
  float weightH = ComputeLineWeights(pix, imgParams.centerLine, imgParams.maxDistArray);
  float weight = weightV*weightH;
  return weight;
}
