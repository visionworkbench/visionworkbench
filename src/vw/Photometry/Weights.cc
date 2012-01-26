// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <iostream>
#include <string>

#include <vector>
#include <sys/types.h>
#include <sys/stat.h>

#include <vw/Core.h>
#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Cartography.h>
#include <vw/Math.h>
using namespace vw;
using namespace vw::cartography;

#include <vw/Photometry/ShapeFromShading.h>
#include <vw/Photometry/Shape.h>
#include <vw/Photometry/Albedo.h>
#include <vw/Photometry/Exposure.h>
#include <vw/Photometry/Reflectance.h>
#include <vw/Photometry/Reconstruct.h>
#include <vw/Photometry/Shadow.h>
#include <vw/Photometry/Index.h>
#include <vw/Photometry/Weights.h>
using namespace vw::photometry;


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

int*
vw::photometry::ComputeImageHCenterLine(std::string input_img_file,
                                       int **r_hMaxDistArray) {

    //compute the center of the image
    DiskImageView<PixelMask<PixelGray<uint8> > >  input_img(input_img_file);

    int k, l;

    int *hCenterLine = new int [input_img.rows()];
    int* hMaxDistArray = new int [input_img.rows()];

    int minVal, maxVal;

    printf("file=%s\n",input_img_file.c_str());

    //initialize  output_img, and numSamples
    for (k = 0 ; k < (int)input_img.rows(); ++k) {

        minVal = input_img.cols();
        maxVal = 0;
        for (l = 0; l < (int)input_img.cols(); ++l) {


           if ( is_valid(input_img(l,k)) ) {

             if (l < minVal){
                 minVal = l;
             }
             if (l > maxVal){
                 maxVal = l;
             }
          }
      }
      hCenterLine[k] = (minVal + maxVal)/2;
      hMaxDistArray[k] = maxVal - minVal;
      if (hMaxDistArray[k] < 0){
          hMaxDistArray[k]=0;
      }
      //printf("cl[%d] = %d, maxDist[%d] = %d\n", k, hCenterLine[k], k, hMaxDistArray[k]);
   }



   *r_hMaxDistArray = hMaxDistArray;
   return hCenterLine;
}


int*
vw::photometry::ComputeImageVCenterLine(std::string input_img_file,
                                       int **r_vMaxDistArray) {

    //compute the center of the image
    DiskImageView<PixelMask<PixelGray<uint8> > >  input_img(input_img_file);

    int k, l;

    int *vCenterLine = new int [input_img.cols()];
    int* vMaxDistArray = new int [input_img.cols()];

    int minVal, maxVal;

    printf("file=%s\n",input_img_file.c_str());

    //initialize  output_img, and numSamples
    for (k = 0 ; k < (int)input_img.cols(); ++k) {

        minVal = input_img.rows();
        maxVal = 0;
        for (l = 0; l < (int)input_img.rows(); ++l) {

          if ( is_valid(input_img(k, l)) ) {

             if (l < minVal){
                 minVal = l;
             }
             if (l > maxVal){
                 maxVal = l;
             }
          }
      }
      vCenterLine[k] = (minVal + maxVal)/2;
      vMaxDistArray[k] = maxVal - minVal;
      if (vMaxDistArray[k] < 0){
          vMaxDistArray[k]=0;
      }
      //printf("cl[%d] = %d, maxDist[%d] = %d\n", k, vCenterLine[k], k, vMaxDistArray[k]);
   }

   *r_vMaxDistArray = vMaxDistArray;
   return vCenterLine;
}


int*
vw::photometry::ComputeDEMHCenterLine(std::string input_DEM_file,int noDataDEMVal,
                                     int **r_hMaxDistArray) {

    //compute the center of the image
    DiskImageView<PixelGray<float> >   input_DEM(input_DEM_file);

    int k, l;

    int *hCenterLine = new int[input_DEM.rows()];
    int *hMaxDistArray = new int[input_DEM.rows()];

    int minVal, maxVal;

    printf("file=%s\n",input_DEM_file.c_str());
    //initialize  output_img, and numSamples
    for (k = 0 ; k < input_DEM.rows(); ++k) {

        minVal = input_DEM.cols();
        maxVal = 0;
        for (l = 0; l < input_DEM.cols(); ++l) {

	  if ( input_DEM(l,k) != noDataDEMVal ) {
             if (l < minVal){
                 minVal = l;
             }
             if (l > maxVal){
                 maxVal = l;
             }
          }
      }
      hCenterLine[k] = (minVal + maxVal)/2;
      hMaxDistArray[k] = maxVal - minVal;
      if (hMaxDistArray[k] < 0){
        hMaxDistArray[k]=0;
      }
      //printf("cl[%d] = %d, maxDist[%d] = %d\n", k, hCenterLine[k], k, hMaxDistArray[k]);
   }



   *r_hMaxDistArray = hMaxDistArray;
   return hCenterLine;
}

int*
vw::photometry::ComputeDEMVCenterLine(std::string input_DEM_file,int noDataDEMVal,
                                     int **r_vMaxDistArray) {

    //compute the center of the image
    DiskImageView<PixelGray<float> >   input_DEM(input_DEM_file);

    int k, l;

    int *vCenterLine = new int[input_DEM.cols()];
    int *vMaxDistArray = new int[input_DEM.cols()];

    int minVal, maxVal;

    printf("file=%s\n",input_DEM_file.c_str());
    //initialize  output_img, and numSamples
    for (k = 0 ; k < input_DEM.cols(); ++k) {

        minVal = input_DEM.rows();
        maxVal = 0;
        for (l = 0; l < input_DEM.rows(); ++l) {

	  if ( input_DEM(k,l) != noDataDEMVal ) {
             if (l < minVal){
                 minVal = l;
             }
             if (l > maxVal){
                 maxVal = l;
             }
          }
      }
      vCenterLine[k] = (minVal + maxVal)/2;
      vMaxDistArray[k] = maxVal - minVal;
      if (vMaxDistArray[k] < 0){
        vMaxDistArray[k]=0;
      }
      //printf("cl[%d] = %d, maxDist[%d] = %d\n", k, vCenterLine[k], k, vMaxDistArray[k]);
   }



   *r_vMaxDistArray = vMaxDistArray;
   return vCenterLine;
}

float
vw::photometry::ComputeLineWeightsH(Vector2 pix,
                                    int *hCenterLine, int *hMaxDistArray){
  int maxDist = hMaxDistArray[(int)pix[1]]/2.0;
  int center = hCenterLine[(int)pix[1]];
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

float
vw::photometry::ComputeLineWeightsV(Vector2 pix,
                                    int *vCenterLine, int *vMaxDistArray) {
  int maxDist = vMaxDistArray[(int)pix[0]]/2.0;
  int center = vCenterLine[(int)pix[0]];
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

float
vw::photometry::ComputeLineWeightsHV(Vector2 pix,
                                     int *hCenterLine, int *hMaxDistArray,
                                     int *vCenterLine, int *vMaxDistArray
                                     )
{
  float weightH = ComputeLineWeightsH(pix, hCenterLine, hMaxDistArray);
  float weightV = ComputeLineWeightsV(pix, vCenterLine, vMaxDistArray);
  float weight = weightH*weightV;
  return weight;
}


// Saves weights to file.
void 
vw::photometry::SaveWeightsParamsToFile(struct ModelParams modelParams)
{

  FILE *fp;
  fp = fopen((char*)(modelParams.weightFilename).c_str(), "w");
  DiskImageView<PixelGray<float> >   dem(modelParams.DEMFilename);
  DiskImageView<PixelMask<PixelGray<uint8> > >  drg(modelParams.inputFilename);

  // Horizontal
  for (int i = 0; i < (int)dem.rows(); i++){
      fprintf(fp, "%d ", modelParams.hCenterLineDEM[i]);
  }
  fprintf(fp, "\n");

  for (int i = 0; i < (int)dem.rows(); i++){
      fprintf(fp, "%d ", modelParams.hMaxDistArrayDEM[i]);
  }
  fprintf (fp, "\n");

  for (int i = 0; i < (int)drg.rows(); i++){
      fprintf(fp, "%d ", modelParams.hCenterLine[i]);
  }
  fprintf(fp, "\n");

  for (int i = 0; i < (int)drg.rows(); i++){
      fprintf(fp, "%d ", modelParams.hMaxDistArray[i]);
  }
  fprintf(fp,"\n");
  
  // Vertical
  
  for (int i = 0; i < (int)dem.cols(); i++){
      fprintf(fp, "%d ", modelParams.vCenterLineDEM[i]);
  }
  fprintf(fp, "\n");

  for (int i = 0; i < (int)dem.cols(); i++){
      fprintf(fp, "%d ", modelParams.vMaxDistArrayDEM[i]);
  }
  fprintf (fp, "\n");

  for (int i = 0; i < (int)drg.cols(); i++){
      fprintf(fp, "%d ", modelParams.vCenterLine[i]);
  }
  fprintf(fp, "\n");

  for (int i = 0; i < (int)drg.cols(); i++){
      fprintf(fp, "%d ", modelParams.vMaxDistArray[i]);
  }
  fprintf(fp,"\n");
  
  fclose(fp);
}

void 
vw::photometry::ReadWeightsParamsFromFile(struct ModelParams *modelParams)
{
  FILE *fp;

  fp = fopen((char*)(modelParams->weightFilename).c_str(), "r");
  if (fp==NULL){
    return;
  }

  DiskImageView<PixelGray<float> >   dem(modelParams->DEMFilename);
  DiskImageView<PixelMask<PixelGray<uint8> > >  drg(modelParams->inputFilename);

  // Horizontal
  
  modelParams->hCenterLineDEM = new int[dem.rows()];
  modelParams->hMaxDistArrayDEM = new int[dem.rows()];
  for (int i = 0; i < (int)dem.rows(); i++){
    fscanf(fp, "%d ", &(modelParams->hCenterLineDEM[i]));
  }
  fscanf(fp, "\n");
  
  for (int i = 0; i < (int)dem.rows(); i++){
    fscanf(fp, "%d ", &(modelParams->hMaxDistArrayDEM[i]));
  }
  fscanf (fp, "\n");
  
  modelParams->hCenterLine = new int[drg.rows()];
  modelParams->hMaxDistArray = new int[drg.rows()];
   
  for (int i = 0; i < (int)drg.rows(); i++){
    fscanf(fp, "%d ", &(modelParams->hCenterLine[i]));
  }
  fscanf(fp, "\n");

  for (int i = 0; i < (int)drg.rows(); i++){
    fscanf(fp, "%d ", &(modelParams->hMaxDistArray[i]));
  }
  fscanf(fp,"\n");

  // Vertical
  
  modelParams->vCenterLineDEM = new int[dem.cols()];
  modelParams->vMaxDistArrayDEM = new int[dem.cols()];
  for (int i = 0; i < (int)dem.cols(); i++){
    fscanf(fp, "%d ", &(modelParams->vCenterLineDEM[i]));
  }
  fscanf(fp, "\n");
  
  for (int i = 0; i < (int)dem.cols(); i++){
    fscanf(fp, "%d ", &(modelParams->vMaxDistArrayDEM[i]));
  }
  fscanf (fp, "\n");
  
  modelParams->vCenterLine = new int[drg.cols()];
  modelParams->vMaxDistArray = new int[drg.cols()];
   
  for (int i = 0; i < (int)drg.cols(); i++){
    fscanf(fp, "%d ", &(modelParams->vCenterLine[i]));
  }
  fscanf(fp, "\n");

  for (int i = 0; i < (int)drg.cols(); i++){
    fscanf(fp, "%d ", &(modelParams->vMaxDistArray[i]));
  }
  fscanf(fp,"\n");
  
  fclose(fp);
  
}
