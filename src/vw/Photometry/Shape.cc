// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <iostream>

#include <math.h>
#include <time.h>

#include <vw/Core.h>
#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Cartography.h>
#include <vw/Math.h>
using namespace vw;
using namespace vw::cartography;

#include <vw/Photometry/Shape.h>
#include <vw/Photometry/Reconstruct.h>
#include <vw/Photometry/Weights.h>
using namespace vw::photometry;

float ComputeGradient_DEM(float /*intensity*/, float T,
                          float albedo, Vector3 s,
                          Vector3 p, Vector3 p_left,
                          Vector3 p_top, Vector3 /*xyz_prior*/) {
  float grad;
  float temp  = (p[2]-p_left[2])*(p[2]-p_left[2]) + (p[2]-p_top[2])*(p[2]-p_top[2]) + 1;
  float temp1 = (p[2]-p_left[2])*s[0]+(p[2]-p_top[2])*s[1] +s[2];
  float temp2 = (p[2]-p_left[2])+(p[2]-p_top[2]);
  grad = T*albedo*((s[0]+s[1])*sqrt(temp)-(temp1*temp2/sqrt(temp)))/temp;

  return grad;
}

float ComputeError_DEM(float intensity, float T, float albedo,
                       float reflectance, Vector3 /*xyz*/, Vector3 xyz_prior) {
  float error;
  error = (intensity-T*albedo*reflectance) + (xyz_prior[2]-xyz_prior[2]);
  return error;
}

//initializes the DEM file by getting the average DEM values in the overlapping areas of consecutive DEM files.
void vw::photometry::InitDEM( ModelParams input_img_params,
                              std::vector<ModelParams> overlap_img_params,
                              GlobalParams globalParams) {

  int i;
  unsigned l, k;

  std::string input_DEM_file = input_img_params.DEMFilename;
  std::string mean_DEM_file = input_img_params.meanDEMFilename;
  std::string var2_DEM_file = input_img_params.var2DEMFilename;

  DiskImageView<PixelGray<float> >  input_DEM_image(input_DEM_file);
  GeoReference input_DEM_geo;
  read_georeference(input_DEM_geo, input_DEM_file);

  ImageView<PixelGray<float> > mean_DEM_image(input_DEM_image.cols(), input_DEM_image.rows());
  ImageView<PixelMask<PixelGray<float> > >var2_DEM_image(input_DEM_image.cols(), input_DEM_image.rows());

  ImageView<PixelGray<int> > numSamples(input_DEM_image.cols(), input_DEM_image.rows());
  ImageView<PixelGray<float> > norm(input_DEM_image.cols(), input_DEM_image.rows());

  printf("corrections to var2\n");
  //initialize  mean_DEM-image, var2_DEM_image and numSamples and norm
  for (k = 0 ; k < (unsigned)input_DEM_image.rows(); ++k) {
    for (l = 0; l < (unsigned)input_DEM_image.cols(); ++l) {

      mean_DEM_image(l, k) = -10000;
      var2_DEM_image(l, k) = 0;
      numSamples(l, k) = 0;
      norm(l,k) = 0;

      Vector2 input_DEM_pix(l,k);

      if ( input_DEM_image(l,k) != -10000 ) {

        if (globalParams.useWeights == 0){
          mean_DEM_image(l, k) = (float)input_DEM_image(l,k);
          var2_DEM_image(l, k) = (float)input_DEM_image(l,k)*(float)input_DEM_image(l,k);
          numSamples(l, k) = 1;
        }
        else{
          float weight = ComputeLineWeights(input_DEM_pix, input_img_params.centerLineDEM, input_img_params.maxDistArrayDEM);
          mean_DEM_image(l, k) = (float)input_DEM_image(l,k)*weight;
          //weight added by Ara 08/28
          var2_DEM_image(l, k) = (float)input_DEM_image(l,k)*(float)input_DEM_image(l,k)*weight;
          //numSamples(l, k) = 1;
          norm(l, k) = weight;
        }
      }
    }
  }

  for (i = 0; i < (int) overlap_img_params.size(); i++){

    printf("DEM = %s\n", overlap_img_params[i].DEMFilename.c_str());
    DiskImageView<PixelGray<float> >  overlap_DEM_image(overlap_img_params[i].DEMFilename);
    GeoReference overlap_DEM_geo;
    read_georeference(overlap_DEM_geo, overlap_img_params[i].DEMFilename);


    ImageViewRef<PixelGray<float> >  interp_overlap_DEM_image = interpolate(edge_extend(overlap_DEM_image.impl(),
                                                                                        ConstantEdgeExtension()),
                                                                            BilinearInterpolation());

    for (k = 0 ; k < (unsigned)input_DEM_image.rows(); ++k) {
      for (l = 0; l < (unsigned)input_DEM_image.cols(); ++l) {

        Vector2 input_DEM_pix(l,k);

        if ( input_DEM_image(l,k) != -10000 ) {

          //check for overlap between the output image and the input DEM image
          Vector2 overlap_dem_pix = overlap_DEM_geo.lonlat_to_pixel(input_DEM_geo.pixel_to_lonlat(input_DEM_pix));
          float x = overlap_dem_pix[0];
          float y = overlap_dem_pix[1];

          //check for valid DEM coordinates
          if ((x>=0) && (x < overlap_DEM_image.cols()) && (y>=0) && (y< overlap_DEM_image.rows())){

            if ( overlap_DEM_image(x, y) != -10000 ) {

              if (globalParams.useWeights == 0){
                mean_DEM_image(l, k) = (float)mean_DEM_image(l, k) + (float)interp_overlap_DEM_image(x, y);
                var2_DEM_image(l, k) = (float)var2_DEM_image(l, k) + (float)interp_overlap_DEM_image(x, y)*(float)interp_overlap_DEM_image(x, y);
                numSamples(l, k) = (int)numSamples(l, k) + 1;
              }
              else{
                float weight = ComputeLineWeights(overlap_dem_pix, overlap_img_params[i].centerLineDEM, overlap_img_params[i].maxDistArrayDEM);
                mean_DEM_image(l, k) = (float)mean_DEM_image(l, k) + (float)interp_overlap_DEM_image(x, y)*weight;
                //weight added by Ara 08/28/
                var2_DEM_image(l, k) = (float)var2_DEM_image(l, k) + (float)interp_overlap_DEM_image(x, y)*(float)interp_overlap_DEM_image(x, y)*weight;
                norm(l, k) = norm(l,k) + weight;
              }
            }
          }
        }
      }
    }
  }


  for (k = 0 ; k < (unsigned)input_DEM_image.rows(); ++k) {
    for (l = 0; l < (unsigned)input_DEM_image.cols(); ++l) {

      //compute variance only where the mean DEM is valid
      if ( input_DEM_image(l,k) != -10000 ) {

        if ((globalParams.useWeights == 0) && (numSamples(l,k)!=0)){
          mean_DEM_image(l, k) = mean_DEM_image(l, k)/numSamples(l,k);
          var2_DEM_image(l, k) = var2_DEM_image(l, k)/numSamples(l, k) - mean_DEM_image(l,k)*mean_DEM_image(l,k);
        }

        if ((globalParams.useWeights == 1) && (norm(l,k)!=0)){
          mean_DEM_image(l, k) = mean_DEM_image(l, k)/norm(l,k);
          var2_DEM_image(l, k) = var2_DEM_image(l, k)/norm(l, k) - mean_DEM_image(l,k)*mean_DEM_image(l,k);
          if (var2_DEM_image(l, k) < 0){ //this should never happen
              var2_DEM_image(l, k) = 0;
              var2_DEM_image(l, k).invalidate();
          }
        }
      }

    }
  }

  write_georeferenced_image(mean_DEM_file,
                            mean_DEM_image,
                            input_DEM_geo, TerminalProgressCallback("{Core}","Processing:"));

  write_georeferenced_image(var2_DEM_file,
                            var2_DEM_image,
                            input_DEM_geo, TerminalProgressCallback("{Core}","Processing:"));

}

/*
//initializes the DEM file by getting the average DEM values in the overlapping areas of consecutive DEM files.
void InitDEM( std::string input_DEM_file, std::string mean_DEM_file, std::string var2_DEM_file, ModelParams input_img_params,
              std::vector<std::string> overlap_DEM_files,  std::vector<ModelParams> overlap_img_params, GlobalParams globalParams)
{

    int i, l, k;

    DiskImageView<PixelGray<float> >  input_DEM_image(input_DEM_file);
    GeoReference input_DEM_geo;
    read_georeference(input_DEM_geo, input_DEM_file);

    ImageView<PixelGray<float> > mean_DEM_image(input_DEM_image.cols(), input_DEM_image.rows());
    ImageView<PixelMask<PixelGray<float> > >var2_DEM_image(input_DEM_image.cols(), input_DEM_image.rows());

    ImageView<PixelGray<int> > numSamples(input_DEM_image.cols(), input_DEM_image.rows());
    ImageView<PixelGray<float> > norm(input_DEM_image.cols(), input_DEM_image.rows());

    float avgStdDevDEM = 0.0;
    //int numTiles = overlap_img_params.size()+1;
    Vector<float, 5> meanDEMOffset;
    Vector<float, 5> numDEMSamples;

    //initialize  mean_DEM-image, var2_DEM_image and numSamples
    for (k = 0 ; k < (int)input_DEM_image.rows(); ++k) {
      for (l = 0; l < (int)input_DEM_image.cols(); ++l) {

           mean_DEM_image(l, k) = -10000;
           numSamples(l, k) = 0;
           Vector2 input_DEM_pix(l,k);

           if ( input_DEM_image(l,k) != -10000 ) {

              if (globalParams.useWeights == 0){
                  mean_DEM_image(l, k) = (float)input_DEM_image(l,k);
                  var2_DEM_image(l, k) = (float)input_DEM_image(l,k)*(float)input_DEM_image(l,k);
                  numSamples(l, k) = 1;

                  meanDEMOffset[0] = meanDEMOffset[0] + (float)input_DEM_image(l,k);
                  numDEMSamples[0] = numDEMSamples[0] + 1;
              }
              else{
                  float weight = ComputeLineWeights(input_DEM_pix, input_img_params.centerLineDEM, input_img_params.maxDistArrayDEM);
                  mean_DEM_image(l, k) = (float)input_DEM_image(l,k)*weight;
                  var2_DEM_image(l, k) = (float)input_DEM_image(l,k)*(float)input_DEM_image(l,k);
                  numSamples(l, k) = 1;
                  norm(l, k) = weight;

                  meanDEMOffset[0] = meanDEMOffset[0] + (float)input_DEM_image(l,k);
                  numDEMSamples[0] = numDEMSamples[0] + 1;
              }

           }

        }
    }


    for (i = 0; i < (int)overlap_DEM_files.size(); i++){

      printf("DEM = %s\n", overlap_DEM_files[i].c_str());

      DiskImageView<PixelGray<float> >  overlap_DEM_image(overlap_DEM_files[i]);
      GeoReference overlap_DEM_geo;
      read_georeference(overlap_DEM_geo, overlap_DEM_files[i]);


      ImageViewRef<PixelGray<float> >  interp_overlap_DEM_image = interpolate(edge_extend(overlap_DEM_image.impl(),
                                                                              ConstantEdgeExtension()),
                                                                              BilinearInterpolation());

      for (k = 0 ; k < input_DEM_image.rows(); ++k) {
        for (l = 0; l < input_DEM_image.cols(); ++l) {

          Vector2 input_DEM_pix(l,k);

          if ( input_DEM_image(l,k) != -10000 ) {

              //check for overlap between the output image and the input DEM image
              Vector2 overlap_dem_pix = overlap_DEM_geo.lonlat_to_pixel(input_DEM_geo.pixel_to_lonlat(input_DEM_pix));
              int x = (int)overlap_dem_pix[0];
              int y = (int)overlap_dem_pix[1];

              //check for valid DEM coordinates
              if ((x>=0) && (x < overlap_DEM_image.cols()) && (y>=0) && (y< overlap_DEM_image.rows())){

                if ( overlap_DEM_image(x, y) != -10000 ) {

                     if (globalParams.useWeights == 0){
                         mean_DEM_image(l, k) = (float)mean_DEM_image(l, k) + (float)interp_overlap_DEM_image(x, y);
                         var2_DEM_image(l, k) = (float)var2_DEM_image(l, k) + (float)interp_overlap_DEM_image(x, y)*(float)interp_overlap_DEM_image(x, y);
                         numSamples(l, k) = (int)numSamples(l, k) + 1;

                         meanDEMOffset[i+1] = meanDEMOffset[i+1] + (float)interp_overlap_DEM_image(x, y);
                         numDEMSamples[i+1] = numDEMSamples[i+1] + 1;
                     }
                     else{
                         float weight = ComputeLineWeights(overlap_dem_pix, overlap_img_params[i].centerLineDEM, overlap_img_params[i].maxDistArrayDEM);
                         mean_DEM_image(l, k) = (float)mean_DEM_image(l, k) + (float)interp_overlap_DEM_image(x, y)*weight;
                         var2_DEM_image(l, k) = (float)var2_DEM_image(l, k) + (float)interp_overlap_DEM_image(x, y)*(float)interp_overlap_DEM_image(x, y);
                         numSamples(l, k) = (int)numSamples(l, k) + 1;
                         norm(l, k) = norm(l,k) + weight;

                         meanDEMOffset[i+1] = meanDEMOffset[i+1] + (float)interp_overlap_DEM_image(x, y);
                         numDEMSamples[i+1] = numDEMSamples[i+1] + 1;
                     }

                  }
              }
          }
        }
      }
    }


    //compute mean and variance
    int totalNumSamples = 0;
    avgStdDevDEM = 0.0;
    for (k = 0 ; k < input_DEM_image.rows(); ++k) {
       for (l = 0; l < input_DEM_image.cols(); ++l) {
         if (( numSamples(l,k)!=0) && (norm(l,k)!=0)){
            if (globalParams.useWeights == 0){
               mean_DEM_image(l, k) = mean_DEM_image(l, k)/numSamples(l,k);
            }
            else{
               mean_DEM_image(l, k) = mean_DEM_image(l, k)/norm(l,k);
            }
            var2_DEM_image(l, k) = var2_DEM_image(l, k)/numSamples(l,k)- mean_DEM_image(l, k)*mean_DEM_image(l,k);
            //printf("var2(%d, %d) =%f\n", l, k,(float)var2_DEM_image(l,k));
            //var2_DEM_image(l, k) = 0.02*(float)var2_DEM_image(l,k);

            //compute the DEM standard deviation
            var2_DEM_image(l, k) = sqrt((float)var2_DEM_image(l,k));
            avgStdDevDEM = avgStdDevDEM + var2_DEM_image(l,k)/numSamples(l,k);
            totalNumSamples++;
         }
       }
    }

    //printf("average DEM error = %f\n", totalVar2/(float)numSamples);

    //write in the previous DEM

    //compute the mean DEM oset from each DEM tile to mean

    //print the average standard dev for DEM and the mean offsets
    for (i = 0; i < 5; i++){
        meanDEMOffset[i] = meanDEMOffset[i]/numDEMSamples[i];
        printf("meanDEMoffset[%d] = %f\n", i, meanDEMOffset[i]);
    }
    printf("avgStdDevDEM = %f\n", avgStdDevDEM);

    write_georeferenced_image(mean_DEM_file,
                              //channel_cast<float>(mean_DEM_image),
                              mean_DEM_image,
                              input_DEM_geo, TerminalProgressCallback("{Core}","Processing:"));

    write_georeferenced_image(var2_DEM_file,
                              //channel_cast<float>(var2_DEM_image),
                              channel_cast<uint8>(clamp(var2_DEM_image,0.0,255.0)),
                              //var2_DEM_image,
                              input_DEM_geo, TerminalProgressCallback("{Core}","Processing:"));


    //write_georeferenced_image(output_file, tm_image, geo1, TerminalProgressCallback());


}
*/

//initializes the DEM file by getting the average DEM values in the overlapping areas of consecutive DEM files.
void DetectDEMOutliers( std::string input_DEM_file,
                        std::string /*mean_DEM_file*/,
                        std::string var2_DEM_file,
                        ModelParams input_img_params,
                        std::vector<std::string> overlap_DEM_files,
                        std::vector<ModelParams> /*overlap_img_params*/,
                        GlobalParams /*globalParams*/) {

    DiskImageView<PixelGray<float> >  input_DEM_image(input_DEM_file);
    GeoReference input_DEM_geo;
    read_georeference(input_DEM_geo, input_DEM_file);

    DiskImageView<PixelGray<float> >  mean_DEM_image(input_DEM_file);

    ImageView<PixelMask<PixelGray<float> > >var2_DEM_image(input_DEM_image.cols(), input_DEM_image.rows());
    ImageView<PixelGray<int> > numSamples(input_DEM_image.cols(), input_DEM_image.rows());

    float avgStdDevDEM;

    float *meanDEMOffset = new float[5];
    //read the meanDEMOffset and avgStdDevDEM from file
    FILE *fp = fopen(input_img_params.infoFilename.c_str(),"r");

    fscanf(fp, "%f %f %f %f %f\n", &meanDEMOffset[0], &meanDEMOffset[1], &meanDEMOffset[2], &meanDEMOffset[3], &avgStdDevDEM);
    fclose(fp);


    //initialize  mean_DEM-image, var2_DEM_image and numSamples
    for (int32 k = 0 ; k < input_DEM_image.rows(); ++k) {
      for (int32 l = 0; l < input_DEM_image.cols(); ++l) {

           numSamples(l, k) = 0;
           Vector2 input_DEM_pix(l,k);

           if ( input_DEM_image(l,k) != -10000 ) {
               var2_DEM_image(l, k) = ((float)input_DEM_image(l,k)-meanDEMOffset[0]) *((float)input_DEM_image(l,k)-meanDEMOffset[0]);
               numSamples(l, k) = 1;
           }

        }
    }

    for (size_t i = 0; i < overlap_DEM_files.size(); i++){

      printf("DEM = %s\n", overlap_DEM_files[i].c_str());

      DiskImageView<PixelGray<float> >  overlap_DEM_image(overlap_DEM_files[i]);
      GeoReference overlap_DEM_geo;
      read_georeference(overlap_DEM_geo, overlap_DEM_files[i]);


      ImageViewRef<PixelGray<float> >  interp_overlap_DEM_image = interpolate(edge_extend(overlap_DEM_image.impl(),
                                                                              ConstantEdgeExtension()),
                                                                              BilinearInterpolation());

      for (int32 k = 0 ; k < input_DEM_image.rows(); ++k) {
        for (int32 l = 0; l < input_DEM_image.cols(); ++l) {

          Vector2 input_DEM_pix(l,k);

          if ( input_DEM_image(l,k) != -10000 ) {

              //check for overlap between the output image and the input DEM image
              Vector2 overlap_dem_pix = overlap_DEM_geo.lonlat_to_pixel(input_DEM_geo.pixel_to_lonlat(input_DEM_pix));
              int x = (int)overlap_dem_pix[0];
              int y = (int)overlap_dem_pix[1];

              //check for valid DEM coordinates
              if ((x>=0) && (x < overlap_DEM_image.cols()) && (y>=0) && (y< overlap_DEM_image.rows())){

                if ( overlap_DEM_image(x, y) != -10000 ) {

                    var2_DEM_image(l, k) = (float)var2_DEM_image(l, k) + ((float)interp_overlap_DEM_image(x, y)-meanDEMOffset[i+1])*
                                                                         ((float)interp_overlap_DEM_image(x, y)-meanDEMOffset[i+1]);
                    numSamples(l, k) = (int)numSamples(l, k) + 1;
                }
             }
          }
        }
      }
    }


    //compute mean and variance
    int totalNumSamples = 0;
    avgStdDevDEM = 0.0;
    for (int32 k = 0 ; k < input_DEM_image.rows(); ++k) {
      for (int32 l = 0; l < input_DEM_image.cols(); ++l) {
         if (numSamples(l,k)!=0){

            var2_DEM_image(l, k) = var2_DEM_image(l, k)/numSamples(l,k) - mean_DEM_image(l, k)*mean_DEM_image(l,k);

            //compute the DEM standard deviation
            var2_DEM_image(l, k) = sqrt((float)var2_DEM_image(l,k));
            avgStdDevDEM = avgStdDevDEM + var2_DEM_image(l,k)/numSamples(l,k);
            totalNumSamples++;
         }
       }
    }

    //printf("average DEM error = %f\n", totalVar2/(float)numSamples);

    //write the DEM
    write_georeferenced_image(var2_DEM_file,
                              //channel_cast<float>(var2_DEM_image),
                              channel_cast<uint8>(clamp(var2_DEM_image,0.0,255.0)),
                              //var2_DEM_image,
                              input_DEM_geo, TerminalProgressCallback("photometry","Processing:"));


    //write_georeferenced_image(output_file, tm_image, geo1, TerminalProgressCallback());
}


//input_files[i], input_files[i-1], output_files[i], output_files[i-1]
//writes the albedo of the current image in the area of overlap with the previous mage
//writes the previous image in the area of overal with the current image
void vw::photometry::ComputeSaveDEM(std::string curr_input_file,
                                    std::string prev_input_file,
                                    std::string prior_DEM_file,
                                    std::string output_DEM_file,
                                    ModelParams currModelParams,
                                    ModelParams prevModelParams) {
    DiskImageView<PixelMask<PixelGray<uint8> > > curr_image(curr_input_file);
    DiskImageView<PixelMask<PixelGray<uint8> > > prev_image(prev_input_file);

    GeoReference prev_geo, curr_geo;
    read_georeference(prev_geo, prev_input_file);
    read_georeference(curr_geo, curr_input_file);

    float prevSunCorrection, currSunCorrection;

    printf("exposure_time = %f, a_rescale = %f, b_rescale = %f\n",
            currModelParams.exposureTime, currModelParams.rescalingParams[0], currModelParams.rescalingParams[1]);

    ImageViewRef<PixelMask<PixelGray<uint8> > >  interp_prev_image = interpolate(edge_extend(prev_image.impl(),
                                                                                 ConstantEdgeExtension()),
                                                                                 BilinearInterpolation());

    //read the prior DEM file
    DiskImageView<PixelGray<float> >  prior_dem_image(prior_DEM_file);
    GeoReference prior_dem_geo;
    read_georeference(prior_dem_geo, prior_DEM_file);

    ImageView<PixelGray<float> > out_dem_image(prior_dem_image.cols(), prior_dem_image.rows());

    for (int k=0; k < (int)curr_image.rows(); ++k) {
      for (int l=0; l < (int)curr_image.cols(); ++l) {

         Vector2 curr_sample_pix(l,k);

         if ( is_valid(curr_image(l,k)) ) {

           Vector2 lon_lat = curr_geo.pixel_to_lonlat(curr_sample_pix);
           Vector2 sample_pix_dem = prior_dem_geo.lonlat_to_pixel(prior_dem_geo.pixel_to_lonlat(curr_sample_pix));

           int x = (int)sample_pix_dem[0];
           int y = (int)sample_pix_dem[1];

           //check for valid DEM coordinates
           if ((x>=0) && (x < out_dem_image.cols()) && (y>=0) && (y< out_dem_image.rows())){

             Vector3 longlat3(lon_lat(0),lon_lat(1),(out_dem_image)(x, y));
             Vector3 xyz = curr_geo.datum().geodetic_to_cartesian(longlat3);

             Vector2 sample_pix_dem_left;
             sample_pix_dem_left(0) = x-1;
             sample_pix_dem_left(1) = y;

             Vector2 sample_pix_dem_top;
             sample_pix_dem_top(0) = x;
             sample_pix_dem_top(1) = y-1;


             //check for valid DEM pixel value and valid left and top coordinates
             if ((sample_pix_dem_left(0) >= 0) && (sample_pix_dem_top(1) >= 0) && (out_dem_image(x,y) != -10000)){

               Vector2 lon_lat_left = prior_dem_geo.pixel_to_lonlat(sample_pix_dem_left);
               Vector3 longlat3_left(lon_lat_left(0),lon_lat_left(1),(out_dem_image)(sample_pix_dem_left(0), sample_pix_dem_left(1)));
               Vector3 xyz_left = curr_geo.datum().geodetic_to_cartesian(longlat3_left);

               Vector2 lon_lat_top = prior_dem_geo.pixel_to_lonlat(sample_pix_dem_top);
               Vector3 longlat3_top(lon_lat_top(0),lon_lat_top(1),(out_dem_image)(sample_pix_dem_top(0), sample_pix_dem_top(1)));
               Vector3 xyz_top = curr_geo.datum().geodetic_to_cartesian(longlat3_top);

               Vector3 normal = computeNormalFrom3DPoints(xyz, xyz_left, xyz_top);

               currSunCorrection = -computeReflectanceFromNormal(currModelParams.sunPosition, xyz,  normal);

               out_dem_image(l, k) = (float)curr_image(l, k)/(currModelParams.exposureTime*currSunCorrection);

               //determine the point in the previous image with the same lon and lat
               Vector2 prev_sample_pix = prev_geo.lonlat_to_pixel(lon_lat);
               prev_sample_pix[0] = floor(prev_sample_pix[0]);
               prev_sample_pix[1] = floor(prev_sample_pix[1]);
               PixelMask<PixelGray<uint8> > prev_image_pixel = interp_prev_image(prev_sample_pix[0], prev_sample_pix[1]);

               if ( is_valid (prev_image_pixel) ) {

                 prevSunCorrection = -computeReflectanceFromNormal(prevModelParams.sunPosition, xyz,  normal);

                 float intensity, albedo;

                 Vector3 xyz_prior;
                 float curr_grad = ComputeGradient_DEM(intensity, currModelParams.exposureTime, albedo, currModelParams.sunPosition, xyz, xyz_left, xyz_top, xyz_prior);
                 float prev_grad = ComputeGradient_DEM(prev_image_pixel, prevModelParams.exposureTime, albedo, prevModelParams.sunPosition, xyz, xyz_left, xyz_top, xyz_prior);
                 float curr_error = ComputeError_DEM(intensity, currModelParams.exposureTime, albedo, currSunCorrection, xyz, xyz_prior);
                 float prev_error = ComputeError_DEM(prev_image_pixel, prevModelParams.exposureTime, albedo, prevSunCorrection, xyz, xyz_prior);
                 float delta;

                 //compute delta
                 delta = (prev_grad*prev_error + curr_grad*curr_error)/(prev_grad*prev_grad + curr_grad*curr_grad);
                 xyz[2] = xyz[2] + delta;

                 out_dem_image(l, k) = xyz[2];

               }

             }
           }
         }
       }
    }

    printf("output_DEM_file = %s\n", output_DEM_file.c_str());

    write_georeferenced_image(output_DEM_file,
                              channel_cast<uint8>(clamp(out_dem_image,0.0,255.0)),
                              curr_geo, TerminalProgressCallback("photometry","Processing:"));


}
