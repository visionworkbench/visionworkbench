// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <iostream>
#include <fstream>

#include <vw/Core.h>
#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Math.h>
#include <vw/Cartography.h>

using namespace vw;
using namespace vw::cartography;

#include <vw/Photometry/Albedo.h>
#include <vw/Photometry/Reconstruct.h>
#include <vw/Photometry/ReconstructError.h>
#include <vw/Photometry/Reflectance.h>
#include <vw/Photometry/Weights.h>

using namespace vw::photometry;


float
vw::photometry::ComputeError(float intensity, float T,
    float albedo, float reflectance) {//Vector3 /*xyz*/, Vector3 /*xyz_prior*/) {

  float error;
  error = (intensity-T*albedo*reflectance);
  // std::cout << "intensity " << intensity << " albedo " << albedo << " reflectance " << reflectance << std::endl;
  return error;
}

void
vw::photometry::ComputeReconstructionErrorMap(ModelParams input_img_params,
    std::vector<ModelParams> overlap_img_params,
    GlobalParams globalParams,
    float *avgError, int *totalNumSamples) {

  int i, l, k;
  std::string input_img_file = input_img_params.inputFilename;
  std::string DEM_file = input_img_params.meanDEMFilename;
  std::string shadow_file = input_img_params.shadowFilename;
  std::string albedo_file = input_img_params.outputFilename;
  std::string error_img_file = input_img_params.errorFilename;

  DiskImageView<PixelMask<PixelGray<uint8> > > input_img(input_img_file);
  GeoReference input_img_geo;
  read_georeference(input_img_geo, input_img_file);

  DiskImageView<PixelGray<float> > input_dem_image(DEM_file);
  GeoReference input_dem_geo;
  read_georeference(input_dem_geo, DEM_file);

  DiskImageView<PixelMask<PixelGray<uint8> > > shadowImage(shadow_file);

  DiskImageView<PixelMask<PixelGray<uint8> > > albedo (albedo_file);

  ImageView<PixelMask<PixelGray<float> > > error_img (input_img.cols(), input_img.rows());

  ImageView<PixelGray<int> > numSamples(input_img.cols(), input_img.rows());

  Vector3 xyz;
  Vector3 xyz_prior;
  int x, y;

  // Wrong way of doing interpolation. See the Stereo module for the right way.
  ImageViewRef<PixelGray<float> > interp_dem_image = interpolate(edge_extend(input_dem_image.impl(),
        ConstantEdgeExtension()),
      BilinearInterpolation());

  //initialize the nominator and denomitor images
  for (k = 0 ; k < input_img.rows(); ++k) {
    for (l = 0; l < input_img.cols(); ++l) {

      Vector2 input_image_pix(l,k);
      numSamples(l,k) = 0;

      //reject invalid pixels and pixels that are in shadow.
      if ( is_valid(input_img(l,k)) && ( shadowImage(l, k) == 0)) {

        //get the corresponding DEM value

        Vector2 lon_lat = input_img_geo.pixel_to_lonlat(input_image_pix);
        Vector2 input_dem_pix = input_dem_geo.lonlat_to_pixel(input_img_geo.pixel_to_lonlat(input_image_pix));

        x = (int)input_dem_pix[0];
        y = (int)input_dem_pix[1];


        //check for valid DEM coordinates
        if ((x>=0) && (x < input_dem_image.cols()) && (y>=0) && (y< input_dem_image.rows())){

          //Vector3 longlat3(lon_lat(0),lon_lat(1),(input_dem_image)(x, y));
          Vector3 longlat3(lon_lat(0),lon_lat(1),(interp_dem_image)(x, y));
          Vector3 xyz = input_img_geo.datum().geodetic_to_cartesian(longlat3);//3D coordinates in the img coordinates


          Vector2 input_img_left_pix;
          input_img_left_pix(0) = l-1;
          input_img_left_pix(1) = k;

          Vector2 input_img_top_pix;
          input_img_top_pix(0) = l;
          input_img_top_pix(1) = k-1;

          //check for valid DEM pixel value and valid left and top coordinates
          //if ((input_dem_left_pix(0) >= 0) && (input_dem_top_pix(1) >= 0) && (input_dem_image(x,y) != -10000)){
          if ((input_img_left_pix(0) >= 0) && (input_img_top_pix(1) >= 0) && (input_dem_image(x,y) != -10000)){

            //determine the 3D coordinates of the pixel left of the current pixel
            Vector2 input_dem_left_pix = input_dem_geo.lonlat_to_pixel(input_img_geo.pixel_to_lonlat(input_img_left_pix));
            Vector2 lon_lat_left = input_img_geo.pixel_to_lonlat(input_img_left_pix);
            Vector3 longlat3_left(lon_lat_left(0),lon_lat_left(1),(interp_dem_image)(input_dem_left_pix(0), input_dem_left_pix(1)));

            //Vector3 xyz_left = input_dem_geo.datum().geodetic_to_cartesian(longlat3_left);
            Vector3 xyz_left = input_img_geo.datum().geodetic_to_cartesian(longlat3_left);

            //determine the 3D coordinates of the pixel top of the current pixel
            Vector2 input_dem_top_pix = input_dem_geo.lonlat_to_pixel(input_img_geo.pixel_to_lonlat(input_img_top_pix));
            Vector2 lon_lat_top = input_img_geo.pixel_to_lonlat(input_img_top_pix);
            Vector3 longlat3_top(lon_lat_top(0),lon_lat_top(1),(interp_dem_image)(input_dem_top_pix(0), input_dem_top_pix(1)));
            Vector3 xyz_top = input_img_geo.datum().geodetic_to_cartesian(longlat3_top);

            //Vector3 normal = computeNormalFrom3DPoints(xyz, xyz_left, xyz_top);
            Vector3 normal = computeNormalFrom3DPointsGeneral(xyz, xyz_left, xyz_top);

            //This part is the only image depedent part - START
            float input_img_reflectance, phaseAngle;
            input_img_reflectance = ComputeReflectance(normal, xyz, input_img_params, globalParams, phaseAngle);

            if (input_img_reflectance > 0){
              float input_img_error = ComputeError((float)input_img(l,k),
                  input_img_params.exposureTime, (float)albedo(l, k), input_img_reflectance);//, xyz, xyz_prior);

              error_img(l, k) = input_img_error*input_img_error;
              numSamples(l,k) = numSamples(l,k)+1;
            }

            //This part is the only image depedent part - END
          }
        }
      }
    }//l
  }//k


  //update from the overlapping images
  for (i = 0; i < (int)overlap_img_params.size(); i++){

    DiskImageView<PixelMask<PixelGray<uint8> > > overlap_img(overlap_img_params[i].inputFilename);
    GeoReference overlap_geo;
    read_georeference(overlap_geo, overlap_img_params[i].inputFilename);

    DiskImageView<PixelMask<PixelGray<uint8> > > overlapShadowImage(/*overlapShadowFileArray[i]*/overlap_img_params[i].shadowFilename);


    ImageViewRef<PixelMask<PixelGray<uint8> > > interp_overlap_img = interpolate(edge_extend(overlap_img.impl(),
          ConstantEdgeExtension()),
        BilinearInterpolation());

    ImageViewRef<PixelMask<PixelGray<uint8> > > interpOverlapShadowImage = interpolate(edge_extend(overlapShadowImage.impl(),
          ConstantEdgeExtension()),
        BilinearInterpolation());

    for (k = 0 ; k < input_img.rows(); ++k) {
      for (l = 0; l < input_img.cols(); ++l) {

        Vector2 input_img_pix(l,k);

        if ( is_valid(input_img(l,k)) ) {

          //get the corresponding DEM value
          Vector2 lon_lat = input_img_geo.pixel_to_lonlat(input_img_pix);
          Vector2 input_dem_pix = input_dem_geo.lonlat_to_pixel(input_img_geo.pixel_to_lonlat(input_img_pix));

          x = (int)input_dem_pix[0];
          y = (int)input_dem_pix[1];

          //check for valid DEM coordinates
          if ((x>=0) && (x < input_dem_image.cols()) && (y>=0) && (y< input_dem_image.rows())){

            //get the top and left DEM value
            Vector3 longlat3(lon_lat(0),lon_lat(1),(interp_dem_image)(x, y));
            Vector3 xyz = input_img_geo.datum().geodetic_to_cartesian(longlat3);

            Vector2 input_img_left_pix;
            input_img_left_pix(0) = l-1;
            input_img_left_pix(1) = k;

            Vector2 input_img_top_pix;
            input_img_top_pix(0) = l;
            input_img_top_pix(1) = k-1;

            //check for valid DEM pixel value and valid left and top coordinates
            if ((input_img_left_pix(0) >= 0) && (input_img_top_pix(1) >= 0) && (input_dem_image(x,y) != -10000)){


              //determine the 3D coordinates of the pixel left of the current pixel
              Vector2 input_dem_left_pix = input_dem_geo.lonlat_to_pixel(input_img_geo.pixel_to_lonlat(input_img_left_pix));
              Vector2 lon_lat_left = input_img_geo.pixel_to_lonlat(input_img_left_pix);
              Vector3 longlat3_left(lon_lat_left(0),lon_lat_left(1),(interp_dem_image)(input_dem_left_pix(0), input_dem_left_pix(1)));
              Vector3 xyz_left = input_img_geo.datum().geodetic_to_cartesian(longlat3_left);


              Vector2 input_dem_top_pix = input_dem_geo.lonlat_to_pixel(input_img_geo.pixel_to_lonlat(input_img_top_pix));
              Vector2 lon_lat_top = input_img_geo.pixel_to_lonlat(input_img_top_pix);
              Vector3 longlat3_top(lon_lat_top(0),lon_lat_top(1),(interp_dem_image)(input_dem_top_pix(0), input_dem_top_pix(1)));
              Vector3 xyz_top = input_img_geo.datum().geodetic_to_cartesian(longlat3_top);

              //Vector3 normal = computeNormalFrom3DPoints(xyz, xyz_left, xyz_top);
              Vector3 normal = computeNormalFrom3DPointsGeneral(xyz, xyz_left, xyz_top);

              //check for overlap between the output image and the input DEM image
              Vector2 overlap_pix = overlap_geo.lonlat_to_pixel(input_img_geo.pixel_to_lonlat(input_img_pix));
              x = (int)overlap_pix[0];
              y = (int)overlap_pix[1];

              //image dependent part of the code - START
              PixelMask<PixelGray<uint8> > overlap_img_pixel = interp_overlap_img(x, y);

              //check for valid overlap_img coordinates
              //remove shadow pixels in the overlap_img.
              if ((x>=0) && (x < overlap_img.cols()) && (y>=0) && (y< overlap_img.rows()) && interpOverlapShadowImage(x, y) == 0){

                if ( is_valid(overlap_img_pixel) ) { //common area between input_img and overlap_img

                  float overlap_img_reflectance, phaseAngle;
                  overlap_img_reflectance = ComputeReflectance(normal, xyz, overlap_img_params[i],
                                                               globalParams, phaseAngle);
                  if (overlap_img_reflectance > 0){
                    float overlap_img_error = ComputeError((float)overlap_img_pixel, overlap_img_params[i].exposureTime, (float)albedo(l,k), overlap_img_reflectance);//, xyz, xyz_prior);

                    error_img(l, k) = error_img(l,k) + overlap_img_error*overlap_img_error;
                    numSamples(l, k) = numSamples(l,k) + 1;
                  }
                }//if
              }//if
              //image dependent part of the code - START
            }
          }
        }
      }// for l
    } // for k
  } //for i

  float l_avgError = 0;
  int l_totalNumSamples = 0;
  //finalize the output image; computes the standard deviation
  for (k = 0 ; k < error_img.rows(); ++k) {
    for (l = 0; l < error_img.cols(); ++l) {
      if ( numSamples(l,k) ) {
        error_img(l, k) = (float)sqrt(error_img(l, k)/numSamples(l, k));
        l_totalNumSamples++;
        l_avgError = l_avgError + error_img(l,k);
      }
    }
  }

  l_avgError = l_avgError/l_totalNumSamples;

  *avgError = l_avgError;
  *totalNumSamples = l_totalNumSamples;
  //write the output (standard deviation of the reconstructed albedo) image
  write_georeferenced_image(error_img_file,
      channel_cast<uint8>(clamp(error_img,0.0,255.0)),
      input_img_geo, TerminalProgressCallback("photometry}","Processing:"));
}
