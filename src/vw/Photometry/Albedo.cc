// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
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

//josh - moved to ReconstructError.cc and changed name to ComputeError
//float ComputeError_Albedo(float intensity, float T,
//                          float albedo, float reflectance,
//                          Vector3 /*xyz*/, Vector3 /*xyz_prior*/) {
//  float error;
//  error = (intensity-T*albedo*reflectance);
//  return error;
//}

float ComputeGradient_Albedo(float T, float reflectance) {
  float grad;
  grad = T*reflectance;

  return grad;
}

void vw::photometry::InitImageMosaic(ModelParams input_img_params,
                                     //std::vector<std::string> overlap_img_files,
                                     std::vector<ModelParams> overlap_img_params,
                                     GlobalParams globalParams) {

    printf("image mosaic initialization\n");

    int i, l, k;

    std::string shadow_file = input_img_params.shadowFilename;
    std::string output_img_file = input_img_params.outputFilename;
    std::string input_img_file = input_img_params.inputFilename;

    DiskImageView<PixelMask<PixelGray<uint8> > >  input_img(input_img_file);
    GeoReference input_img_geo;
    read_georeference(input_img_geo, input_img_file);

    DiskImageView<PixelMask<PixelGray<uint8> > >  shadowImage(shadow_file);

    ImageView<PixelMask<PixelGray<float> > > output_img (input_img.cols(), input_img.rows());

    printf("temp mem allocation-START\n");
    ImageView<PixelGray<int> > numSamples(input_img.cols(), input_img.rows());
    printf("numSamples allocation\n");
    ImageView<PixelGray<float> > norm(input_img.cols(), input_img.rows());
    printf("temp mem allocation-END\n");

    int x,y;
    //initialize  output_img, and numSamples
    for (k = 0 ; k < input_img.rows(); ++k) {
        for (l = 0; l < input_img.cols(); ++l) {

           numSamples(l, k) = 0;
           Vector2 input_image_pix(l,k);

           if ( is_valid(input_img(l,k)) && ( shadowImage(l, k) == 0) ) {

              //compute the local reflectance
              //Vector2 lon_lat = input_img_geo.pixel_to_lonlat(input_image_pix);
              float input_img_reflectance;
              input_img_reflectance = 1.0;
              if (globalParams.useWeights == 0){
                  output_img(l, k) = (float)input_img(l,k)/(input_img_params.exposureTime*input_img_reflectance);
                  numSamples(l, k) = 1;
              }
              else{
                  float weight = ComputeLineWeights(input_image_pix, input_img_params.centerLine, input_img_params.maxDistArray);
                  output_img(l, k) = ((float)input_img(l,k)*weight)/(input_img_params.exposureTime*input_img_reflectance);
                  norm(l, k) = weight;
                  numSamples(l, k) = 1;
              }
           }
       }
    }

    //update the initial image mosaic
    for (i = 0; i < (int)overlap_img_params.size(); i++){

      printf("overlap_img = %s\n", overlap_img_params[i].inputFilename.c_str());

      DiskImageView<PixelMask<PixelGray<uint8> > >  overlap_img(overlap_img_params[i].inputFilename);
      GeoReference overlap_geo;
      read_georeference(overlap_geo, overlap_img_params[i].inputFilename);

      ImageViewRef<PixelMask<PixelGray<uint8> > >  interp_overlap_img = interpolate(edge_extend(overlap_img.impl(),
                                                                                    ConstantEdgeExtension()),
                                                                                    BilinearInterpolation());
      DiskImageView<PixelMask<PixelGray<uint8> > >  overlapShadowImage(overlap_img_params[i].shadowFilename);
      ImageViewRef<PixelMask<PixelGray<uint8> > >  interpOverlapShadowImage = interpolate(edge_extend(overlapShadowImage.impl(),
                                                                                           ConstantEdgeExtension()),
                                                                                           BilinearInterpolation());

      for (k = 0 ; k < input_img.rows(); ++k) {
        for (l = 0; l < input_img.cols(); ++l) {

          Vector2 input_img_pix(l,k);

          if ( is_valid(input_img(l,k)) ) {

              //get the corresponding DEM value
              Vector2 lon_lat = input_img_geo.pixel_to_lonlat(input_img_pix);

              //check for overlap between the output image and the input DEM image
              Vector2 overlap_pix = overlap_geo.lonlat_to_pixel(input_img_geo.pixel_to_lonlat(input_img_pix));
              x = (int)overlap_pix[0];
              y = (int)overlap_pix[1];

              //image dependent part of the code  - START
              PixelMask<PixelGray<uint8> > overlap_img_pixel = interp_overlap_img(x, y);

              //check for valid overlap_img coordinates
              //TO DO: remove shadow pixels in the overlap_img.
              if ((x>=0) && (x < overlap_img.cols()) && (y>=0) && (y < overlap_img.rows()) && (interpOverlapShadowImage(x, y) == 0)){

                    if ( is_valid(overlap_img_pixel) ) { //overlaping area between input_img and overlap_img

                        float overlap_img_reflectance = 1.0;
                        if (globalParams.useWeights == 0){
                              output_img(l, k) = (float)output_img(l, k) + (float)overlap_img_pixel/(overlap_img_params[i].exposureTime*overlap_img_reflectance);
                              numSamples(l, k) = numSamples(l,k) + 1;
                          }
                          else{
                             float weight = ComputeLineWeights(overlap_pix, overlap_img_params[i].centerLine, overlap_img_params[i].maxDistArray);
                             output_img(l, k) = (float)output_img(l, k) + ((float)overlap_img_pixel*weight)/(overlap_img_params[i].exposureTime*overlap_img_reflectance);
                             numSamples(l, k) = numSamples(l,k) + 1;
                             norm(l,k) = norm(l,k) + weight;
                          }

                    }//if ( is_valid(overlap_img_pixel) )
              }//if
          }
        }
      }
    }

    //compute the average image mosaic value
    for (k = 0 ; k < input_img.rows(); ++k) {
       for (l = 0; l < input_img.cols(); ++l) {

         if ( (is_valid(input_img(l,k))) && (numSamples(l, k)!=0) ) {
              if (globalParams.useWeights == 0){
                  output_img(l, k) = output_img(l, k)/numSamples(l,k);
              }
              else{
                  output_img(l, k) = output_img(l, k)/norm(l,k);
              }
           }
       }
    }

    /*
    //TODO: compute the image variance (standard deviation)
    for (k = 0 ; k < input_img.rows(); ++k) {
       for (l = 0; l < input_img.cols(); ++l) {

         if ( (is_valid(input_img(l,k))) && (numSamples(l, k)!=0) ) {
              output_img(l, k) = output_img(l, k)/numSamples(l,k);
           }
       }
    }
    */

    //write in the previous DEM
    write_georeferenced_image(output_img_file,
                              channel_cast<uint8>(clamp(output_img,0.0,255.0)),
                              input_img_geo, TerminalProgressCallback("photometry","Processing:"));

}

void vw::photometry::InitImageMosaicByBlocks(ModelParams input_img_params,
                                             std::vector<ModelParams> overlap_img_params,
                                             GlobalParams globalParams) {

    printf("image mosaic by block initialization\n");

    std::string input_img_file = input_img_params.inputFilename;
    std::string shadow_file = input_img_params.shadowFilename;
    std::string output_img_file = input_img_params.outputFilename;

    int horBlockSize = 500;
    int verBlockSize = 500;

    int i, l, k, lb, kb;
    int x,y;

    DiskImageView<PixelMask<PixelGray<uint8> > >  input_img(input_img_file);
    GeoReference input_img_geo;
    read_georeference(input_img_geo, input_img_file);

    DiskImageView<PixelMask<PixelGray<uint8> > >  shadowImage(shadow_file);

    ImageView<PixelMask<PixelGray<float> > > output_img (input_img.cols(), input_img.rows());

    int numHorBlocks = input_img.cols()/horBlockSize + 1;
    int numVerBlocks = input_img.rows()/verBlockSize + 1;

    printf("numHorBlocks = %d, numVerBlocks = %d\n", numHorBlocks, numVerBlocks);

    ImageView<PixelGray<int> > numSamples(horBlockSize, verBlockSize);
    ImageView<PixelGray<float> > norm(horBlockSize, verBlockSize);

    for (kb = 0 ; kb < numVerBlocks; ++kb) {
       for (lb = 0; lb < numHorBlocks; ++lb) {

         printf("kb = %d, lb=%d\n", kb, lb);

          //initialize  output_img, numSamples and norm
         for (k = 0 ; k < verBlockSize; ++k) {
           for (l = 0; l < horBlockSize; ++l) {

              int ii = kb*horBlockSize+k;
              int jj = lb*verBlockSize+l;

              if ((ii < input_img.rows()) && (jj < input_img.cols())){

                 numSamples(l, k) = 0;

                 Vector2 input_image_pix(jj,ii);

                 if ( is_valid(input_img(jj,ii)) && ( shadowImage(l, k) == 0)) {

                   float input_img_reflectance = 1.0;
                   if (globalParams.useWeights == 0){
                      output_img(jj, ii) = (float)input_img(jj, ii)/(input_img_params.exposureTime*input_img_reflectance);
                      numSamples(l, k) = 1;
                   }
                   else{
                      float weight = ComputeLineWeights(input_image_pix, input_img_params.centerLine, input_img_params.maxDistArray);
                      output_img(jj,ii) = ((float)input_img(jj,ii)*weight)/(input_img_params.exposureTime*input_img_reflectance);
                      norm(l, k) = weight;
                      numSamples(l, k) = 1;
                   }
                }
             }
           } //l
         } //k

         printf ("done with initialization block index %d %d\n", kb, lb);

         //update the initial image mosaic
         for (i = 0; i < (int)overlap_img_params.size(); i++){
           printf("overlap_img = %s\n", overlap_img_params[i].inputFilename.c_str());

           DiskImageView<PixelMask<PixelGray<uint8> > >  overlap_img(overlap_img_params[i].inputFilename);
           GeoReference overlap_geo;
           read_georeference(overlap_geo, overlap_img_params[i].inputFilename);

           ImageViewRef<PixelMask<PixelGray<uint8> > >  interp_overlap_img = interpolate(edge_extend(overlap_img.impl(),
                                                                                         ConstantEdgeExtension()),
                                                                                         BilinearInterpolation());

           DiskImageView<PixelMask<PixelGray<uint8> > >  overlapShadowImage(overlap_img_params[i].shadowFilename);
           ImageViewRef<PixelMask<PixelGray<uint8> > >  interpOverlapShadowImage = interpolate(edge_extend(overlapShadowImage.impl(),
                                                                                               ConstantEdgeExtension()),
                                                                                               BilinearInterpolation());


           for (k = 0 ; k < verBlockSize; ++k) {
             for (l = 0; l < horBlockSize; ++l) {

               int ii = kb*horBlockSize+k;
               int jj = lb*verBlockSize+l;

               Vector2 input_img_pix (jj,ii);

               if ((ii < input_img.rows()) && (jj < input_img.cols())){

                 if ( is_valid(input_img(jj,ii)) ) {

                   //get the corresponding DEM value
                   Vector2 lon_lat = input_img_geo.pixel_to_lonlat(input_img_pix);

                   //check for overlap between the output image and the input DEM image
                   Vector2 overlap_pix = overlap_geo.lonlat_to_pixel(input_img_geo.pixel_to_lonlat(input_img_pix));
                   x = (int)overlap_pix[0];
                   y = (int)overlap_pix[1];


                   PixelMask<PixelGray<uint8> > overlap_img_pixel = interp_overlap_img(x, y);

                   //check for valid overlap_img coordinates
                   //TO DO: remove shadow pixels in the overlap_img.
                   if ((x>=0) && (x < overlap_img.cols()) && (y>=0) && (y < overlap_img.rows()) && (interpOverlapShadowImage(x, y) == 0)){

                     if ( is_valid(overlap_img_pixel) ) { //overlaping area between input_img and overlap_img

                        float overlap_img_reflectance = 1.0;
                        if (globalParams.useWeights == 0){
                            output_img(jj,ii) = (float)output_img(jj,ii) + (float)overlap_img_pixel/(overlap_img_params[i].exposureTime*overlap_img_reflectance);
                            numSamples(l, k) = numSamples(l,k) + 1;
                         }
                         else{
                             float weight = ComputeLineWeights(overlap_pix, overlap_img_params[i].centerLine, overlap_img_params[i].maxDistArray);
                             output_img(jj,ii) = (float)output_img(jj,ii) + ((float)overlap_img_pixel*weight)/(overlap_img_params[i].exposureTime*overlap_img_reflectance);
                             numSamples(l, k) = numSamples(l,k) + 1;
                             norm(l,k) = norm(l,k) + weight;
                         }
                     }//if ( is_valid(overlap_img_pixel) )
                   }//if
                 }
               }
             }
           }
         }
         printf ("done with update block index %d %d\n", kb, lb);

         //compute the estimated image mosaic value
         for (k = 0 ; k < verBlockSize; ++k) {
            for (l = 0; l < horBlockSize; ++l) {

               int ii = kb*horBlockSize+k;
               int jj = lb*verBlockSize+l;
               if ((ii < input_img.rows()) && (jj < input_img.cols())){

                 if ( (is_valid(input_img(jj,ii))) && (numSamples(l, k)!=0) ) {
                   if (globalParams.useWeights == 0){
                     output_img(jj,ii) = output_img(jj,ii)/numSamples(l,k);
                   }
                   else{
                     output_img(jj,ii) = output_img(jj,ii)/norm(l,k);
                   }
                 }
               }
           }
         }
         printf("done computed the final init value\n");
         /*
         //TODO: compute the image variance (standard deviation)
         for (k = 0 ; k < input_img.rows(); ++k) {
            for (l = 0; l < input_img.cols(); ++l) {
               if ( (is_valid(input_img(l,k))) && (numSamples(l, k)!=0) ) {
                  output_img(l, k) = output_img(l, k)/numSamples(l,k);
               }
            }
         }
         */

       } //lb
    }//kb



    //write in the previous DEM
    write_georeferenced_image(output_img_file,
                              channel_cast<uint8>(clamp(output_img,0.0,255.0)),
                              input_img_geo, TerminalProgressCallback("photometry","Processing:"));

}

//updates the image mosaic
//author: Ara Nefian
void vw::photometry::UpdateImageMosaic(ModelParams input_img_params,
                                       std::vector<ModelParams> overlap_img_params,
                                       GlobalParams globalParams) {
    int i, l, k;

    std::string input_img_file = input_img_params.inputFilename;
    std::string shadow_file = input_img_params.shadowFilename;
    std::string output_img_file = input_img_params.outputFilename;

    DiskImageView<PixelMask<PixelGray<uint8> > >  input_img(input_img_file);
    GeoReference input_img_geo;
    read_georeference(input_img_geo, input_img_file);


    DiskImageView<PixelMask<PixelGray<uint8> > >  shadowImage(shadow_file);

    DiskImageView<PixelMask<PixelGray<uint8> > > output_img_r(output_img_file);


    ImageView<PixelMask<PixelGray<float> > > output_img (output_img_r.cols(), output_img_r.rows());

    ImageView<PixelGray<float> > nominator(input_img.cols(), input_img.rows());
    ImageView<PixelGray<float> > denominator(input_img.cols(), input_img.rows());

    Vector3 xyz;
    Vector3 xyz_prior;

    //initialize the nominator and denomitor images
    for (k = 0 ; k < input_img.rows(); ++k) {
        for (l = 0; l < input_img.cols(); ++l) {

           nominator(l, k) = 0;
           denominator(l, k) = 0;

           Vector2 input_image_pix(l,k);

           //reject invalid pixels and pixels that are in shadow.
           if ( is_valid(input_img(l,k)) && ( shadowImage(l, k) == 0)) {

              Vector2 lon_lat = input_img_geo.pixel_to_lonlat(input_image_pix);


              float input_img_reflectance = 1.0;

              float input_img_error = ComputeError((float)input_img(l,k), input_img_params.exposureTime,
                                                                 (float)output_img_r(l, k),  input_img_reflectance);

//              float input_img_error = ComputeError_Albedo((float)input_img(l,k), input_img_params.exposureTime,
//                                                                 (float)output_img_r(l, k),  input_img_reflectance, xyz, xyz_prior);

              float input_albedo_grad = ComputeGradient_Albedo(input_img_params.exposureTime, input_img_reflectance);

              if (globalParams.useWeights == 0){
                  nominator(l, k) = input_albedo_grad*input_img_error;
                  denominator(l, k) = input_albedo_grad*input_albedo_grad;
              }
              else{
                  float weight = ComputeLineWeights(input_image_pix, input_img_params.centerLine, input_img_params.maxDistArray);
                  nominator(l, k)   = input_albedo_grad*input_img_error*weight;
                  denominator(l, k) = input_albedo_grad*input_albedo_grad*weight;
              }
              output_img(l, k) = 0;//(float)(output_img_r(l, k));
              //This part is the only image depedent part - END
           }
        }
    }


    //update from the overlapping images
    for (i = 0; i < (int)overlap_img_params.size(); i++){

      printf("overlap_img = %s\n", overlap_img_params[i].inputFilename.c_str());

      DiskImageView<PixelMask<PixelGray<uint8> > >  overlap_img(overlap_img_params[i].inputFilename);
      GeoReference overlap_geo;
      read_georeference(overlap_geo, overlap_img_params[i].inputFilename);


      DiskImageView<PixelMask<PixelGray<uint8> > >  overlapShadowImage(overlap_img_params[i].shadowFilename);

      ImageViewRef<PixelMask<PixelGray<uint8> > >  interp_overlap_img = interpolate(edge_extend(overlap_img.impl(),
                                                                                    ConstantEdgeExtension()),
                                                                                    BilinearInterpolation());

      ImageViewRef<PixelMask<PixelGray<uint8> > >  interpOverlapShadowImage = interpolate(edge_extend(overlapShadowImage.impl(),
                                                                                           ConstantEdgeExtension()),
                                                                                           BilinearInterpolation());

      for (k = 0 ; k < input_img.rows(); ++k) {
        for (l = 0; l < input_img.cols(); ++l) {

         Vector2 input_img_pix(l,k);

         if ( is_valid(input_img(l,k)) ) {

              //determine the corresponding pixel in the overlapping image
              Vector2 overlap_pix = overlap_geo.lonlat_to_pixel(input_img_geo.pixel_to_lonlat(input_img_pix));
              int x = (int)overlap_pix[0];
              int y = (int)overlap_pix[1];

              PixelMask<PixelGray<uint8> > overlap_img_pixel = interp_overlap_img(x, y);

              //check for valid overlap_img coordinates and non-shadow pixels
              if ((x>=0) && (x < overlap_img.cols()) && (y>=0) && (y< overlap_img.rows()) && interpOverlapShadowImage(x, y) == 0){

                 if (is_valid(overlap_img_pixel)) { //common area between input_img and overlap_img

                     float overlap_img_reflectance = 1.0;
                     Vector3 xyz;
                     Vector3 xyz_prior;

                     float overlap_img_error = ComputeError((float)overlap_img_pixel, overlap_img_params[i].exposureTime,
                                                                   (float)output_img_r(l, k), overlap_img_reflectance);

//                     float overlap_img_error = ComputeError_Albedo((float)overlap_img_pixel, overlap_img_params[i].exposureTime,
//                                                                   (float)output_img_r(l, k), overlap_img_reflectance, xyz, xyz_prior);

                     float overlap_albedo_grad = ComputeGradient_Albedo(overlap_img_params[i].exposureTime, overlap_img_reflectance);
                     if (globalParams.useWeights == 0){
                         nominator(l, k) = nominator(l, k) + overlap_albedo_grad*overlap_img_error;
                         denominator(l, k) = denominator(l, k) + overlap_albedo_grad*overlap_albedo_grad;
                     }
                     else{
                         float weight = ComputeLineWeights(overlap_pix, overlap_img_params[i].centerLine, overlap_img_params[i].maxDistArray);
                         nominator(l, k)   = nominator(l,k) + overlap_albedo_grad*overlap_img_error*weight;
                         denominator(l, k) = denominator(l,k) + overlap_albedo_grad*overlap_albedo_grad*weight;
                     }
                 }//if
             }//if
           }
        }// for l
      } // for k
    } //for i


    //finalize the output image
    for (k = 0 ; k < output_img.rows(); ++k) {
       for (l = 0; l < output_img.cols(); ++l) {

           if ( is_valid(output_img(l,k)) ) {
             if ((float)denominator(l, k) != 0){
                float delta = (float)nominator(l, k)/(float)denominator(l, k);
                //printf("k = %d, l = %d, output_img = %f, delta = %f\n", k, l, (float)output_img(l,k), delta);
                output_img(l, k) = (float)output_img_r(l, k) + delta;
             }
           }
       }
    }

    //write the output (albedo) image
    write_georeferenced_image(output_img_file,
                              channel_cast<uint8>(clamp(output_img,0.0,255.0)),
                              input_img_geo, TerminalProgressCallback("photometry","Processing:"));

}

//-------------------------------------------------------------------------------
//Below are the functions for albedo reconstruction
//-------------------------------------------------------------------------------

//initializes the albedo mosaic
//TO DO: build the version that does the blockwise processing to deal with large scale images
void
vw::photometry::InitAlbedoMosaic(ModelParams input_img_params,
                                 std::vector<ModelParams> overlap_img_params,
                                 GlobalParams globalParams) {
    int horBlockSize = 500;
    int verBlockSize = 500;
    
    int i, l, k;
    std::string input_img_file = input_img_params.inputFilename;
    std::string reflectance_file = input_img_params.reliefFilename;
    std::string shadow_file = input_img_params.shadowFilename;
    std::string output_img_file = input_img_params.outputFilename;


    DiskImageView<PixelMask<PixelGray<uint8> > >  input_img(input_img_file);
    GeoReference input_img_geo;
    read_georeference(input_img_geo, input_img_file);

    DiskImageView<PixelMask<PixelGray<uint8> > >  shadowImage(shadow_file);

    ImageView<PixelMask<PixelGray<float> > > output_img (input_img.cols(), input_img.rows());
    ImageView<PixelGray<int> > numSamples(input_img.cols(), input_img.rows());
    ImageView<PixelGray<float> > norm(input_img.cols(), input_img.rows());

 
    DiskImageView<PixelMask<PixelGray<float> > > reflectance_image(reflectance_file);
    GeoReference reflectance_geo;
    read_georeference(reflectance_geo, reflectance_file);

    ImageViewRef<PixelMask<PixelGray<float> > >  interp_reflectance_image = interpolate(edge_extend(reflectance_image.impl(),
											ConstantEdgeExtension()),
									    BilinearInterpolation());


    int numHorBlocks = input_img.cols()/horBlockSize + 1;
    int numVerBlocks = input_img.rows()/verBlockSize + 1;

    int x,y;
    //initialize  output_img, and numSamples
    for (k = 0 ; k < input_img.rows(); ++k) {
        for (l = 0; l < input_img.cols(); ++l) {

           numSamples(l, k) = 0;
           Vector2 input_image_pix(l,k);

           if ( is_valid(input_img(l,k)) && ( shadowImage(l, k) == 0)) { //valid image point

              //compute the local reflectance
            
              Vector2 lon_lat = input_img_geo.pixel_to_lonlat(input_image_pix);
              Vector2 reflectance_pix = reflectance_geo.lonlat_to_pixel(input_img_geo.pixel_to_lonlat(input_image_pix));

              float x = reflectance_pix[0];
              float y = reflectance_pix[1];
             
              //check for valid reflectance coordinates
              if ((x>=0) && (x < reflectance_image.cols()) && (y>=0) && (y< reflectance_image.rows())){ //valid reflectance point  
               
		if (is_valid(interp_reflectance_image(x,y))){
                 
                  float input_img_reflectance = interp_reflectance_image(x,y);
                  if (input_img_reflectance != 0.0){
                      if (globalParams.useWeights == 0){
                          output_img(l, k) = (float)input_img(l,k)/(input_img_params.exposureTime*input_img_reflectance);
                          numSamples(l, k) = 1;
                      }
                      else{
                         float weight = ComputeLineWeights(input_image_pix, input_img_params.centerLine, input_img_params.maxDistArray);
                         output_img(l, k) = ((float)input_img(l,k)*weight)/(input_img_params.exposureTime*input_img_reflectance);
                         norm(l, k) = weight;
                         numSamples(l, k) = 1;
                      }
                  }
                }
             }
          }
       }
    }

    for (i = 0; i < (int)overlap_img_params.size(); i++){
      
      printf("overlap_img = %s\n", overlap_img_params[i].inputFilename.c_str());

      DiskImageView<PixelMask<PixelGray<uint8> > >  overlap_img(overlap_img_params[i].inputFilename);
      GeoReference overlap_geo;
      read_georeference(overlap_geo, overlap_img_params[i].inputFilename);

      ImageViewRef<PixelMask<PixelGray<uint8> > >  interp_overlap_img = interpolate(edge_extend(overlap_img.impl(),
                                                                                    ConstantEdgeExtension()),
                                                                                    BilinearInterpolation());


      
      DiskImageView<PixelMask<PixelGray<uint8> > >  overlapShadowImage(overlap_img_params[i].shadowFilename);

      ImageViewRef<PixelMask<PixelGray<uint8> > >  interpOverlapShadowImage = interpolate(edge_extend(overlapShadowImage.impl(),
												      ConstantEdgeExtension()),
											  BilinearInterpolation());
      
      for (k = 0 ; k < input_img.rows(); ++k) {
        for (l = 0; l < input_img.cols(); ++l) {

          Vector2 input_img_pix(l,k);

          if ( is_valid(input_img(l,k)) ) {

              //get the corresponding reflectance value
              Vector2 lon_lat = input_img_geo.pixel_to_lonlat(input_img_pix);
	     
              Vector2 reflectance_pix = reflectance_geo.lonlat_to_pixel(input_img_geo.pixel_to_lonlat(input_img_pix));
              float x = reflectance_pix[0];
              float y = reflectance_pix[1];

              //check for valid reflectance_image coordinates
              if ((x>=0) && (x < reflectance_image.cols()) && (y>=0) && (y< reflectance_image.rows())){
	
                if (is_valid(interp_reflectance_image(x,y))){
                 
                  //check for overlap between the output image and the input DEM image
                  Vector2 overlap_pix = overlap_geo.lonlat_to_pixel(input_img_geo.pixel_to_lonlat(input_img_pix));
                  float overlap_x = overlap_pix[0];
                  float overlap_y = overlap_pix[1];

                  //image dependent part of the code  - START
                  PixelMask<PixelGray<uint8> > overlap_img_pixel = interp_overlap_img(overlap_x, overlap_y);

                  //check for valid overlap_img coordinates
                  //TO DO: remove shadow pixels in the overlap_img.
                  if ((overlap_x>=0) && (overlap_x < overlap_img.cols()) && (overlap_y >= 0) && (overlap_y< overlap_img.rows()) && (interpOverlapShadowImage(x, y) == 0)){

                    if ( is_valid(overlap_img_pixel) ) { //common area between input_img and overlap_img

                      float overlap_img_reflectance = interp_reflectance_image(x,y);

                      if (overlap_img_reflectance != 0.0){
                          if (globalParams.useWeights == 0){
                              output_img(l, k) = (float)output_img(l, k) + (float)overlap_img_pixel/(overlap_img_params[i].exposureTime*overlap_img_reflectance);
                              numSamples(l, k) = numSamples(l,k) + 1;
                          }
                          else{
                             float weight = ComputeLineWeights(overlap_pix, overlap_img_params[i].centerLine, overlap_img_params[i].maxDistArray);
                             output_img(l, k) = (float)output_img(l, k) + ((float)overlap_img_pixel*weight)/(overlap_img_params[i].exposureTime*overlap_img_reflectance);
                             numSamples(l, k) = numSamples(l,k) + 1;
                             norm(l,k) = norm(l,k) + weight;
                          }
                      }

                    }//if
                  }//if

                  //image dependent part of the code  - END
               }
             }
          }
        }
      }
    }

    //compute the mean albedo value
    int numValid = 0;
    for (k = 0 ; k < input_img.rows(); ++k) {
       for (l = 0; l < input_img.cols(); ++l) {

         //output_img(l,k).invalidate();

         if ( (is_valid(input_img(l,k))) && (numSamples(l, k)!=0) ) {

              output_img(l,k).validate();

              if (globalParams.useWeights == 0){
                  output_img(l, k) = output_img(l, k)/numSamples(l,k);
              }
              else{
		output_img(l, k) = output_img(l, k)/norm(l,k);
              }
	      numValid++;
         }

      }
    }

    printf("numValid = %d, total = %d\n", numValid, input_img.rows()*input_img.cols());
    
    //TODO: compute the albedo variance (standard deviation)
   
   
    //write in the albedo image
    write_georeferenced_image(output_img_file,
                              channel_cast<uint8>(clamp(output_img,0.0,255.0)),
                              input_img_geo, TerminalProgressCallback("{Core}","Processing:"));
    
}



//writes the current albedo of the current image in the area of overlap with the previous mage
//writes the previous albedo in the area of overlap with the current image
void
vw::photometry::UpdateAlbedoMosaic(ModelParams input_img_params,
                                   std::vector<ModelParams> overlap_img_params,
                                   GlobalParams globalParams) {
    int i, l, k;

    std::string input_img_file = input_img_params.inputFilename;
    std::string DEM_file = input_img_params.meanDEMFilename;
    std::string shadow_file = input_img_params.shadowFilename;
    std::string output_img_file = input_img_params.reliefFilename;

  
    DiskImageView<PixelMask<PixelGray<uint8> > >  input_img(input_img_file);
    GeoReference input_img_geo;
    read_georeference(input_img_geo, input_img_file);

    DiskImageView<PixelGray<float> >  input_dem_image(DEM_file);
    GeoReference input_dem_geo;
    read_georeference(input_dem_geo, DEM_file);

    //TO DO: read the reflectance image instead.
 
    DiskImageView<PixelMask<PixelGray<uint8> > >  shadowImage(shadow_file);

    DiskImageView<PixelMask<PixelGray<uint8> > > output_img_r(output_img_file);


    ImageView<PixelMask<PixelGray<float> > > output_img (output_img_r.cols(), output_img_r.rows());

    ImageView<PixelGray<float> > nominator(input_img.cols(), input_img.rows());
    ImageView<PixelGray<float> > denominator(input_img.cols(), input_img.rows());

    Vector3 xyz;
    Vector3 xyz_prior;
    int x, y;

    ImageViewRef<PixelGray<float> >  interp_dem_image = interpolate(edge_extend(input_dem_image.impl(),
                                                                                ConstantEdgeExtension()),
                                                                                BilinearInterpolation());

    //initialize the nominator and denomitor images
    for (k = 0 ; k < input_img.rows(); ++k) {
        for (l = 0; l < input_img.cols(); ++l) {

           nominator(l, k) = 0;
           denominator(l, k) = 0;

           Vector2 input_image_pix(l,k);

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
                if ((input_img_left_pix(0) >= 0) && (input_img_top_pix(1) >= 0) && (input_dem_image(x,y) != globalParams.noDEMDataValue)){

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
                  float input_img_reflectance;
                  input_img_reflectance = ComputeReflectance(normal, xyz, input_img_params, globalParams);

                  if (input_img_reflectance > 0){
                     float input_img_error = ComputeError((float)input_img(l,k), input_img_params.exposureTime,
                                                                 (float)output_img_r(l, k),  input_img_reflectance);

//                     float input_img_error = ComputeError_Albedo((float)input_img(l,k), input_img_params.exposureTime,
//                                                                 (float)output_img_r(l, k),  input_img_reflectance, xyz, xyz_prior);

                     float input_albedo_grad = ComputeGradient_Albedo(input_img_params.exposureTime, input_img_reflectance);


                     if (globalParams.useWeights == 0){
                         nominator(l, k) = input_albedo_grad*input_img_error;
                         denominator(l, k) = input_albedo_grad*input_albedo_grad;
                         }
                     else{
                         float weight = ComputeLineWeights(input_image_pix, input_img_params.centerLine, input_img_params.maxDistArray);
                         nominator(l, k)   = input_albedo_grad*input_img_error*weight;
                         denominator(l, k) = input_albedo_grad*input_albedo_grad*weight;
                     }

                     output_img(l, k) = 0;//(float)(output_img_r(l, k));
                  }

                  //This part is the only image depedent part - END
                }
              }
           }
        }
    }


    //update from the overlapping images
    for (i = 0; i < (int)overlap_img_params.size(); i++){

      printf("overlap_img = %s\n", overlap_img_params[i].inputFilename.c_str());

      DiskImageView<PixelMask<PixelGray<uint8> > >  overlap_img(overlap_img_params[i].inputFilename);
      GeoReference overlap_geo;
     
      read_georeference(overlap_geo, overlap_img_params[i].inputFilename);

      DiskImageView<PixelMask<PixelGray<uint8> > >  overlapShadowImage(overlap_img_params[i].shadowFilename);

      ImageViewRef<PixelMask<PixelGray<uint8> > >  interp_overlap_img = interpolate(edge_extend(overlap_img.impl(),
                                                                                    ConstantEdgeExtension()),
                                                                                    BilinearInterpolation());

      ImageViewRef<PixelMask<PixelGray<uint8> > >  interpOverlapShadowImage = interpolate(edge_extend(overlapShadowImage.impl(),
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
                if ((input_img_left_pix(0) >= 0) && (input_img_top_pix(1) >= 0) && (input_dem_image(x,y) != globalParams.noDEMDataValue)){


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

                  //image dependent part of the code  - START
                  PixelMask<PixelGray<uint8> > overlap_img_pixel = interp_overlap_img(x, y);

                  //check for valid overlap_img coordinates
                  //TO DO: remove shadow pixels in the overlap_img.
                  if ((x>=0) && (x < overlap_img.cols()) && (y>=0) && (y< overlap_img.rows()) && interpOverlapShadowImage(x, y) == 0){

                    if ( is_valid(overlap_img_pixel) ) { //common area between input_img and overlap_img

                      float overlap_img_reflectance;
                      overlap_img_reflectance = ComputeReflectance(normal, xyz, overlap_img_params[i], globalParams);
                      if (overlap_img_reflectance > 0){
                         float overlap_img_error = ComputeError((float)overlap_img_pixel, overlap_img_params[i].exposureTime,
                                                                       (float)output_img_r(l, k), overlap_img_reflectance);
//                         float overlap_img_error = ComputeError_Albedo((float)overlap_img_pixel, overlap_img_params[i].exposureTime,
//                                                                       (float)output_img_r(l, k), overlap_img_reflectance, xyz, xyz_prior);

                         float overlap_albedo_grad = ComputeGradient_Albedo(overlap_img_params[i].exposureTime, overlap_img_reflectance);
                         if (globalParams.useWeights == 0){
                             nominator(l, k) = nominator(l, k) + overlap_albedo_grad*overlap_img_error;
                             denominator(l, k) = denominator(l, k) + overlap_albedo_grad*overlap_albedo_grad;
                         }
                         else{

                           //float weight = ComputeWeights(overlap_pix, overlap_img_params[i].center2D, overlap_img_params[i].maxDistance);
                            float weight = ComputeLineWeights(overlap_pix, overlap_img_params[i].centerLine, overlap_img_params[i].maxDistArray);
                            nominator(l, k)   = nominator(l,k) + overlap_albedo_grad*overlap_img_error*weight;
                            denominator(l, k) = denominator(l,k) + overlap_albedo_grad*overlap_albedo_grad*weight;
                         }
                      }
                    }//if
                  }//if

                  //image dependent part of the code  - START
               }
            }
          }
        }// for l
      } // for k
    } //for i


    //finalize the output image
    for (k = 0 ; k < output_img.rows(); ++k) {
       for (l = 0; l < output_img.cols(); ++l) {
           
           output_img(l,k).invalidate();
           
           if ( is_valid(output_img(l,k)) ) {
             if ((float)denominator(l, k) != 0){
                float delta = (float)nominator(l, k)/(float)denominator(l, k);
                //printf("k = %d, l = %d, output_img = %f, delta = %f\n", k, l, (float)output_img(l,k), delta);
                output_img(l,k) = (float)output_img_r(l, k) + delta;
                output_img(l,k).validate();
             }
           }

       }
    }



    //write the output (albedo) image
     write_georeferenced_image(output_img_file,
                              channel_cast<uint8>(clamp(output_img,0.0,255.0)),
                              input_img_geo, TerminalProgressCallback("photometry","Processing:"));


}
//input_files[i], input_files[i-1], output_files[i], output_files[i-1]
//writes the current albedo of the current image in the area of overlap with the previous mage
//writes the previous albedo in the area of overlap with the current image


//josh - moved to ReconstructError.cc where all error related functions will be defined
//void
//vw::photometry::ComputeReconstructionErrorMap(ModelParams input_img_params,
//                                              std::vector<ModelParams> overlap_img_params,
//                                              GlobalParams globalParams,
//                                              float *avgError, int *totalNumSamples) {
//
//    int i, l, k;
//
//    std::string input_img_file = input_img_params.inputFilename;
//    std::string DEM_file = input_img_params.meanDEMFilename;
//    std::string shadow_file = input_img_params.shadowFilename;
//    std::string albedo_file = input_img_params.outputFilename;
//    std::string error_img_file = input_img_params.errorFilename;
//
//    DiskImageView<PixelMask<PixelGray<uint8> > >  input_img(input_img_file);
//    GeoReference input_img_geo;
//    read_georeference(input_img_geo, input_img_file);
//
//    DiskImageView<PixelGray<float> >  input_dem_image(DEM_file);
//    GeoReference input_dem_geo;
//    read_georeference(input_dem_geo, DEM_file);
//
//    DiskImageView<PixelMask<PixelGray<uint8> > >  shadowImage(shadow_file);
//
//    DiskImageView<PixelMask<PixelGray<uint8> > > albedo (albedo_file);
//
//    ImageView<PixelMask<PixelGray<float> > > error_img (input_img.cols(), input_img.rows());
//
//    ImageView<PixelGray<int> > numSamples(input_img.cols(), input_img.rows());
//
//    Vector3 xyz;
//    Vector3 xyz_prior;
//    int x, y;
//
//    ImageViewRef<PixelGray<float> >  interp_dem_image = interpolate(edge_extend(input_dem_image.impl(),
//                                                                                ConstantEdgeExtension()),
//                                                                                BilinearInterpolation());
//
//    //initialize the nominator and denomitor images
//    for (k = 0 ; k < input_img.rows(); ++k) {
//        for (l = 0; l < input_img.cols(); ++l) {
//
//           Vector2 input_image_pix(l,k);
//           numSamples(l,k) = 0;
//
//           //reject invalid pixels and pixels that are in shadow.
//           if ( is_valid(input_img(l,k)) && ( shadowImage(l, k) == 0)) {
//
//              //get the corresponding DEM value
//
//              Vector2 lon_lat = input_img_geo.pixel_to_lonlat(input_image_pix);
//              Vector2 input_dem_pix = input_dem_geo.lonlat_to_pixel(input_img_geo.pixel_to_lonlat(input_image_pix));
//
//              x = (int)input_dem_pix[0];
//              y = (int)input_dem_pix[1];
//
//
//              //check for valid DEM coordinates
//              if ((x>=0) && (x < input_dem_image.cols()) && (y>=0) && (y< input_dem_image.rows())){
//
//                //Vector3 longlat3(lon_lat(0),lon_lat(1),(input_dem_image)(x, y));
//                Vector3 longlat3(lon_lat(0),lon_lat(1),(interp_dem_image)(x, y));
//                Vector3 xyz = input_img_geo.datum().geodetic_to_cartesian(longlat3);//3D coordinates in the img coordinates
//
//
//                Vector2 input_img_left_pix;
//                input_img_left_pix(0) = l-1;
//                input_img_left_pix(1) = k;
//
//                Vector2 input_img_top_pix;
//                input_img_top_pix(0) = l;
//                input_img_top_pix(1) = k-1;
//
//                //check for valid DEM pixel value and valid left and top coordinates
//                //if ((input_dem_left_pix(0) >= 0) && (input_dem_top_pix(1) >= 0) && (input_dem_image(x,y) != -10000)){
//                if ((input_img_left_pix(0) >= 0) && (input_img_top_pix(1) >= 0) && (input_dem_image(x,y) != -10000)){
//
//                  //determine the 3D coordinates of the pixel left of the current pixel
//                  Vector2 input_dem_left_pix = input_dem_geo.lonlat_to_pixel(input_img_geo.pixel_to_lonlat(input_img_left_pix));
//                  Vector2 lon_lat_left = input_img_geo.pixel_to_lonlat(input_img_left_pix);
//                  Vector3 longlat3_left(lon_lat_left(0),lon_lat_left(1),(interp_dem_image)(input_dem_left_pix(0), input_dem_left_pix(1)));
//
//                  //Vector3 xyz_left = input_dem_geo.datum().geodetic_to_cartesian(longlat3_left);
//                  Vector3 xyz_left = input_img_geo.datum().geodetic_to_cartesian(longlat3_left);
//
//                  //determine the 3D coordinates of the pixel top of the current pixel
//                  Vector2 input_dem_top_pix = input_dem_geo.lonlat_to_pixel(input_img_geo.pixel_to_lonlat(input_img_top_pix));
//                  Vector2 lon_lat_top = input_img_geo.pixel_to_lonlat(input_img_top_pix);
//                  Vector3 longlat3_top(lon_lat_top(0),lon_lat_top(1),(interp_dem_image)(input_dem_top_pix(0), input_dem_top_pix(1)));
//                  Vector3 xyz_top = input_img_geo.datum().geodetic_to_cartesian(longlat3_top);
//
//                  //Vector3 normal = computeNormalFrom3DPoints(xyz, xyz_left, xyz_top);
//                  Vector3 normal = computeNormalFrom3DPointsGeneral(xyz, xyz_left, xyz_top);
//
//                  //This part is the only image depedent part - START
//                  float input_img_reflectance;
//                  input_img_reflectance = ComputeReflectance(normal, xyz, input_img_params, globalParams);
//
//                  if (input_img_reflectance > 0){
//                     float input_img_error = ComputeError_Albedo((float)input_img(l,k), input_img_params.exposureTime,
//                                                                 (float)albedo(l, k),  input_img_reflectance, xyz, xyz_prior);
//
//                     error_img(l, k) = input_img_error*input_img_error;
//                     numSamples(l,k) = numSamples(l,k)+1;
//                  }
//
//                  //This part is the only image depedent part - END
//                }
//              }
//           }
//        }
//    }
//
//
//    //update from the overlapping images
//    for (i = 0; i < (int)overlap_img_params.size(); i++){
//
//      DiskImageView<PixelMask<PixelGray<uint8> > >  overlap_img(overlap_img_params[i].inputFilename);
//      GeoReference overlap_geo;
//      read_georeference(overlap_geo, overlap_img_params[i].inputFilename);
//
//      DiskImageView<PixelMask<PixelGray<uint8> > >  overlapShadowImage(/*overlapShadowFileArray[i]*/overlap_img_params[i].shadowFilename);
//
//
//      ImageViewRef<PixelMask<PixelGray<uint8> > >  interp_overlap_img = interpolate(edge_extend(overlap_img.impl(),
//                                                                                    ConstantEdgeExtension()),
//                                                                                    BilinearInterpolation());
//
//      ImageViewRef<PixelMask<PixelGray<uint8> > >  interpOverlapShadowImage = interpolate(edge_extend(overlapShadowImage.impl(),
//                                                                                           ConstantEdgeExtension()),
//                                                                                           BilinearInterpolation());
//
//      for (k = 0 ; k < input_img.rows(); ++k) {
//        for (l = 0; l < input_img.cols(); ++l) {
//
//         Vector2 input_img_pix(l,k);
//
//         if ( is_valid(input_img(l,k)) ) {
//
//              //get the corresponding DEM value
//              Vector2 lon_lat = input_img_geo.pixel_to_lonlat(input_img_pix);
//              Vector2 input_dem_pix = input_dem_geo.lonlat_to_pixel(input_img_geo.pixel_to_lonlat(input_img_pix));
//
//              x = (int)input_dem_pix[0];
//              y = (int)input_dem_pix[1];
//
//              //check for valid DEM coordinates
//              if ((x>=0) && (x < input_dem_image.cols()) && (y>=0) && (y< input_dem_image.rows())){
//
//                //get the top and left DEM value
//                Vector3 longlat3(lon_lat(0),lon_lat(1),(interp_dem_image)(x, y));
//                Vector3 xyz = input_img_geo.datum().geodetic_to_cartesian(longlat3);
//
//                Vector2 input_img_left_pix;
//                input_img_left_pix(0) = l-1;
//                input_img_left_pix(1) = k;
//
//                Vector2 input_img_top_pix;
//                input_img_top_pix(0) = l;
//                input_img_top_pix(1) = k-1;
//
//                //check for valid DEM pixel value and valid left and top coordinates
//                if ((input_img_left_pix(0) >= 0) && (input_img_top_pix(1) >= 0) && (input_dem_image(x,y) != -10000)){
//
//
//                  //determine the 3D coordinates of the pixel left of the current pixel
//                  Vector2 input_dem_left_pix = input_dem_geo.lonlat_to_pixel(input_img_geo.pixel_to_lonlat(input_img_left_pix));
//                  Vector2 lon_lat_left = input_img_geo.pixel_to_lonlat(input_img_left_pix);
//                  Vector3 longlat3_left(lon_lat_left(0),lon_lat_left(1),(interp_dem_image)(input_dem_left_pix(0), input_dem_left_pix(1)));
//                  Vector3 xyz_left = input_img_geo.datum().geodetic_to_cartesian(longlat3_left);
//
//
//                  Vector2 input_dem_top_pix = input_dem_geo.lonlat_to_pixel(input_img_geo.pixel_to_lonlat(input_img_top_pix));
//                  Vector2 lon_lat_top = input_img_geo.pixel_to_lonlat(input_img_top_pix);
//                  Vector3 longlat3_top(lon_lat_top(0),lon_lat_top(1),(interp_dem_image)(input_dem_top_pix(0), input_dem_top_pix(1)));
//                  Vector3 xyz_top = input_img_geo.datum().geodetic_to_cartesian(longlat3_top);
//
//                  //Vector3 normal = computeNormalFrom3DPoints(xyz, xyz_left, xyz_top);
//                  Vector3 normal = computeNormalFrom3DPointsGeneral(xyz, xyz_left, xyz_top);
//
//                  //check for overlap between the output image and the input DEM image
//                  Vector2 overlap_pix = overlap_geo.lonlat_to_pixel(input_img_geo.pixel_to_lonlat(input_img_pix));
//                  x = (int)overlap_pix[0];
//                  y = (int)overlap_pix[1];
//
//                  //image dependent part of the code  - START
//                  PixelMask<PixelGray<uint8> > overlap_img_pixel = interp_overlap_img(x, y);
//
//                  //check for valid overlap_img coordinates
//                  //remove shadow pixels in the overlap_img.
//                  if ((x>=0) && (x < overlap_img.cols()) && (y>=0) && (y< overlap_img.rows()) && interpOverlapShadowImage(x, y) == 0){
//
//                    if ( is_valid(overlap_img_pixel) ) { //common area between input_img and overlap_img
//
//                      float overlap_img_reflectance;
//                      overlap_img_reflectance = ComputeReflectance(normal, xyz, overlap_img_params[i], globalParams);
//                      if (overlap_img_reflectance > 0){
//                         float overlap_img_error = ComputeError_Albedo((float)overlap_img_pixel, overlap_img_params[i].exposureTime,
//                                                                       (float)albedo(l, k), overlap_img_reflectance, xyz, xyz_prior);
//                         error_img(l, k) = error_img(l,k) + overlap_img_error*overlap_img_error;
//                         numSamples(l, k) = numSamples(l,k) + 1;
//                      }
//                    }//if
//                  }//if
//                 //image dependent part of the code  - START
//               }
//            }
//          }
//        }// for l
//      } // for k
//    } //for i
//
//    float l_avgError = 0;
//    int l_totalNumSamples = 0;
//    //finalize the output image; computes the standard deviation
//    for (k = 0 ; k < error_img.rows(); ++k) {
//       for (l = 0; l < error_img.cols(); ++l) {
//         if ( numSamples(l,k) ) {
//             error_img(l, k) = (float)sqrt(error_img(l, k)/numSamples(l, k));
//             l_totalNumSamples++;
//             l_avgError = l_avgError + error_img(l,k);
//          }
//       }
//    }
//
//    l_avgError = l_avgError/l_totalNumSamples;
//
//    *avgError = l_avgError;
//    *totalNumSamples = l_totalNumSamples;
//    //write the output (standard deviation of the reconstructed albedo) image
//    write_georeferenced_image(error_img_file,
//                              channel_cast<uint8>(clamp(error_img,0.0,255.0)),
//                              input_img_geo, TerminalProgressCallback("photometry}","Processing:"));
//
//
//}
