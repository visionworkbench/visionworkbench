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
#include <vw/Photometry/Misc.h>

using namespace vw::photometry;
using namespace std;

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
                  float weight = ComputeLineWeightsHV(input_image_pix, input_img_params);
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

      // This is the wrong way of doing interpolation. The the Stereo module for  the right way.
      ImageViewRef<PixelMask<PixelGray<uint8> > >  interp_overlap_img = interpolate(edge_extend(overlap_img.impl(),
                                                                                    ConstantEdgeExtension()),
                                                                                    BilinearInterpolation());
      DiskImageView<PixelMask<PixelGray<uint8> > >  overlap_shadow_image(overlap_img_params[i].shadowFilename);
      ImageViewRef<PixelMask<PixelGray<uint8> > >  interp_overlap_shadow_image = interpolate(edge_extend(overlap_shadow_image.impl(),
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
              if ((x>=0) && (x < overlap_img.cols()) && (y>=0) && (y < overlap_img.rows()) && (interp_overlap_shadow_image(x, y) == 0)){

                    if ( is_valid(overlap_img_pixel) ) { //overlaping area between input_img and overlap_img

                        float overlap_img_reflectance = 1.0;
                        if (globalParams.useWeights == 0){
                              output_img(l, k) = (float)output_img(l, k) + (float)overlap_img_pixel/(overlap_img_params[i].exposureTime*overlap_img_reflectance);
                              numSamples(l, k) = numSamples(l,k) + 1;
                          }
                          else{
                             float weight = ComputeLineWeightsHV(overlap_pix, overlap_img_params[i]);
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
                      float weight = ComputeLineWeightsHV(input_image_pix, input_img_params);
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
           // This is the wrong way of doing interpolation
           ImageViewRef<PixelMask<PixelGray<uint8> > >  interp_overlap_img = interpolate(edge_extend(overlap_img.impl(),
                                                                                         ConstantEdgeExtension()),
                                                                                         BilinearInterpolation());

           DiskImageView<PixelMask<PixelGray<uint8> > >  overlap_shadow_image(overlap_img_params[i].shadowFilename);
           // This is the wrong way of doing interpolation
           ImageViewRef<PixelMask<PixelGray<uint8> > >  interp_overlap_shadow_image = interpolate(edge_extend(overlap_shadow_image.impl(),
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
                   if ((x>=0) && (x < overlap_img.cols()) && (y>=0) && (y < overlap_img.rows()) && (interp_overlap_shadow_image(x, y) == 0)){

                     if ( is_valid(overlap_img_pixel) ) { //overlaping area between input_img and overlap_img

                        float overlap_img_reflectance = 1.0;
                        if (globalParams.useWeights == 0){
                            output_img(jj,ii) = (float)output_img(jj,ii) + (float)overlap_img_pixel/(overlap_img_params[i].exposureTime*overlap_img_reflectance);
                            numSamples(l, k) = numSamples(l,k) + 1;
                         }
                         else{
                             float weight = ComputeLineWeightsHV(overlap_pix, overlap_img_params[i]);
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
                  float weight = ComputeLineWeightsHV(input_image_pix, input_img_params);
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


      DiskImageView<PixelMask<PixelGray<uint8> > >  overlap_shadow_image(overlap_img_params[i].shadowFilename);
      // This is the wrong way of doing interpolation
      ImageViewRef<PixelMask<PixelGray<uint8> > >  interp_overlap_img = interpolate(edge_extend(overlap_img.impl(),
                                                                                    ConstantEdgeExtension()),
                                                                                    BilinearInterpolation());
      // This is the wrong way of doing interpolation
      ImageViewRef<PixelMask<PixelGray<uint8> > >  interp_overlap_shadow_image = interpolate(edge_extend(overlap_shadow_image.impl(),
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
              if ((x>=0) && (x < overlap_img.cols()) && (y>=0) && (y< overlap_img.rows()) && interp_overlap_shadow_image(x, y) == 0){

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
                         float weight = ComputeLineWeightsHV(overlap_pix, overlap_img_params[i]);
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

namespace {
  
  template <class realType>
  void cropImageAndGeoRefInPlace(// Inputs
                                 Vector2 begLonLat, Vector2 endLonLat,
                                 ImageView< PixelMask<PixelGray<realType> > > & Image,
                                 cartography::GeoReference & geoRef){

    // Crop an image to the region with corners given by begLonLat,
    // endLonLat.  Adjust the GeoReference accordingly.  The point
    // begLonLat is the upper-left region corner, and endLonLat is the
    // lower-right region corner.
    
    // This is not a general purpose routine, it gets tweaked as
    // needed for the purpose of saving albedo, and that's why it does
    // not belong in a shared file.
    
    Vector2 begPixel = geoRef.lonlat_to_pixel(begLonLat);
    Vector2 endPixel = geoRef.lonlat_to_pixel(endLonLat);

    int begCol = std::max(0, (int)ceil(begPixel(0)));
    int begRow = std::max(0, (int)ceil(begPixel(1)));
    // We add 1 below to keep the bottom/right pixels which
    // intersect the boundary of the desired box.
    int endCol = std::min( Image.cols(), (int)round(endPixel(0))+1 );
    int endRow = std::min( Image.rows(), (int)round(endPixel(1))+1 );

    int numCols = std::max(endCol - begCol, 0);
    int numRows = std::max(endRow - begRow, 0);

    Image = crop(Image, BBox2i(begCol, begRow, numCols, numRows));
    geoRef = vw::cartography::crop(geoRef, begCol, begRow);

    return;
  }
  
}    

//-------------------------------------------------------------------------------
//Below are the functions for albedo reconstruction
//-------------------------------------------------------------------------------

double
vw::photometry::actOnTile(bool isLastIter, bool computeErrors,
                          std::string blankTileFile,
                          std::string DEMTileFile,   std::string albedoTileFile,
                          std::string errorTileFile, std::string weightsSumFile,
                          std::vector<ModelParams> & overlap_img_params,
                          GlobalParams globalParams,
                          phaseCoeffsData & PCD
                          ){

  // Perform one of the several calculations for the given tile:

  // 1. Sum the weights at each pixel in the tile (skip weights at which
  //    the image is below threshold or the reflectance is invalid)
  // 2. Initialize the albedo tile
  // 3. Update the albedo tile
  // 4. update the phase coefficients
  // 5. Compute the albedo errors.
  
  // We use double precision for the numerical computations, even
  // though some inputs are float. This is more important when the
  // albedo is updated, and less so when it is initialized.

  // If we initialize the albedo for the given tile, read all the
  // images which overlap with it, and do a weighted average involving
  // the images, weights, reflectance, and exposure.

  bool useReflectance = (globalParams.reflectanceType != NO_REFL);
  if ( globalParams.updateTilePhaseCoeffs && !useReflectance ){
    std::cout << "ERROR: Cannot update the phase coefficients if "
              << "we don't use perform the reflectance computation" << std::endl;
    exit(1);
  }
  if ( (globalParams.updateTilePhaseCoeffs || globalParams.computeWeightsSum) && isLastIter){
    std::cout << "ERROR: Cannot update sum the weights/compute the phase coefficients "
              << "at the last iteration" << std::endl;
    exit(1);
  }

  // The cost function is the weighted sum of squares of errors over
  // all the pixels in the given tile (excluding the padded region at
  // the boundary).
  double costFunVal = 0.0;

  DiskImageView<PixelMask<PixelGray<uint8> > >  blankTile(blankTileFile);
  std::cout << "Reading " << blankTileFile << std::endl;
  ImageView<PixelMask<PixelGray<float> > > albedoTile(blankTile.cols(), blankTile.rows());
  GeoReference albedoTile_geo;
  read_georeference(albedoTile_geo, blankTileFile);

  ImageView<float> weightsSum;
  if (globalParams.useNormalizedWeights && !globalParams.computeWeightsSum){
    // The weights were already summed up. Read them.
    weightsSum = copy(DiskImageView<float>(weightsSumFile));
  }

  ImageView<Vector3> dem_xyz, surface_normal;
  if (useReflectance){
    // Use a block to de-allocate DEMTile as soon as it is no longer needed.
    DiskImageView<PixelGray<float> > DEMTile(DEMTileFile);
    std::cout << "Reading file: "<< DEMTileFile << std::endl;
    GeoReference DEMGeo;
    read_georeference(DEMGeo, DEMTileFile);
    float noDEMDataValue;
    if ( !readNoDEMDataVal(DEMTileFile, noDEMDataValue)){
      std::cerr << "ERROR: Could not read the NoData Value from " << DEMTileFile << std::endl;
      exit(1);
    }
      
    vw::photometry::computeXYZandSurfaceNormal(DEMTile.impl(), DEMGeo, noDEMDataValue,
                                               dem_xyz, surface_normal
                                               );
  }
    
  ImageView<float> norm(albedoTile.cols(), albedoTile.rows());
  GeoReference weightsSum_geo;
  if (globalParams.useNormalizedWeights && globalParams.computeWeightsSum){
    read_georeference(weightsSum_geo, blankTileFile);
  }

  //initialize  albedoTile
  for (int k = 0 ; k < albedoTile.rows(); ++k) {
    for (int l = 0; l < albedoTile.cols(); ++l) {
      albedoTile(l, k) = 0;
      albedoTile(l, k).invalidate(); // a pixel is invalid unless proven otherwise
      norm(l, k) = 0;
    }
  }

  // If we do an update or compute the errors, we need to read in the albedo before the update.
  ImageView<PixelMask<PixelGray<float> > > inputAlbedoTile;
  bool willReadInputTile = (!globalParams.computeWeightsSum && !globalParams.initAlbedo);  
  if (willReadInputTile){
    std::ifstream aFile(albedoTileFile.c_str());
    if (!aFile) return 0;
    std::cout << "Reading " << albedoTileFile << std::endl;
    // Copy to a float, to be able to do accurate calculations
    inputAlbedoTile = copy(DiskImageView<PixelMask<PixelGray<uint8> > >(albedoTileFile));

    // Sanity check: if the user is not careful, the padding may have been stripped
    // from the albedo tile by now, which of course would cause a size mis-match.
    if (useReflectance){
      DiskImageView<PixelGray<float> > DEMTile(DEMTileFile);
      if (DEMTile.rows() != inputAlbedoTile.rows() || DEMTile.cols() != inputAlbedoTile.cols()){
        std::cout << "ERROR: We expect the DEM and albedo tiles to have the same number "
                  << "of rows/columns." << std::endl;
        exit(1);
      }
    }
  }
  
  ImageView<PixelMask<PixelGray<float> > > errorTile;
  GeoReference errorTile_geo;
  if (computeErrors){
    read_georeference(errorTile_geo, blankTileFile);
    errorTile.set_size(blankTile.cols(), blankTile.rows());
    for (int k = 0 ; k < errorTile.rows(); ++k) {
      for (int l = 0; l < errorTile.cols(); ++l) {
        errorTile(l, k) = 0;
        errorTile(l, k).invalidate(); // a pixel is invalid unless proven otherwise
      }
    }
  }
  
  // The lon-lat coordinates of tile corners without padding.
  double min_tile_x, max_tile_x, min_tile_y, max_tile_y;
  getTileCornersWithoutPadding(// Inputs
                               albedoTile.cols(), albedoTile.rows(), albedoTile_geo,  
                               globalParams.tileSize, globalParams.pixelPadding,  
                               // Outputs
                               min_tile_x, max_tile_x,  
                               min_tile_y, max_tile_y
                               );

  // These are need for finding the components of the phase coefficients per tile
  PCD.phaseCoeffA1_num = 0.0;
  PCD.phaseCoeffA1_den = 0.0;
  PCD.phaseCoeffA2_num = 0.0;
  PCD.phaseCoeffA2_den = 0.0;

  for (int i = 0; i < (int)overlap_img_params.size(); i++){

    printf("overlap_img = %s with %s\n", albedoTileFile.c_str(), overlap_img_params[i].inputFilename.c_str());
    //system("echo date1 is $(date)");
    
    // Read the weights only if really needed
    bool useWeights = (globalParams.useWeights != 0);
    if (useWeights){
      bool useTiles = true;
      ReadWeightsParamsFromFile(useTiles, &overlap_img_params[i]);
    }

    // The georeference for the entire overlap image
    GeoReference overlap_geo_orig;
    read_georeference(overlap_geo_orig, overlap_img_params[i].inputFilename);

    // Skip the current tile if there are no valid image pixels inside
    // of it. This optimization is necessary, but is a big of a hack
    // since it won't apply if we don't use weights (which is rare).
    if (useWeights){
      bool hasGoodPixels = false;
      for (int k = 0; k < albedoTile.rows(); ++k) {
        if (hasGoodPixels) break;          
        for (int l = 0; l < albedoTile.cols(); ++l) {
          // We need the pixel coordinates in the ENTIRE image, as opposed to its coordinates
          // in the subimage, for the purpose of computing the weight.
          Vector2 albedoTile_pix(l,k);
          Vector2 lon_lat = albedoTile_geo.pixel_to_lonlat(albedoTile_pix);
          Vector2 overlap_pix_orig = overlap_geo_orig.lonlat_to_pixel(lon_lat);
          double weight = ComputeLineWeightsHV(overlap_pix_orig, overlap_img_params[i]);
          if (weight != 0.0){
            hasGoodPixels = true;
            break;
          }
        }
      }
      if (!hasGoodPixels) continue;
    }
    
    // The input DRG must be uint8
    enforceUint8Img(overlap_img_params[i].inputFilename);

    // Read only the portion of the overlap image which intersects
    // the current tile to save on memory.  Create the appropriate
    // georeference for it.
    ImageView<PixelMask<PixelGray<float> > > overlap_img;
    GeoReference overlap_geo;
    Vector2 begTileLonLat = albedoTile_geo.pixel_to_lonlat(Vector2(0, 0));
    Vector2 endTileLonLat = albedoTile_geo.pixel_to_lonlat(Vector2(albedoTile.cols()-1, albedoTile.rows()-1));
    bool success = getSubImageWithMargin< PixelMask<PixelGray<uint8> >, PixelMask<PixelGray<float> > >
      (// Inputs
       begTileLonLat, endTileLonLat, overlap_img_params[i].inputFilename,  
       // Outputs
       overlap_img, overlap_geo
       );
    if (!success) continue;

    InterpolationView<EdgeExtensionView<ImageView<PixelMask<PixelGray<float> > >, ConstantEdgeExtension>, BilinearInterpolation>
      interp_overlap_img = interpolate(overlap_img, BilinearInterpolation(), ConstantEdgeExtension());

    ImageView<PixelMask<PixelGray<float> > > overlapReflectance;
    bool savePhaseAngle = globalParams.updateTilePhaseCoeffs; // Need the phase angle only when we update the phase coeffs
    ImageView<PixelMask<PixelGray<float> > > phaseAngle;
    if (useReflectance){
      computeReflectanceAux(dem_xyz, surface_normal,  
                            overlap_img_params[i],
                            globalParams,  
                            overlapReflectance, // output
                            savePhaseAngle,
                            phaseAngle          // output
                            );
    }else{
      // The reflectance is set to 1.
      overlapReflectance.set_size(albedoTile.cols(), albedoTile.rows());
      for (int row = 0; row < (int)albedoTile.rows(); row++){
        for (int col = 0; col < (int)albedoTile.cols(); col++){
          overlapReflectance(col, row) = 1.0;
          overlapReflectance(col, row).validate();
        }
      }
    }

    for (int k = 0; k < albedoTile.rows(); ++k) {
      for (int l = 0; l < albedoTile.cols(); ++l) {

        Vector2 albedoTile_pix(l,k);

        Vector2 lon_lat = albedoTile_geo.pixel_to_lonlat(albedoTile_pix);
              
        if (globalParams.useNormalizedWeights && !globalParams.computeWeightsSum && weightsSum(l, k) == 0) continue;
        
        double weight = 1.0;
        if (useWeights){
          // We need the pixel coordinates in the ENTIRE image, as opposed to its coordinates
          // in the subimage, for the purpose of computing the weight.
          Vector2 overlap_pix_orig = overlap_geo_orig.lonlat_to_pixel(lon_lat);
          weight = ComputeLineWeightsHV(overlap_pix_orig, overlap_img_params[i]);
        }
        if (weight == 0.0) continue;

        // Normalize to ensure that all weights for all images at given pixel add up to 1
        if (globalParams.useNormalizedWeights && !globalParams.computeWeightsSum) weight /= weightsSum(l, k);

        //check for overlap between the output image and the input DEM image
        Vector2 overlap_pix = overlap_geo.lonlat_to_pixel(lon_lat);
        double overlap_x    = overlap_pix[0];
        double overlap_y    = overlap_pix[1];

        //image dependent part of the code  - START
        
        // Check for valid overlap_img pixels
        int j0 = (int)floor(overlap_x), j1 = (int)ceil(overlap_x);
        int i0 = (int)floor(overlap_y), i1 = (int)ceil(overlap_y);
        // Note: Below must do overlap_x <= overlap_img.cols()-1
        // rather than overlap_x < overlap_img.cols().
        bool isValidOverlap = ( (overlap_x >= 0) && (overlap_x <= overlap_img.cols()-1 )        &&
                                (overlap_y >= 0) && (overlap_y <= overlap_img.rows()-1 )        &&
                                is_valid(overlap_img(j0, i0))  && is_valid(overlap_img(j0, i1)) &&
                                is_valid(overlap_img(j1, i0))  && is_valid(overlap_img(j1, i1)) 
                                );
        if (!isValidOverlap) continue;
        
        // We have a valid pixel.
        albedoTile(l,k).validate();
        if (computeErrors) errorTile(l,k).validate();
        
        double R = (double)overlapReflectance(l, k);
        if (R == 0.0) continue; // no reflectance data

        double exposureRefl = overlap_img_params[i].exposureTime*R;
        
        // Check if the image is above the shadow threshold
        double t = getShadowThresh(globalParams, exposureRefl);
        bool isBlack = !((double)overlap_img(j0, i0) >= t &&
                         (double)overlap_img(j0, i1) >= t &&
                         (double)overlap_img(j1, i0) >= t &&
                         (double)overlap_img(j1, i1) >= t
                         );
        if (isBlack) continue; // No image data

        // See an explanation of the logic below in reconstruct.cc.
        if (globalParams.forceMosaic){
          overlap_img_params[i].exposureTime = globalParams.TRConst;
          R = 1.0;
          exposureRefl = overlap_img_params[i].exposureTime*R;
        }
        
        PixelMask<PixelGray<float> > overlap_img_pixel = interp_overlap_img(overlap_x, overlap_y);
        
        if (globalParams.computeWeightsSum){

          if (globalParams.useNormalizedWeights) norm(l,k) = norm(l,k) + weight;

        }else if (globalParams.initAlbedo){
          
          //New averaging
          albedoTile(l, k) = (double)albedoTile(l, k) + ((double)overlap_img_pixel)*exposureRefl*weight;
          norm(l,k)        = norm(l,k) + exposureRefl*exposureRefl*weight;
          // Old averaging
          //albedoTile(l, k) = (float)albedoTile(l, k) + ((float)overlap_img_pixel*weight)/exposureRefl;
          //norm(l,k) = norm(l,k) + weight;
          
        }else if (willReadInputTile && is_valid(inputAlbedoTile(l, k))){

          double lx = lon_lat(0), ly = lon_lat(1);
          bool isInTile = (min_tile_x <= lx && lx < max_tile_x && min_tile_y <= ly && ly < max_tile_y);

          // Update albedo or update phase coefficients or compute errors
          double diff = (double)overlap_img_pixel - ((double)inputAlbedoTile(l, k))*exposureRefl;
          if (globalParams.updateAlbedo){
            // Update the albedo
            albedoTile(l, k) = (double)albedoTile(l, k) + diff*exposureRefl*weight;
            norm(l,k)        = norm(l,k) + exposureRefl*exposureRefl*weight;
          }else if (globalParams.updateTilePhaseCoeffs && isInTile){
            // Update the phase coefficients components for the current tile
            double alpha   = phaseAngle(l, k);
            double e       = exp(-globalParams.phaseCoeffA1*alpha);
            // The derivative of A*T*R in respect to the phase coeffs,
            // where R = (exp(-A1*alpha) + A2)*lambertian
            double derivA2 = ((double)inputAlbedoTile(l, k))*exposureRefl/(e + globalParams.phaseCoeffA2);
            double derivA1 = derivA2*(-alpha*e);
            PCD.phaseCoeffA1_num += diff*derivA1*weight;
            PCD.phaseCoeffA1_den += derivA1*derivA1*weight;
            PCD.phaseCoeffA2_num += diff*derivA2*weight;
            PCD.phaseCoeffA2_den += derivA2*derivA2*weight;
          }else if (computeErrors){
            // Compute the errors
            errorTile(l, k) = (double)errorTile(l, k) + (diff/exposureRefl)*(diff/exposureRefl)*weight;
            norm(l,k)       = norm(l,k) + weight;
          }
                      
          if (isInTile){
            // Accumulate the cost function only for pixels in the tile proper rather than
            // in the padded region
            costFunVal += diff*diff*weight; // Weighted sum of squares
          }
        }

      }
                  
      //image dependent part of the code  - END
    }

    //std::cout << "sun and spacecraft: " << overlap_img_params[i].inputFilename  << ' '
    //          << overlap_img_params[i].sunPosition  << ' ' << overlap_img_params[i].spacecraftPosition
    //          << std::endl;
    //system("echo albedo top is $(top -u $(whoami) -b -n 1|grep lt-reconstruct)");
  }

  // Nothing else to do if all we care about is updating the phase coefficients
  if (globalParams.updateTilePhaseCoeffs) return costFunVal;

  if (globalParams.useNormalizedWeights && globalParams.computeWeightsSum){
    std::cout << "Writing: " << weightsSumFile << std::endl;
    write_georeferenced_image(weightsSumFile,
                              norm,
                              weightsSum_geo, TerminalProgressCallback("{Core}","Processing:"));
    
    return 0;
  }
  
  if (computeErrors){
    // If we computed the errors, we did not update the albedo, as such, the albedo
    // must be set to the input albedo.
    albedoTile = copy(inputAlbedoTile);
  }
    
  //compute the mean albedo value
  int numValid = 0;
  for (int k = 0 ; k < albedoTile.rows(); ++k) {
    for (int l = 0; l < albedoTile.cols(); ++l) {

      // When updating the albedo, if the input albedo is invalid, then the updated
      // albedo is also invalid.
      if (willReadInputTile && globalParams.updateAlbedo && !is_valid(inputAlbedoTile(l, k))){
        albedoTile(l,k).invalidate();
      }
      if (willReadInputTile && computeErrors && !is_valid(inputAlbedoTile(l, k))){
        errorTile(l,k).invalidate();
      }
      
      if ( (float)albedoTile(l,k) == 0 ) albedoTile(l,k).invalidate(); // temporary!!!
        
      if ( (is_valid(albedoTile(l,k))) && (norm(l, k) != 0) ) {
          
        if (globalParams.initAlbedo){
          // Init tile
          albedoTile(l, k) = albedoTile(l, k)/norm(l,k);
        }else if (willReadInputTile && globalParams.updateAlbedo){
          // Update the tile
          albedoTile(l, k) = inputAlbedoTile(l, k) + albedoTile(l, k)/norm(l,k);
        }else{
          // Compute the errors
          errorTile(l, k) = errorTile(l, k)/norm(l,k);
        }
        numValid++;
      }
        
    }
  }

  printf("numValid = %d, total = %d\n", numValid, albedoTile.rows()*albedoTile.cols());
  
  if (isLastIter){

    // The true albedo is in fact A0*( exp(-A1*alpha) + A2 ), for alpha = 0.
    // The quantity A0 is what we computed so far, so we must multiply by the remaining factor.
    if (globalParams.postScaleAlbedo){
      for (int k = 0; k < albedoTile.rows(); ++k) {
        for (int l = 0; l < albedoTile.cols(); ++l) {
          if ( is_valid(albedoTile(l, k)) ) albedoTile(l, k)
                                              = (double)albedoTile(l, k)*(1 + globalParams.phaseCoeffA2);
        }
      }
    }
    
    // Remove the temporary work padding if this is the last albedo update iteration.
    cropImageAndGeoRefInPlace(Vector2(min_tile_x, max_tile_y),
                              Vector2(max_tile_x, min_tile_y),
                              albedoTile, albedoTile_geo // inputs-outputs
                              );

  }
    
  // Write the albedo tile. Note that we write the albedo tile even if we are in the mode in which
  // we compute the error. That because, if the error computation is the last iteration, this is the time
  // at which to crop the albedo tile itself (we could not have cropped it at a previous step,
  // because we cannot compute the error with a cropped albedo tile).
  // After cropping the image, if it is made up of only empty pixels, then don't save it.
  numValid = 0;
  for (int k = 0; k < albedoTile.rows(); ++k) {
    for (int l = 0; l < albedoTile.cols(); ++l) {
      if ( (double)albedoTile(l, k) > 0 ) numValid++;
    }
  }
  if (numValid > 0){
  std::cout << "Writing: " << albedoTileFile << std::endl;
  write_georeferenced_image(albedoTileFile,
                            channel_cast<uint8>(clamp(albedoTile, 0.0, 255.0)),
                            albedoTile_geo, TerminalProgressCallback("{Core}","Processing:"));
  }
  
  if (numValid > 0 && computeErrors){
    
    if (isLastIter){
      // Strip the padding before saving the error to disk
      cropImageAndGeoRefInPlace(Vector2(min_tile_x, max_tile_y),
                                Vector2(max_tile_x, min_tile_y),
                                errorTile, errorTile_geo // inputs-outputs
                                );
    }
    
    std::cout << "Writing: " << errorTileFile << std::endl;
    // The scaling below is pretty arbitrary
    write_georeferenced_image(errorTileFile,
                              channel_cast<uint8>(clamp(errorTile/10.0, 0.0, 255.0)),
                              errorTile_geo, TerminalProgressCallback("{Core}","Processing:"));
  }
  
  // Return the cost function value
  return costFunVal;
}

void vw::photometry::AppendCostFunToFile(double costFunVal, std::string fileName)
{
  // Append the current cost function to the file.
  FILE *fp;
  fp = fopen((char*)fileName.c_str(), "a");
  std::cout << "Writing " << fileName << std::endl;
  fprintf(fp, "%0.16g\n", costFunVal);
  fclose(fp);
}


//initializes the albedo mosaic
//TO DO: build the version that does the blockwise processing to deal with large scale images
void
vw::photometry::InitAlbedoMosaic(ModelParams input_img_params,
                                 std::vector<ModelParams> overlap_img_params,
                                 GlobalParams globalParams) {
    
    int i, l, k;
    std::string input_img_file = input_img_params.inputFilename;
    std::string reflectance_file = input_img_params.reliefFilename;
    std::string shadow_file = input_img_params.shadowFilename;
    std::string output_img_file = input_img_params.outputFilename;

    // Read the uint8 data from disk. Copy it to float, to be able to interpolate accurately.
    ImageView<PixelMask<PixelGray<float> > > input_img
      = copy(DiskImageView<PixelMask<PixelGray<uint8> > >(input_img_file));
    GeoReference input_img_geo;
    read_georeference(input_img_geo, input_img_file);

    DiskImageView<PixelMask<PixelGray<uint8> > > shadowImage(shadow_file);
    
    DiskImageView<PixelMask<PixelGray<float> > > reflectance_image(reflectance_file);
    GeoReference reflectance_geo;
    read_georeference(reflectance_geo, reflectance_file);
    InterpolationView<EdgeExtensionView<DiskImageView<PixelMask<PixelGray<float> > >, ConstantEdgeExtension>, BilinearInterpolation>
      interp_reflectance_image = interpolate(reflectance_image.impl(), BilinearInterpolation(), ConstantEdgeExtension());
    // Wrong below
    //ImageViewRef<PixelMask<PixelGray<float> > >  interp_reflectance_image
    //  = interpolate(edge_extend(reflectance_image.impl(), ConstantEdgeExtension()),
    //                BilinearInterpolation());

    ImageView<PixelMask<PixelGray<float> > > output_img (input_img.cols(), input_img.rows());
    ImageView<PixelGray<int> > numSamples(input_img.cols(), input_img.rows());
    ImageView<PixelGray<float> > norm(input_img.cols(), input_img.rows());
    //ImageView<PixelGray<float> > weights_img(input_img.cols(), input_img.rows());
 
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
                          float weight       = ComputeLineWeightsHV(input_image_pix, input_img_params);
                          float exposureRefl = input_img_params.exposureTime*input_img_reflectance;

                          // New averaging
                          output_img(l, k)  = ((float)input_img(l,k))*exposureRefl*weight;
                          norm(l, k)        = exposureRefl*exposureRefl*weight;

                          // Old averaging
                          //output_img(l, k) = ((float)input_img(l,k)*weight)/exposureRefl;
                          //norm(l, k) = weight;
                          
                          numSamples(l, k)  = 1;
                          //weights_img(l, k) = weight;

                          
                      }
                  }
                }
              }
           }
        }
    }

    //for (k = 0 ; k < input_img.rows(); ++k) {
    //  for (l = 0; l < input_img.cols(); ++l) {
    //    float val = weights_img(l, k)[0];
    //    weights_img(l, k) = std::min(255.0, val*155.0);
    //  }
    //}

    for (i = 0; i < (int)overlap_img_params.size(); i++){
      
      printf("overlap_img = %s\n", overlap_img_params[i].inputFilename.c_str());

      // Read the images from disk and copy them into float, to get higher precision
      // when interpolating.
      ImageView<PixelMask<PixelGray<float> > >  overlap_img
        = copy(DiskImageView<PixelMask<PixelGray<uint8> > >(overlap_img_params[i].inputFilename));
      InterpolationView<EdgeExtensionView<ImageView<PixelMask<PixelGray<float> > >, ConstantEdgeExtension>, BilinearInterpolation>
        interp_overlap_img = interpolate(overlap_img, BilinearInterpolation(), ConstantEdgeExtension());
      //DiskImageView<PixelMask<PixelGray<uint8> > >  overlap_img(overlap_img_params[i].inputFilename);
      //ImageViewRef<PixelMask<PixelGray<uint8> > >  interp_overlap_img = interpolate(edge_extend(overlap_img.impl(),
      //                                                                              ConstantEdgeExtension()),
      //                                                                              BilinearInterpolation());
      
      GeoReference overlap_geo;
      read_georeference(overlap_geo, overlap_img_params[i].inputFilename);

    std::string overlap_reflectance_file = overlap_img_params[i].reliefFilename;
    DiskImageView<PixelMask<PixelGray<float> > >   overlap_reflectance_image(overlap_reflectance_file);
    InterpolationView<EdgeExtensionView<DiskImageView<PixelMask<PixelGray<float> > >, ConstantEdgeExtension>, BilinearInterpolation>
      interp_overlap_reflectance_image
      = interpolate(overlap_reflectance_image, BilinearInterpolation(), ConstantEdgeExtension());
    //ImageViewRef<PixelMask<PixelGray<float> > >  interp_overlap_reflectance_image
    //  = interpolate(edge_extend(overlap_reflectance_image.impl(),
    // ConstantEdgeExtension()), BilinearInterpolation());
    GeoReference overlap_reflectance_geo;
    read_georeference(overlap_reflectance_geo, overlap_reflectance_file);

    DiskImageView<PixelMask<PixelGray<uint8> > > overlap_shadow_image(overlap_img_params[i].shadowFilename);
    InterpolationView<EdgeExtensionView<DiskImageView<PixelMask<PixelGray<uint8> > >, ConstantEdgeExtension>, BilinearInterpolation>
    interp_overlap_shadow_image = interpolate(overlap_shadow_image, BilinearInterpolation(), ConstantEdgeExtension());
    //ImageViewRef<PixelMask<PixelGray<uint8> > >  interp_overlap_shadow_image
    //  = interpolate(edge_extend(overlap_shadow_image.impl(), ConstantEdgeExtension()), BilinearInterpolation());
    
      for (k = 0 ; k < input_img.rows(); ++k) {
        for (l = 0; l < input_img.cols(); ++l) {
          
          Vector2 input_img_pix(l,k);
          
          if ( is_valid(input_img(l,k)) ) {

              //get the corresponding overlap_reflectance value
              Vector2 lon_lat = input_img_geo.pixel_to_lonlat(input_img_pix);
	     
              Vector2 overlap_reflectance_pix = overlap_reflectance_geo.lonlat_to_pixel(input_img_geo.pixel_to_lonlat(input_img_pix));
              float x = overlap_reflectance_pix[0];
              float y = overlap_reflectance_pix[1];

              //check for valid overlap_reflectance_image coordinates
              if ((x>=0) && (x < overlap_reflectance_image.cols()) && (y>=0) && (y< overlap_reflectance_image.rows())){
	
                if (is_valid(interp_overlap_reflectance_image(x,y))){
                 
                  //check for overlap between the output image and the input DEM image
                  Vector2 overlap_pix = overlap_geo.lonlat_to_pixel(input_img_geo.pixel_to_lonlat(input_img_pix));
                  float overlap_x = overlap_pix[0];
                  float overlap_y = overlap_pix[1];

                  //image dependent part of the code  - START
                  PixelMask<PixelGray<float> > overlap_img_pixel = interp_overlap_img(overlap_x, overlap_y);

                  //check for valid overlap_img coordinates
                  //TO DO: remove shadow pixels in the overlap_img.
                  if ((overlap_x>=0) && (overlap_x < overlap_img.cols()) && (overlap_y >= 0) && (overlap_y< overlap_img.rows()) && (interp_overlap_shadow_image(x, y) == 0)){

                    if ( is_valid(overlap_img_pixel) ) { //common area between input_img and overlap_img

                      float overlap_img_reflectance = interp_overlap_reflectance_image(x,y);

                      if (overlap_img_reflectance != 0.0){
                          if (globalParams.useWeights == 0){
                              output_img(l, k) = (float)output_img(l, k) + (float)overlap_img_pixel/(overlap_img_params[i].exposureTime*overlap_img_reflectance);
                              numSamples(l, k) = numSamples(l,k) + 1;
                          }
                          else{
                            float weight       = ComputeLineWeightsHV(overlap_pix, overlap_img_params[i]);
                            float exposureRefl = overlap_img_params[i].exposureTime*overlap_img_reflectance;

                            //New averaging
                            output_img(l, k) = (float)output_img(l, k) + (float)overlap_img_pixel*exposureRefl*weight;
                            norm(l,k)        = norm(l,k) + exposureRefl*exposureRefl*weight;
                            
                            // Old averaging
                            //output_img(l, k) = (float)output_img(l, k)
                            // + ((float)overlap_img_pixel*weight)/ exposureRefl;
                            //norm(l,k) = norm(l,k) + weight;
                            
                            numSamples(l, k) = numSamples(l,k) + 1;


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
    std::cout << "Writing: " << output_img_file << std::endl;
    write_georeferenced_image(output_img_file,
                              channel_cast<uint8>(clamp(output_img,0.0,255.0)),
                              input_img_geo, TerminalProgressCallback("{Core}","Processing:"));
    //write_georeferenced_image(output_img_file,
    //                          output_img,
    //                          input_img_geo, TerminalProgressCallback("{Core}","Processing:"));
    //std::string weights_file = output_img_file, str2 = "DRG";
    //weights_file.replace(weights_file.find(str2),str2.length(),"wt");
    //std::cout << "Writing the weights to "  <<  weights_file << std::endl;
    //write_georeferenced_image(weights_file,
    //                          channel_cast<uint8>(clamp(weights_img,0.0,255.0)),
    //                          input_img_geo, TerminalProgressCallback("{Core}","Processing:"));
    

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

    // This is the wrong way of doing interpolation
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
                  float input_img_reflectance, phaseAngle;
                  input_img_reflectance = ComputeReflectance(normal, xyz, input_img_params, globalParams, phaseAngle);

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
                         float weight = ComputeLineWeightsHV(input_image_pix, input_img_params);
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

      DiskImageView<PixelMask<PixelGray<uint8> > >  overlap_shadow_image(overlap_img_params[i].shadowFilename);

      // This is the wrong way of doing interpolation
      ImageViewRef<PixelMask<PixelGray<uint8> > >  interp_overlap_img = interpolate(edge_extend(overlap_img.impl(),
                                                                                    ConstantEdgeExtension()),
                                                                                    BilinearInterpolation());

      // This is the wrong way of doing interpolation
      ImageViewRef<PixelMask<PixelGray<uint8> > >  interp_overlap_shadow_image = interpolate(edge_extend(overlap_shadow_image.impl(),
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
                  if ((x>=0) && (x < overlap_img.cols()) && (y>=0) && (y< overlap_img.rows()) && interp_overlap_shadow_image(x, y) == 0){

                    if ( is_valid(overlap_img_pixel) ) { //common area between input_img and overlap_img

                      float overlap_img_reflectance, phaseAngle;
                      overlap_img_reflectance = ComputeReflectance(normal, xyz, overlap_img_params[i],
                                                                   globalParams, phaseAngle);
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
                            float weight = ComputeLineWeightsHV(overlap_pix, overlap_img_params[i]);
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
