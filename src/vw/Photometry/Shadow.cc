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

#include <math.h>
#include <vw/Photometry/Reconstruct.h>


//void ComputeSaveShadowMap(std::string origfile, std::string shadowMapFile, GlobalParams globalParams)
void ComputeSaveShadowMap( ModelParams input_img_params, GlobalParams globalParams)
{

  string shadowMapFile = input_img_params.shadowFilename;
  string origfile = input_img_params.inputFilename;
  DiskImageView<PixelMask<PixelGray<uint8> > > originalImage(origfile);
  ImageView<PixelMask<PixelGray<uint8> > > shadowImage(originalImage.cols(), originalImage.rows());


  GeoReference originalGeo;
  read_georeference(originalGeo, origfile);

  for (int k=0; k < (int)originalImage.rows(); ++k) {
    for (int l=0; l < (int)originalImage.cols(); ++l) {


         Vector2 sample_pix(l,k);

         if ( is_valid(originalImage(l,k)) ) {
           shadowImage(l, k) = 0;
           if (originalImage(l, k) < globalParams.shadowThresh){
               shadowImage(l, k) = 255;
           }
         }

      }
    }

    write_georeferenced_image(shadowMapFile,
                              channel_cast<uint8>(shadowImage),
                              originalGeo, TerminalProgressCallback("{Core}","Processing:"));
}


//input_img_file is the original image
//output_img_file is the brightness compensated image file with invalid values for shadow
//this is also the filename of the output image where shadows are added
//
void AddShadows(std::string input_img_file,  std::string output_img_file, std::string shadow_file)
{
    DiskImageView<PixelMask<PixelGray<uint8> > >  input_img(input_img_file);
    GeoReference input_img_geo;
    read_georeference(input_img_geo, input_img_file);

    DiskImageView<PixelMask<PixelGray<uint8> > >  output_img(output_img_file);

    DiskImageView<PixelMask<PixelGray<uint8> > >  shadowImage(shadow_file);

    ImageView<PixelMask<PixelGray<uint8> > > r_img (input_img.cols(), input_img.rows());
    int l,k;
    //initialize  output_img, and numSamples
    for (k = 0 ; k < input_img.rows(); ++k) {
        for (l = 0; l < input_img.cols(); ++l) {
          if ( (is_valid(input_img(l,k))) && (shadowImage(l,k) == 255)  ){
              r_img(l,k) = (uint8)(input_img(l,k));
          }
          else{
              r_img(l,k) = (uint8)(output_img(l,k));
          }
        }
    }

    //write in the previous DEM
    write_georeferenced_image(output_img_file,
                              channel_cast<uint8>(r_img),
                              input_img_geo, TerminalProgressCallback("{Core}","Processing:"));

}
