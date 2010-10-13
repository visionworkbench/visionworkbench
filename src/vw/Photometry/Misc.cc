// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <iostream>
#include <string>
#include <complex>
#include <vector>
#include <limits>
#include <time.h>
#include <unistd.h>
#include <proj_api.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include <boost/filesystem/operations.hpp>
namespace fs = boost::filesystem;

#include <vw/Core.h>
#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Cartography.h>
#include <vw/Math.h>
using namespace vw;
using namespace vw::cartography;

#include <vw/Photometry/Reconstruct.h>
#include <vw/Photometry/Misc.h>
#include <vw/Photometry/Shape.h>
using namespace vw::photometry;


//upsamples a geo referenced tiff image by four- used in sfs

void upsample_image(std::string output_file, std::string input_file, int upsampleFactor) {
  GeoReference geo;
  read_georeference(geo, input_file);
  DiskImageView<PixelGray<float> >   image(input_file);

  int cols = (image.cols())*upsampleFactor, rows = (image.rows())*upsampleFactor;
  ImageView<PixelGray<float> >  tm_image(cols, rows);

  ImageViewRef<PixelGray<float> >   interp = interpolate(edge_extend(image.impl(),
                                                                     ConstantEdgeExtension()),
                                                         BilinearInterpolation());

  int x, y;

  for (x=0; x<cols; ++x){
    for (y=0; y<rows; ++y){
      //if ( is_valid(image(2*x,2*y)) || is_valid(image(2*x+1,2*y)) || is_valid(image(2*x,2*y+1)) || is_valid(image(2*x+1,2*y+1)) ) {
      //if (is_valid(image(2*x,2*y)) && is_valid(image(2*x+1,2*y)) && is_valid(image(2*x,2*y+1)) && is_valid(image(2*x+1,2*y+1)) ) {
      float xx = x/upsampleFactor;
      float yy = y/upsampleFactor;

      if ( interp(x, y) != -10000 ){
           tm_image(x,y) = interp(x, y);
      }
      else{
           tm_image(x,y) = -10000;
      }
    }
  }

  Matrix<double> H = geo.transform();
  H(0,0) /= upsampleFactor;
  H(1,1) /= upsampleFactor;
  geo.set_transform(H);

  write_georeferenced_image(output_file, tm_image, geo, TerminalProgressCallback("photometry","Processing:"));
}


//subsamples a geo referenced tiff image by two
void subsample_image(std::string output_file, std::string input_file) {
  GeoReference geo;
  read_georeference(geo, input_file);
  DiskImageView<PixelMask<PixelGray<uint8> > >  image(input_file);
  int cols = (image.cols()+1)/2, rows = (image.rows()+1)/2;
  ImageView<PixelMask<PixelGray<uint8> > > tm_image(cols, rows);


  ImageViewRef<PixelMask<PixelGray<uint8> > >  interp = interpolate(edge_extend(image.impl(),
                                                                                ConstantEdgeExtension()),
                                                                    BilinearInterpolation());

  int x, y;

  for (x=0; x<cols; ++x){
    for (y=0; y<rows; ++y){
      //if ( is_valid(image(2*x,2*y)) || is_valid(image(2*x+1,2*y)) || is_valid(image(2*x,2*y+1)) || is_valid(image(2*x+1,2*y+1)) ) {
      //if (is_valid(image(2*x,2*y)) && is_valid(image(2*x+1,2*y)) && is_valid(image(2*x,2*y+1)) && is_valid(image(2*x+1,2*y+1)) ) {
      if ( is_valid(interp(2*x+0.5, 2*y+0.5)) ){
        tm_image(x,y) = interp(2*x+0.5, 2*y+0.5);
      }
      else{
        tm_image(x,y).invalidate();
      }
    }
  }

  Matrix<double> H = geo.transform();
  H(0,0) *= 2;
  H(1,1) *= 2;
  geo.set_transform(H);

  write_georeferenced_image(output_file, tm_image, geo, TerminalProgressCallback("photometry","Processing:"));
}


// Given two images and two georeferences, this function picks a set
// of matching pixel samples between the two images.  It rejects
// pixels that are not valid, and it should probably also reject
// pixels that are near saturation (though it does not yet!).
template <class ViewT>
std::vector<Vector4> sample_images(ImageViewBase<ViewT> const& image1,
                                   ImageViewBase<ViewT> const& image2,
                                   GeoReference const& geo1,
                                   GeoReference const& geo2,
                                   int num_samples,
                                   std::string const& DEM_file,
                                   std::vector<Vector3> *normalArray,
                                   std::vector<Vector3> *xyzArray ) {
  int sample = 0;
  int numtries = 0;
  std::vector<Vector4> result;

  // Random numbers
  srandom((unsigned int) clock());

  ImageViewRef<typename ViewT::pixel_type> interp_image1 = interpolate(edge_extend(image1.impl(),
                                                                       ConstantEdgeExtension()),
                                                                       BilinearInterpolation());
  ImageViewRef<typename ViewT::pixel_type> interp_image2 = interpolate(edge_extend(image2.impl(),
                                                                       ConstantEdgeExtension()),
                                                                       BilinearInterpolation());

  // This block of code samples the images comprehensively, adding a
  // sample pair for every valid pixel in interp_image1.
  // for (unsigned j=0; j < interp_image1.rows(); ++j) {
  //   for (unsigned i=0; i < interp_image1.cols(); ++i) {
  //     Vector2 sample_pix1(i,j);
  //     Vector2 sample_pix2 = geo2.lonlat_to_pixel(geo1.pixel_to_lonlat(sample_pix1));

  //     // Check to see whether these pixels are valid
  //     typename ViewT::pixel_type pix1 = interp_image1(sample_pix1[0], sample_pix1[1]);
  //     typename ViewT::pixel_type pix2 = interp_image2(sample_pix2[0], sample_pix2[1]);
  //     if ( is_valid(pix1) && is_valid(pix2) &&
  //          pix1[0] > 10 && pix1[0] < 245 &&
  //          pix2[0] > 10 && pix2[0] < 245 ) {
  //       result.push_back(Vector2(pix1[0],pix2[0]));
  //       ++sample;
  //       //        std::cout << result[result.size()-1][0] << " " << result[result.size()-1][1] << "\n";
  //     }
  //     ++numtries;
  //   }
  // }


  //added by Ara to support DEMs - START
  DiskImageView<PixelGray<float> >  dem_image(DEM_file);
  GeoReference GR;
  read_georeference(GR, DEM_file);
  //added by Ara to support DEMs - END

  // This block of code samples the images randomly, gathering up to
  // num_samples samples from the images.
  while (sample < num_samples && numtries < num_samples*10) {

    Vector2 sample_pix1(float(random())/RAND_MAX * image1.impl().cols(),
                        float(random())/RAND_MAX * image1.impl().rows());
    Vector2 sample_pix2 = geo2.lonlat_to_pixel(geo1.pixel_to_lonlat(sample_pix1));

    Vector2 sample_pix_dem = GR.lonlat_to_pixel(geo1.pixel_to_lonlat(sample_pix1));

    Vector2 lonlat = geo1.pixel_to_lonlat(sample_pix1);

    // Check to see whether these pixels are valid
    typename ViewT::pixel_type pix1 = interp_image1(sample_pix1[0], sample_pix1[1]);
    typename ViewT::pixel_type pix2 = interp_image2(sample_pix2[0], sample_pix2[1]);
    if ( is_valid(pix1) && is_valid(pix2) ) {

       //result.push_back(Vector4(pix1[0],pix2[0],lonlat[0],lonlat[1]));

       int x = (int)sample_pix_dem[0];
       int y = (int)sample_pix_dem[1];

       if (x < 0){
           x = 0;
       }
       if (x > dem_image.cols()-1){
           x = dem_image.cols()-1;
       }
       if (y < 0){
           y = 0;
       }
       if (y > dem_image.rows()-1){
           y = dem_image.rows()-1;
       }

       Vector3 longlat3(lonlat(0),lonlat(1),(dem_image)(x, y));
       Vector3 xyz = geo1.datum().geodetic_to_cartesian(longlat3);

       Vector2 sample_pix_dem_left;
       sample_pix_dem_left(0) = x-1;
       if (sample_pix_dem_left(0) < 0){
          sample_pix_dem_left(0) = 0;
          //break;
       }
       sample_pix_dem_left(1) = y;
       lonlat = GR.pixel_to_lonlat(sample_pix_dem_left);

       Vector3 longlat3_left(lonlat(0),lonlat(1),(dem_image)(sample_pix_dem_left(0), sample_pix_dem_left(1)));
       Vector3 xyz_left = geo1.datum().geodetic_to_cartesian(longlat3_left);

       Vector2 sample_pix_dem_top;
       sample_pix_dem_top(0) = x;
       sample_pix_dem_top(1) = y-1;
       if (sample_pix_dem_top(1) < 0){
         sample_pix_dem_top(1) = 0;
         //break;
       }

       lonlat = GR.pixel_to_lonlat(sample_pix_dem_top);
       Vector3 longlat3_top(lonlat(0),lonlat(1),(dem_image)(sample_pix_dem_top(0), sample_pix_dem_top(1)));
       Vector3 xyz_top = geo1.datum().geodetic_to_cartesian(longlat3_top);



       Vector3 normal = computeNormalFrom3DPoints(xyz, xyz_left, xyz_top);

       //printf("normal:%f %f %f\n", normal(0), normal(1), normal(2));
       //printf("%f %f %f\n", loc_longlat3(0), loc_longlat3(1), loc_longlat3(2));

       result.push_back(Vector4(pix1[0],pix2[0],lonlat[0],lonlat[1]));

       normalArray->push_back(normal);
       xyzArray->push_back(xyz);

       ++sample;
       //      std::cout << result[result.size()-1][0] << " " << result[result.size()-1][1] << "\n";
    }
    ++numtries;
  }
  return result;
}
/*
/// Erases a file suffix if one exists and returns the base string
static std::string prefix_from_filename(std::string const& filename) {
        std::string result = filename;
        int index = result.rfind(".");
        if (index != -1)
                result.erase(index, result.size());
        return result;
}

/// Erases a file suffix if one exists and returns the base string less3 characters
static std::string prefix_less3_from_filename(std::string const& filename) {
  std::string result = filename;
  int index = result.rfind(".");
  if (index != -1)
    result.erase(index-3, result.size()+3);
  return result;
}

/// Erases a file suffix if one exists and returns the base string less3 characters
static std::string sufix_from_filename(std::string const& filename) {
  std::string result = filename;
  int index = result.rfind("/");
  if (index != -1)
    result.erase(0, index);
  return result;
}
*/

//reads the tiff DEM into a 3D coordinate
//pos is a Vector2 of pixel coordinates, GR is georeference
template <class ViewT>
Vector3 pixel_to_cart (Vector2 pos, ImageViewBase<ViewT> const& img,  GeoReference GR) {
    Vector2 loc_longlat2=GR.point_to_lonlat(GR.pixel_to_point(pos));
    Vector3 loc_longlat3(loc_longlat2(0),loc_longlat2(1),img((int)pos[0],(int)pos[1]));
    Vector3 loc_cartesian=GR.datum().geodetic_to_cartesian(loc_longlat3);
    return loc_cartesian;
}


// Create the output, index, and radiance file names

std::vector<std::string> parse_command_arguments(int argc, char *argv[] ) {
        int num_matches;
        std::vector<std::string> input_files;

        po::options_description general_options("Options");
        general_options.add_options()
        ("help", "Display this help message")
        ("num-matches,m", po::value<int>(&num_matches)->default_value(1000), "Number of points to match for linear regression.");

        po::options_description hidden_options("");
        hidden_options.add_options()
        ("input-files", po::value<std::vector<std::string> >(&input_files));

        po::options_description options("Allowed Options");
        options.add(general_options).add(hidden_options);

        po::positional_options_description p;
        p.add("input-files", -1);

        po::variables_map vm;
        po::store( po::command_line_parser( argc, argv ).options(options).positional(p).run(), vm );
        po::notify( vm );

        std::ostringstream usage;
        usage << "Description: tonematches several images" << std::endl << std::endl;
        usage << "Usage: histeq [options] <filename1> <filename2> ..." << std::endl << std::endl;
        usage << general_options << std::endl;

        if( vm.count("help") ) {
                std::cerr << usage.str() << std::endl;
                exit(1);
        }

        if( vm.count("input-files")<1 ) {
                std::cerr << "Error: Must specify at least one input file!" << std::endl << std::endl;
                std::cerr << usage.str();
                exit(1);
        }

        return input_files;
}



