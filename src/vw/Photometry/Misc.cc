// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
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

#include <boost/tokenizer.hpp>
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
using namespace std;

#include <vw/Photometry/Reconstruct.h>
#include <vw/Photometry/Misc.h>
#include <vw/Photometry/Shape.h>
using namespace vw::photometry;

// Upsample a geo-referenced tiff image. Care is taken to deal properly with invalid pixels.
void vw::photometry::upsample_uint8_image(std::string output_file, std::string input_file, int upsampleFactor){

  // Note: This is function is not quite correct. I noticed that the
  // original and upsampled image are not exactly on top of each
  // other. The georeference may need a bit of adjustment. I think
  // there is a good example in orthoproject.cc in Stereo Pipeline
  // about how to transform an image correctly.
  
  GeoReference geo;
  read_georeference(geo, input_file);
  DiskImageView<PixelMask<PixelGray<uint8> > > img(input_file);

  int cols = (img.cols())*upsampleFactor, rows = (img.rows())*upsampleFactor;
  ImageView<PixelMask<PixelGray<uint8> > > up_img(cols, rows);
  
  InterpolationView<EdgeExtensionView<DiskImageView<PixelMask<PixelGray<uint8> > >, ConstantEdgeExtension>, BilinearInterpolation>
    interp_img = interpolate(img, BilinearInterpolation(), ConstantEdgeExtension());
  
  for (int y=0; y<rows; ++y){
    for (int x=0; x<cols; ++x){

      double xx = (double)x/upsampleFactor;
      double yy = (double)y/upsampleFactor;

      if ( is_valid(img( floor(xx), floor(yy))) &&
           is_valid(img( floor(xx), ceil (yy))) &&
           is_valid(img( ceil (xx), floor(yy))) &&
           is_valid(img( ceil (xx), ceil (yy)))
           ){
        up_img(x, y) = interp_img(xx, yy);
        up_img(x, y).validate();
      }else{
        up_img(x, y) = 0;
        up_img(x, y).invalidate();
      }
    }
  }
  
  Matrix<double> H = geo.transform();
  H(0,0) /= upsampleFactor;
  H(1,1) /= upsampleFactor;
  geo.set_transform(H);

  std::cout << "Writing: " << output_file << std::endl;
  write_georeferenced_image(output_file, up_img, geo, TerminalProgressCallback("photometry","Processing:"));
  return;
}

//upsamples a geo referenced tiff image by four- used in sfs
void upsample_image(std::string output_file, std::string input_file, int upsampleFactor) {
  GeoReference geo;
  read_georeference(geo, input_file);
  DiskImageView<PixelGray<float> >   image(input_file);

  int cols = (image.cols())*upsampleFactor, rows = (image.rows())*upsampleFactor;
  ImageView<PixelGray<float> >  tm_image(cols, rows);

  // This is the wrong way of doing interpolation. See the Stereo module for the right way.
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


  // Wrong way of doing interpolation
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
                                   cartography::GeoReference const& geo1,
                                   cartography::GeoReference const& geo2,
                                   int num_samples,
                                   std::string const& DEM_file,
                                   std::vector<Vector3> *normalArray,
                                   std::vector<Vector3> *xyzArray ) {
  int sample = 0;
  int numtries = 0;
  std::vector<Vector4> result;

  // Random numbers
  srandom((unsigned int) clock());

  // Wrong way of doing interpolation
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

//reads the tiff DEM into a 3D coordinate
//pos is a Vector2 of pixel coordinates, GR is georeference
template <class ViewT>
Vector3 pixel_to_cart (Vector2 pos, ImageViewBase<ViewT> const& img,  cartography::GeoReference GR) {
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

void vw::photometry::getTileCornersWithoutPadding(// Inputs
                                                  int numCols, int numRows,
                                                  cartography::GeoReference const& geoRef,
                                                  double tileSize, int pixelPadding,
                                                  // Outputs
                                                  double & min_x, double & max_x,
                                                  double & min_y, double & max_y
                                                  ){
  
  // Given a tile which we know is padded by pixelPadding on each side, find the lon-lat coordinates
  // of the tile corners without the padding.
  
  // The functions applyPaddingToTileCorners() and getTileCornersWithoutPadding() are intimately
  // related.
  
  // The north-west and south-east tile corners
  Vector2 tile_NW = geoRef.pixel_to_lonlat(Vector2(pixelPadding, pixelPadding));
  Vector2 tile_SE = geoRef.pixel_to_lonlat(Vector2(numCols - 1 - pixelPadding, numRows - 1 - pixelPadding));

  // Snap to the corners of the tile proper, ignoring the half-pixel discrepancy
  // which I still need to understand. 
  min_x = tileSize*round(tile_NW(0)/tileSize);
  max_x = tileSize*round(tile_SE(0)/tileSize);
  min_y = tileSize*round(tile_SE(1)/tileSize);
  max_y = tileSize*round(tile_NW(1)/tileSize);
  
  return;
}
  

void vw::photometry::applyPaddingToTileCorners(// Inputs
                                               cartography::GeoReference const& geoRef,
                                               int pixelPadding,
                                               double min_x, double max_x,
                                               double min_y, double max_y,
                                               // Outputs
                                               double & min_x_padded, double & max_x_padded,
                                               double & min_y_padded, double & max_y_padded){
  
  // Given a tile, put a padding of pixelPadding pixels on each side. Return the lon-lat
  // coordinates of the obtained tile.
  
  // The functions applyPaddingToTileCorners() and getTileCornersWithoutPadding() are intimately
  // related.
  
  // Upper left corner lon lat
  Vector2 A = Vector2(min_x, max_y);
  
  // Right and down by pixelPadding
  Vector2 B = geoRef.pixel_to_lonlat(geoRef.lonlat_to_pixel(A) + Vector2(pixelPadding, pixelPadding));
  Vector2 D = B - A;
  
  // Careful with the signs below.
  min_x_padded = min_x - D(0); max_x_padded = max_x + D(0);
  min_y_padded = min_y + D(1); max_y_padded = max_y - D(1);

  return;
}

void vw::photometry::readDEMTilesIntersectingBox(// Inputs
                                                 double noDEMDataValue,
                                                 Vector2 boxNW, Vector2 boxSE,
                                                 std::vector<std::string> const& DEMTiles,
                                                 // Outputs
                                                 ImageView<PixelGray<float> > & combinedDEM,
                                                 cartography::GeoReference    & combinedDEM_geo){
  
  // Given a set of int16 DEM tiles and a box, get all the pixels from
  // all the tiles which are contained within the box. The newly
  // created image will be float.

  // The box is specified by the North-West and South-East corners.

  // First thing initialize the outputs
  combinedDEM.set_size(0, 0);
  combinedDEM_geo = GeoReference();

  ImageView<int> count;
  
  bool isFirstImage = true;
  for (int i = 0; i < (int)DEMTiles.size(); i++){
    
    std::string DEMTileFile = DEMTiles[i];
    GeoReference DEMTile_geo;

    ImageFormat img_fmt;
    {
      boost::scoped_ptr<SrcImageResource> img_rsrc( DiskImageResource::open(DEMTileFile) );
      img_fmt = img_rsrc->format();
    }
    
    // Get just the portion of the tile which overlaps with the current box
    // Note that we accept float or int16 input DEM tfiles.
    ImageView<PixelGray<float> > DEMTile;
    //std::cout << "Will read: " << DEMTileFile  << std::endl;
    bool success;
    
    if (img_fmt.channel_type == VW_CHANNEL_INT16){
      success = getSubImageWithMargin< PixelGray<int16>, PixelGray<float> > 
        (boxNW, boxSE, DEMTileFile, // Inputs
         DEMTile, DEMTile_geo       // Outputs
         );
    }else{
      success = getSubImageWithMargin< PixelGray<float>, PixelGray<float> > 
        (boxNW, boxSE, DEMTileFile, // Inputs
         DEMTile, DEMTile_geo       // Outputs
         );
    }
    if (!success) continue;
    
    if (isFirstImage){
      isFirstImage = false;
      // The first iteration. The right time to initialize combinedDEM and create its GeoReference.
      Vector2 begPixel = DEMTile_geo.lonlat_to_pixel(boxNW);
      Vector2 endPixel = DEMTile_geo.lonlat_to_pixel(boxSE);

      // Make the image a bit larger than necessary to help with
      // bilinear interpolation below. A padding of 1 pixel would
      // probably be enough here.
      int extra = 2;
      begPixel(0) = floor(begPixel(0)) - extra; begPixel(1) = floor(begPixel(1)) - extra;
      endPixel(0) = ceil(endPixel(0))  + extra; endPixel(1) = ceil(endPixel(1))  + extra;
      int numCols = (int)round(endPixel(0) - begPixel(0));
      int numRows = (int)round(endPixel(1) - begPixel(1));
      combinedDEM.set_size(numCols, numRows);
      count.set_size(numCols, numRows);
      for (int row = 0; row < combinedDEM.rows(); row++){
        for (int col = 0; col < combinedDEM.cols(); col++){
          combinedDEM(col, row) = noDEMDataValue;
          count(col, row) = 0;
        }
      }
      
      // In combinedDEM_geo, the (0, 0) pixel will where
      // begPixel is in DEMTile_geo.
      combinedDEM_geo = vw::cartography::crop(DEMTile_geo, begPixel(0), begPixel(1));
    }
    
    for (int col = 0; col < DEMTile.cols(); col++){
      for (int row = 0; row < DEMTile.rows(); row++){
        Vector2 pix = combinedDEM_geo.lonlat_to_pixel(DEMTile_geo.pixel_to_lonlat(Vector2(col, row)));
        int lCol = (int)round(pix(0));
        int lRow = (int)round(pix(1));
        if (0 <= lCol && lCol < combinedDEM.cols() &&
            0 <= lRow && lRow < combinedDEM.rows() &&
            DEMTile(col, row) != noDEMDataValue
            ){
          if (count(lCol, lRow) == 0) combinedDEM(lCol, lRow)  = DEMTile(col, row);
          else                        combinedDEM(lCol, lRow) += DEMTile(col, row);
          count(lCol, lRow)++;
        }
      }
    }
  } // Done visiting the overlapping tiles

  // Average the obtained pixels
  for (int row = 0; row < combinedDEM.rows(); row++){
    for (int col = 0; col < combinedDEM.cols(); col++){
      if (count(col, row) > 1) combinedDEM(col, row) /= count(col, row);
    }
  }
  
  return;
}

void vw::photometry::listTifsInDir(const std::string & dirName,
                                   std::vector<std::string> & tifsInDir
                                   ){
  
  tifsInDir.clear();
  
  fs::path dir(dirName);
  if ( !fs::exists(dirName) || !fs::is_directory(dirName ) ) return;

  fs::directory_iterator end_iter; // default construction yields past-the-end
  for  ( fs::directory_iterator dir_iter(dirName); dir_iter != end_iter; ++dir_iter)
    {
      if (! fs::is_regular_file(dir_iter->status()) ) continue;
      std::string fileName = (*dir_iter).string();
      int len = fileName.size();
      if (len >= 4 && fileName.substr(len - 4, 4) == ".tif"){
        //std::cout << "Now adding " << fileName << std::endl;
        tifsInDir.push_back( fileName );
      }
    }

  // Sort the files in lexicographic order
  std::sort(tifsInDir.begin(), tifsInDir.end());
  
  return;
}

void vw::photometry::writeSunAndSpacecraftPosition(std::string prefix,
                                                   std::string sunFile, std::string spacecraftFile,
                                                   Vector3 sunPosition, Vector3 spacecraftPosition){

  // convert from meters to kilometers
  double tokm = 0.001;
  sunPosition        *= tokm;
  spacecraftPosition *= tokm;
  
  std::ofstream sunF(sunFile.c_str());
  sunF.precision(16);
  sunF << prefix << " "
       << sunPosition(0) << " "
       << sunPosition(1) << " "
       << sunPosition(2)<< std::endl;
  sunF.close();
  
  std::ofstream spacecraftF(spacecraftFile.c_str());
  spacecraftF.precision(16);
  spacecraftF << prefix << " "
              << spacecraftPosition(0) << " "
              << spacecraftPosition(1) << " "
              << spacecraftPosition(2)<< std::endl;
  spacecraftF.close();
  
  return;
}

std::string vw::photometry::getFirstElevenCharsFromFileName(std::string fileName){

  // Out of path/to/AS15-M-1723_1724-DEM.tif extract the 11 characters
  // "AS15-M-1723".
  int index = fileName.rfind("/");
  if (index != -1) fileName.erase(0, index + 1);
  
  return fileName.substr(0, 11);
}

void vw::photometry::indexFilesByKey(std::string dirName, std::map<std::string, std::string> & index){

  // Create the map: AS15-M-1723 -> path/to/AS15-M-1723_1724-DEM.tif
  // for all the files in the given directory.

  index.clear();
  
  fs::path dir(dirName);
  if ( !fs::exists(dirName) || !fs::is_directory(dirName ) ){
    std::cerr << "Cannot find directory: " << dirName << std::endl;
    exit(1);
  }
  
  fs::directory_iterator end_iter; // default construction yields past-the-end
  for (fs::directory_iterator dir_iter(dirName); dir_iter != end_iter; ++dir_iter){
    if (!fs::is_regular_file(dir_iter->status())) continue;
    std::string fileName = (*dir_iter).string();
    index[ getFirstElevenCharsFromFileName(fileName) ] = fileName;
  }

  return;
}

void vw::photometry::enforceUint8Img(std::string imgName){

  ImageFormat img_fmt;
  {
    boost::scoped_ptr<SrcImageResource> img_rsrc( DiskImageResource::open(imgName) );
    img_fmt = img_rsrc->format();
  }
  if (img_fmt.channel_type != VW_CHANNEL_UINT8){
    std::cout << "ERROR: The input DRG images must be uint8." << std::endl;
    exit(1); // We check this exit status later
  }

  return;
}

bool vw::photometry::readNoDEMDataVal(std::string DEMFile, float & noDEMDataValue){

  boost::scoped_ptr<SrcImageResource> rsrc( DiskImageResource::open(DEMFile) );
  if ( rsrc->has_nodata_read() ){
    noDEMDataValue = rsrc->nodata_read();
    return true;
  }
  
  return false;
}


void vw::photometry::maskPixels(std::string imgFile, std::string maskFile, double shadowThresh, std::string outDir){
  
  // Any pixels in imgFile, which are below shadowThresh in maskFile, will be set to black.
  
  fs::create_directories(outDir);

  std::cout << "Reading " << imgFile << std::endl;
  DiskImageView<PixelMask<PixelGray<uint8> > > inputImg(imgFile);
  ImageViewRef<PixelMask<PixelGray<uint8> > > inputImgRef = inputImg;
  GeoReference imgGeo;
  read_georeference(imgGeo, imgFile);

  std::cout << "Reading " << maskFile << std::endl;
  DiskImageView<PixelMask<PixelGray<uint8> > >  maskImg(maskFile);
  ImageViewRef<PixelMask<PixelGray<uint8> > > maskImgRef = maskImg;
  GeoReference maskGeo;
  read_georeference(maskGeo, maskFile);

  std::string outFile = outDir + suffix_from_filename(imgFile);

  std::cout << "Writing: " << outFile << std::endl;
  write_georeferenced_image(outFile,
                            mask_image(inputImgRef, maskImgRef, shadowThresh, imgGeo, maskGeo),
                            imgGeo, TerminalProgressCallback("asp",""));
  
  return;
}

void vw::photometry::ReadPhaseCoeffsFromFile(std::string phaseDir, GlobalParams& settings)
{

  // Read the latest values of the phase coefficients from the file

  settings.phaseCoeffsFileName = phaseDir + "/phaseCoeffs.txt";
  
  std::ifstream fp(settings.phaseCoeffsFileName.c_str());
  if (fp){
    float a1, a2;
    while( fp >> a1 >> a2){
      settings.phaseCoeffA1 = a1;
      settings.phaseCoeffA2 = a2;
    }
  }
  fp.close();

  return;
}

void vw::photometry::AppendPhaseCoeffsToFile(const GlobalParams& settings){

  // Append the current phase coefficients to the file. This way when
  // we do multiple albedo iterations we keep all the current and
  // previous values of phase coefficients.
  
  if (settings.phaseCoeffsFileName == ""){
    std::cerr << "The phase coefficients file was not set yet." << std::endl;
    exit(1);
  }
  
  FILE *fp;
  fp = fopen(settings.phaseCoeffsFileName.c_str(), "a");
  std::cout << "Writing " << settings.phaseCoeffsFileName << std::endl;
  fprintf(fp, "%f %f\n", settings.phaseCoeffA1, settings.phaseCoeffA2);
  fclose(fp);
}

float vw::photometry::getShadowThresh(const GlobalParams& settings, float exposureRefl){

  float t = -1.0; // shadow threshold
  if (settings.shadowRemovalType == CONSTANT_THRESHOLD_SHADOW_REMOVAL)
    t = settings.shadowThresh;
  else if (settings.shadowRemovalType == LAMBERTIAN_THRESHOLD_SHADOW_REMOVAL ||
           settings.shadowRemovalType == LUNAR_LAMBERTIAN_THRESHOLD_SHADOW_REMOVAL)
    t = settings.shadowThresh/exposureRefl;

  return t;
}

void vw::photometry::resampleImage(std::string initFilename, std::string outputFilename, int factor){

  DiskImageView<float>  initImg(initFilename);
  GeoReference initGeo;
  read_georeference(initGeo, initFilename);

  InterpolationView<EdgeExtensionView<EdgeExtensionView<DiskImageView<float>, ConstantEdgeExtension>,
    ConstantEdgeExtension>, BilinearInterpolation> interpInitImg
    = interpolate(edge_extend(initImg.impl(),ConstantEdgeExtension()), BilinearInterpolation());

  ImageView<float> outImg(initImg.rows()/factor, initImg.cols()/factor);
    
  //create the outputGeo - START
  GeoReference outputGeo = initGeo;

  Matrix<double> init_H;
  init_H = initGeo.transform();
  cout<<"init_H="<<init_H<<endl;
    
  Matrix<double> output_H;
  output_H = initGeo.transform();
  //lon = H(0,0)*i + 0*j + H(0,2)
  //lat = 0*i + H(1,1)*j + H(1,2)
   
  output_H(0,2) = init_H(0,2);
  output_H(1,2) = init_H(1,2);
  output_H(0,0) = factor*init_H(0,0);
  output_H(1,1) = factor*init_H(1,1);

  outputGeo.set_transform(output_H);
  //create the outputGeo - END

  for (int j = 0; j < outImg.rows(); j++){
    for (int i = 0; i < outImg.cols(); i++){

      Vector2 outputPix;
      outputPix(0) = i;
      outputPix(1) = j;
      Vector2 outLonLat = outputGeo.pixel_to_lonlat(outputPix);
	 
      Vector2 initPix;
      initPix = initGeo.lonlat_to_pixel(outLonLat);
       
      outImg.impl()(i,j) = interpInitImg.impl()(initPix(0), initPix(1));
    }
  }
       
 

  //write the corrected file
  write_georeferenced_image(outputFilename,
                            outImg,
                            outputGeo, TerminalProgressCallback("photometry","Processing:"));
    

}  

bool vw::photometry::boxesOverlap(const Vector4 & box1Corners, const Vector4 & box2Corners){

  int lonOverlap = 0;
  int latOverlap = 0; 
  
  if (box1Corners(0) > box1Corners(1) || box2Corners(0) > box2Corners(1))
    {
      std::cout << "ERROR: Must never happen: " << __FILE__ << " at line " << __LINE__ << std::endl;
      exit(1);
    }
  
  if ( std::max(box1Corners(0), box2Corners(0)) < std::min(box1Corners(1), box2Corners(1)) )
    {
      lonOverlap = 1;
    }
  
  if (box1Corners(2) > box1Corners(3) || box2Corners(2) > box2Corners(3))
    {
      std::cout << "ERROR: Must never happen: " << __FILE__ << " at line " << __LINE__ << std::endl;
      exit(1);
    }
  
  if ( std::max(box1Corners(2), box2Corners(2)) < std::min(box1Corners(3), box2Corners(3)) )
    {
      latOverlap = 1;
    }
  
  return (lonOverlap == 1 && latOverlap == 1);
         
}
      

Vector4 vw::photometry::ComputeGeoBoundary(cartography::GeoReference Geo, int width, int height){

  // Get the lonlat coordinates of the four pixels corners of the image.
  
  Vector4 corners;
  Vector2 leftTopPixel(0,0);
  Vector2 leftTopLonLat = Geo.pixel_to_lonlat(leftTopPixel);

  Vector2 rightBottomPixel(width-1, height-1);
  Vector2 rightBottomLonLat = Geo.pixel_to_lonlat(rightBottomPixel);
  
  float minLon = leftTopLonLat(0);
  float minLat = leftTopLonLat(1);
  float maxLon = rightBottomLonLat(0);
  float maxLat = rightBottomLonLat(1);

  if (maxLat<minLat){
    float temp = minLat;
    minLat = maxLat;
    maxLat = temp;    
  }

  if (maxLon<minLon){
    float temp = minLon;
    minLon = maxLon;
    maxLon = temp;    
  }

  corners(0) = minLon;
  corners(1) = maxLon;
  corners(2) = minLat;
  corners(3) = maxLat;

  return corners;
}

Vector4 vw::photometry::getImageCorners(std::string imageFile){

  // Get the four corners of an image, that is the lon-lat coordinates of
  // the pixels in the image corners.

  // Note: Below we assume that the image is uint8. In fact, for the
  // purpose of calculation of corners the type of the image being
  // read does not matter.
  DiskImageView<PixelMask<PixelGray<uint8> > >  image(imageFile);
  
  GeoReference imageGeo;
  read_georeference(imageGeo, imageFile);
  Vector4 imageCorners = ComputeGeoBoundary(imageGeo, image.cols(), image.rows());
  return imageCorners;
}

void vw::photometry::listTifsInDirOverlappingWithBox(const std::string & dirName,
                                                     Vector4 & boxCorners,
                                                     const std::string & outputListName){
  
  std::vector<std::string> tifsInDir;
  listTifsInDir(dirName, tifsInDir);

  ofstream fh(outputListName.c_str());
  if (!fh) {
    std::cerr << "ERROR: listTifsInDirOverlappingWithBox: can't open " << outputListName
              << " for writing" << std::endl;
    exit(1);
  }

  fh.precision(20);
  for (int fileIter = 0; fileIter < (int)tifsInDir.size(); fileIter++){
    const std::string & currFile = tifsInDir[fileIter];
    Vector4 currCorners = getImageCorners(currFile);
    if (!boxesOverlap(currCorners, boxCorners)) continue;
    fh << 1 << " " << currFile << " " << currCorners(0) << " " << currCorners(1) << " " <<
      currCorners(2) << " " << currCorners(3) << std::endl;
  }
  fh.close();

  return;
}

void vw::photometry::createAlbedoTilesOverlappingWithDRG(double tileSize, int pixelPadding,
                                                         std::string imageFile, Vector4 const& simulationBox,
                                                         std::vector<ImageRecord> const& drgRecords,
                                                         std::string blankTilesList,  std::string blankTilesDir,
                                                         std::string DEMTilesList,    std::string meanDEMDir,
                                                         std::string albedoTilesList, std::string albedoDir
                                                         ){

  // Create all the tiles which overlap with all DRG images which in
  // turn overlap with the simulation box.

  // The georeference of tiles will be obtained from the georeference
  // of an input image (any one of those images would work as well as
  // any other).

  //  Write to disk both the tiles themselves and their list.

  // Note that the tiles have a padding, so they are a few pixels larger than what
  // they should be. We need that in order to be able to compute the normals
  // for  DEM, and also for SfS. The padding will be removed when at the end
  // of all computations we save the final albedo.

  if (tileSize <= 0.0 || pixelPadding < 0){
    std::cout << "ERROR: Must have positive tile size and non-negative pixel padding." << std::endl;
    exit(1);
  }

  // To do: it is more intuitive if one iterates from north to south than from south to north.

  // Find all tiles overlapping with given DRGs. Use a set to avoid duplicates.
  std::set< std::pair<double, double> > Tiles;
  for (int j = 0; j < (int)drgRecords.size(); j++){
    const ImageRecord& rec = drgRecords[j];
    Vector4 currCorners = Vector4(rec.west, rec.east, rec.south, rec.north);
    for (double min_y = tileSize*floor(rec.south/tileSize); min_y < rec.north; min_y += tileSize){
      for (double min_x = tileSize*floor(rec.west/tileSize); min_x < rec.east; min_x += tileSize){
        Tiles.insert(std::make_pair(min_y, min_x));
      }
    }
  }

  // The input DRG must be uint8
  enforceUint8Img(imageFile);
  
  GeoReference geo; read_georeference(geo, imageFile);
  ofstream fht(blankTilesList.c_str());  fht.precision(20);
  ofstream fhd(DEMTilesList.c_str());    fhd.precision(20);
  ofstream fha(albedoTilesList.c_str()); fha.precision(20);

  for (std::set< std::pair<double, double> >::iterator it = Tiles.begin(); it != Tiles.end(); it++){
    
    std::pair<double, double> Tile = *it;

    // Tile corners coordinates without padding
    double min_y = Tile.first;
    double min_x = Tile.second;
    double max_y = min_y + tileSize;
    double max_x = min_x + tileSize;
    
    // Tile corners coordinates with padding
    double min_x_padded, max_x_padded, min_y_padded, max_y_padded;
    applyPaddingToTileCorners(// Inputs
                              geo, pixelPadding, min_x, max_x,  min_y, max_y,  
                              // Outputs
                              min_x_padded, max_x_padded, min_y_padded, max_y_padded
                              );
  
    // Set the upper-left corner in the tile
    Matrix3x3 T = geo.transform();
    T(0,2) = min_x_padded;
    T(1,2) = max_y_padded;
    geo.set_transform(T);

    // Determine the size of the tile
    Vector2 pixUL = geo.lonlat_to_pixel(Vector2(min_x_padded, max_y_padded));
    // Note: The value we get for  pixUL is (-0.5, -0.5).
    // I was expecting (0, 0). (Oleg)
    Vector2 pixLR = geo.lonlat_to_pixel(Vector2(max_x_padded, min_y_padded));

    // To do: Below nrows and ncols may need to be interchanged.
    int nrows = (int)round(pixLR(0) - pixUL(0));
    int ncols = (int)round(pixLR(1) - pixUL(1));

    double uE = min_x, uN = max_y; // uppper-left corner without padding
    std::string sN = "N", sE = "E";
    if (uE < 0){ uE = -uE; sE = "W";}
    if (uN < 0){ uN = -uN; sN = "S";}
    ostringstream os;
    os << blankTilesDir << "/tile_" << uE << sE << uN << sN << ".tif";
    std::string blankTileFile = os.str();
    std::string DEMTileFile    = meanDEMDir + suffix_from_filename(blankTileFile);
    std::string albedoTileFile = albedoDir  + suffix_from_filename(blankTileFile);
    
    // The blank tiles themselves have no information, they are just
    // templates which we will later cycle through and create DRG
    // and tiles at each pixel.
    ImageView<PixelMask<PixelGray<float> > > blankTile(nrows, ncols);
      
    std::cout << "Writing " << blankTileFile << std::endl;
    write_georeferenced_image(blankTileFile,
                              channel_cast<uint8>(clamp(blankTile, 0.0,255.0)),
                              geo, TerminalProgressCallback("{Core}","Processing:"));

    // The actual corners of the tile may differ slightly than what we intended
    // due to rounding. Compute the actual corners now that the tile was created.
    Vector4 C = ComputeGeoBoundary(geo, blankTile.cols(), blankTile.rows());

    fht << 1 << " " << blankTileFile << " "
        << C(0) << " " << C(1) << " " << C(2) << " " << C(3) << std::endl;

    fha << 1 << " " << albedoTileFile << " "
        << C(0) << " " << C(1) << " " << C(2) << " " << C(3) << std::endl;

    fhd << 1 << " " << DEMTileFile << " "
        << C(0) << " " << C(1) << " " << C(2) << " " << C(3) << std::endl;

  } // End iterating over tiles
    
  fht.close();
  fhd.close();
  fha.close();
  
}

std::vector<int> vw::photometry::GetInputIndices( std::vector<std::string> inputFiles, std::vector<std::string> DRGFiles){

  std::vector<int>  inputIndices;
  for (int j = 0; j < (int)inputFiles.size(); j++){
    for (int i = 0; i < (int)DRGFiles.size(); i++){
      if (DRGFiles[i].compare(inputFiles[j])==0){
        inputIndices.push_back(i);
      }
    }
  }

  return inputIndices;
}

//this function determines the image overlap for the general case
//it takes into consideration any set of overlapping images.
std::vector<int> vw::photometry::makeOverlapList(const std::vector<ModelParams>& drgFiles,
                                                 const std::string& currFile) {

  std::vector<int> overlapIndices; overlapIndices.clear();
  Vector4 currCorners = getImageCorners(currFile);

  //std::cout << "file " << currFile << " overlaps with ";
  for (unsigned int i = 0; i < drgFiles.size(); i++){

    const ModelParams& params = drgFiles[i];
    //std::cout << params.inputFilename << " ";

    Vector4 corners;
    if (params.corners(3) == ImageRecord::defaultCoord) {
      std::cerr << "ERROR: Missing the bounding box information for image: " << params.inputFilename
                << std::endl;
      cerr << "This should have been specified in the list of images." << endl;
      exit(1);
      corners = getImageCorners(params.inputFilename);
      //std::cout << "Reading from disk: " << corners << std::endl;
    } else {
      corners = params.corners;
      //std::cout << "Cached from list: " << corners << std::endl;
    }

    if (boxesOverlap(corners, currCorners) && currFile != params.inputFilename){
      overlapIndices.push_back(i);
      //std::cout << params.inputFilename << " ";
    }
  }
  
  //std::cout << std::endl;
  return overlapIndices;
}

std::vector<int> vw::photometry::makeOverlapList(const std::vector<ImageRecord>& drgRecords,
                                                 const std::string& currFile) {

  // To do: Merge this function with the one above it and together with other
  // overlap logic seen in this file.
  std::vector<int> overlapIndices;
  Vector4 corners = getImageCorners(currFile);

  for (int j = 0; j < (int)drgRecords.size(); j++){
    const ImageRecord& rec = drgRecords[j];
    Vector4 currCorners = Vector4(rec.west, rec.east, rec.south, rec.north);
    if (! boxesOverlap(currCorners, corners)){
      continue;
    }
    overlapIndices.push_back(j);
  }

  return overlapIndices;
}

std::vector<std::vector<int> > vw::photometry::makeOverlapList(const std::vector<std::string>& inputFiles,
                                                               const std::vector<ModelParams>& DRGFiles){
  std::vector<std::vector<int> > overlapIndices;
  overlapIndices.resize(inputFiles.size());
  
  for (unsigned int i = 0; i < inputFiles.size(); ++i) {  
    overlapIndices[i] = makeOverlapList(DRGFiles, inputFiles[i]);
  }
  return overlapIndices;  
}

void vw::photometry::printOverlapList(std::vector<std::vector<int> > overlapIndices){
  
  for (int i=0; i < (int)overlapIndices.size(); i++){
    printf("%d: ", i);
    for (int j = 0; j < (int)overlapIndices[i].size(); j++){
      printf("%d ", overlapIndices[i][j]);
    }
    printf("\n");
  }
}


Vector4 vw::photometry::parseSimBox(std::string simulationBoxStr){

  // Parse the string "13:49:-12:28" to extract the vector of
  // numbers 13, 49, -12, 28 (lonMin, lonMax, latMin, latMax).
  
  typedef boost::tokenizer<boost::char_separator<char> >  tokenizer;
  boost::char_separator<char> colon(":");
  tokenizer tokens(simulationBoxStr, colon);

  Vector4 simulationBox;
  int count = 0;
  for (tokenizer::iterator tok_iter = tokens.begin(); tok_iter != tokens.end(); ++tok_iter){
    std::string tok = *tok_iter;
    if (tok == "") continue;
    simulationBox(count) = atof(tok.c_str());
    count++;
    if (count >= 4) break;
  }

  // If parsing did not succeed, then fail
  if (count < 4){
    cerr << "ERROR: Could not extract the simulation box from the string: " << simulationBoxStr << endl;
    exit(1);
  }

  if (simulationBox(0) >= simulationBox(1) || simulationBox(2) >= simulationBox(3)){
    std::cerr << "ERROR: Invalid simulationBox: " << simulationBox << std::endl;
    simulationBox = Vector4(0, 0, 0, 0);
  }

  // If we simulate the full sphere, we need to go beyond [-180, 180], since
  // the images can have pixels outside of this range.
  if (simulationBox(0) <= -180.0) simulationBox(0) = std::min(-360.0, simulationBox(0));
  if (simulationBox(1) >=  180.0) simulationBox(1) = std::max( 360.0, simulationBox(1));
  if (simulationBox(2) <= -180.0) simulationBox(2) = std::min(-360.0, simulationBox(2));
  if (simulationBox(3) >=  180.0) simulationBox(3) = std::max( 360.0, simulationBox(3));
  
  return simulationBox;
}

void vw::photometry::extractSimBox(char * line, Vector4 & simulationBox){

  // Out of the string "SIMULATION_BOX            6 : 10 : -10 : -9 "
  // extract the value  "6 : 10 : -10 : -9", then parse it to extract
  // the individual numbers in a vector.
  
  istringstream is(line);
  std::string token, boxStr;

  // First token
  if ( !(is >> token) || token != "SIMULATION_BOX"){
    return;
  }

  // Subsequent tokens
  boxStr = "";
  while(is >> token){
    boxStr += token + " ";
  }

  simulationBox = parseSimBox(boxStr);
  
  return;
}

int vw::photometry::ReadConfigFile(char *config_filename, struct GlobalParams & settings){
  int MAX_LENGTH = 5000;
  char line[MAX_LENGTH];
  char inName[MAX_LENGTH];
  char inVal[MAX_LENGTH];
  char buffer[MAX_LENGTH];
  char *commentPos;
  ifstream configFile(config_filename);
  int ret;

  // Default values
  settings.drgDir                = "";
  settings.demDir                = "";
  settings.sunPosFile            = "";
  settings.spacecraftPosFile     = "";
  settings.initialSetup          = 0;
  settings.tileSize              = -1.0; // invalid on purpose, must be set later
  settings.useTiles              = 1;    // 1 or 0, this variable will go away 
  settings.pixelPadding          = 5;    // By how much to pad the albedo tiles
  // Default simulation box, simulate the full sphere. We need to go
  // beyond [-180, 180], since the images can have pixels
  // outside of this range.
  settings.simulationBox         = Vector4(-360.0, 360.0, -360.0, 360.0);
  settings.reflectanceType       = -1; // invalid on purpose, must be set later
  settings.saveReflectance       = 0;
  settings.slopeType             = 0;
  settings.initDEM               = 0;
  settings.initExposure          = 0;
  settings.initAlbedo            = 0;
  settings.shadowRemovalType     = CONSTANT_THRESHOLD_SHADOW_REMOVAL;
  settings.shadowThresh          = -1; // invalid on purpose, must be set later
  settings.exposureInfoFilename  = "";
  //settings.TRConst             = 1.24;
  settings.TRConst               = 1.0;
  settings.updateAlbedo          = 0;
  settings.updateExposure        = 0;
  settings.updateHeight          = 0;
  // Two parameters used in the formula for the reflectance
  // Initialize the phase coefficients. They will be updated later.
  settings.phaseCoeffA1          = 1.4;
  settings.phaseCoeffA2          = 0.5;
  settings.updatePhaseCoeffs     = 0;
  settings.updateTilePhaseCoeffs = 0;
  settings.phaseCoeffsFileName   = "";
  settings.useWeights            = 0;
  settings.saveWeights           = 0;
  settings.useNormalizedWeights  = 0;
  settings.postScaleAlbedo       = 0;
  settings.computeWeightsSum     = 0;
  settings.maxNumIter            = 0;
  settings.computeErrors         = 0;
  settings.noDEMDataValue        = 0;

#define CHECK_VAR(name, fmt, assignTo)          \
  if (0 == strcmp(inName, name)) {              \
    sscanf(inVal, fmt, &(settings.assignTo));   \
  }
  
#define CHECK_STR(name, fmt, assignTo)                          \
  if (0 == strcmp(inName, name)) {                              \
    sscanf(inVal, fmt, buffer); settings.assignTo = buffer;     \
  }
  
  if (!configFile.is_open()){
    std::cout << "ERROR: Config file " << config_filename << " not found."<< std::endl;
    exit(1);
  }
  
  printf("CONFIG FILE FOUND\n");
  
  while (!configFile.eof()) {
    configFile.getline(line, MAX_LENGTH);
    
    // truncate comments
    commentPos = strchr(line, '#');
    if (NULL != commentPos) {
      *commentPos = '\0';
    }

    ret = sscanf(line, "%s %s\n", inName, inVal);
    if (ret < 2) continue;

    // Files/directories
    CHECK_STR("DRG_DIR",                  "%s", drgDir);
    CHECK_STR("DEM_DIR",                  "%s", demDir);
    CHECK_STR("SUN_POSITION_FILE",        "%s", sunPosFile);
    CHECK_STR("SPACECRAFT_POSITION_FILE", "%s", spacecraftPosFile);

    // Constants
    //CHECK_VAR("USE_TILES",              "%d", useTiles); // hidden from user, true by default
    CHECK_VAR("TILE_SIZE",                "%f", tileSize);
    //CHECK_VAR("PIXEL_PADDING",          "%d", pixelPadding); // internal variable
    extractSimBox(line, settings.simulationBox);
    CHECK_VAR("REFLECTANCE_TYPE",         "%d", reflectanceType);
    CHECK_VAR("SHADOW_THRESH",            "%f", shadowThresh);
    CHECK_VAR("SHADOW_REMOVAL_TYPE",      "%d", shadowRemovalType);
    CHECK_VAR("TR_CONST",                 "%f", TRConst);
    CHECK_VAR("PHASE_COEFF_A1",           "%f", phaseCoeffA1);
    CHECK_VAR("PHASE_COEFF_A2",           "%f", phaseCoeffA2);
    CHECK_VAR("MAX_NUM_ITER",             "%d", maxNumIter);
    CHECK_VAR("NO_DEM_DATA_VAL",          "%d", noDEMDataValue);

    // Actions
    CHECK_VAR("USE_WEIGHTS",              "%d", useWeights);
    CHECK_VAR("USE_NORMALIZED_WEIGHTS",   "%d", useNormalizedWeights);
    CHECK_VAR("POST_SCALE_ALBEDO",        "%d", postScaleAlbedo);
    //CHECK_VAR("UPDATE_HEIGHT",          "%d", updateHeight); // handled via cmd-line option
    //CHECK_VAR("COMPUTE_ERRORS",         "%d", computeErrors); // handled via cmd-line option

    // Parameters controlling the flow (optional)
    CHECK_VAR("INITIAL_SETUP",            "%d", initialSetup);
    CHECK_VAR("SAVE_WEIGHTS",             "%d", saveWeights);
    CHECK_VAR("COMPUTE_WEIGHTS_SUM",      "%d", computeWeightsSum);
    CHECK_VAR("INIT_DEM",                 "%d", initDEM);
    CHECK_VAR("INIT_EXPOSURE",            "%d", initExposure);
    CHECK_VAR("INIT_ALBEDO",              "%d", initAlbedo);
    CHECK_VAR("UPDATE_EXPOSURE",          "%d", updateExposure);
    CHECK_VAR("UPDATE_TILE_PHASE_COEFFS", "%d", updateTilePhaseCoeffs);
    CHECK_VAR("UPDATE_ALBEDO",            "%d", updateAlbedo);
    CHECK_VAR("UPDATE_PHASE_COEFFS",      "%d", updatePhaseCoeffs);
    
    // Potentially obsolete
    CHECK_VAR("SAVE_REFLECTANCE",         "%d", saveReflectance);
    CHECK_VAR("SLOPE_TYPE",               "%d", slopeType);
  }

  configFile.close();

  // Validation

  if (settings.reflectanceType != NO_REFL &&
      settings.reflectanceType != LAMBERT &&
      settings.reflectanceType != LUNAR_LAMBERT){
    std::cerr << "Expecting one of the following types of reflectance: "
              << "no reflectance ("                  << NO_REFL       << "), "
              << "Lambertian reflectance ("          << LAMBERT       << "), "
              << "or Lunar-Lambertian reflectance (" << LUNAR_LAMBERT << ")" << std::endl;
    exit(1);
  }

  // When we want to do adaptive threshold shadow removal, we need to compute the reflectance.
  // If we do just mosaic, that reflectance can be computed using either the Lambertian or
  // Lunar-Lambertian models. If we do albedo, the reflectance used for adaptive shadow removal
  // must use the same model as the reflectance used for albedo.
  if  (settings.reflectanceType   == LAMBERT &&
       settings.shadowRemovalType == LUNAR_LAMBERTIAN_THRESHOLD_SHADOW_REMOVAL){
    std::cerr << "ERROR: The reflectance type is Lambertian. Then the value of SHADOW_REMOVAL_TYPE cannot be "
              << "Lunar-Lambertian (" << LUNAR_LAMBERTIAN_THRESHOLD_SHADOW_REMOVAL << ")" << endl;
    exit(1);
  }
  if  (settings.reflectanceType   == LUNAR_LAMBERT &&
       settings.shadowRemovalType == LAMBERTIAN_THRESHOLD_SHADOW_REMOVAL){
    std::cerr << "ERROR: The reflectance type is Lunar-Lambertian. Then the value of SHADOW_REMOVAL_TYPE cannot be "
              << "Lambertian (" << LAMBERTIAN_THRESHOLD_SHADOW_REMOVAL << ")" << endl;
    exit(1);
  }

  // If the user wants to do mosaic with adaptive threshold, infer
  // what kind of reflectance we need but never forget that at the end
  // what we want is mosaic, not albedo, so don't use that reflectance
  // in the final computation. This logic must happen after
  // validating the shadowRemovalType field.
  settings.forceMosaic = 0;
  if (settings.reflectanceType == NO_REFL){
    if (settings.shadowRemovalType == LAMBERTIAN_THRESHOLD_SHADOW_REMOVAL){
      settings.reflectanceType = LAMBERT;
      settings.forceMosaic = 1;
    }
    if (settings.shadowRemovalType == LUNAR_LAMBERTIAN_THRESHOLD_SHADOW_REMOVAL){
      settings.reflectanceType = LUNAR_LAMBERT;
      settings.forceMosaic = 1;
    }
  }
  
  if ( !fs::exists(settings.drgDir) ){
    std::cerr << "ERROR: Directory " << settings.drgDir << " does not exist." << std::endl;
    exit(1);
  }
  if ( settings.reflectanceType != NO_REFL && !fs::exists(settings.demDir) ){
    std::cerr << "ERROR: Directory " << settings.demDir << " does not exist." << std::endl;
    exit(1);
  }
  if (settings.reflectanceType != NO_REFL && !fs::exists(settings.sunPosFile) ){
    std::cerr << "ERROR: File " << settings.sunPosFile << " does not exist." << std::endl;
    exit(1);
  }
  if ( settings.reflectanceType != NO_REFL && !fs::exists(settings.spacecraftPosFile) ){
    std::cerr << "ERROR: File " << settings.spacecraftPosFile << " does not exist." << std::endl;
    exit(1);
  }
  if (settings.useTiles != 0 && settings.tileSize <= 0.0){
    std::cerr << "ERROR: The tile size must be positive." << std::endl;
    exit(1);
  }
  if (settings.useTiles != 0 && settings.pixelPadding < 0){
    std::cerr << "ERROR: The pixel padding must be non-negative." << std::endl;
    exit(1);
  }
  if (settings.shadowThresh < 0.0 || settings.shadowThresh > 255.0){
    std::cerr << "ERROR: The shadow threshold must be between 0 and 255." << std::endl;
    exit(1);
  }
    
  return 0;
}

void vw::photometry::PrintGlobalParams(GlobalParams& settings){

  
  // Files/directories
  printf("DRG_DIR                      %s\n", settings.drgDir.c_str());
  printf("DEM_DIR                      %s\n", settings.demDir.c_str());
  printf("SUN_POSITION_FILE            %s\n", settings.sunPosFile.c_str());
  printf("SPACECRAFT_POSITION_FILE     %s\n", settings.spacecraftPosFile.c_str());

  // Constants
  //printf("USE_TILES                  %d\n", settings.useTiles);  // hidden from user, true by default
  printf("TILE_SIZE                    %f\n", settings.tileSize);
  //printf("PIXEL_PADDING              %d\n", settings.pixelPadding); // internal variable
  Vector4 s = settings.simulationBox;
  printf("SIMULATION_BOX               %f:%f:%f:%f\n", s[0], s[1], s[2], s[3]);
  printf("REFLECTANCE_TYPE             %d\n", settings.reflectanceType);
  printf("SHADOW_THRESH                %f\n", settings.shadowThresh);
  printf("SHADOW_REMOVAL_TYPE          %d\n", settings.shadowRemovalType);
  printf("TR_CONST                     %f\n", settings.TRConst);
  printf("PHASE_COEFF_A1               %f\n", settings.phaseCoeffA1);
  printf("PHASE_COEFF_A2               %f\n", settings.phaseCoeffA2);
  printf("MAX_NUM_ITER                 %d\n", settings.maxNumIter);
  printf("NO_DEM_DATA_VAL              %d\n", settings.noDEMDataValue);

  // Actions
  printf("USE_WEIGHTS                  %d\n", settings.useWeights);
  printf("USE_NORMALIZED_WEIGHTS       %d\n", settings.useNormalizedWeights);
  printf("POST_SCALE_ALBEDO            %d\n", settings.postScaleAlbedo);
  //printf("UPDATE_HEIGHT              %d\n", settings.updateHeight);
  //printf("COMPUTE_ERRORS             %d\n", settings.computeErrors); // handled via cmd-line option

  // Parameters controlling the flow
  printf("INITIAL_SETUP                %d\n", settings.initialSetup);
  printf("SAVE_WEIGHTS                 %d\n", settings.saveWeights);
  printf("COMPUTE_WEIGHTS_SUM          %d\n", settings.computeWeightsSum);
  printf("INIT_DEM                     %d\n", settings.initDEM);
  printf("INIT_EXPOSURE                %d\n", settings.initExposure);
  printf("INIT_ALBEDO                  %d\n", settings.initAlbedo);
  printf("UPDATE_EXPOSURE              %d\n", settings.updateExposure);
  printf("UPDATE_TILE_PHASE_COEFFS     %d\n", settings.updateTilePhaseCoeffs);
  printf("UPDATE_PHASE_COEFFS          %d\n", settings.updatePhaseCoeffs);
  printf("UPDATE_ALBEDO                %d\n", settings.updateAlbedo);

  // Potentially obsolete
  printf("SAVE_REFLECTANCE             %d\n", settings.saveReflectance);
  printf("SLOPE_TYPE                   %d\n", settings.slopeType);
  
  return;
}

bool vw::photometry::readImagesFile(std::vector<ImageRecord>& images,
                                    const std::string& imagesListName){
  std::ifstream imagesList(imagesListName.c_str());
  if (!imagesList) {
    std::cerr << "ERROR: readImagesFile: can't open " << imagesListName
              << " for reading: " << strerror(errno) << std::endl;
    return false;
  }

  images.clear();

  std::string line;
  int lineNo = 1;
  while (std::getline(imagesList, line)) {
    char path[1024];
    int useImage;
    ImageRecord rec;

    // ignore blank lines
    if (line.size() == 0) continue;
    // ignore comment lines
    if ('#' == line[0]) continue;

    if (6 != sscanf(line.c_str(), "%d %1023s %lf %lf %lf %lf",
                    &useImage,
                    path,
                    &rec.west,
                    &rec.east,
                    &rec.south,
                    &rec.north
                    )
        ) {
      std::cerr << "ERROR: readImagesFile: " << imagesListName
                << ": line " << lineNo
                << ": expected '%d %s %f %f %f %f' format" << std::endl;
      return false;
    }
    rec.useImage = useImage;
    rec.path = path;
    
    images.push_back(rec);

    lineNo++;
  }

  return true;
}

void vw::photometry::list_DRG_in_box_and_all_DEM(bool useTiles, bool useReflectance,
                                                 std::string allDRGIndex, std::string allDEMIndex,
                                                 Vector4 simulationBox, 
                                                 std::string DRGDir,  std::string DEMDir, 
                                                 std::string DRGInBoxList
                                                 ){

  // Create the lists of ALL DRG and DEM images in DRGDir and
  // DEMDir, if these lists don't exist already.
  
  // Create the list of all DRG files intersecting the current simulationBox.

  Vector4 bigBox = Vector4(-360.0, 360.0, -360.0, 360.0);

  // Create the index of all DRG images if it does not exist already.
  std::vector<ImageRecord> imageRecords;
  if (!readImagesFile(imageRecords, allDRGIndex)){
    std::cout << "WILL create the file " << allDRGIndex << std::endl;
    listTifsInDirOverlappingWithBox(DRGDir, bigBox, allDRGIndex);
    if (!readImagesFile(imageRecords, allDRGIndex)) exit(1); // Second attempt at reading
  }

  // Create the list of all DRG files intersecting the current box.
  ofstream fh(DRGInBoxList.c_str());
  if (!fh){
    std::cerr << "ERROR: list_DRG_in_box_and_all_DEM: can't open " << DRGInBoxList
              << " for writing" << std::endl;
    exit(1);
  }
  fh.precision(20);
  for (int j = 0; j < (int)imageRecords.size(); j++){
    const ImageRecord& rec = imageRecords[j];
    Vector4 currCorners = Vector4(rec.west, rec.east, rec.south, rec.north);
    if (! boxesOverlap(currCorners, simulationBox)) continue;
    fh << rec.useImage << " " << rec.path << " "
       << currCorners(0) << " " << currCorners(1) << " "
       << currCorners(2) << " " << currCorners(3) << std::endl;
  }
  fh.close();

  if (!useTiles || !useReflectance){
    return;
  }
  
  // Create the index of all DEM tiles if it does not exist already.
  std::vector<ImageRecord> DEMTilesRecords;
  if (!readImagesFile(DEMTilesRecords, allDEMIndex)){
    listTifsInDirOverlappingWithBox(DEMDir, bigBox, allDEMIndex);
  }

  return;
}
