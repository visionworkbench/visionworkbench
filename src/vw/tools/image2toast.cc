// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NASA Vision Workbench is licensed under the Apache License,
//  Version 2.0 (the "License"); you may not use this file except in
//  compliance with the License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__


#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <vw/Core/Log.h>
#include <vw/Core/ProgressCallback.h>
#include <vw/Math/BBox.h>
#include <vw/Math/Matrix.h>
#include <vw/Math/Vector.h>
#include <vw/Image/AlgorithmFunctions.h>
#include <vw/Image/EdgeExtension.h>
#include <vw/Image/ImageIO.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Image/Interpolation.h>
#include <vw/Image/Manipulation.h>
#include <vw/Image/Transform.h>
#include <vw/FileIO/DiskImageView.h>
#include <vw/FileIO/DiskImageResourceJPEG.h>
#include <vw/FileIO/DiskImageResourcePNG.h>
#include <vw/FileIO/DiskImageResourceGDAL.h>
#include <vw/Cartography/GeoReference.h>
#include <vw/Cartography/ToastTransform.h>
#include <vw/Mosaic/ImageComposite.h>
#include <vw/Mosaic/QuadTreeGenerator.h>
#include <vw/Mosaic/ToastQuadTreeConfig.h>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

using namespace vw;
using namespace vw::cartography;
using namespace vw::mosaic;

int main(int argc, char **argv) {

  std::string output_file_name;
  std::string output_file_type;
  int tile_size;
  float jpeg_quality;
  int png_compression;
  std::vector<std::string> image_files;

  po::options_description general_options("Turns georeferenced image(s) into a TOAST quadtree.\n\nGeneral Options");
  general_options.add_options()
    ("output-name,o", po::value<std::string>(&output_file_name), "Specify the base output directory")
    ("file-type", po::value<std::string>(&output_file_type)->default_value("png"), "Output file type")
    ("tile-size", po::value<int>(&tile_size)->default_value(256), "Tile size, in pixels")
    ("jpeg-quality", po::value<float>(&jpeg_quality)->default_value(0.75), "JPEG quality factor (0.0 to 1.0)")
    ("png-compression", po::value<int>(&png_compression)->default_value(3), "PNG compression level (0 to 9)")
    ("help,h", "Display this help message");

  po::options_description hidden_options("");
  hidden_options.add_options()
    ("input-file", po::value<std::vector<std::string> >(&image_files));

  po::options_description options("Allowed Options");
  options.add(general_options).add(hidden_options);

  po::positional_options_description p;
  p.add("input-file", -1);

  std::ostringstream usage;
  usage << "Usage: toast [options] <filename>..." <<std::endl << std::endl;
  usage << general_options << std::endl;

  po::variables_map vm;
  try {
    po::store( po::command_line_parser( argc, argv ).options(options).positional(p).run(), vm );
    po::notify( vm );
  } catch (const po::error& e) {
    std::cout << "An error occured while parsing command line arguments.\n";
    std::cout << "\t" << e.what() << "\n\n";
    std::cout << usage.str();
    return 1;
  }

  if( vm.count("help") ) {
    std::cout << usage.str();
    return 0;
  }

  if( vm.count("input-file") < 1 ) {
    std::cerr << "Error: must specify at least one input file!" << std::endl << std::endl;
    std::cout << usage.str();
    return 1;
  }

  if( output_file_name.empty() )
    output_file_name = fs::path(image_files[0]).replace_extension("toast").string();

  if( tile_size <= 0 || tile_size != pow(2.,floor(log((double)tile_size)/log(2.))) ) {
    std::cerr << "Error: The tile size must be a power of two!  (You specified: " << tile_size << ")" << std::endl << std::endl;
    std::cout << usage.str();
    return 1;
  }

  DiskImageResourceJPEG::set_default_quality( jpeg_quality );
  DiskImageResourcePNG::set_default_compression_level( png_compression );

  std::vector<DiskImageView<PixelRGBA<uint8> > > images;
  std::vector<GeoReference> georefs;
  int max_level = 0;

  vw_out() << "Scanning input files...." << std::endl;
  for( unsigned i=0; i<image_files.size(); ++i ) {
    DiskImageView<PixelRGBA<uint8> > image(image_files[i]);
    images.push_back(image);

    GeoReference georef;
    DiskImageResourceGDAL diskrsrc(image_files[i]);
    read_georeference( georef, diskrsrc );
    if( georef.transform() == identity_matrix<3>() ) {
      std::cout << "No georeferencing info found for " << image_files[i] << ".  Assuming global plate carree." << std::endl;
      Matrix3x3 M;
      M(0,0) = 360.0 / image.cols();
      M(0,2) = -180.0;
      M(1,1) = -180.0 / image.rows();
      M(1,2) = 90.0;
      M(2,2) = 1;
      georef.set_transform( M );
    }
    georefs.push_back(georef);

    Vector2 p0 = georef.pixel_to_lonlat(Vector2(image.cols()/2,image.rows()/2));
    Vector2 p1 = georef.pixel_to_lonlat(Vector2(image.cols()/2+1,image.rows()/2));
    Vector2 p2 = georef.pixel_to_lonlat(Vector2(image.cols()/2,image.rows()/2+1));
    double delta = sqrt(pow(p1.y()-p0.y(),2)+pow(p2.y()-p0.y(),2));
    int level = (int)round(log(180/delta/(tile_size-1))/log(2.));
    if( level > max_level ) max_level = level;
  }

  int32 resolution = (1<<max_level)*(tile_size-1)+1;
  std::cout << "Using " << max_level+1 << " levels. (Total resolution = " << resolution << " pixels.)" << std::endl;

  ImageComposite<PixelRGBA<uint8> > composite;
  composite.set_draft_mode(true);
  for( unsigned i=0; i<image_files.size(); ++i ) {
    GeoReference georef = georefs[i];
    DiskImageView<PixelRGBA<uint8> > image = images[i];

    bool global = georef.proj4_str()=="+proj=longlat" &&
      fabs(georef.lonlat_to_pixel(Vector2(-180,0)).x()) < 1 &&
      fabs(georef.lonlat_to_pixel(Vector2(180,0)).x() - image.cols()) < 1 &&
      fabs(georef.lonlat_to_pixel(Vector2(0,90)).y()) < 1 &&
      fabs(georef.lonlat_to_pixel(Vector2(0,-90)).y() - image.rows()) < 1;

    ToastTransform toast_tx( georef, resolution );
    ImageViewRef<PixelRGBA<uint8> > toast_image;

    if( global ) {
      vw_out() << "\t--> Detected global overlay.  Using cylindrical edge extension to hide the seam.\n";
      composite.insert(transform(image,toast_tx,resolution,resolution,CylindricalEdgeExtension(),BicubicInterpolation()), 0, 0);
    }
    else {
      composite.insert(transform(image,toast_tx,resolution,resolution,ZeroEdgeExtension(),BicubicInterpolation()), 0, 0);
    }
  }

  composite.prepare();
  vw_out() << "Composite dimensions: " << composite.cols() << " x " << composite.rows() << "\n";

  vw_out() << "Generating tile set at " << output_file_name << "." << std::endl;
  QuadTreeGenerator qtree( composite, output_file_name );
  ToastQuadTreeConfig tqtc;
  tqtc.configure( qtree, composite );
  qtree.set_file_type( output_file_type );
  qtree.generate( TerminalProgressCallback( "tools.image2toast","") );

  return 0;
}
