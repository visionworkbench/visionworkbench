// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Cartography.h>
#include <vw/Mosaic.h>
#include <vw/Plate/PlateFile.h>
#include <vw/Plate/ToastPlateManager.h>

using namespace vw;
using namespace vw::platefile;
using namespace vw::mosaic;
using namespace vw::cartography;

#include <boost/program_options.hpp>
namespace po = boost::program_options;

// Erases a file suffix if one exists and returns the base string
static std::string prefix_from_filename(std::string const& filename) {
  std::string result = filename;
  int index = result.rfind(".");
  if (index != -1) 
    result.erase(index, result.size());
  return result;
}

int main( int argc, char *argv[] ) {

  std::string output_file_name;
  std::string output_file_type;
  int tile_size;
  float jpeg_quality;
  int png_compression;
  unsigned cache_size;
  unsigned num_threads;
  std::vector<std::string> image_files;

  po::options_description general_options("Turns georeferenced image(s) into a TOAST quadtree.\n\nGeneral Options");
  general_options.add_options()
    ("output-name,o", po::value<std::string>(&output_file_name), "Specify the base output directory")
    ("file-type", po::value<std::string>(&output_file_type)->default_value("png"), "Output file type")
    ("tile-size", po::value<int>(&tile_size)->default_value(256), "Tile size, in pixels")
    ("jpeg-quality", po::value<float>(&jpeg_quality)->default_value(0.75), "JPEG quality factor (0.0 to 1.0)")
    ("png-compression", po::value<int>(&png_compression)->default_value(3), "PNG compression level (0 to 9)")
    ("cache", po::value<unsigned>(&cache_size)->default_value(1024), "Soure data cache size, in megabytes")
    ("num-threads,t", po::value<unsigned>(&num_threads)->default_value(1), "Set number of threads (set to 0 to use system default")
    ("help", "Display this help message");

  po::options_description hidden_options("");
  hidden_options.add_options()
    ("input-file", po::value<std::vector<std::string> >(&image_files));

  po::options_description options("Allowed Options");
  options.add(general_options).add(hidden_options);

  po::positional_options_description p;
  p.add("input-file", -1);

  po::variables_map vm;
  po::store( po::command_line_parser( argc, argv ).options(options).positional(p).run(), vm );
  po::notify( vm );

  std::ostringstream usage;
  usage << "Usage: " << argv[0] << " [options] <filename>..." <<std::endl << std::endl;
  usage << general_options << std::endl;

  if( vm.count("help") ) {
    std::cout << usage.str();
    return 0;
  }

  if( vm.count("input-file") < 1 ) {
    std::cerr << "Error: must specify at least one input file!" << std::endl << std::endl;
    std::cout << usage.str();
    return 1;
  }

  if( output_file_name == "" )
    output_file_name = prefix_from_filename(image_files[0]) + ".plate";

  if( tile_size <= 0 || tile_size != pow(2,int(log(tile_size)/log(2))) ) {
    std::cerr << "Error: The tile size must be a power of two!  (You specified: " << tile_size << ")" << std::endl << std::endl;
    std::cout << usage.str();
    return 1;
  }

  DiskImageResourceJPEG::set_default_quality( jpeg_quality );
  DiskImageResourcePNG::set_default_compression_level( png_compression );
  vw_settings().set_system_cache_size( cache_size*1024*1024 );

  std::vector<DiskImageView<PixelRGBA<uint8> > > images;
  std::vector<GeoReference> georefs;
  int max_level = 0;

  std::cout << "Scanning input files...." << std::endl;
  for( unsigned i=0; i<image_files.size(); ++i ) {
    DiskImageView<PixelRGBA<uint8> > image(image_files[i]);
    images.push_back(image);

    GeoReference georef;
    read_georeference( georef, DiskImageResourceGDAL(image_files[i]) );
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
    int level = (int)round(log(360/delta/(tile_size-1))/log(2));
    if( level > max_level ) max_level = level;
  }
  
  // This is the right dimension if we create tiles that overlap by
  // one pixel on one side.  We do this in the TOAST projection to
  // correctly render tiles in a 3D texturing environment.
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
      vw_out(0) << "\t--> Detected global overlay.  Using cylindrical edge extension to hide the seam.\n";

      composite.insert(transform(image,toast_tx,resolution,resolution,CylindricalEdgeExtension(),BicubicInterpolation()), 0, 0);
    }
    else {
      composite.insert(transform(image,toast_tx,resolution,resolution,ZeroEdgeExtension(),BicubicInterpolation()), 0, 0);
    }
  }
  composite.prepare();

  // -------------------------- PLATE FILE GENERATION --------------------------------

  // Create the plate file
  boost::shared_ptr<PlateFile> platefile = boost::shared_ptr<PlateFile>( new PlateFile(output_file_name, tile_size, output_file_type));
  ToastPlateManager pm(platefile, num_threads);

  std::cout << "\nWriting data to plate file: " << output_file_type << "\n";
  std::cout << "\t--> Building full-resolution tiles\n";
  pm.insert(composite, max_level);
  
  std::cout << "\t--> Building mipmap levels\n";
  pm.mipmap();
}
