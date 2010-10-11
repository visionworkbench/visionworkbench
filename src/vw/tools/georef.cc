// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iostream>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <vw/Core/Cache.h>
#include <vw/Core/ProgressCallback.h>
#include <vw/Math/Matrix.h>
#include <vw/Image/Palette.h>
#include <vw/Image/Transform.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/FileIO/DiskImageResourceGDAL.h>
#include <vw/FileIO/DiskImageView.h>
#include <vw/Cartography/GeoReference.h>

using namespace vw;
using namespace vw::math;
using namespace vw::cartography;

int main( int argc, char *argv[] ) {

  std::string input_filename, output_filename, copy_filename, tfw_filename;
  double north_lat=90.0, south_lat=-90.0;
  double east_lon=180.0, west_lon=-180.0;
  double proj_lat=0, proj_lon=0, proj_scale=1;
  unsigned utm_zone;
  double nudge_x=0, nudge_y=0;

  po::options_description general_options("General Options");
  general_options.add_options()
    ("output-file,o", po::value<std::string>(&output_filename)->default_value("output.tif"), "Specify the base output filename")
    ("help,h", "Display this help message");

  po::options_description projection_options("Projection Options");
  projection_options.add_options()
    ("copy", po::value<std::string>(&copy_filename), "Copy the projection from the given file")
    ("tfw", po::value<std::string>(&tfw_filename), "Create a .tfw sidecar file with the given filename rather than a full copy of the image file")
    ("north", po::value<double>(&north_lat), "The northernmost latitude in degrees")
    ("south", po::value<double>(&south_lat), "The southernmost latitude in degrees")
    ("east", po::value<double>(&east_lon), "The easternmost longitude in degrees")
    ("west", po::value<double>(&west_lon), "The westernmost longitude in degrees")
    ("sinusoidal", "Assume a sinusoidal projection")
    ("mercator", "Assume a Mercator projection")
    ("transverse-mercator", "Assume a transverse Mercator projection")
    ("orthographic", "Assume an orthographic projection")
    ("stereographic", "Assume a stereographic projection")
    ("lambert-azimuthal", "Assume a Lambert azimuthal projection")
    ("utm", po::value<unsigned>(&utm_zone), "Assume UTM projection with the given zone")
    ("proj-lat", po::value<double>(&proj_lat), "The center of projection latitude (if applicable)")
    ("proj-lon", po::value<double>(&proj_lon), "The center of projection longitude (if applicable)")
    ("proj-scale", po::value<double>(&proj_scale), "The projection scale (if applicable)")
    ("nudge-x", po::value<double>(&nudge_x), "Nudge the image, in projected coordinates")
    ("nudge-y", po::value<double>(&nudge_y), "Nudge the image, in projected coordinates")
    ("pixel-as-point", "Encode that the pixel location (0,0) is the center of the upper left hand pixel (the default, if you specify nothing, is to set the upper left hand corner of the upper left pixel as (0,0) (i.e. PixelAsArea).");

  po::options_description hidden_options("");
  hidden_options.add_options()
    ("input-file", po::value<std::string>(&input_filename));

  po::options_description options("Allowed Options");
  options.add(general_options).add(projection_options).add(hidden_options);

  po::positional_options_description p;
  p.add("input-file", -1);

  std::ostringstream usage;
  usage << "Description: Specify planetary coordinates for an image" << std::endl << std::endl;
  usage << "Usage: " << argv[0] << " [options] <filename>..." << std::endl << std::endl;
  usage << general_options << std::endl;
  usage << projection_options << std::endl;

  po::variables_map vm;
  try {
    po::store( po::command_line_parser( argc, argv ).options(options).positional(p).run(), vm );
    po::notify( vm );
  } catch(po::error &e) {
    std::cout << "An error occured while parsing command line arguments.\n";
    std::cout << "\t" << e.what() << "\n\n";
    std::cout << usage.str();
    return 1;
  }

  if( vm.count("help") ) {
    std::cout << usage.str();
    return 1;
  }

  if( vm.count("input-file") < 1 ) {
    std::cout << "Error: Must specify at least one input file!" << std::endl << std::endl;
    std::cout << usage.str();
    return 1;
  }

  GeoReference output_georef;
  output_georef.set_well_known_geogcs("WGS84");

  // Read in georeference info and compute total resolution
  bool manual = vm.count("north") || vm.count("south") || vm.count("east") || vm.count("west");

  DiskImageResourceGDAL file_resource( input_filename );
  GeoReference georef;
  if( vm.count("copy") ) {
    read_georeference( georef, copy_filename );
  } else {
    read_georeference( georef, file_resource );
  }

  if ( georef.proj4_str() == "" ) georef.set_well_known_geogcs("WGS84");
  if( manual || georef.transform() == identity_matrix<3>() ) {
    if( manual ) {
      vw_out() << "Using manual Plate Carree coordinates: ";
      georef = GeoReference(georef.datum());
    } else {
      vw_out() << "No georeferencing info found.  Assuming Plate Carree WGS84: ";
      georef = GeoReference();
      georef.set_well_known_geogcs("WGS84");
    }
    vw_out() << east_lon << " E to " << west_lon << " W, "
             << south_lat << " S to " << north_lat << " N." << std::endl;

    Matrix3x3 m;
    m(0,0) = (east_lon - west_lon) / file_resource.cols();
    m(0,2) = west_lon;
    m(1,1) = (south_lat - north_lat) / file_resource.rows();
    m(1,2) = north_lat;
    m(2,2) = 1;
    georef.set_transform( m );
    manual = true;
  }
  else if( vm.count("sinusoidal") ) georef.set_sinusoidal(proj_lon);
  else if( vm.count("mercator") ) georef.set_mercator(proj_lat,proj_lon,proj_scale);
  else if( vm.count("transverse-mercator") ) georef.set_transverse_mercator(proj_lat,proj_lon,proj_scale);
  else if( vm.count("orthographic") ) georef.set_orthographic(proj_lat,proj_lon);
  else if( vm.count("stereographic") ) georef.set_stereographic(proj_lat,proj_lon,proj_scale);
  else if( vm.count("lambert-azimuthal") ) georef.set_lambert_azimuthal(proj_lat,proj_lon);
  else if( vm.count("utm") ) georef.set_UTM( utm_zone );

  if (vm.count("pixel-as-point"))
    georef.set_pixel_interpretation(GeoReference::PixelAsPoint);
  else // Default: PixelAsArea
    georef.set_pixel_interpretation(GeoReference::PixelAsArea);

  if( vm.count("nudge-x") || vm.count("nudge-y") ) {
    Matrix3x3 m = georef.transform();
    m(0,2) += nudge_x;
    m(1,2) += nudge_y;
    georef.set_transform( m );
  }

  vw_out() << "Writing file with Proj4 String: " << georef.proj4_str() << "\n";

  // Our file readers do a better job that GDAL's of coping with large
  // image files in some cases, so we make a fresh DiskImageView rather
  // than making an ImageResourceView of the existing GDAL resource.

  if (vm.count("tfw")) {
    std::ofstream tfw_file(tfw_filename.c_str());
    tfw_file << georef.transform()(0,0) << "\n";
    tfw_file << georef.transform()(0,1) << "\n";
    tfw_file << georef.transform()(1,0) << "\n";
    tfw_file << georef.transform()(1,1) << "\n";
    if (vm.count("pixel-as-point")) {
      tfw_file << georef.transform()(0,2) << "\n";
      tfw_file << georef.transform()(1,2) << "\n";
    } else {
      tfw_file << georef.transform()(0,2) + 0.5*georef.transform()(0,0) << "\n";
      tfw_file << georef.transform()(1,2) + 0.5*georef.transform()(1,1) << "\n";
    }
    tfw_file.close();
  } else {
    TerminalProgressCallback bar( "tools.georef", "Writing:" );
    switch( file_resource.channel_type() ) {
    case VW_CHANNEL_INT16:
      switch( file_resource.pixel_format() ) {
      case VW_PIXEL_SCALAR: {
        DiskImageView<int16> input_image( input_filename );
        write_georeferenced_image( output_filename, input_image, georef, bar );
        break;
      }
      case VW_PIXEL_GRAY: {
        DiskImageView<PixelGray<int16> > input_image( input_filename );
        write_georeferenced_image( output_filename, input_image, georef, bar );
        break;
      }
      case VW_PIXEL_GRAYA: {
        DiskImageView<PixelGrayA<int16> > input_image( input_filename );
        write_georeferenced_image( output_filename, input_image, georef, bar );
        break;
      }
      case VW_PIXEL_RGB: {
        DiskImageView<PixelRGB<int16> > input_image( input_filename );
        write_georeferenced_image( output_filename, input_image, georef, bar );
        break;
      }
      case VW_PIXEL_RGBA: {
        DiskImageView<PixelRGBA<int16> > input_image( input_filename );
        write_georeferenced_image( output_filename, input_image, georef, bar );
        break;
      }
      default: {
        vw_throw( NoImplErr() << "Unsupported pixel format: " << file_resource.pixel_format() );
        break;
      }
      }
      break;
    case VW_CHANNEL_UINT16:
      switch( file_resource.pixel_format() ) {
      case VW_PIXEL_SCALAR: {
        DiskImageView<uint16> input_image( input_filename );
        write_georeferenced_image( output_filename, input_image, georef, bar );
        break;
      }
      case VW_PIXEL_GRAY: {
        DiskImageView<PixelGray<uint16> > input_image( input_filename );
        write_georeferenced_image( output_filename, input_image, georef, bar );
        break;
      }
      case VW_PIXEL_GRAYA: {
        DiskImageView<PixelGrayA<uint16> > input_image( input_filename );
        write_georeferenced_image( output_filename, input_image, georef, bar );
        break;
      }
      case VW_PIXEL_RGB: {
        DiskImageView<PixelRGB<uint16> > input_image( input_filename );
        write_georeferenced_image( output_filename, input_image, georef, bar );
        break;
      }
      case VW_PIXEL_RGBA: {
        DiskImageView<PixelRGBA<uint16> > input_image( input_filename );
        write_georeferenced_image( output_filename, input_image, georef, bar );
        break;
      }
      default: {
        vw_throw( NoImplErr() << "Unsupported pixel format: " << file_resource.pixel_format() );
        break;
      }
      }
      break;
    case VW_CHANNEL_FLOAT32:
      switch( file_resource.pixel_format() ) {
      case VW_PIXEL_SCALAR: {
        DiskImageView<float32> input_image( input_filename );
        write_georeferenced_image( output_filename, input_image, georef, bar );
        break;
      }
      case VW_PIXEL_GRAY: {
        DiskImageView<PixelGray<float32> > input_image( input_filename );
        write_georeferenced_image( output_filename, input_image, georef, bar );
        break;
      }
      case VW_PIXEL_GRAYA: {
        DiskImageView<PixelGrayA<float32> > input_image( input_filename );
        write_georeferenced_image( output_filename, input_image, georef, bar );
        break;
      }
      case VW_PIXEL_RGB: {
        DiskImageView<PixelRGB<float32> > input_image( input_filename );
        write_georeferenced_image( output_filename, input_image, georef, bar );
        break;
      }
      case VW_PIXEL_RGBA: {
        DiskImageView<PixelRGBA<float32> > input_image( input_filename );
        write_georeferenced_image( output_filename, input_image, georef, bar );
        break;
      }
      default: {
        vw_throw( NoImplErr() << "Unsupported pixel format: " << file_resource.pixel_format() );
        break;
      }
      }
      break;
    default:
      switch( file_resource.pixel_format() ) {
      case VW_PIXEL_SCALAR: {
        DiskImageView<uint8> input_image( input_filename );
        write_georeferenced_image( output_filename, input_image, georef, bar );
        break;
      }
      case VW_PIXEL_GRAY: {
        DiskImageView<PixelGray<uint8> > input_image( input_filename );
        write_georeferenced_image( output_filename, input_image, georef, bar );
        break;
      }
      case VW_PIXEL_GRAYA: {
        DiskImageView<PixelGrayA<uint8> > input_image( input_filename );
        write_georeferenced_image( output_filename, input_image, georef, bar );
        break;
      }
      case VW_PIXEL_RGB: {
        DiskImageView<PixelRGB<uint8> > input_image( input_filename );
        write_georeferenced_image( output_filename, input_image, georef, bar );
        break;
      }
      case VW_PIXEL_RGBA: {
        DiskImageView<PixelRGBA<uint8> > input_image( input_filename );
        write_georeferenced_image( output_filename, input_image, georef, bar );
        break;
      }
      default: {
        vw_throw( NoImplErr() << "Unsupported pixel format: " << file_resource.pixel_format() );
        break;
      }
      }
      break;
    }
  }
  return 0;
}
