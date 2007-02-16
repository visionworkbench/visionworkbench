// __BEGIN_LICENSE__
// 
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
// 
// Copyright 2006 Carnegie Mellon University. All rights reserved.
// 
// This software is distributed under the NASA Open Source Agreement
// (NOSA), version 1.3.  The NOSA has been approved by the Open Source
// Initiative.  See the file COPYING at the top of the distribution
// directory tree for the complete NOSA document.
// 
// THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
// KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
// LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
// SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
// A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
// THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
// DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
// 
// __END_LICENSE__

#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <stdlib.h>
#include <iostream>
#include <fstream>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <vw/Core/Cache.h>
#include <vw/Math/Matrix.h>
#include <vw/Image/Transform.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/FileIO/DiskImageResourceJPEG.h>
#include <vw/FileIO/DiskImageView.h>
#include <vw/Cartography/GeoReference.h>
#include <vw/Cartography/GeoTransform.h>
#include <vw/Cartography/FileIO.h>
#include <vw/Mosaic/ImageComposite.h>
#include <vw/Mosaic/KMLQuadTreeGenerator.h>
using namespace vw;
using namespace vw::math;
using namespace vw::cartography;
using namespace vw::mosaic;

int main( int argc, char *argv[] ) {

  std::vector<std::string> image_files;
  std::string output_file_name;
  std::string output_file_type;
  float north_lat, south_lat;
  float east_lon, west_lon;
  unsigned utm_zone;
  int patch_size, patch_overlap;
  float jpeg_quality;
  unsigned cache_size;
  int max_lod_pixels;

  po::options_description desc("Options");
  desc.add_options()
    ("help", "Display this help message")
    ("input-file", po::value<std::vector<std::string> >(&image_files), "Explicitly specify the input file")
    ("output-name,o", po::value<std::string>(&output_file_name)->default_value("output"), "Specify the base output filename")
    ("file-type", po::value<std::string>(&output_file_type)->default_value("jpg"), "Output file type")
    ("north", po::value<float>(&north_lat)->default_value(90.0), "The northernmost latitude in degrees")
    ("south", po::value<float>(&south_lat)->default_value(-90.0), "The southernmost latitude in degrees")
    ("east", po::value<float>(&east_lon)->default_value(180.0), "The easternmost latitude in degrees")
    ("west", po::value<float>(&west_lon)->default_value(-180.0), "The westernmost latitude in degrees")
    ("utm", po::value<unsigned>(&utm_zone), "Specify a UTM zone")
    ("size", po::value<int>(&patch_size)->default_value(256), "Patch size, in pixels")
    ("overlap", po::value<int>(&patch_overlap)->default_value(0), "Patch overlap, in pixels (must be even)")
    ("jpeg-quality", po::value<float>(&jpeg_quality)->default_value(0.75), "JPEG quality factor (0.0 to 1.0)")
    ("cache", po::value<unsigned>(&cache_size)->default_value(1024), "Cache size, in megabytes")
    ("max-lod-pixels", po::value<int>(&max_lod_pixels)->default_value(-1), "Max LoD in pixels, or -1 for none")
    ("verbose", "Verbose output");
  po::positional_options_description p;
  p.add("input-file", -1);

  po::variables_map vm;
  po::store( po::command_line_parser( argc, argv ).options(desc).positional(p).run(), vm );
  po::notify( vm );

  if( vm.count("help") ) {
    std::cout << desc << std::endl;
    return 1;
  }

  if( vm.count("input-file") < 1 ) {
    std::cout << "Error: Must specify exactly one input file!" << std::endl;
    std::cout << desc << std::endl;
    return 1;
  }

  if( patch_size <= 0 ) {
    std::cerr << "Error: The patch size must be a positive number!  (You specified " << patch_size << ".)" << std::endl;
    std::cout << desc << std::endl;
    return 1;
  }
    
  if( patch_overlap<0 || patch_overlap>=patch_size || patch_overlap%2==1 ) {
    std::cerr << "Error: The patch overlap must be an even number nonnegative number" << std::endl;
    std::cerr << "smaller than the patch size!  (You specified " << patch_overlap << ".)" << std::endl;
    std::cout << desc << std::endl;
    return 1;
  }
    
  if( vm.count("verbose") ) {
    set_debug_level(VerboseDebugMessage);
  }

  DiskImageResourceJPEG::set_default_quality( jpeg_quality );
  Cache::system_cache().resize( cache_size*1024*1024 );

  GeoReference output_georef;
  output_georef.set_well_known_geogcs("WGS84");
  int total_resolution = 1024;

  // Read in georeference info and compute total resolution
  bool manual = false;
  std::vector<GeoReference> georeferences;
  for( unsigned i=0; i<image_files.size(); ++i ) {
    DiskImageResourceGDAL file_resource( image_files[i] );
    GeoReference input_georef;
    file_resource.read_georeference( input_georef );
    // This is sort of a klugey way to check for a georef
    if( input_georef.proj4_str() == "" && input_georef.transform() == identity_matrix<3>() ) {
      if( image_files.size() == 1 ) {
        vw_out(InfoMessage) << "No georeferencing info found.  Assuming Plate Carree WGS84: " 
                            << east_lon << " to " << west_lon << " E, " << south_lat << " to " << north_lat << " N." << std::endl;
        input_georef.set_well_known_geogcs("WGS84");
        Matrix3x3 m;
        m(0,0) = (east_lon - west_lon) / file_resource.cols();
        m(0,2) = west_lon + 0.5*m(0,0);
        m(1,1) = (south_lat - north_lat) / file_resource.rows();
        m(1,2) = north_lat + 0.5*m(1,1);
        m(2,2) = 1;
        input_georef.set_transform( m );
        manual = true;
      }
      else {
        vw_out(ErrorMessage) << "Error: No georeferencing info found for input file \"" << image_files[i] << "\"!" << std::endl;
        vw_out(ErrorMessage) << "(Manually-specified bounds are only allowed for single image files.)" << std::endl;
        exit(1);
      }
    }
    else if( vm.count("utm") ) input_georef.set_UTM( utm_zone );
    georeferences.push_back( input_georef );
    
    GeoTransform geotx( input_georef, output_georef );
    Vector2 center_pixel( file_resource.cols()/2, file_resource.rows()/2 );
    int resolution = GlobalKMLTransform::compute_resolution( geotx, center_pixel );
    if( resolution > total_resolution ) total_resolution = resolution;
  }

  // Configure the composite
  ImageComposite<PixelRGBA<uint8> > composite;
  composite.set_draft_mode( true );
  GlobalKMLTransform kmltx( total_resolution );

  // Add the transformed input files to the composite
  for( unsigned i=0; i<image_files.size(); ++i ) {
    std::cout << "Adding file " << image_files[i] << std::endl;
    GeoTransform geotx( georeferences[i], output_georef );
    DiskImageView<PixelRGBA<uint8> > source( image_files[i] );
    BBox2i bbox = compute_transformed_bbox_fast_int( source, compose(kmltx,geotx) );
    // Constant edge extension is better for transformations that 
    // preserve the rectangularity of the image.  At the moment we 
    // only do this for manual transforms, alas.
    if( manual ) {
      // If the image is being super-sampled the computed bounding 
      // box may be missing a pixel at the edges relative to what 
      // you might expect, which can create visible artifacts if 
      // it happens at the boundaries of the coordinate system.
      if( west_lon == -180 ) bbox.min().x() = 0;
      if( east_lon == 180 ) bbox.max().x() = total_resolution;
      if( north_lat == 90 ) bbox.min().y() = total_resolution/4;
      if( south_lat == -90 ) bbox.max().y() = 3*total_resolution/4;
      composite.insert( crop( transform( source, compose(kmltx,geotx), ConstantEdgeExtension() ), bbox ),
                        bbox.min().x(), bbox.min().y() );
    }
    else {
      composite.insert( crop( transform( source, compose(kmltx,geotx) ), bbox ),
                        bbox.min().x(), bbox.min().y() );
    }
  }

  // Prepare the composite
  std::cout << "Preparing composite..." << std::endl;
  BBox2i total_bbox(0,0,total_resolution,total_resolution);
  composite.prepare( total_bbox );
  BBox2i data_bbox = composite.bbox();
  data_bbox.crop( total_bbox );
  
  // Prepare the quadtree
  KMLQuadTreeGenerator<PixelRGBA<uint8> > quadtree( output_file_name, composite, BBox2(-180,-180,360,360) );
  quadtree.set_max_lod_pixels(max_lod_pixels);
  quadtree.set_crop_bbox( data_bbox );
  quadtree.set_crop_images( true );
  quadtree.set_output_image_file_type( output_file_type );

  // Generate the composite
  std::cout << "Compositing..." << std::endl;
  quadtree.generate();

  return 0;
}
