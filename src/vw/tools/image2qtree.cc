// __BEGIN_LICENSE__
//
// Copyright (C) 2008 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
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

/// \file image2geotree.cc
///
/// This program takes a georeferenced image as its input, and outputs 
/// a quadtree for that image that is viewable in various terrain display
/// programs, such as Google Earth. Currently, the program supports output
/// in KML and TMS format.

#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <map>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <vw/Core/Cache.h>
#include <vw/Core/ProgressCallback.h>
#include <vw/Math/Matrix.h>
#include <vw/Image/Transform.h>
#include <vw/Image/Palette.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/FileIO/DiskImageResourceJPEG.h>
#include <vw/FileIO/DiskImageResourceGDAL.h>
#include <vw/FileIO/DiskImageView.h>
#include <vw/Cartography/GeoReference.h>
#include <vw/Cartography/GeoTransform.h>
#include <vw/Cartography/FileIO.h>
#include <vw/Mosaic/ImageComposite.h>
#include <vw/Mosaic/KMLQuadTreeGenerator.h>
#include <vw/Mosaic/TMSQuadTreeGenerator.h>
#include <vw/Mosaic/UniviewQuadTreeGenerator.h>
using namespace vw;
using namespace vw::math;
using namespace vw::cartography;
using namespace vw::mosaic;

// Global variables
std::vector<std::string> image_files;
std::string output_file_name;
std::string output_file_type;
std::string output_metadata;
std::string module_name;
double north_lat=90.0, south_lat=-90.0;
double east_lon=180.0, west_lon=-180.0;
double proj_lat=0, proj_lon=0, proj_scale=1;
unsigned utm_zone;
int patch_size;
int patch_overlap;
float jpeg_quality;
unsigned cache_size;
std::string palette_file;
std::string channel_type_str;
float palette_scale=1.0, palette_offset=0.0;
int draw_order_offset; // KML only.
int max_lod_pixels; // KML only.
float pixel_scale=1.0, pixel_offset=0.0;
double lcc_parallel1, lcc_parallel2;
int aspect_ratio=1;
bool verbose=false;
bool quiet=false;
bool terrain=false;
// For image stretching.
float lo_value = ScalarTypeLimits<float>::highest();
float hi_value = ScalarTypeLimits<float>::lowest();

// Function pointers for computing resolution.
std::map<std::string, int (*)(const GeoTransform&, const Vector2&)> str_to_resolution_fn_map;

// Fill the maps for converting input strings to function pointers.
static void fill_input_maps() {
  str_to_resolution_fn_map[std::string("none")] = NULL;
  str_to_resolution_fn_map[std::string("kml")] = &vw::cartography::output::kml::compute_resolution;
  str_to_resolution_fn_map[std::string("tms")] = &vw::cartography::output::tms::compute_resolution;
  str_to_resolution_fn_map[std::string("uniview")] = &vw::cartography::output::tms::compute_resolution;
}

static void get_normalize_vals(std::string filename, DiskImageResourceGDAL &file_resource) {
  float no_data_value = file_resource.get_no_data_value(0);
  DiskImageView<PixelRGBA<float> > min_max_file(filename);
  float new_lo, new_hi;
  min_max_channel_values(min_max_file, new_lo, new_hi, no_data_value);
  lo_value = std::min(new_lo, lo_value);
  hi_value = std::max(new_hi, hi_value);
  std::cout << "Pixel range for \"" << filename << "\": [" << new_lo << " " << new_hi << "]    Output dynamic range: [" << lo_value << " " << hi_value << "]" << std::endl;
}

template <class PixelT>
void do_normal_mosaic(po::variables_map const& vm, const ProgressCallback *progress) {
    DiskImageView<PixelT> img(image_files[0]);
    ImageQuadTreeGenerator<PixelT> quadtree(output_file_name, img);
    quadtree.set_patch_size( patch_size );
    quadtree.set_patch_overlap( patch_overlap );
    quadtree.set_output_image_file_type( output_file_type );

    quadtree.generate( *progress );
}

template <class ChannelT, class PixelT>
void do_mosaic(po::variables_map const& vm) {
  TerminalProgressCallback tpc;
  const ProgressCallback *progress = &tpc;
  if(verbose) {
    set_debug_level(VerboseDebugMessage);
    progress = &ProgressCallback::dummy_instance();
  } else if(quiet) {
    set_debug_level(WarningMessage);
  }

  // If we're not outputting any special sort of mosaic (just a regular old
  // quadtree, no georeferencing, no metadata), we use a different
  // function.
  if(output_metadata == "none") {
    if(image_files.size() != 1) {
      std::cerr << "Error: can only have 1 image as input when not creating a geo-referenced quadtree." << std::endl;
      std::cerr << "       (Use the `blend' program to create a quadtree with multiple images.)" << std::endl;
      return;
    }
    do_normal_mosaic<PixelT>(vm, progress);
    return;
  }

  DiskImageResourceJPEG::set_default_quality( jpeg_quality );
  Cache::system_cache().resize( cache_size*1024*1024 );

  // Read in georeference info and compute total resolution.
  GeoReference output_georef;
  output_georef.set_well_known_geogcs("WGS84");
  int total_resolution = 1024;
  std::vector<GeoReference> georeferences;

  for(int i=0; i < image_files.size(); i++) {
    std::cout << "Adding file " << image_files[i] << std::endl;
    DiskImageResourceGDAL file_resource( image_files[i] );

    if( vm.count("normalize") ) get_normalize_vals(image_files[i], file_resource);

    GeoReference input_georef;
    read_georeference( input_georef, file_resource );

    if(vm.count("force-wgs84"))
      input_georef.set_well_known_geogcs("WGS84");
    if(input_georef.proj4_str() == "")
      input_georef.set_well_known_geogcs("WGS84");
    if(input_georef.transform() == identity_matrix<3>() ) {
      vw_out(ErrorMessage) << "ERror: No georeferencing info found for input file \"" << image_files[i] << "\"!" << std::endl;
      exit(1);
    }

    georeferences.push_back( input_georef );

    // Right now, we only need a WGS84 output geoereference to compute 
    // the resolution. The rest of the output info will get set later.
    GeoTransform geotx( input_georef, output_georef );
    // Calculate the best resolution at 5 different points in the image, 
    // as occasionally there's a singularity at the center pixel that 
    // makes it extremely tiny (such as in pole-centered images).
    const int cols = file_resource.cols();
    const int rows = file_resource.rows();
    Vector2 res_pixel[5];
    res_pixel[0] = Vector2( cols/2, rows/2 );
    res_pixel[1] = Vector2( cols/2 + cols/4, rows/2 );
    res_pixel[2] = Vector2( cols/2 - cols/4, rows/2 );
    res_pixel[3] = Vector2( cols/2, rows/2 + rows/4 );
    res_pixel[4] = Vector2 (cols/2, rows/2 - rows/4 );
    int resolution;
    for(int i=0; i < 5; i++) {
      resolution = str_to_resolution_fn_map[output_metadata](geotx, res_pixel[i]);
      if( resolution > total_resolution ) total_resolution = resolution;
    }
  }
  // Now that we have the best resolution, we can get our output_georef.
  int xresolution = total_resolution / aspect_ratio, yresolution = total_resolution;
  if(output_metadata == "kml") {
    output_georef = output::kml::get_output_georeference(xresolution,yresolution);
  } else if(output_metadata == "tms" || output_metadata == "uniview") {
    output_georef = output::tms::get_output_georeference(total_resolution);
  }

  // Configure the composite.
  ImageComposite<PixelT> composite;

  // Add the transformed image files to the composite.
  for(int i=0; i < image_files.size(); i++) {
    GeoTransform geotx( georeferences[i], output_georef );
    ImageViewRef<PixelT> source = DiskImageView<PixelT>( image_files[i] );

    // Do various modifications to the input image here.
    if( pixel_scale != 1.0 || pixel_offset != 0.0 )
      source = channel_cast_rescale<ChannelT>( DiskImageView<PixelT>( image_files[i] ) * pixel_scale + pixel_offset );

    if( vm.count("normalize") )
      source = pixel_cast<PixelT>(channel_cast_rescale<ChannelT>( normalize_retain_alpha(DiskImageView<PixelRGBA<float> >( image_files[i] ), lo_value, hi_value, 0.0, 1.0) ) );

    if( vm.count("palette-file") ) {
      DiskImageView<float> disk_image( image_files[i] );
      if( vm.count("palette-scale") || vm.count("palette-offset") ) {
        source = per_pixel_filter( disk_image*palette_scale+palette_offset, PaletteFilter<PixelT>(palette_file) );
      } else {
        source = per_pixel_filter( disk_image, PaletteFilter<PixelT>(palette_file) );
      }
    }

    BBox2i bbox = geotx.forward_bbox( BBox2i(0,0,source.cols(),source.rows()) );
    source = crop( transform( source, geotx ), bbox );
    composite.insert( source, bbox.min().x(), bbox.min().y() );
    // Images that wrap the date line must be added to the composite on both sides.
    if( bbox.max().x() > total_resolution ) {
      composite.insert( source, bbox.min().x()-total_resolution, bbox.min().y() );
    }
  }

  BBox2i bbox = composite.bbox();
  BBox2i data_bbox;
  BBox2i total_bbox;
  BBox2 ll_bbox;
  // Now we differ a bit based on our output.
  if(output_metadata == "kml") {
    // Compute a tighter Google Earth coordinate system aligned bounding box.
    bbox.crop( BBox2i(0,0,xresolution,yresolution) );
    int dim = 2 << (int)(log( (std::max)(bbox.width(),bbox.height()) )/log(2));
    if( dim > total_resolution ) dim = total_resolution;
    total_bbox = BBox2i( (bbox.min().x()/dim)*dim, (bbox.min().y()/dim)*dim, dim, dim );
    if( ! total_bbox.contains( bbox ) ) {
      if( total_bbox.max().x() == xresolution ) total_bbox.min().x() -= dim;
      else total_bbox.max().x() += dim;
      if( total_bbox.max().y() == yresolution ) total_bbox.min().y() -= dim;
      else total_bbox.max().y() += dim;
    }

    ll_bbox = BBox2( -180.0 + (360.0*total_bbox.min().x())/xresolution,
                      180.0 - (360.0*total_bbox.max().y())/yresolution,
                      (360.0*total_bbox.width())/xresolution,
                      (360.0*total_bbox.height())/yresolution );
  } else if(output_metadata == "tms" || output_metadata == "uniview") {
    total_bbox = composite.bbox();
    total_bbox.grow( BBox2i(0,0,total_resolution,total_resolution) );
    total_bbox.crop( BBox2i(0,0,total_resolution,total_resolution) );

    Vector2 invmin = output_georef.pixel_to_lonlat(total_bbox.min());
    Vector2 invmax = output_georef.pixel_to_lonlat(total_bbox.max());
    ll_bbox.min().x() = invmin[0];
    ll_bbox.max().y() = invmin[1];
    ll_bbox.max().x() = invmax[0];
    ll_bbox.min().y() = invmax[1];
  }
  std::cerr << "DEBUG: ll_bbox = " << ll_bbox << std::endl;

  // Prepare the composite.
  if( vm.count("composite-multiband") ) {
    std::cout << "Preparing composite..." << std::endl;
    composite.prepare( total_bbox, *progress );
  } else {
    composite.set_draft_mode( true );
    composite.prepare( total_bbox );
  }

  // Data bbox.
  if(output_metadata == "kml") {
    data_bbox = composite.bbox();
    data_bbox.crop( BBox2i(0,0,total_bbox.width(),total_bbox.height()) );
  } else if(output_metadata == "tms" || output_metadata == "uniview") {
    data_bbox = BBox2i((int)std::floor(double(bbox.min().x())/patch_size)*patch_size,
                       (int)std::floor(double(bbox.min().y())/patch_size)*patch_size,
                       (int)std::ceil(double(bbox.width())/patch_size)*patch_size,
                       (int)std::ceil(double(bbox.height())/patch_size)*patch_size);
    data_bbox.crop(total_bbox);
  }

  ImageQuadTreeGenerator<PixelT> *quadtree;

 // KML specific things. 
  if( output_metadata == "kml" ) {
    KMLQuadTreeGenerator<PixelT> *tmp;
    tmp = new KMLQuadTreeGenerator<PixelT>( output_file_name, composite, ll_bbox );
    quadtree = tmp;

    tmp->set_max_lod_pixels( max_lod_pixels );
    tmp->set_draw_order_offset( draw_order_offset );

  // TMS specific things.
  } else if( output_metadata == "tms" ) {
    TMSQuadTreeGenerator<PixelT> *tmp;
    tmp = new TMSQuadTreeGenerator<PixelT>( output_file_name, composite );
    quadtree = tmp;

  // Uniview specific things.
  } else if( output_metadata == "uniview" ) {
    UniviewQuadTreeGenerator<PixelT> *tmp;
    tmp = new UniviewQuadTreeGenerator<PixelT>( output_file_name, composite, terrain );
    quadtree = tmp;

    // Have to delay the uniview-specific things output config file until
    // the generation is over, due to the way variables are set in the 
    // QuadTreeGenerator.
  } else {
    // Unreachable
    vw_throw(LogicErr() << "Unreachable statement reached: bad value for output_metadata (value was " << output_metadata << ")");
  }
    
  quadtree->set_crop_bbox(data_bbox);
  if( vm.count("crop") ) quadtree->set_crop_images(true);
  quadtree->set_patch_size(patch_size);
  quadtree->set_patch_overlap(patch_overlap);
  quadtree->set_output_image_file_type(output_file_type);


  // Generate the composite.
  vw_out(InfoMessage) << "Generating " << output_metadata << " overlay..." << std::endl;
  quadtree->generate(*progress);

  // Now output the Uniview config file. Had to delay it to here because 
  // the quadtree generator does not calculate the tree levels until 
  // the generation function is called.
  if(output_metadata == "uniview") {
    std::string config_filename = output_file_name + ".conf";
    std::ofstream conf( config_filename.c_str() );
    if(terrain) {
      conf << "// Terrain\n";
      conf << "HeightmapCacheLocation=modules/" << module_name << "/Offlinedatasets/" << output_file_name << "/Terrain/\n";
      conf << "HeightmapCallstring=Generated by the NASA Vision Workbench image2qtree tool.\n";
      conf << "HeightmapFormat=" << quadtree->get_output_image_file_type() << '\n';
      conf << "NrHeightmapLevels=" << quadtree->get_tree_levels()-1 << '\n';
      conf << "NrLevelsPerHeightmap=1\n";
    } else {
      conf << "[Offlinedataset]\n";
      conf << "NrRows=1\n";
      conf << "NrColumns=2\n";
      conf << "Bbox= -180 -90 180 90\n";
      conf << "DatasetTitle=" << output_file_name << "\n";
      conf << "Tessellation=19\n\n";

      conf << "// Texture\n";
      conf << "TextureCacheLocation=modules/" << module_name << "/Offlinedatasets/" << output_file_name << "/Texture/\n";
      conf << "TextureCallstring=Generated by the NASA Vision Workbench image2qtree tool.\n";
      conf << "TextureFormat=" << quadtree->get_output_image_file_type() << "\n";
      conf << "TextureLevels= " << quadtree->get_tree_levels()-1 << "\n";
      conf << "TextureSize= " << patch_size << "\n\n";
    }
    conf.close();
    std::cout << "Note: You must merge the texture and terrain config files into a single file (Terrain info should go below texture info.)" << std::endl;
    std::cout << "Both output sets should be in the same directory, with the texture in a subdirectory named Texture and the terrain in a subdirectory named Terrain." << std::endl;
  }
}

int main(int argc, char **argv) {
  fill_input_maps();

  po::options_description general_options("General Options");
  general_options.add_options()
    ("output-name,o", po::value<std::string>(&output_file_name)->default_value("output"), "Specify the base output directory")
    ("quiet,q", "Quiet output")
    ("verbose,v", "Verbose output")
    ("cache", po::value<unsigned>(&cache_size)->default_value(1024), "Cache size, in megabytes")
    ("help", "Display this help message");

  po::options_description input_options("Input Options");
  input_options.add_options()
    ("force-wgs84", "Use WGS84 as the input images' geographic coordinate systems, even if they're not (old behavior)")
    ("pixel-scale", po::value<float>(&pixel_scale)->default_value(1.0), "Scale factor to apply to pixels")
    ("pixel-offset", po::value<float>(&pixel_offset)->default_value(0.0), "Offset to apply to pixels")
    ("normalize", "Normalize input images so that their full dynamic range falls in between [0,255].");

  po::options_description output_options("Output Options");
  output_options.add_options()
    ("output-metadata,m", po::value<std::string>(&output_metadata)->default_value("none"), "Specify the output metadata type. One of [kml, tms, uniview, none]")
    ("file-type", po::value<std::string>(&output_file_type)->default_value("png"), "Output file type")
    ("channel-type", po::value<std::string>(&channel_type_str)->default_value("uint8"), "Output (and input) channel type. One of [uint8, uint16, int16, float]")
    ("module-name", po::value<std::string>(&module_name)->default_value("marsds"), "Uniview module name (Uniview only). The module where the output will be placed")
    ("terrain", "Outputs image files suitable for a Uniview terrain view. Implies output format as PNG, channel type uint16. Uniview only")
    ("jpeg-quality", po::value<float>(&jpeg_quality)->default_value(0.75), "JPEG quality factor (0.0 to 1.0)")
    ("palette-file", po::value<std::string>(&palette_file), "Apply a palette from the given file")
    ("palette-scale", po::value<float>(&palette_scale), "Apply a scale factor before applying the palette")
    ("palette-offset", po::value<float>(&palette_offset), "Apply an offset before applying the palette")
    ("patch-size", po::value<int>(&patch_size)->default_value(256), "Patch size, in pixels")
    ("patch-overlap", po::value<int>(&patch_overlap)->default_value(0), "Patch overlap, in pixels")
    ("max-lod-pixels", po::value<int>(&max_lod_pixels)->default_value(1024), "Max LoD in pixels, or -1 for none (kml only)")
    ("draw-order-offset", po::value<int>(&draw_order_offset)->default_value(0), "Offset for the <drawOrder> tag for this overlay (kml only)")
    ("composite-multiband", "Composite images using multi-band blending")
    ("aspect-ratio", po::value<int>(&aspect_ratio)->default_value(1), "Pixel aspect ratio (for polar overlays; should be a power of two)");

  po::options_description hidden_options("");
  hidden_options.add_options()
    ("input-file", po::value<std::vector<std::string> >(&image_files));

  po::options_description options("Allowed Options");
  options.add(general_options).add(input_options).add(output_options).add(hidden_options);

  po::positional_options_description p;
  p.add("input-file", -1);

  po::variables_map vm;
  po::store( po::command_line_parser( argc, argv ).options(options).positional(p).run(), vm );
  po::notify( vm );

  std::ostringstream command_line;
  for(int i=0; i < argc; i++) {
    command_line << argv[i];
    if(i < argc-1) command_line << ' ';
  }

  std::ostringstream usage;
  usage << "Usage: image2qtree [options] <filename>..." <<std::endl << std::endl;
  usage << general_options << std::endl;
  usage << input_options << std::endl;
  usage << output_options << std::endl;

  if( vm.count("help") ) {
    std::cout << usage.str();
    return 0;
  }

  if( vm.count("input-file") < 1 ) {
    std::cerr << "Error: must specify at least one input file!" << std::endl << std::endl;
    std::cout << usage.str();
    return 1;
  }

  if( patch_size <= 0 ) {
    std::cerr << "Error: The patch size must be a positive number!  (You specified: " << patch_size << ".)" << std::endl << std::endl;
    std::cout << usage.str();
    return 1;
  }

  if(str_to_resolution_fn_map.find(output_metadata) == str_to_resolution_fn_map.end()) {
    std::cerr << "Error: Output metadata must be one of [none, kml, tms, uniview]!  (You specified: " << output_metadata << ".)" << std::endl << std::endl;
    std::cout << usage.str();
    return 1;
  }

  // Set a few booleans based on input values.
  if(vm.count("verbose")) verbose = true;
  if(vm.count("quiet")) quiet = true;
  if(verbose && quiet) {
    std::cerr << "Error: Cannot be verbose and quiet at the same time." << std::endl << std::endl;
    std::cout << usage.str();
    return 1;
  }

  if(vm.count("terrain")) {
    terrain = true;
    output_file_type = std::string("png");
  }

  // Get the right pixel/channel type, and call the mosaic.
  DiskImageResource *first_resource = DiskImageResource::open(image_files[0]);
  ChannelTypeEnum channel_type = first_resource->channel_type();
  PixelFormatEnum pixel_format = first_resource->pixel_format();
  delete first_resource;
  if(vm.count("channel-type")) {
    if(channel_type_str == "uint8") channel_type = VW_CHANNEL_UINT8;
    else if(channel_type_str == "int16") channel_type = VW_CHANNEL_INT16;
    else if(channel_type_str == "uint16") channel_type = VW_CHANNEL_UINT16;
    else if(channel_type_str == "float") channel_type = VW_CHANNEL_FLOAT32;
    else {
      std::cerr << "Error: Channel type must be one of [uint8, uint16, int16, float]!  (You specified: " << channel_type_str << ".)" << std::endl << std::endl;
      std::cout << usage.str();
      return 1;
    }
  }

  // Convert non-alpha channel images into images with one for the
  // composite.
  switch(pixel_format) {
    case VW_PIXEL_GRAY:
    case VW_PIXEL_GRAYA:
      switch(channel_type) {
        case VW_CHANNEL_UINT8:  do_mosaic<uint8, PixelGrayA<uint8> >(vm); break;
        case VW_CHANNEL_INT16:  do_mosaic<int16, PixelGrayA<int16> >(vm); break;
        case VW_CHANNEL_UINT16: do_mosaic<uint16, PixelGrayA<uint16> >(vm); break;
        default:                do_mosaic<float32, PixelGrayA<float32> >(vm); break;
      }
      break;
    case VW_PIXEL_RGB:
    case VW_PIXEL_RGBA:
    default:
      switch(channel_type) {
        case VW_CHANNEL_UINT8:  do_mosaic<uint8, PixelRGBA<uint8> >(vm); break;
        case VW_CHANNEL_INT16:  do_mosaic<int16, PixelRGBA<int16> >(vm); break;
        case VW_CHANNEL_UINT16: do_mosaic<uint16, PixelRGBA<uint16> >(vm); break;
        default:                do_mosaic<float32, PixelRGBA<float32> >(vm); break;
      }
      break;
  }

  return 0;
}
