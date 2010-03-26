// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file image2geotree.cc
///
/// This program takes a georeferenced image as its input, and outputs 
/// a quadtree for that image that is viewable in various terrain display
/// programs, such as Google Earth. Currently, the program supports output
/// in KML, TMS, Uniview, and Google Maps formats.

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
#include <boost/algorithm/string.hpp>
#include <boost/filesystem/path.hpp>
namespace po = boost::program_options;

#include <vw/Core/Cache.h>
#include <vw/Core/ProgressCallback.h>
#include <vw/Math/Matrix.h>
#include <vw/Image/Transform.h>
#include <vw/Image/Palette.h>
#include <vw/Image/MaskViews.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/FileIO/DiskImageResourceJPEG.h>
#include <vw/FileIO/DiskImageResourcePNG.h>
#include <vw/FileIO/DiskImageResourceGDAL.h>
#include <vw/FileIO/DiskImageView.h>
#include <vw/Cartography/GeoReference.h>
#include <vw/Cartography/GeoTransform.h>
#include <vw/Cartography/FileIO.h>
#include <vw/Mosaic.h>
using namespace vw;
using namespace vw::math;
using namespace vw::cartography;
using namespace vw::mosaic;

// Global variables
std::vector<std::string> image_files;
std::string output_file_name = "";
std::string output_file_type;
std::string output_metadata;
std::string module_name;
double nudge_x=0, nudge_y=0;
double north_lat=90.0, south_lat=-90.0;
double east_lon=180.0, west_lon=-180.0;
double proj_lat=0, proj_lon=0, proj_scale=1;
int utm_zone;
int tile_size;
float jpeg_quality;
int png_compression;
std::string palette_file;
std::string channel_type_str;
float palette_scale=1.0, palette_offset=0.0;
int draw_order_offset; // KML only.
int max_lod_pixels; // KML only.
float pixel_scale=1.0, pixel_offset=0.0;
double lcc_parallel1, lcc_parallel2;
int aspect_ratio=1;
int global_resolution=0;
bool terrain=false;
float nodata=0;

// For image stretching.
float lo_value = ScalarTypeLimits<float>::highest();
float hi_value = ScalarTypeLimits<float>::lowest();

// Function pointers for computing resolution.
std::map<std::string, vw::int32 (*)(const GeoTransform&, const Vector2&)> str_to_resolution_fn_map;

// Fill the maps for converting input strings to function pointers.
static void fill_input_maps() {
  str_to_resolution_fn_map[std::string("none")]     = NULL;
  str_to_resolution_fn_map[std::string("kml")]      = &vw::cartography::output::kml::compute_resolution;
  str_to_resolution_fn_map[std::string("tms")]      = &vw::cartography::output::tms::compute_resolution;
  str_to_resolution_fn_map[std::string("uniview")]  = &vw::cartography::output::tms::compute_resolution;
  str_to_resolution_fn_map[std::string("gmap")]     = &vw::cartography::output::tms::compute_resolution;
  str_to_resolution_fn_map[std::string("celestia")] = &vw::cartography::output::tms::compute_resolution;
  str_to_resolution_fn_map[std::string("gigapan")]  = &vw::cartography::output::tms::compute_resolution;
  str_to_resolution_fn_map[std::string("gigapan-noproj")] = NULL;
}

static void get_normalize_vals(std::string filename, DiskImageResourceGDAL &file_resource) {
  
  DiskImageView<PixelRGBA<float> > min_max_file(filename);
  float new_lo, new_hi;
  if ( file_resource.has_nodata_value() ) {
    PixelRGBA<float> no_data_value( file_resource.nodata_value() );
    min_max_channel_values( create_mask(min_max_file,no_data_value), new_lo, new_hi );
  } else {
    min_max_channel_values( create_mask(min_max_file), new_lo, new_hi );
  }
  lo_value = std::min(new_lo, lo_value);
  hi_value = std::max(new_hi, hi_value);
  std::cout << "Pixel range for \"" << filename << "\": [" << new_lo << " " << new_hi << "]    Output dynamic range: [" << lo_value << " " << hi_value << "]" << std::endl;
}

template <class PixelT>
void do_normal_mosaic(po::variables_map const& /*vm*/, const ProgressCallback *progress) {
    DiskImageView<PixelT> img(image_files[0]);
    QuadTreeGenerator quadtree(img, output_file_name);
    quadtree.set_tile_size( tile_size );
    quadtree.set_file_type( output_file_type );

    if (output_metadata == "gigapan-noproj") {
      GigapanQuadTreeConfig config;
      config.configure( quadtree );
    }

    quadtree.generate( *progress );
}

template <class PixelT>
void do_mosaic(po::variables_map const& vm, const ProgressCallback *progress)
{
  typedef typename PixelChannelType<PixelT>::type ChannelT;

  // If we're not outputting any special sort of mosaic (just a regular old
  // quadtree, no georeferencing, no metadata), we use a different
  // function.
  if(output_metadata == "none" || output_metadata == "gigapan-noproj") {
    if(image_files.size() != 1) {
      std::cerr << "Error: can only have 1 image as input when not creating a geo-referenced quadtree." << std::endl;
      std::cerr << "       (Use the `blend' program to create a quadtree with multiple images.)" << std::endl;
      return;
    }
    do_normal_mosaic<PixelT>(vm, progress);
    return;
  }

  DiskImageResourceJPEG::set_default_quality( jpeg_quality );
  DiskImageResourcePNG::set_default_compression_level( png_compression );

  // Read in georeference info and compute total resolution.
  GeoReference output_georef;
  int total_resolution = 1024;
  std::vector<GeoReference> georeferences;

  for(unsigned i=0; i < image_files.size(); i++) {
    std::cout << "Adding file " << image_files[i] << std::endl;
    DiskImageResourceGDAL file_resource( image_files[i] );

    if( vm.count("normalize") ) get_normalize_vals(image_files[i], file_resource);

    GeoReference input_georef;
    read_georeference( input_georef, file_resource );

    if(vm.count("force-lunar-datum")) {
      const double LUNAR_RADIUS = 1737400;
      vw_out() << "\t--> Using standard lunar spherical datum: "
                << LUNAR_RADIUS << "\n";
      cartography::Datum datum("D_MOON",
                               "MOON",
                               "Reference Meridian",
                               LUNAR_RADIUS,
                               LUNAR_RADIUS,
                               0.0);
      input_georef.set_datum(datum);
    } else if(vm.count("force-mars-datum")) {
      const double MOLA_PEDR_EQUATORIAL_RADIUS = 3396000.0;
      vw_out() << "\t--> Using standard MOLA spherical datum: "
                << MOLA_PEDR_EQUATORIAL_RADIUS << "\n";
      cartography::Datum datum("D_MARS",
                               "MARS",
                               "Reference Meridian",
                               MOLA_PEDR_EQUATORIAL_RADIUS,
                               MOLA_PEDR_EQUATORIAL_RADIUS,
                               0.0);
      input_georef.set_datum(datum);
    }

    if(vm.count("force-wgs84")) {
      input_georef.set_well_known_geogcs("WGS84");
    }
    if(input_georef.proj4_str() == "") {
      input_georef.set_well_known_geogcs("WGS84");
    }
    if(i==0) {
      output_georef.set_datum( input_georef.datum() );
    }

    bool manual = vm.count("north") || vm.count("south") || vm.count("east") || vm.count("west");
    if( manual || input_georef.transform() == identity_matrix<3>() ) {
      if( image_files.size() == 1 ) {
        vw_out() << "No georeferencing info found.  Assuming Plate Carree WGS84: "
                 << east_lon << " to " << west_lon << " E, " << south_lat
                 << " to " << north_lat << " N." << std::endl;
        input_georef = GeoReference();
        input_georef.set_well_known_geogcs("WGS84");
        Matrix3x3 m;
        m(0,0) = (east_lon - west_lon) / file_resource.cols();
        m(0,2) = west_lon;
        m(1,1) = (south_lat - north_lat) / file_resource.rows();
        m(1,2) = north_lat;
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
    else if( vm.count("sinusoidal") ) input_georef.set_sinusoidal(proj_lon);
    else if( vm.count("mercator") ) input_georef.set_mercator(proj_lat,proj_lon,proj_scale);
    else if( vm.count("transverse-mercator") ) input_georef.set_transverse_mercator(proj_lat,proj_lon,proj_scale);
    else if( vm.count("orthographic") ) input_georef.set_orthographic(proj_lat,proj_lon);
    else if( vm.count("stereographic") ) input_georef.set_stereographic(proj_lat,proj_lon,proj_scale);
    else if( vm.count("lambert-azimuthal") ) input_georef.set_lambert_azimuthal(proj_lat,proj_lon);
    else if( vm.count("lambert-conformal-conic") ) input_georef.set_lambert_azimuthal(lcc_parallel1, lcc_parallel2, proj_lat, proj_lon);
    else if( vm.count("utm") ) {
      input_georef.set_UTM( abs(utm_zone), utm_zone > 0 );
    }

    if( vm.count("nudge-x") || vm.count("nudge-y") ) {
      Matrix3x3 m = input_georef.transform();
      m(0,2) += nudge_x;
      m(1,2) += nudge_y;
      input_georef.set_transform( m );
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

  if( global_resolution) total_resolution = global_resolution;

  // Now that we have the best resolution, we can get our output_georef.
  int xresolution = total_resolution / aspect_ratio, yresolution = total_resolution;

  if(output_metadata == "kml") {
    output_georef = output::kml::get_output_georeference(xresolution,yresolution);
    output_georef.set_datum( georeferences[0].datum() );
  } else if(output_metadata != "none") {
    output_georef = output::tms::get_output_georeference(total_resolution);
    output_georef.set_datum( georeferences[0].datum() );
  }

  // Configure the composite.
  ImageComposite<PixelT> composite;

  // Add the transformed image files to the composite.
  for(unsigned i=0; i < image_files.size(); i++) {
    GeoTransform geotx( georeferences[i], output_georef );
    ImageViewRef<PixelT> source = DiskImageView<PixelT>( image_files[i] );

    if( vm.count("nodata") )
      source = mask_to_alpha(create_mask(pixel_cast<typename PixelWithoutAlpha<PixelT>::type >(source),ChannelT(nodata)));

    bool global = boost::trim_copy(georeferences[i].proj4_str())=="+proj=longlat" &&
      fabs(georeferences[i].lonlat_to_pixel(Vector2(-180,0)).x()) < 1 &&
      fabs(georeferences[i].lonlat_to_pixel(Vector2(180,0)).x() - source.cols()) < 1 &&
      fabs(georeferences[i].lonlat_to_pixel(Vector2(0,90)).y()) < 1 &&
      fabs(georeferences[i].lonlat_to_pixel(Vector2(0,-90)).y() - source.rows()) < 1;

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
    if (global) {
      vw_out() << "\t--> Detected global overlay.  Using cylindrical edge extension to hide the seam.\n";
      source = crop( transform( source, geotx, source.cols(), source.rows(), CylindricalEdgeExtension() ), bbox );
    }
    else
      source = crop( transform( source, geotx ), bbox );
    // Images that wrap the date line must be added to the composite on both sides.
    if( bbox.max().x() > total_resolution ) {
      composite.insert( source, bbox.min().x()-total_resolution, bbox.min().y() );
    }
    // Images that are in the 180-360 range *only* go on the other side.
    if( bbox.min().x() < xresolution ) {
      composite.insert( source, bbox.min().x(), bbox.min().y() );
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
    int dim = 2 << (int)(log( (double)(std::max)(bbox.width(),bbox.height()) )/log(2.));
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
  } else if (output_metadata == "gigapan") {
    total_bbox = bbox; 
    bbox.crop( BBox2i(0,0,xresolution,yresolution) );
    ll_bbox = BBox2( -180.0 + (360.0*total_bbox.min().x())/xresolution,
                      180.0 - (360.0*total_bbox.max().y())/yresolution,
                      (360.0*total_bbox.width())/xresolution,
                      (360.0*total_bbox.height())/yresolution );
  }  else if(output_metadata != "none") {
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

  // Prepare the composite.
  if( vm.count("composite-multiband") ) {
    std::cout << "Preparing composite..." << std::endl;
    composite.prepare( total_bbox, *progress );
  } else {
    composite.set_draft_mode( true );
    composite.prepare( total_bbox );
  }

  // Data bbox.
  if(output_metadata == "kml" || output_metadata == "gigapan") {
    data_bbox = composite.bbox();
    data_bbox.crop( BBox2i(0,0,total_bbox.width(),total_bbox.height()) );
  } else if(output_metadata != "none") {
    data_bbox = BBox2i((int)std::floor(double(bbox.min().x())/tile_size)*tile_size,
                       (int)std::floor(double(bbox.min().y())/tile_size)*tile_size,
                       (int)std::ceil(double(bbox.width())/tile_size)*tile_size,
                       (int)std::ceil(double(bbox.height())/tile_size)*tile_size);
    data_bbox.crop(total_bbox);
  }

  QuadTreeGenerator quadtree( composite, output_file_name );

  // KML specific things. 
  if( output_metadata == "kml" ) {
    KMLQuadTreeConfig config;
    config.set_longlat_bbox( ll_bbox );
    config.set_max_lod_pixels( max_lod_pixels );
    config.set_draw_order_offset( draw_order_offset );
    config.configure( quadtree );

  // TMS specific things.
  } else if( output_metadata == "tms" ) {
    TMSQuadTreeConfig config;
    config.configure( quadtree );

  // Uniview specific things.
  } else if( output_metadata == "uniview" ) {
    UniviewQuadTreeConfig config( terrain );
    config.configure( quadtree );

  // Google Maps specific things.
  } else if( output_metadata == "gmap" ) {
    GMapQuadTreeConfig config;
    config.configure( quadtree );

  // Celestia specific things.
  } else if ( output_metadata == "celestia" ) {
    CelestiaQuadTreeConfig config;
    config.configure( quadtree );

  // Gigapan specific things
  } else if ( output_metadata == "gigapan" ) {
    GigapanQuadTreeConfig config;
    config.set_longlat_bbox( ll_bbox );
    config.configure( quadtree );
  
  // Unreachable
  } else {
    vw_throw(LogicErr() << "Unreachable statement reached: bad value for output_metadata (value was " << output_metadata << ")");
  }

  quadtree.set_crop_bbox(data_bbox);
  if( vm.count("crop") ) quadtree.set_crop_images(true);
  quadtree.set_tile_size(tile_size);
  quadtree.set_file_type(output_file_type);

  // Generate the composite.
  vw_out() << "Generating " << output_metadata << " overlay..." << std::endl;
  quadtree.generate(*progress);

  // This should really get moved into a metadata function for 
  // UniviewQuadTreeConfig.
  if(output_metadata == "uniview") {
    std::string config_filename = output_file_name + ".conf";
    std::ofstream conf( config_filename.c_str() );
    if(terrain) {
      conf << "// Terrain\n";
      conf << "HeightmapCacheLocation=modules/" << module_name << "/Offlinedatasets/" << output_file_name << "/Terrain/\n";
      conf << "HeightmapCallstring=Generated by the NASA Vision Workbench image2qtree tool.\n";
      conf << "HeightmapFormat=" << quadtree.get_file_type() << '\n';
      conf << "NrHeightmapLevels=" << quadtree.get_tree_levels()-1 << '\n';
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
      conf << "TextureFormat=" << quadtree.get_file_type() << "\n";
      conf << "TextureLevels= " << quadtree.get_tree_levels()-1 << "\n";
      conf << "TextureSize= " << tile_size << "\n\n";
    }
    conf.close();
    std::cout << "Note: You must merge the texture and terrain config files into a single file (Terrain info should go below texture info.)" << std::endl;
    std::cout << "Both output sets should be in the same directory, with the texture in a subdirectory named Texture and the terrain in a subdirectory named Terrain." << std::endl;
  } else if (output_metadata == "celestia") {
    std::string fn = output_file_name + ".ctx";
    std::ofstream ctx( fn.c_str() );
    ctx << "VirtualTexture\n";
    ctx << "{\n";
    ctx << "        ImageDirectory \"" << output_file_name << "\"\n";
    ctx << "        BaseSplit 0\n";
    ctx << "        TileSize " << (tile_size >> 1) << "\n";
    ctx << "        TileType \"" << output_file_type << "\"\n";
    ctx << "}\n";
    ctx.close();

    fn = output_file_name + ".ssc";
    std::ofstream ssc( fn.c_str() );

    ssc << "AltSurface \"" << output_file_name << "\" \"" << module_name << "\"\n";
    ssc << "{\n";
    ssc << "    Texture \"" << output_file_name << ".ctx" << "\"\n";
    ssc << "}\n";
    ssc.close();
    std::cout << "Place " << output_file_name << ".ssc" << " in Celestia's extras dir" << std::endl;
    std::cout << "Place " << output_file_name << ".ctx" << " and the output dir ("
                          << output_file_name << ") in extras/textures/hires" << std::endl;
  } 
}

int main(int argc, char **argv) {
  fill_input_maps();

  po::options_description general_options("Description: Turns georeferenced image(s) into a quadtree with geographical metadata\n\nGeneral Options");
  general_options.add_options()
    ("output-name,o", po::value<std::string>(&output_file_name), "Specify the base output directory")
    ("help,h", "Display this help message");

  po::options_description input_options("Input Options");
  input_options.add_options()
    ("force-wgs84", "Use WGS84 as the input images' geographic coordinate systems, even if they're not (old behavior)")
    ("force-lunar-datum", "Use the lunar spherical datum for the input images' geographic coordinate systems, even if they are not encoded to do so.")
    ("force-mars-datum", "Use the Mars spherical datum for the input images' geographic coordinate systems, even if they are not encoded to do so.")
    ("pixel-scale", po::value<float>(&pixel_scale)->default_value(1.0), "Scale factor to apply to pixels")
    ("pixel-offset", po::value<float>(&pixel_offset)->default_value(0.0), "Offset to apply to pixels")
    ("normalize", "Normalize input images so that their full dynamic range falls in between [0,255].")
    ("nodata",po::value<float>(&nodata),"Set the input's nodata value so that it will be transparent in output");

  po::options_description output_options("Output Options");
  output_options.add_options()
    ("output-metadata,m", po::value<std::string>(&output_metadata)->default_value("kml"), "Specify the output metadata type. One of [kml, tms, uniview, gmap, celestia, none]")
    ("file-type", po::value<std::string>(&output_file_type)->default_value("png"), "Output file type.  (Choose \'auto\' to generate jpgs in opaque areas and png images where there is transparency.)")
    ("channel-type", po::value<std::string>(&channel_type_str)->default_value("uint8"), "Output (and input) channel type. One of [uint8, uint16, int16, float]")
    ("module-name", po::value<std::string>(&module_name)->default_value("marsds"), "The module where the output will be placed. Ex: marsds for Uniview, or Sol/Mars for Celestia")
    ("terrain", "Outputs image files suitable for a Uniview terrain view. Implies output format as PNG, channel type uint16. Uniview only")
    ("jpeg-quality", po::value<float>(&jpeg_quality)->default_value(0.75), "JPEG quality factor (0.0 to 1.0)")
    ("png-compression", po::value<int>(&png_compression)->default_value(3), "PNG compression level (0 to 9)")
    ("palette-file", po::value<std::string>(&palette_file), "Apply a palette from the given file")
    ("palette-scale", po::value<float>(&palette_scale), "Apply a scale factor before applying the palette")
    ("palette-offset", po::value<float>(&palette_offset), "Apply an offset before applying the palette")
    ("tile-size", po::value<int>(&tile_size)->default_value(256), "Tile size, in pixels")
    ("max-lod-pixels", po::value<int>(&max_lod_pixels)->default_value(1024), "Max LoD in pixels, or -1 for none (kml only)")
    ("draw-order-offset", po::value<int>(&draw_order_offset)->default_value(0), "Offset for the <drawOrder> tag for this overlay (kml only)")
    ("composite-multiband", "Composite images using multi-band blending")
    ("aspect-ratio", po::value<int>(&aspect_ratio)->default_value(1), "Pixel aspect ratio (for polar overlays; should be a power of two)")
    ("global-resolution", po::value<int>(&global_resolution)->default_value(0), "Override the global pixel resolution; should be a power of two");

  po::options_description projection_options("Projection Options");
  projection_options.add_options()
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
    ("lambert-conformal-conic", "Assume a Lambert Conformal Conic projection")
    ("utm", po::value(&utm_zone), "Assume UTM projection with the given zone (+ for North, - for South)")
    ("proj-lat", po::value<double>(&proj_lat), "The center of projection latitude (if applicable)")
    ("proj-lon", po::value<double>(&proj_lon), "The center of projection longitude (if applicable)")
    ("proj-scale", po::value<double>(&proj_scale), "The projection scale (if applicable)")
    ("std-parallel1", po::value<double>(&lcc_parallel1), "Standard parallels for Lambert Conformal Conic projection")
    ("std-parallel2", po::value<double>(&lcc_parallel2), "Standard parallels for Lambert Conformal Conic projection")
    ("nudge-x", po::value<double>(&nudge_x), "Nudge the image, in projected coordinates")
    ("nudge-y", po::value<double>(&nudge_y), "Nudge the image, in projected coordinates");

  po::options_description hidden_options("");
  hidden_options.add_options()
    ("input-file", po::value<std::vector<std::string> >(&image_files));

  po::options_description options("Allowed Options");
  options.add(general_options).add(input_options).add(output_options).add(projection_options).add(hidden_options);

  po::positional_options_description p;
  p.add("input-file", -1);

  std::ostringstream usage;
  usage << "Usage: image2qtree [options] <filename>..." <<std::endl << std::endl;
  usage << general_options << std::endl;
  usage << input_options << std::endl;
  usage << output_options << std::endl;
  usage << projection_options << std::endl;

  po::variables_map vm;
  try {
    po::store( po::command_line_parser( argc, argv ).options(options).positional(p).run(), vm );
    po::notify( vm );
  } catch (po::error &e) {
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

  if( output_file_name == "" )
    output_file_name = fs::path(image_files[0]).replace_extension().string();

  if( tile_size <= 0 ) {
    std::cerr << "Error: The tile size must be a positive number!  (You specified: " << tile_size << ")." << std::endl << std::endl;
    std::cout << usage.str();
    return 1;
  }

  if(str_to_resolution_fn_map.find(output_metadata) == str_to_resolution_fn_map.end()) {
    std::cerr << "Error: Output metadata must be one of [none, kml, tms, uniview, gmap, celestia]!  (You specified: " << output_metadata << ")." << std::endl << std::endl;
    std::cout << usage.str();
    return 1;
  }

  TerminalProgressCallback tpc( "tools.image2qtree", "");
  const ProgressCallback *progress = &tpc;

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
    channel_type = channel_name_to_enum(channel_type_str);
    switch (channel_type) {
      case VW_CHANNEL_UINT8:  case VW_CHANNEL_INT16:
      case VW_CHANNEL_UINT16: case VW_CHANNEL_FLOAT32:
        break;
      default:
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
        case VW_CHANNEL_UINT8:  do_mosaic<PixelGrayA<uint8>   >(vm, progress); break;
        case VW_CHANNEL_INT16:  do_mosaic<PixelGrayA<int16>   >(vm, progress); break;
        case VW_CHANNEL_UINT16: do_mosaic<PixelGrayA<uint16>  >(vm, progress); break;
        default:                do_mosaic<PixelGrayA<float32> >(vm, progress); break;
      }
      break;
    case VW_PIXEL_RGB:
    case VW_PIXEL_RGBA:
    default:
      switch(channel_type) {
        case VW_CHANNEL_UINT8:  do_mosaic<PixelRGBA<uint8>   >(vm, progress); break;
        case VW_CHANNEL_INT16:  do_mosaic<PixelRGBA<int16>   >(vm, progress); break;
        case VW_CHANNEL_UINT16: do_mosaic<PixelRGBA<uint16>  >(vm, progress); break;
        default:                do_mosaic<PixelRGBA<float32> >(vm, progress); break;
      }
      break;
  }

  return 0;
}
