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


/// \file image2qtree.cc
///
/// This program takes a georeferenced image as its input, and outputs
/// a quadtree for that image that is viewable in various terrain display
/// programs, such as Google Earth. Currently, the program supports output
/// in KML, TMS, Uniview, and Google Maps formats.

#include <vw/tools/image2qtree.h>
#include <vw/FileIO/GdalWriteOptions.h>
#include <vw/FileIO/GdalWriteOptionsDesc.h>
#include <vw/Cartography/GeoReferenceUtils.h>
#include <boost/program_options.hpp>

using namespace vw;
using namespace vw::math;
using namespace vw::cartography;
using namespace vw::mosaic;
using std::string;
using vw::tools::Usage;
using std::cout;
using std::cerr;
using std::endl;

namespace po = boost::program_options;

// Call the correct version of the do_mosaic function based on pixel type.
#define SWITCH_ON_CHANNEL_TYPE( PIXELTYPE )                             \
  switch (fmt.channel_type) {                                           \
  case VW_CHANNEL_UINT8:   do_mosaic_##PIXELTYPE##_uint8(opt,progress); break; \
  case VW_CHANNEL_INT16:   do_mosaic_##PIXELTYPE##_int16(opt,progress); break; \
  case VW_CHANNEL_UINT16:  do_mosaic_##PIXELTYPE##_uint16(opt,progress); break; \
  default:                 do_mosaic_##PIXELTYPE##_float32(opt,progress); break; \
  }

int32 compute_resolution(const std::string& p, const GeoTransform& t, const Vector2& v) {
  if (p == "KML")      return vw::cartography::output::kml::compute_resolution(t,v);
  if (p == "TMS")      return vw::cartography::output::tms::compute_resolution(t,v);
  if (p == "UNIVIEW")  return vw::cartography::output::tms::compute_resolution(t,v);
  if (p == "GMAP")     return vw::cartography::output::tms::compute_resolution(t,v);
  if (p == "CELESTIA") return vw::cartography::output::tms::compute_resolution(t,v);
  if (p == "GIGAPAN")  return vw::cartography::output::tms::compute_resolution(t,v);
  vw_throw(LogicErr() << "Asked to compute resolution for unknown quadtree type.");
}

void get_normalize_vals(boost::shared_ptr<DiskImageResource> file,
                               const Options& opt) {
  DiskImageView<PixelRGB<float> > min_max_file(file);
  float new_lo, new_hi;
  if ( opt.nodata_set ) {
    PixelRGB<float> no_data_value( opt.nodata );
    min_max_channel_values( create_mask(min_max_file,no_data_value), new_lo, new_hi );
  } else if ( file->has_nodata_read() ) {
    PixelRGB<float> no_data_value( file->nodata_read() );
    min_max_channel_values( create_mask(min_max_file,no_data_value), new_lo, new_hi );
  } else {
    min_max_channel_values( min_max_file, new_lo, new_hi );
  }
  lo_value = std::min(new_lo, lo_value);
  hi_value = std::max(new_hi, hi_value);
  cout << "Pixel range for \"" << file->filename() << ": [" << new_lo << " " << new_hi << "]    Output dynamic range: [" << lo_value << " " << hi_value << "]" << endl;
}

std::vector<GeoReference>
vw::load_image_georeferences( const Options& opt, int& total_resolution ) {
  std::vector<GeoReference> georeferences;
  georeferences.reserve( opt.input_files.size() );

  BOOST_FOREACH(const std::string filename, opt.input_files) {
    boost::shared_ptr<DiskImageResource> file( DiskImageResourcePtr(filename) );
    vw_out() << "Adding file " << file->filename() << std::endl;

    if( opt.normalize ) get_normalize_vals(file, opt);

    GeoReference input_georef = vw::make_input_georef(file, opt);
    georeferences.push_back( input_georef );

    GeoReference output_georef(input_georef.datum());

    // Right now, we only need a WGS84 output geoereference to compute
    // the resolution. The rest of the output info will get set later.
    GeoTransform geotx( input_georef, output_georef );

    // Calculate the best resolution at 5 different points in the image,
    // as occasionally there's a singularity at the center pixel that
    // makes it extremely tiny (such as in pole-centered images).
    const int cols = file->cols();
    const int rows = file->rows();
    Vector2 res_pixel[5];
    res_pixel[0] = Vector2( cols/2, rows/2 );
    res_pixel[1] = Vector2( cols/2 + cols/4, rows/2 );
    res_pixel[2] = Vector2( cols/2 - cols/4, rows/2 );
    res_pixel[3] = Vector2( cols/2, rows/2 + rows/4 );
    res_pixel[4] = Vector2 (cols/2, rows/2 - rows/4 );
    int resolution;
    for(int i=0; i < 5; i++) {
      resolution = compute_resolution(opt.mode, geotx, res_pixel[i]);
      if( resolution > total_resolution ) total_resolution = resolution;
    }
  }

  if(opt.global_resolution > 0) {
    vw_out(VerboseDebugMessage) << "Overriding calculated resolution " << total_resolution << " with " << opt.global_resolution << std::endl;
    total_resolution = opt.global_resolution;
  }

  return georeferences;
}

GeoReference vw::make_input_georef(boost::shared_ptr<DiskImageResource> file,
                                   const Options& opt) {
  GeoReference input_georef;
  bool fail_read_georef = false;
  try {
    fail_read_georef = !read_georeference( input_georef, *file );
  } catch ( const InputErr& e ) {
    vw_out(ErrorMessage) << "Input " << file->filename() << " has malformed georeferencing information.\n";
    fail_read_georef = true;
  }

  if (opt.datum.type == "WGS84") input_georef.set_well_known_geogcs("WGS84" ); 
  if (opt.datum.type == "LUNAR") input_georef.set_well_known_geogcs("D_MOON");
  if (opt.datum.type == "MARS" ) input_georef.set_well_known_geogcs("D_MARS");
  if (opt.datum.type == "SPHERE") {
      cartography::Datum datum("USER SUPPLIED DATUM", "SPHERICAL DATUM", "Reference Meridian",
          opt.datum.sphere_radius, opt.datum.sphere_radius, 0.0);
      input_georef.set_datum(datum);
  }

  if( opt.manual ) {
    Matrix3x3 m;
    m(0,0) = double(opt.east - opt.west) / file->cols();
    m(0,2) = opt.west;
    m(1,1) = double(opt.south - opt.north) / file->rows();
    m(1,2) = opt.north;
    m(2,2) = 1;
    input_georef.set_transform( m );
  } else if ( fail_read_georef ) {
    vw_out(ErrorMessage) << "Missing input georeference. Please provide --north --south --east and --west.\n";
    exit(1);
  }

  // If the user passed in georef parameters, process them here
  if (opt.proj.type == "LAMBERT_AZIMUTHAL") input_georef.set_lambert_azimuthal(opt.proj.lat,opt.proj.lon);
  if (opt.proj.type == "LAMBERT_CONFORMAL_CONIC") 
    input_georef.set_lambert_conformal(opt.proj.p1, opt.proj.p2, opt.proj.lat, opt.proj.lon);
  if (opt.proj.type == "MERCATOR")    input_georef.set_mercator(opt.proj.lat,opt.proj.lon,opt.proj.scale);
  if (opt.proj.type == "ORTHOGRAPHIC")            input_georef.set_orthographic(opt.proj.lat,opt.proj.lon);
  if (opt.proj.type == "PLATE_CARREE")            input_georef.set_geographic();
  if (opt.proj.type == "SINUSOIDAL")              input_georef.set_sinusoidal(opt.proj.lon);
  if (opt.proj.type == "STEREOGRAPHIC")           
    input_georef.set_stereographic(opt.proj.lat,opt.proj.lon,opt.proj.scale);
  if (opt.proj.type == "TRANSVERSE_MERCATOR")     
    input_georef.set_transverse_mercator(opt.proj.lat,opt.proj.lon,opt.proj.scale);
  if (opt.proj.type == "UTM") input_georef.set_UTM( abs(opt.proj.utm_zone), opt.proj.utm_zone > 0 );

  // Handle nudge arguments (x and y shifts in projected coordinates)
  if( opt.nudge_x || opt.nudge_y ) {
    Matrix3x3 m = input_georef.transform();
    m(0,2) += opt.nudge_x;
    m(1,2) += opt.nudge_y;
    input_georef.set_transform( m );
  }

  return input_georef;
}

int handle_options(int argc, char *argv[], Options& opt) {
  po::options_description general_options("Description: Turns georeferenced image(s) into a quadtree with geographical metadata\n\nGeneral options");
  general_options.add_options()
    ("output-name,o", po::value(&opt.output_file_name), "Specify the base output directory");

  general_options.add(vw::GdalWriteOptionsDescription(opt));

const std::string channel_options_str    = "DEFAULT, UINT8, UINT16, INT16, FLOAT";
const std::string mode_options_str       = "NONE, KML, TMS, UNIVIEW, GMAP, CELESTIA, GIGAPAN";
const std::string datum_options_str      = "NONE, WGS84, LUNAR, MARS, SPHERE";
const std::string projection_options_str = "DEFAULT, NONE, SINUSOIDAL, MERCATOR, TRANSVERSE_MERCATOR, ORTHOGRAPHIC, STEREOGRAPHIC, LAMBERT_AZIMUTHAL, LAMBERT_CONFORMAL_CONIC, UTM, PLATE_CARREE";

  po::options_description input_options("Input options");
  string datum_desc = string("Override input datum ["            ) + datum_options_str      + "]";
  string mode_desc  = string("Specify the output metadata type [") + mode_options_str       + "]";
  string proj_desc  = string("Projection type ["                 ) + projection_options_str + "]";
  string chan_desc  = string("Output channel type ["             ) + channel_options_str    + "]";

  input_options.add_options()
    ("force-datum" , po::value(&opt.datum.type)                       , datum_desc.c_str())
    ("datum-radius", po::value(&opt.datum.sphere_radius)              , "Radius to use for --force-datum SPHERE")
    ("pixel-scale" , po::value(&opt.pixel_scale)->default_value(1.0)  , "Scale factor to apply to pixels")
    ("pixel-offset", po::value(&opt.pixel_offset)->default_value(0.0) , "Offset to apply to pixels")
    ("normalize"   , po::bool_switch(&opt.normalize)                  , "Normalize input images so that their full dynamic range falls in between [0,255].")
    ("nodata"      , po::value(&opt.nodata)                           , "Set the input's nodata value so that it will be transparent in output");

  po::options_description output_options("Output options");
  output_options.add_options()
    ("mode,m"           , po::value(&opt.mode)->default_value("KML"), mode_desc.c_str())
    ("file-type"        , po::value(&opt.output_file_type)                       , "Output file type.  (Choose \'auto\' to generate jpgs in opaque areas and png images where there is transparency.)")
    ("channel-type"     , po::value(&opt.channel_type)->default_value("DEFAULT"), chan_desc.c_str())
    ("module-name"      , po::value(&opt.module_name)                            , "The module where the output will be placed. Ex: marsds for Uniview,  or Sol/Mars for Celestia")
    ("terrain"          , po::bool_switch(&opt.terrain)                          , "Outputs image files suitable for a Uniview terrain view. Implies output format as PNG, channel type uint16. Uniview only")
    ("jpeg-quality"     , po::value(&opt.jpeg_quality)                           , "JPEG quality factor (0.0 to 1.0)")
    ("png-compression"  , po::value(&opt.png_compression)                        , "PNG compression level (0 to 9)")
    ("tile-size"        , po::value(&opt.tile_size)                              , "Tile size in pixels")
    ("max-lod-pixels"   , po::value(&opt.kml.max_lod_pixels)->default_value(1024), "Max LoD in pixels, or -1 for none (kml only)")
    ("draw-order-offset", po::value(&opt.kml.draw_order_offset)->default_value(0), "Offset for the <drawOrder> tag for this overlay (kml only)")
    ("multiband"        , po::bool_switch(&opt.multiband)                        , "Composite images using multi-band blending")
    ("aspect-ratio"     , po::value(&opt.aspect_ratio)                           , "Pixel aspect ratio (for polar overlays; should be a power of two)")
    ("global-resolution", po::value(&opt.global_resolution)                      , "Override the global pixel resolution; should be a power of two");

  po::options_description projection_options("Input projection options");
  projection_options.add_options()
    ("north"      ,  po::value(&opt.north)         , "The northernmost latitude in projection units")
    ("south"      ,  po::value(&opt.south)         , "The southernmost latitude in projection units")
    ("east"       ,  po::value(&opt.east)          , "The easternmost longitude in projection units")
    ("west"       ,  po::value(&opt.west)          , "The westernmost longitude in projection units")
    ("global"     ,  po::bool_switch(&opt.global)  , "Override image size to global (in lonlat)")
    ("projection" ,  po::value(&opt.proj.type)->default_value("DEFAULT")     , proj_desc.c_str())
    ("utm-zone"   ,  po::value(&opt.proj.utm_zone) , "Set zone for --projection UTM (negative for south)")
    ("proj-lat"   ,  po::value(&opt.proj.lat)      , "The center of projection latitude")
    ("proj-lon"   ,  po::value(&opt.proj.lon)      , "The center of projection longitude")
    ("proj-scale" ,  po::value(&opt.proj.scale)    , "The projection scale")
    ("p1"         ,  po::value(&opt.proj.p1)       , "parallel for Lambert Conformal Conic projection")
    ("p2"         ,  po::value(&opt.proj.p2)       , "parallel for Lambert Conformal Conic projection")
    ("nudge-x"    ,  po::value(&opt.nudge_x)       , "Nudge the image, in projected coordinates")
    ("nudge-y"    ,  po::value(&opt.nudge_y)       , "Nudge the image, in projected coordinates");

  po::options_description hidden_options("");
  hidden_options.add_options()
    ("input-file", po::value(&opt.input_files));

  po::options_description options("Allowed options");
  options.add(general_options).add(input_options).add(projection_options).add(output_options).add(hidden_options);

  po::positional_options_description p;
  p.add("input-file", -1);

  std::ostringstream usage;
  usage << "Usage: image2qtree [options] <filename>..." <<endl << endl;
  usage << general_options    << endl;
  usage << input_options      << endl;
  usage << output_options     << endl;
  usage << projection_options << endl;

  try {
    namespace ps = po::command_line_style;
    int style = ps::unix_style & ~ps::allow_guessing;
    po::variables_map vm;
    po::store( po::command_line_parser( argc, argv ).style(style).options(options).positional(p).run(), vm );
    po::notify( vm );
    
    if (!vm.count("input-nodata-value"))
      opt.nodata_set = false;
    else
      opt.nodata_set = true;
    
    opt.validate();
  } catch (const po::error& e) {
    cerr << usage.str() << endl
         << "Failed to parse command line arguments:" << endl
         << "\t" << e.what() << endl;
    return false;
  } catch (const Usage& e) {
    const char* msg = e.what();
      cerr << usage.str() << endl;
    if (strlen(msg) > 0) {
           cerr << endl
           << "Invalid argument:" << endl
           << "\t" << msg << endl;
    }
    return false;
  }

  opt.setVwSettingsFromOpt();
  
  // Make sure all these string parameters are upper case
  boost::to_upper(opt.mode);
  boost::to_upper(opt.datum.type);
  boost::to_upper(opt.channel_type);
  boost::to_upper(opt.proj.type);
  return true;
}

/// This function picks the correct function call to do the work.
void run(const Options& opt) {
  TerminalProgressCallback tpc( "image2qtree", "");
  const ProgressCallback *progress = &tpc;

  // Get the right pixel/channel type, and call the mosaic.
  ImageFormat fmt = vw::image_format(opt.input_files[0]);

  if(opt.channel_type != "DEFAULT")
    fmt.channel_type = channel_name_to_enum(opt.channel_type);
  if (fmt.channel_type == VW_CHANNEL_UNKNOWN)
    vw_throw(LogicErr() << "Unrecognized channel type: " << opt.channel_type);

  // Convert non-alpha channel images into images with one for the composite.
  switch(fmt.pixel_format) {
  case VW_PIXEL_GRAY:
  case VW_PIXEL_GRAYA:
    SWITCH_ON_CHANNEL_TYPE(PixelGrayA); break;
  case VW_PIXEL_RGB:
  case VW_PIXEL_RGBA:
  default:
    SWITCH_ON_CHANNEL_TYPE(PixelRGBA); break;
  }
}

/// Main, input parsing and error handling, calls run()
int main(int argc, char **argv) {
  Options opt;
  if (!handle_options(argc, argv, opt))
    return 1;

  try {
    run(opt);
  } catch ( const ArgumentErr& e ) {
    vw_out() << e.what() << std::endl;
    return 1;
  } catch ( const Exception& e ) {
    std::cerr << "\n\nVW Error: " << e.what() << std::endl;
    return 1;
  } catch ( const std::bad_alloc& e ) {
    std::cerr << "\n\nError: Ran out of Memory!" << std::endl;
    return 1;
  } catch ( const std::exception& e ) {
    std::cerr << "\n\nError: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}
