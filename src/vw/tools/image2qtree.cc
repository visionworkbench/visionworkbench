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

#include <vw/tools/Common.h>

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <map>

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/foreach.hpp>
namespace po = boost::program_options;

#include <vw/Core/Cache.h>
#include <vw/Core/ProgressCallback.h>
#include <vw/Math/Matrix.h>
#include <vw/Image/Transform.h>
#include <vw/Image/MaskViews.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/FileIO/DiskImageResourceJPEG.h>
#include <vw/FileIO/DiskImageResourcePNG.h>
#include <vw/FileIO/DiskImageResourceGDAL.h>
#include <vw/FileIO/DiskImageView.h>
#include <vw/Cartography/GeoReference.h>
#include <vw/Cartography/GeoTransform.h>
#include <vw/Mosaic.h>
using namespace vw;
using namespace vw::math;
using namespace vw::cartography;
using namespace vw::mosaic;
using std::string;
using vw::tools::Usage;
using vw::tools::Tristate;
using std::cout;
using std::cerr;
using std::endl;

VW_DEFINE_ENUM(Channel, 5, (NONE, UINT8, UINT16, INT16, FLOAT));
VW_DEFINE_ENUM(Mode, 8, (NONE, KML, TMS, UNIVIEW, GMAP, CELESTIA, GIGAPAN, GIGAPAN_NOPROJ));
VW_DEFINE_ENUM(DatumOverride, 5, (NONE, WGS84, LUNAR, MARS, SPHERE))
VW_DEFINE_ENUM(Projection, 10, (
    NONE,
    SINUSOIDAL,
    MERCATOR,
    TRANSVERSE_MERCATOR,
    ORTHOGRAPHIC,
    STEREOGRAPHIC,
    LAMBERT_AZIMUTHAL,
    LAMBERT_CONFORMAL_CONIC,
    UTM,
    PLATE_CARREE))

struct Options {
  Options() :
    output_file_type("png"),
    module_name("", true),
    tile_size(256),
    jpeg_quality(0, true),
    png_compression(0, true),
    pixel_scale(1),
    pixel_offset(0),
    aspect_ratio(1),
    global_resolution(0, true),
    nodata(0, true),
    north(0, true), south(0, true),
    east(0, true), west(0, true),
    channel_type(Channel::NONE),
    multiband(false),
    help(false),
    normalize(false),
    terrain(false),
    manual(false),
    global(false)
    {}

  std::vector<string> input_files;

  string output_file_name;
  Tristate<string> output_file_type;
  Tristate<string> module_name;
  Tristate<double> nudge_x, nudge_y;
  Tristate<uint32> tile_size;
  Tristate<float>  jpeg_quality;
  Tristate<uint32> png_compression;
  Tristate<float>  pixel_scale, pixel_offset;
  Tristate<int32>  aspect_ratio;
  Tristate<uint32> global_resolution;
  Tristate<float>  nodata;
  Tristate<float>  north, south, east, west;

  Channel channel_type;
  Mode mode;

  bool multiband;
  bool help;
  bool normalize;
  bool terrain;
  bool manual;
  bool global;

  struct {
    uint32 draw_order_offset;
    uint32 max_lod_pixels;
  } kml;

  struct proj_{
    Projection type;
    Tristate<double> lat, lon, scale /*=1*/;
    Tristate<double> p1, p2;
    Tristate<int32> utm_zone;
    proj_() :
      type(Projection::NONE),
      scale(1),
      utm_zone(0, true) {}
  } proj;

  struct datum_ {
    DatumOverride type;
    Tristate<float> sphere_radius;
    datum_() : type(DatumOverride::NONE), sphere_radius(0, true) {}
  } datum;

  void validate() {
    VW_ASSERT(!help, Usage());
    VW_ASSERT(input_files.size() > 0, Usage() << "Need at least one input image");

    if (datum.type == DatumOverride::SPHERE)
      VW_ASSERT(datum.sphere_radius.set(), Usage() << "Sphere datum override requires a radius");

    if(output_file_name.empty())
      output_file_name = fs::path(input_files[0]).replace_extension().string();

    if (global || north.set() || south.set() || east.set() || west.set()) {
      VW_ASSERT(input_files.size() == 1,
          Usage() << "Cannot override georeference information on multiple images");
      VW_ASSERT(global || (north.set() && south.set() && east.set() && west.set()),
          Usage() << "If you provide one, you must provide all of: --north --south --east --west");
      if (global) {
        north = 90; south = -90; east = 180; west = -180;
      }
      manual = true;
    }

    switch (mode) {
      case Mode::NONE:
      case Mode::GIGAPAN_NOPROJ:
        VW_ASSERT(input_files.size() == 1, Usage() << "Non-georeferenced images cannot be composed");
        break;
      case Mode::CELESTIA:
      case Mode::UNIVIEW:
        VW_ASSERT(module_name.set(), Usage() << "Uniview and Celestia require --module-name");
        break;
      default:
        /* nothing */
        break;
    }

    if (jpeg_quality.set())
      DiskImageResourceJPEG::set_default_quality( jpeg_quality );
    if (png_compression.set())
      DiskImageResourcePNG::set_default_compression_level( png_compression );
  }
};

// For image stretching.
static float lo_value = ScalarTypeLimits<float>::highest();
static float hi_value = ScalarTypeLimits<float>::lowest();

vw::int32 compute_resolution(const Mode& p, const GeoTransform& t, const Vector2& v) {
  switch(p.value()) {
    case Mode::KML:      return vw::cartography::output::kml::compute_resolution(t,v);
    case Mode::TMS:      return vw::cartography::output::tms::compute_resolution(t,v);
    case Mode::UNIVIEW:  return vw::cartography::output::tms::compute_resolution(t,v);
    case Mode::GMAP:     return vw::cartography::output::tms::compute_resolution(t,v);
    case Mode::CELESTIA: return vw::cartography::output::tms::compute_resolution(t,v);
    case Mode::GIGAPAN:  return vw::cartography::output::tms::compute_resolution(t,v);
    default: vw_throw(LogicErr() << "Asked to compute resolution for unknown profile " << p.string());
  }
}

static void get_normalize_vals(boost::shared_ptr<DiskImageResource> file,
                               const Options& opt) {

  DiskImageView<PixelRGB<float> > min_max_file(file);
  float new_lo, new_hi;
  if ( opt.nodata.set() ) {
    PixelRGB<float> no_data_value( opt.nodata.value() );
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

template <class PixelT>
void do_normal_mosaic(const Options& opt, const ProgressCallback *progress) {
    DiskImageView<PixelT> img(opt.input_files[0]);
    QuadTreeGenerator quadtree(img, opt.output_file_name);
    quadtree.set_tile_size( opt.tile_size );
    quadtree.set_file_type( opt.output_file_type );

    if (opt.mode == Mode::GIGAPAN_NOPROJ) {
      GigapanQuadTreeConfig config;
      config.configure( quadtree );
    }

    quadtree.generate( *progress );
}

GeoReference make_input_georef(boost::shared_ptr<DiskImageResource> file,
                               const Options& opt) {
  GeoReference input_georef;
  bool fail_read_georef = false;
  try {
    fail_read_georef = !read_georeference( input_georef, *file );
  } catch ( InputErr const& e ) {
    vw_out(ErrorMessage) << "Input " << file->filename() << " has malformed georeferencing information.\n";
    fail_read_georef = true;
  }

  switch(opt.datum.type) {
    case DatumOverride::WGS84: input_georef.set_well_known_geogcs("WGS84");  break;
    case DatumOverride::LUNAR: input_georef.set_well_known_geogcs("D_MOON"); break;
    case DatumOverride::MARS:  input_georef.set_well_known_geogcs("D_MARS"); break;
    case DatumOverride::SPHERE: {
      cartography::Datum datum("USER SUPPLIED DATUM", "SPHERICAL DATUM", "Reference Meridian",
          opt.datum.sphere_radius, opt.datum.sphere_radius, 0.0);
      input_georef.set_datum(datum);
      break;
    }
    case DatumOverride::NONE: break;
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

  switch (opt.proj.type) {
    case Projection::LAMBERT_AZIMUTHAL:       input_georef.set_lambert_azimuthal(opt.proj.lat,opt.proj.lon); break;
    case Projection::LAMBERT_CONFORMAL_CONIC: input_georef.set_lambert_conformal(opt.proj.p1, opt.proj.p2, opt.proj.lat, opt.proj.lon); break;
    case Projection::MERCATOR:                input_georef.set_mercator(opt.proj.lat,opt.proj.lon,opt.proj.scale); break;
    case Projection::ORTHOGRAPHIC:            input_georef.set_orthographic(opt.proj.lat,opt.proj.lon); break;
    case Projection::PLATE_CARREE:            input_georef.set_geographic(); break;
    case Projection::SINUSOIDAL:              input_georef.set_sinusoidal(opt.proj.lon); break;
    case Projection::STEREOGRAPHIC:           input_georef.set_stereographic(opt.proj.lat,opt.proj.lon,opt.proj.scale); break;
    case Projection::TRANSVERSE_MERCATOR:     input_georef.set_transverse_mercator(opt.proj.lat,opt.proj.lon,opt.proj.scale); break;
    case Projection::UTM:                     input_georef.set_UTM( abs(opt.proj.utm_zone), opt.proj.utm_zone > 0 ); break;
    case Projection::NONE: break;
  }

  if( opt.nudge_x || opt.nudge_y ) {
    Matrix3x3 m = input_georef.transform();
    m(0,2) += opt.nudge_x;
    m(1,2) += opt.nudge_y;
    input_georef.set_transform( m );
  }

  return input_georef;
}

template <class PixelT>
void do_mosaic(const Options& opt, const ProgressCallback *progress)
{
  typedef typename PixelChannelType<PixelT>::type ChannelT;

  // If we're not outputting any special sort of mosaic (just a regular old
  // quadtree, no georeferencing, no metadata), we use a different
  // function.
  if(opt.mode == Mode::NONE || opt.mode == Mode::GIGAPAN_NOPROJ) {
    do_normal_mosaic<PixelT>(opt, progress);
    return;
  }

  // Read in georeference info and compute total resolution.
  int total_resolution = 1024;
  std::vector<GeoReference> georeferences;

  BOOST_FOREACH(const string filename, opt.input_files) {
    boost::shared_ptr<DiskImageResource> file( DiskImageResource::open(filename) );
    cout << "Adding file " << file->filename() << endl;

    if( opt.normalize ) get_normalize_vals(file, opt);

    GeoReference input_georef = make_input_georef(file, opt);
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

  if(opt.global_resolution.set()) {
    vw_out(VerboseDebugMessage) << "Overriding calculated resolution " << total_resolution << " with " << opt.global_resolution.value() << endl;
    total_resolution = opt.global_resolution;
  }

  boost::shared_ptr<QuadTreeConfig> config = QuadTreeConfig::make(opt.mode.string());

  // Now that we have the best resolution, we can get our output_georef.
  int xresolution = total_resolution / opt.aspect_ratio, yresolution = total_resolution;

  GeoReference output_georef = config->output_georef(xresolution, yresolution);
  vw_out(VerboseDebugMessage, "tool") << "Output Georef:\n" << output_georef << endl;

  // Configure the composite.
  ImageComposite<PixelT> composite;

  // Add the transformed image files to the composite.
  for(unsigned i=0; i < opt.input_files.size(); i++) {
    const std::string& filename = opt.input_files[i];
    const GeoReference& input_ref = georeferences[i];

    boost::shared_ptr<DiskImageResource> file( DiskImageResource::open(filename) );
    GeoTransform geotx( input_ref, output_georef );
    ImageViewRef<PixelT> source = DiskImageView<PixelT>( file );

    if ( opt.nodata.set() ) {
      vw_out(VerboseDebugMessage, "tool") << "Using nodata value: "
                                          << opt.nodata.value() << "\n";
      source = mask_to_alpha(create_mask(pixel_cast<typename PixelWithoutAlpha<PixelT>::type >(source),ChannelT(opt.nodata.value())));
    } else if ( file->has_nodata_read() ) {
      vw_out(VerboseDebugMessage, "tool") << "Using nodata value: "
                                          << file->nodata_read() << "\n";
      source = mask_to_alpha(create_mask(pixel_cast<typename PixelWithoutAlpha<PixelT>::type >(source),ChannelT(file->nodata_read())));
    }

    bool global = boost::trim_copy(input_ref.proj4_str())=="+proj=longlat" &&
      fabs(input_ref.lonlat_to_pixel(Vector2(-180,0)).x()) < 1 &&
      fabs(input_ref.lonlat_to_pixel(Vector2(180,0)).x() - source.cols()) < 1 &&
      fabs(input_ref.lonlat_to_pixel(Vector2(0,90)).y()) < 1 &&
      fabs(input_ref.lonlat_to_pixel(Vector2(0,-90)).y() - source.rows()) < 1;

    // Do various modifications to the input image here.
    if( opt.pixel_scale.set() || opt.pixel_offset.set() ) {
      vw_out(VerboseDebugMessage, "tool") << "Apply input scaling: "
                           << opt.pixel_scale.value() << " offset: "
                           << opt.pixel_offset.value() << "\n";
      source = channel_cast_rescale<ChannelT>( source * opt.pixel_scale.value() + opt.pixel_offset.value() );
    }

    if( opt.normalize ) {
      vw_out(VerboseDebugMessage, "tool") << "Apply normalizing: ["
                                          << lo_value << ", "
                                          << hi_value << "]\n";
      typedef ChannelRange<ChannelT> range_type;
      source = normalize_retain_alpha(source, lo_value, hi_value,
                                      range_type::min(), range_type::max());
    }

    BBox2i bbox = geotx.forward_bbox( BBox2i(0,0,source.cols(),source.rows()) );
    if (global) {
      vw_out() << "\t--> Detected global overlay.  Using cylindrical edge extension to hide the seam.\n";
      source = crop( transform( source, geotx, source.cols(), source.rows(), CylindricalEdgeExtension() ), bbox );
    }
    else
      source = crop( transform( source, geotx ), bbox );

    // Images that wrap the date line must be added to the composite
    // on both sides.
    if( bbox.max().x() > total_resolution ) {
      composite.insert( source, bbox.min().x()-total_resolution, bbox.min().y() );
    }
    // Images that are in the 180-360 range *only* go on the other side.
    if( bbox.min().x() < xresolution ) {
      composite.insert( source, bbox.min().x(), bbox.min().y() );
    }
  }

  // This box represents the entire input data set, in pixels, in the output
  // projection space. This should NOT include the extra data used to hide
  // seams and such.
  BBox2i total_bbox = composite.bbox();
  total_bbox.crop(BBox2i(0,0,xresolution,yresolution));

  if(opt.mode == Mode::KML) {
    BBox2i bbox = total_bbox;
    // Compute a tighter Google Earth coordinate system aligned bounding box.
    int dim = 2 << (int)(log( (double)(std::max)(bbox.width(),bbox.height()) )/log(2.));
    if( dim > total_resolution ) dim = total_resolution;
    total_bbox = BBox2i( (bbox.min().x()/dim)*dim, (bbox.min().y()/dim)*dim, dim, dim );
    if( ! total_bbox.contains( bbox ) ) {
      if( total_bbox.max().x() == xresolution ) total_bbox.min().x() -= dim;
      else total_bbox.max().x() += dim;
      if( total_bbox.max().y() == yresolution ) total_bbox.min().y() -= dim;
      else total_bbox.max().y() += dim;
    }
  }

  // Prepare the composite.
  if(!opt.multiband)
    composite.set_draft_mode( true );
  composite.prepare( total_bbox, *progress );

  QuadTreeGenerator quadtree( composite, opt.output_file_name );

  // This whole bit here is terrible. This functionality should be moved into
  // the Config base class somehow.
  if( opt.mode == Mode::KML ) {
    KMLQuadTreeConfig *c2 = dynamic_cast<KMLQuadTreeConfig*>(config.get());
    BBox2 ll_bbox( -180.0 + (360.0*total_bbox.min().x())/xresolution,
                    180.0 - (360.0*total_bbox.max().y())/yresolution,
                    (360.0*total_bbox.width())/xresolution,
                    (360.0*total_bbox.height())/yresolution );

    c2->set_longlat_bbox( ll_bbox );
    c2->set_max_lod_pixels( opt.kml.max_lod_pixels );
    c2->set_draw_order_offset( opt.kml.draw_order_offset );
  } else if( opt.mode == Mode::CELESTIA ) {
    CelestiaQuadTreeConfig *c2 = dynamic_cast<CelestiaQuadTreeConfig*>(config.get());
    c2->set_module(opt.module_name.value());
  } else if( opt.mode == Mode::UNIVIEW ) {
    UniviewQuadTreeConfig *c2 = dynamic_cast<UniviewQuadTreeConfig*>(config.get());
    c2->set_terrain(opt.terrain);
    c2->set_module(opt.module_name.value());
  } else if ( opt.mode == Mode::GIGAPAN ) {
    GigapanQuadTreeConfig *c2 = dynamic_cast<GigapanQuadTreeConfig*>(config.get());
    BBox2 ll_bbox( -180.0 + (360.0*total_bbox.min().x())/xresolution,
                    180.0 - (360.0*total_bbox.max().y())/yresolution,
                    (360.0*total_bbox.width())/xresolution,
                    (360.0*total_bbox.height())/yresolution );
    c2->set_longlat_bbox( ll_bbox );
  }

  config->configure(quadtree);

  if (opt.tile_size.set())
    quadtree.set_tile_size(opt.tile_size);
  if (opt.output_file_type.set())
    quadtree.set_file_type(opt.output_file_type);

  // This box represents the input data, shifted such that total_bbox.min() is
  // the origin, and cropped to the size of the output resolution.
  BBox2i data_bbox = composite.bbox();
  data_bbox.crop( BBox2i(0,0,total_bbox.width(),total_bbox.height()));

  quadtree.set_crop_bbox(data_bbox);

  // Generate the composite.
  vw_out() << "Generating " << opt.mode.string() << " overlay..." << endl;
  quadtree.generate(*progress);
}



int handle_options(int argc, char *argv[], Options& opt) {
  po::options_description general_options("Description: Turns georeferenced image(s) into a quadtree with geographical metadata\n\nGeneral Options");
  general_options.add_options()
    ("output-name,o", po::value(&opt.output_file_name), "Specify the base output directory")
    ("help,h", po::bool_switch(&opt.help), "Display this help message");

  po::options_description input_options("Input Options");
  string datum_desc = string("Override input datum [") + DatumOverride::list() + "]";
  string mode_desc  = string("Specify the output metadata type [") + Mode::list() + "]";
  string proj_desc  = string("Projection type [") + Projection::list() + "]";
  string chan_desc  = string("Output channel type [") + Channel::list() + "]";

  input_options.add_options()
    ("force-datum" , po::value(&opt.datum.type)                       , datum_desc.c_str())
    ("datum-radius", po::value(&opt.datum.sphere_radius)              , "Radius to use for --force-datum SPHERE")
    ("pixel-scale" , po::value(&opt.pixel_scale)->default_value(1.0)  , "Scale factor to apply to pixels")
    ("pixel-offset", po::value(&opt.pixel_offset)->default_value(0.0) , "Offset to apply to pixels")
    ("normalize"   , po::bool_switch(&opt.normalize)                  , "Normalize input images so that their full dynamic range falls in between [0,255].")
    ("nodata"      , po::value(&opt.nodata)                           , "Set the input's nodata value so that it will be transparent in output");

  po::options_description output_options("Output Options");
  output_options.add_options()
    ("mode,m"           , po::value(&opt.mode)->default_value(Mode("kml"))       , mode_desc.c_str())
    ("file-type"        , po::value(&opt.output_file_type)                       , "Output file type.  (Choose \'auto\' to generate jpgs in opaque areas and png images where there is transparency.)")
    ("channel-type"     , po::value(&opt.channel_type)                           , chan_desc.c_str())
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

  po::options_description projection_options("Input Projection Options");
  projection_options.add_options()
    ("north"      ,  po::value(&opt.north)         , "The northernmost latitude in projection units")
    ("south"      ,  po::value(&opt.south)         , "The southernmost latitude in projection units")
    ("east"       ,  po::value(&opt.east)          , "The easternmost longitude in projection units")
    ("west"       ,  po::value(&opt.west)          , "The westernmost longitude in projection units")
    ("global"     ,  po::bool_switch(&opt.global)  , "Override image size to global (in lonlat)")
    ("projection" ,  po::value(&opt.proj.type)     , proj_desc.c_str())
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

  po::options_description options("Allowed Options");
  options.add(general_options).add(input_options).add(projection_options).add(output_options).add(hidden_options);

  po::positional_options_description p;
  p.add("input-file", -1);

  std::ostringstream usage;
  usage << "Usage: image2qtree [options] <filename>..." <<endl << endl;
  usage << general_options << endl;
  usage << input_options << endl;
  usage << output_options << endl;
  usage << projection_options << endl;

  try {
    namespace ps = po::command_line_style;
    int style = ps::unix_style & ~ps::allow_guessing;
    po::variables_map vm;
    po::store( po::command_line_parser( argc, argv ).style(style).options(options).positional(p).run(), vm );
    po::notify( vm );
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
  return true;
}

int run(const Options& opt) {
  TerminalProgressCallback tpc( "tools.image2qtree", "");
  const ProgressCallback *progress = &tpc;

  // Get the right pixel/channel type, and call the mosaic.
  ImageFormat fmt = tools::taste_image(opt.input_files[0]);

  if(opt.channel_type != Channel::NONE)
    fmt.channel_type = channel_name_to_enum(opt.channel_type.string());

  // Convert non-alpha channel images into images with one for the
  // composite.
  switch(fmt.pixel_format) {
    case VW_PIXEL_GRAY:
    case VW_PIXEL_GRAYA:
      switch(fmt.channel_type) {
        case VW_CHANNEL_UINT8:  do_mosaic<PixelGrayA<uint8>   >(opt, progress); break;
        case VW_CHANNEL_INT16:  do_mosaic<PixelGrayA<int16>   >(opt, progress); break;
        case VW_CHANNEL_UINT16: do_mosaic<PixelGrayA<uint16>  >(opt, progress); break;
        default:                do_mosaic<PixelGrayA<float32> >(opt, progress); break;
      }
      break;
    case VW_PIXEL_RGB:
    case VW_PIXEL_RGBA:
    default:
      switch(fmt.channel_type) {
        case VW_CHANNEL_UINT8:  do_mosaic<PixelRGBA<uint8>   >(opt, progress); break;
        case VW_CHANNEL_INT16:  do_mosaic<PixelRGBA<int16>   >(opt, progress); break;
        case VW_CHANNEL_UINT16: do_mosaic<PixelRGBA<uint16>  >(opt, progress); break;
        default:                do_mosaic<PixelRGBA<float32> >(opt, progress); break;
      }
      break;
  }

  return 0;
}

int main(int argc, char **argv) {
  Options opt;
  if (!handle_options(argc, argv, opt))
    return 1;
  return run(opt);
}
