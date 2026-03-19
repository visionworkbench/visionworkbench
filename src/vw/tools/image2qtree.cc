// __BEGIN_LICENSE__
//  Copyright (c) 2006-2026, United States Government as represented by the
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
/// a quadtree for that image that is viewable in Google Earth via KML.

#include <vw/Core/ProgressCallback.h>
#include <vw/Core/Log.h>
#include <vw/Math/BBox.h>
#include <vw/Math/Vector.h>
#include <vw/Image/Transform.h>
#include <vw/Image/MaskViews.h>
#include <vw/Image/Statistics.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/FileIO/DiskImageResourceJPEG.h>
#include <vw/FileIO/DiskImageResourcePNG.h>
#include <vw/FileIO/DiskImageView.h>
#include <vw/FileIO/GdalWriteOptions.h>
#include <vw/FileIO/GdalWriteOptionsDesc.h>
#include <vw/Cartography/GeoReference.h>
#include <vw/Cartography/GeoReferenceUtils.h>
#include <vw/Cartography/GeoTransform.h>
#include <vw/Mosaic/ImageComposite.h>
#include <vw/Mosaic/KMLQuadTreeConfig.h>
#include <vw/Mosaic/QuadTreeConfig.h>
#include <vw/Mosaic/QuadTreeGenerator.h>
#include <vw/tools/Common.h>

#include <boost/algorithm/string/trim.hpp>
#include <boost/smart_ptr/shared_ptr.hpp>
#include <boost/foreach.hpp>
#include <boost/program_options.hpp>

#include <fstream>

namespace fs = boost::filesystem;
namespace po = boost::program_options;

using namespace vw;
using namespace vw::math;
using namespace vw::cartography;
using namespace vw::mosaic;
using std::string;
using vw::tools::Usage;
using std::cout;
using std::cerr;
using std::endl;

// The single pixel type used throughout. This is a preview tool,
// so uint8 RGBA is sufficient for all inputs.
typedef PixelRGBA<uint8> PixelT;
typedef uint8 ChannelT;

template <class ImageT, class TransformT>
TransformView<ImageT, TransformT>
transform_only(ImageViewBase<ImageT> const& v,
               TransformT const& transform_func) {
  return TransformView<ImageT, TransformT>(v.impl(), transform_func);
}

struct Options: GdalWriteOptions {

  Options():
    output_file_type(""),
    module_name(""),
    nudge_x(0), nudge_y(0),
    tile_size(0),
    jpeg_quality(-9999),
    png_compression(99999),
    pixel_scale(0),
    pixel_offset(0),
    aspect_ratio(1),
    global_resolution(0),
    nodata(0),
    nodata_set(false),
    north(0), south(0),
    east(0), west(0),
    multiband(false),
    help(false),
    normalize(false),
    manual(false),
    global(false)
  {}

  std::vector<std::string> input_files;

  std::string output_file_name;
  std::string output_file_type;
  std::string module_name;
  double      nudge_x, nudge_y;
  uint32      tile_size;
  float       jpeg_quality;
  uint32      png_compression;
  float       pixel_scale, pixel_offset;
  int32       aspect_ratio;
  uint32      global_resolution;
  float       nodata;
  bool        nodata_set;
  float       north, south, east, west;

  std::string mode; // Quadtree type

  bool multiband;
  bool help;
  bool normalize;
  bool manual;
  bool global;

  struct {
    uint32 draw_order_offset;
    uint32 max_lod_pixels;
  } kml;

  struct proj_ {
    std::string type;
    double    lat, lon, scale;
    double    p1, p2;
    int32     utm_zone;
    proj_():
      type("DEFAULT"),
      scale(1),
      utm_zone(0) {}
  } proj;

  struct datum_ {
    std::string type;
    float sphere_radius;
    datum_(): type("NONE"), sphere_radius(0) {}
  } datum;

  void validate() {
    VW_ASSERT(!help, Usage());
    VW_ASSERT(input_files.size() > 0,
              Usage() << "Need at least one input image");

    if (datum.type == "SPHERE")
      VW_ASSERT(datum.sphere_radius > 0,
                Usage() << "Sphere datum override requires a radius");

    if (output_file_name.empty())
      output_file_name =
        fs::path(input_files[0]).replace_extension().string();

    // Handle options for manually specifying projection bounds
    if (global || north != 0 || south != 0 || east != 0 || west != 0) {
      VW_ASSERT(input_files.size() == 1,
                Usage() << "Cannot override georeference information "
                "on multiple images");
      VW_ASSERT(global || ((north != south) && (east != west)),
                Usage() << "If you provide one, you must provide all "
                "of: --north --south --east --west");
      if (global) {
        north = 90; south = -90; east = 180; west = -180;
      }
      manual = true;
    }

    if (mode == "NONE")
      VW_ASSERT(input_files.size() == 1,
                Usage() << "Non-georeferenced images cannot be composed");
    if (proj.type == "NONE")
      VW_ASSERT(input_files.size() == 1,
                Usage() << "Non-georeferenced images cannot be composed");
    if (jpeg_quality > 0)
      DiskImageResourceJPEG::set_default_quality(jpeg_quality);
    if (png_compression != 99999)
      DiskImageResourcePNG::set_default_compression_level(png_compression);
  }
};

// For image stretching.
static float lo_value = ScalarTypeLimits<float>::highest();
static float hi_value = ScalarTypeLimits<float>::lowest();

// Small local utilities
namespace local {
  namespace kml {
    // Returns the number of pixels per planetary circumference,
    // rounding up to a power of two.
    template <class TransformT>
    inline int32 compute_resolution(TransformT const& tx,
                                    Vector2 const& pixel) {
      Vector2 pos      = tx.forward(pixel);
      Vector2 x_vector = tx.forward(pixel + Vector2(1, 0)) - pos;
      Vector2 y_vector = tx.forward(pixel + Vector2(0, 1)) - pos;
      double degrees_per_pixel =
        (std::min)(norm_2(x_vector), norm_2(y_vector));
      double pixels_per_circumference = 360.0 / degrees_per_pixel;
      int scale_exponent =
        (int)ceil(log(pixels_per_circumference) / log(2.0));
      if (scale_exponent >= 31) scale_exponent = 30;
      return 1 << scale_exponent;
    }
  } // namespace local::kml
} // namespace local

int32 compute_resolution(const std::string& p, const GeoTransform& t,
                         const Vector2& v) {
  if (p == "KML")
    return local::kml::compute_resolution(t, v);
  vw_throw(LogicErr() << "Asked to compute resolution for unknown "
           "quadtree type.");
}

void get_normalize_vals(boost::shared_ptr<DiskImageResource> file,
                        const Options& opt) {
  DiskImageView<PixelRGB<float>> min_max_file(file);
  float new_lo, new_hi;
  if (opt.nodata_set) {
    PixelRGB<float> no_data_value(opt.nodata);
    min_max_channel_values(create_mask(min_max_file, no_data_value),
                           new_lo, new_hi);
  } else if (file->has_nodata_read()) {
    PixelRGB<float> no_data_value(file->nodata_read());
    min_max_channel_values(create_mask(min_max_file, no_data_value),
                           new_lo, new_hi);
  } else {
    min_max_channel_values(min_max_file, new_lo, new_hi);
  }
  lo_value = std::min(new_lo, lo_value);
  hi_value = std::max(new_hi, hi_value);
  vw_out() << "Pixel range for \"" << file->filename() << ": ["
           << new_lo << " " << new_hi
           << "]    Output dynamic range: ["
           << lo_value << " " << hi_value << "]\n";
}

GeoReference make_input_georef(boost::shared_ptr<DiskImageResource> file,
                               const Options& opt) {
  GeoReference input_georef;
  bool fail_read_georef = false;
  try {
    fail_read_georef = !read_georeference(input_georef, *file);
  } catch (const InputErr& e) {
    vw_out(ErrorMessage) << "Input " << file->filename()
                         << " has malformed georeferencing information.\n";
    fail_read_georef = true;
  }

  if (opt.datum.type == "WGS84")
    input_georef.set_well_known_geogcs("WGS84");
  if (opt.datum.type == "LUNAR")
    input_georef.set_well_known_geogcs("D_MOON");
  if (opt.datum.type == "MARS")
    input_georef.set_well_known_geogcs("D_MARS");
  if (opt.datum.type == "SPHERE") {
    Datum datum("USER SUPPLIED DATUM", "SPHERICAL DATUM",
                "Reference Meridian",
                opt.datum.sphere_radius, opt.datum.sphere_radius, 0.0);
    input_georef.set_datum(datum);
  }

  if (opt.manual) {
    Matrix3x3 m;
    m(0, 0) = double(opt.east - opt.west) / file->cols();
    m(0, 2) = opt.west;
    m(1, 1) = double(opt.south - opt.north) / file->rows();
    m(1, 2) = opt.north;
    m(2, 2) = 1;
    input_georef.set_transform(m);
  } else if (fail_read_georef) {
    vw_out(ErrorMessage) << "Missing input georeference. Please provide "
                         << "--north --south --east and --west.\n";
    exit(1);
  }

  // If the user passed in georef parameters, process them here
  if (opt.proj.type == "LAMBERT_AZIMUTHAL")
    input_georef.set_lambert_azimuthal(opt.proj.lat, opt.proj.lon);
  if (opt.proj.type == "LAMBERT_CONFORMAL_CONIC")
    input_georef.set_lambert_conformal(opt.proj.p1, opt.proj.p2,
                                       opt.proj.lat, opt.proj.lon);
  if (opt.proj.type == "MERCATOR")
    input_georef.set_mercator(opt.proj.lat, opt.proj.lon, opt.proj.scale);
  if (opt.proj.type == "ORTHOGRAPHIC")
    input_georef.set_orthographic(opt.proj.lat, opt.proj.lon);
  if (opt.proj.type == "PLATE_CARREE")
    input_georef.set_geographic();
  if (opt.proj.type == "SINUSOIDAL")
    input_georef.set_sinusoidal(opt.proj.lon);
  if (opt.proj.type == "STEREOGRAPHIC")
    input_georef.set_stereographic(opt.proj.lat, opt.proj.lon,
                                    opt.proj.scale);
  if (opt.proj.type == "TRANSVERSE_MERCATOR")
    input_georef.set_transverse_mercator(opt.proj.lat, opt.proj.lon,
                                          opt.proj.scale);
  if (opt.proj.type == "UTM")
    input_georef.set_UTM(abs(opt.proj.utm_zone), opt.proj.utm_zone > 0);

  // Handle nudge arguments (x and y shifts in projected coordinates)
  if (opt.nudge_x || opt.nudge_y) {
    Matrix3x3 m = input_georef.transform();
    m(0, 2) += opt.nudge_x;
    m(1, 2) += opt.nudge_y;
    input_georef.set_transform(m);
  }

  return input_georef;
}

std::vector<GeoReference>
load_image_georeferences(const Options& opt, int& total_resolution) {
  std::vector<GeoReference> georeferences;
  georeferences.reserve(opt.input_files.size());

  BOOST_FOREACH(const std::string filename, opt.input_files) {
    boost::shared_ptr<DiskImageResource> file(DiskImageResourcePtr(filename));
    vw_out() << "Adding file " << file->filename() << "\n";

    if (opt.normalize) get_normalize_vals(file, opt);

    GeoReference input_georef = make_input_georef(file, opt);
    georeferences.push_back(input_georef);

    GeoReference output_georef(input_georef.datum());

    // Right now, we only need a WGS84 output geoereference to compute
    // the resolution. The rest of the output info will get set later.
    GeoTransform geotx(input_georef, output_georef);

    // Calculate the best resolution at 5 different points in the image,
    // as occasionally there's a singularity at the center pixel that
    // makes it extremely tiny (such as in pole-centered images).
    const int cols = file->cols();
    const int rows = file->rows();
    Vector2 res_pixel[5];
    res_pixel[0] = Vector2(cols / 2, rows / 2);
    res_pixel[1] = Vector2(cols / 2 + cols / 4, rows / 2);
    res_pixel[2] = Vector2(cols / 2 - cols / 4, rows / 2);
    res_pixel[3] = Vector2(cols / 2, rows / 2 + rows / 4);
    res_pixel[4] = Vector2(cols / 2, rows / 2 - rows / 4);
    int resolution = 0;
    for (int i = 0; i < 5; i++) {
      resolution = compute_resolution(opt.mode, geotx, res_pixel[i]);
      if (resolution > total_resolution) total_resolution = resolution;
    }
  }

  if (opt.global_resolution > 0) {
    vw_out(VerboseDebugMessage) << "Overriding calculated resolution "
      << total_resolution << " with " << opt.global_resolution << "\n";
    total_resolution = opt.global_resolution;
  }

  return georeferences;
}

void do_normal_mosaic(const Options& opt,
                      const ProgressCallback *progress) {
  DiskImageView<PixelT> img(opt.input_files[0]);
  QuadTreeGenerator quadtree(img, opt.output_file_name);
  quadtree.set_tile_size(256);
  quadtree.set_file_type("png");

  if (opt.mode != "NONE") {
    boost::shared_ptr<QuadTreeConfig> config =
      QuadTreeConfig::make(opt.mode);
    config->configure(quadtree);
  }

  vw_out() << "Generating overlay...\n";
  vw_out() << "Writing: " << opt.output_file_name << "\n";

  quadtree.generate(*progress);
}

void do_mosaic(const Options& opt, const ProgressCallback *progress) {

  // If we're not outputting any special sort of mosaic (just a regular old
  // quadtree, no georeferencing, no metadata), use a different function.
  if (opt.mode == "NONE" || opt.proj.type == "NONE") {
    do_normal_mosaic(opt, progress);
    return;
  }

  // Read in georeference info and compute total resolution.
  int total_resolution = 1024;
  std::vector<GeoReference> georeferences =
    load_image_georeferences(opt, total_resolution);

  boost::shared_ptr<QuadTreeConfig> config =
    QuadTreeConfig::make(opt.mode);

  // Now that we have the best resolution, we can get our output_georef.
  int xresolution = total_resolution / opt.aspect_ratio;
  int yresolution = total_resolution;

  GeoReference output_georef =
    config->output_georef(xresolution, yresolution);
  vw_out(VerboseDebugMessage, "tool") << "Output Georef:\n"
                                      << output_georef << "\n";

  // Configure the composite.
  ImageComposite<PixelT> composite;

  // Add the transformed image files to the composite.
  for (size_t i = 0; i < opt.input_files.size(); i++) {
    const std::string& filename = opt.input_files[i];
    const GeoReference& input_georef = georeferences[i];

    boost::shared_ptr<DiskImageResource> file(
      DiskImageResource::open(filename));
    GeoTransform geotx(input_georef, output_georef);

    ImageViewRef<PixelT> source = DiskImageView<PixelT>(file);

    // Handle nodata values/mask
    if (opt.nodata_set) {
      vw_out() << "Using nodata value: " << opt.nodata << "\n";
      source = mask_to_alpha(create_mask(
        pixel_cast<PixelWithoutAlpha<PixelT>::type>(source),
        ChannelT(opt.nodata)));
    } else if (file->has_nodata_read()) {
      vw_out() << "Using nodata value: " << file->nodata_read() << "\n";
      source = mask_to_alpha(create_mask(
        pixel_cast<PixelWithoutAlpha<PixelT>::type>(source),
        ChannelT(file->nodata_read())));
    }

    bool is_global =
      ((boost::trim_copy(input_georef.proj4_str()) == "+proj=longlat") &&
       (fabs(input_georef.lonlat_to_pixel(Vector2(-180, 0)).x()) < 1) &&
       (fabs(input_georef.lonlat_to_pixel(Vector2(180, 0)).x() -
             source.cols()) < 1) &&
       (fabs(input_georef.lonlat_to_pixel(Vector2(0, 90)).y()) < 1) &&
       (fabs(input_georef.lonlat_to_pixel(Vector2(0, -90)).y() -
             source.rows()) < 1));

    // Do various modifications to the input image here.
    if ((opt.pixel_scale != 1.0) || (opt.pixel_offset != 0.0)) {
      vw_out() << "Apply input scaling: " << opt.pixel_scale
               << " offset: " << opt.pixel_offset << "\n";
      source = channel_cast<ChannelT>(
        source * opt.pixel_scale + opt.pixel_offset);
    }

    // Normalize pixel intensity if desired
    if (opt.normalize) {
      vw_out() << "Apply normalizing: [" << lo_value << ", "
               << hi_value << "]\n";
      typedef ChannelRange<ChannelT> range_type;
      source = normalize_retain_alpha(source, lo_value, hi_value,
                                      range_type::min(),
                                      range_type::max());
    }

    BBox2i bbox =
      geotx.forward_bbox(BBox2i(0, 0, source.cols(), source.rows()));
    if (is_global) {
      vw_out() << "\t--> Detected global overlay. Using cylindrical "
               << "edge extension to hide the seam.\n";
      source = crop(transform(source, geotx, source.cols(),
                              source.rows(),
                              CylindricalEdgeExtension()), bbox);
    } else {
      if (norm_2(geotx.reverse(geotx.forward(Vector2()))) >
          0.01 * norm_2(Vector2(source.cols(), source.rows()))) {
        BBox2i correction(
          -geotx.reverse(geotx.forward(Vector2()))[0], 0,
          source.cols(), source.rows());
        source = crop(transform_only(
          crop(interpolate(source), correction), geotx), bbox);
      } else {
        source = transform(source, geotx, bbox);
      }
    }

    // Images that wrap the date line must be added to the composite
    // on both sides.
    if (bbox.max().x() > total_resolution)
      composite.insert(source, bbox.min().x() - total_resolution,
                       bbox.min().y());
    // Images that are in the 180-360 range *only* go on the other side.
    if (bbox.min().x() < xresolution)
      composite.insert(source, bbox.min().x(), bbox.min().y());
  } // End loop through input files

  // This box represents the entire input data set, in pixels, in the
  // output projection space.
  BBox2i total_bbox = composite.bbox();
  total_bbox.crop(BBox2i(0, 0, xresolution, yresolution));

  VW_ASSERT(total_bbox.width() > 0 && total_bbox.height() > 0,
            LogicErr() << "Total bbox is empty. Georeference "
            "calculation is probably incorrect.");

  if (opt.mode == "KML") {
    BBox2i bbox = total_bbox;
    // Compute a tighter Google Earth coordinate system aligned bbox.
    int dim = 2 << (int)(log((double)(std::max)(bbox.width(),
                                                bbox.height())) /
                         log(2.));
    if (dim > total_resolution)
      dim = total_resolution;
    total_bbox = BBox2i((bbox.min().x() / dim) * dim,
                        (bbox.min().y() / dim) * dim, dim, dim);
    if (!total_bbox.contains(bbox)) {
      if (total_bbox.max().x() == xresolution)
        total_bbox.min().x() -= dim;
      else
        total_bbox.max().x() += dim;
      if (total_bbox.max().y() == yresolution)
        total_bbox.min().y() -= dim;
      else
        total_bbox.max().y() += dim;
    }
  }

  // Prepare the composite.
  if (!opt.multiband)
    composite.set_draft_mode(true);
  composite.prepare(total_bbox, *progress);
  VW_ASSERT(composite.rows() > 0 && composite.cols() > 0,
            LogicErr() << "Composite image is empty. Georeference "
            "calculation is probably incorrect.");

  QuadTreeGenerator quadtree(composite, opt.output_file_name);

  if (opt.mode == "KML") {
    KMLQuadTreeConfig *c2 =
      dynamic_cast<KMLQuadTreeConfig*>(config.get());
    BBox2 ll_bbox(
      -180.0 + (360.0 * total_bbox.min().x()) / xresolution,
       180.0 - (360.0 * total_bbox.max().y()) / yresolution,
                (360.0 * total_bbox.width())  / xresolution,
                (360.0 * total_bbox.height()) / yresolution);
    c2->set_longlat_bbox(ll_bbox);
    c2->set_max_lod_pixels(opt.kml.max_lod_pixels);
    c2->set_draw_order_offset(opt.kml.draw_order_offset);
  }

  config->configure(quadtree);

  if (opt.tile_size > 0)
    quadtree.set_tile_size(opt.tile_size);
  if (!opt.output_file_type.empty())
    quadtree.set_file_type(opt.output_file_type);

  BBox2i data_bbox = composite.bbox();
  data_bbox.crop(BBox2i(0, 0, total_bbox.width(), total_bbox.height()));
  quadtree.set_crop_bbox(data_bbox);

  vw_out() << "Generating overlay...\n";
  vw_out() << "Writing: " << opt.output_file_name << "\n";

  quadtree.generate(*progress);
}

int handle_options(int argc, char *argv[], Options& opt) {
  po::options_description general_options(
    "Description: Turns georeferenced image(s) into a quadtree "
    "with geographical metadata\n\nGeneral options");
  general_options.add_options()
    ("output-name,o", po::value(&opt.output_file_name),
     "Specify the base output directory");

  general_options.add(GdalWriteOptionsDescription(opt));

  const std::string mode_options_str = "NONE, KML";
  const std::string datum_options_str = "NONE, WGS84, LUNAR, MARS, SPHERE";
  const std::string projection_options_str =
    "DEFAULT, NONE, SINUSOIDAL, MERCATOR, TRANSVERSE_MERCATOR, "
    "ORTHOGRAPHIC, STEREOGRAPHIC, LAMBERT_AZIMUTHAL, "
    "LAMBERT_CONFORMAL_CONIC, UTM, PLATE_CARREE";

  po::options_description input_options("Input options");
  string datum_desc =
    string("Override input datum [") + datum_options_str + "]";
  string mode_desc =
    string("Specify the output metadata type [") + mode_options_str + "]";
  string proj_desc =
    string("Projection type [") + projection_options_str + "]";
  input_options.add_options()
    ("force-datum", po::value(&opt.datum.type), datum_desc.c_str())
    ("datum-radius", po::value(&opt.datum.sphere_radius),
     "Radius to use for --force-datum SPHERE")
    ("pixel-scale", po::value(&opt.pixel_scale)->default_value(1.0),
     "Scale factor to apply to pixels")
    ("pixel-offset", po::value(&opt.pixel_offset)->default_value(0.0),
     "Offset to apply to pixels")
    ("normalize", po::bool_switch(&opt.normalize),
     "Normalize input images so that their full dynamic range "
     "falls in between [0,255].")
    ("nodata", po::value(&opt.nodata),
     "Set the input's nodata value so that it will be transparent "
     "in output");

  po::options_description output_options("Output options");
  output_options.add_options()
    ("mode,m", po::value(&opt.mode)->default_value("KML"),
     mode_desc.c_str())
    ("file-type", po::value(&opt.output_file_type),
     "Output file type. Choose 'auto' to generate jpgs in opaque "
     "areas and png images where there is transparency.")
    ("module-name", po::value(&opt.module_name),
     "The module where the output will be placed.")
    ("jpeg-quality", po::value(&opt.jpeg_quality),
     "JPEG quality factor (0.0 to 1.0)")
    ("png-compression", po::value(&opt.png_compression),
     "PNG compression level (0 to 9)")
    ("tile-size", po::value(&opt.tile_size), "Tile size in pixels")
    ("max-lod-pixels",
     po::value(&opt.kml.max_lod_pixels)->default_value(1024),
     "Max LoD in pixels, or -1 for none (kml only)")
    ("draw-order-offset",
     po::value(&opt.kml.draw_order_offset)->default_value(0),
     "Offset for the <drawOrder> tag for this overlay (kml only)")
    ("multiband", po::bool_switch(&opt.multiband),
     "Composite images using multi-band blending")
    ("aspect-ratio", po::value(&opt.aspect_ratio),
     "Pixel aspect ratio (for polar overlays; should be a power of two)")
    ("global-resolution", po::value(&opt.global_resolution),
     "Override the global pixel resolution; should be a power of two");

  po::options_description projection_options("Input projection options");
  projection_options.add_options()
    ("north", po::value(&opt.north),
     "The northernmost latitude in projection units")
    ("south", po::value(&opt.south),
     "The southernmost latitude in projection units")
    ("east", po::value(&opt.east),
     "The easternmost longitude in projection units")
    ("west", po::value(&opt.west),
     "The westernmost longitude in projection units")
    ("global", po::bool_switch(&opt.global),
     "Override image size to global (in lonlat)")
    ("projection",
     po::value(&opt.proj.type)->default_value("DEFAULT"),
     proj_desc.c_str())
    ("utm-zone", po::value(&opt.proj.utm_zone),
     "Set zone for --projection UTM (negative for south)")
    ("proj-lat", po::value(&opt.proj.lat),
     "The center of projection latitude")
    ("proj-lon", po::value(&opt.proj.lon),
     "The center of projection longitude")
    ("proj-scale", po::value(&opt.proj.scale),
     "The projection scale")
    ("p1", po::value(&opt.proj.p1),
     "parallel for Lambert Conformal Conic projection")
    ("p2", po::value(&opt.proj.p2),
     "parallel for Lambert Conformal Conic projection")
    ("nudge-x", po::value(&opt.nudge_x),
     "Nudge the image, in projected coordinates")
    ("nudge-y", po::value(&opt.nudge_y),
     "Nudge the image, in projected coordinates");

  po::options_description hidden_options("");
  hidden_options.add_options()
    ("input-file", po::value(&opt.input_files));

  po::options_description options("Allowed options");
  options.add(general_options).add(input_options)
    .add(projection_options).add(output_options).add(hidden_options);

  po::positional_options_description p;
  p.add("input-file", -1);

  std::ostringstream usage;
  usage << "Usage: image2qtree [options] <filename>..." << endl << endl;
  usage << general_options    << endl;
  usage << input_options      << endl;
  usage << output_options     << endl;
  usage << projection_options << endl;

  try {
    namespace ps = po::command_line_style;
    int style = ps::unix_style & ~ps::allow_guessing;
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv)
              .style(style).options(options).positional(p).run(), vm);
    po::notify(vm);

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
    if (strlen(msg) > 0)
      cerr << endl << "Invalid argument:" << endl
           << "\t" << msg << endl;
    return false;
  }

  opt.setVwSettingsFromOpt();

  // Make sure all these string parameters are upper case
  boost::to_upper(opt.mode);
  boost::to_upper(opt.datum.type);
  boost::to_upper(opt.proj.type);
  return true;
}

int main(int argc, char **argv) {
  Options opt;
  if (!handle_options(argc, argv, opt))
    return 1;

  try {
    TerminalProgressCallback tpc("image2qtree", "");
    do_mosaic(opt, &tpc);
  } catch (const ArgumentErr& e) {
    vw_out() << e.what() << "\n";
    return 1;
  } catch (const Exception& e) {
    std::cerr << "\n\nVW Error: " << e.what() << "\n";
    return 1;
  } catch (const std::bad_alloc& e) {
    std::cerr << "\n\nError: Ran out of Memory!\n";
    return 1;
  } catch (const std::exception& e) {
    std::cerr << "\n\nError: " << e.what() << "\n";
    return 1;
  }

  return 0;
}
