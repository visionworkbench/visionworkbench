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


/// \file image2plate.cc
///
/// This utility adds a georefernced image to a plate file.
///
/// By default, this tool will add an image to the mosaic, along with
/// mipmapped tiles for that image, to its own layer.  Layers are
/// indexed by transaction IDs.
///

#include <vw/Plate/PlateFile.h>
#include <vw/Plate/PlateCarreePlateManager.h>
#include <vw/Plate/PolarStereoPlateManager.h>
#include <vw/Plate/ToastPlateManager.h>
#include <vw/Cartography/GeoReference.h>
#include <vw/FileIO.h>
#include <vw/Image/MaskViews.h>
#include <vw/Image/Filter.h>

using namespace vw;
using namespace vw::platefile;
using namespace vw::cartography;

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>
namespace fs = boost::filesystem;

using boost::optional;

VW_DEFINE_EXCEPTION(Usage, Exception);

template <typename T>
std::istream& operator>>(std::istream& i, optional<T>& val) {
  T x;
  i >> x;
  val = x;
  return i;
}

struct Options {
  optional<Url> url;
  optional<TransactionOrNeg> transaction_id;
  std::string filetype;
  optional<double> nodata_value;
  std::string mode;
  size_t tile_size;
  optional<float> jpeg_quality;
  optional<unsigned> png_compression;
  size_t cache_size;
  bool terrain;
  double nudge_x;
  double nudge_y;
  bool force_lunar;
  bool force_mars;
  optional<double> force_spherical;
  bool force_float;
  bool debug;
  bool help;
  std::vector<std::string> image_files;
  optional<float> north, south, east, west;
  bool manual;
  bool global;

  Options() :
    filetype("png"),
    tile_size(256),
    cache_size(0),
    terrain(false),
    nudge_x(0),
    nudge_y(0),
    force_lunar(false),
    force_mars(false),
    force_float(false),
    debug(false),
    help(false),
    manual(false),
    global(false)
    {}


  void validate() {
    VW_ASSERT(!help, Usage());
    VW_ASSERT(image_files.size() > 0, Usage() << "Need at least one input image");
    VW_ASSERT(url, Usage() << "Output url required");
    if (global || north || south || east || west) {
      VW_ASSERT(image_files.size() == 1,
          Usage() << "Cannot override georeference information on multiple images");
      VW_ASSERT(global || (north && south && east && west),
          Usage() << "If you provide one, you must provide all of: --north --south --east --west");
      if (global) {
        north = 90; south = -90; east = 180; west = -180;
      }
      manual = true;
    }
    if (transaction_id) {
      if (transaction_id.get() < 1)
        vw_throw(Usage() << "Transaction ID must be > 0 (got " << transaction_id.get() << ")");
    } else
      transaction_id = -1;

    VW_ASSERT(bool(force_lunar) + bool(force_mars) + bool(force_spherical) < 2,
              Usage() << "Cannot force more than one of: lunar, mars, spherical");

    if (transaction_id && image_files.size() != 1)
      vw_throw(Usage() << "You cannot override the transaction-id while processing multiple images");

    if (transaction_id && transaction_id.get() < 1u)
      vw_throw(Usage() << "you must specify a positive transaction-id.");

    VW_ASSERT(mode == "toast" || mode == "equi" || mode == "polar",
        Usage() << "Unknown mode: " << mode);

    if (cache_size)
      vw_settings().set_system_cache_size( cache_size*1024*1024 );

    if (jpeg_quality)
      DiskImageResourceJPEG::set_default_quality( jpeg_quality.get() );
    if (png_compression)
      DiskImageResourcePNG::set_default_compression_level( png_compression.get() );
  }
};

// --------------------------------------------------------------------------
//                                DO_MOSAIC
// --------------------------------------------------------------------------

template <class ViewT>
void do_mosaic(boost::shared_ptr<PlateFile> platefile,
               ImageViewBase<ViewT> const& view,
               std::string filename, GeoReference const& georef, const Options& opt) {
  typedef typename ViewT::pixel_type PixelT;

  typedef PlateManager<PixelT> PM;

  boost::scoped_ptr<PM> pm(PM::make(opt.mode, platefile));

  pm->insert(view.impl(), filename, opt.transaction_id.get(), georef,
             opt.terrain, opt.debug, TerminalProgressCallback( "plate.tools.image2plate", "\t    Processing") );
}

// --------------------------------------------------------------------------
//                                    MAIN
// --------------------------------------------------------------------------

bool handle_options( int argc, char *argv[], Options& opt ) {
  po::options_description general_options("Insert image into a platefile.\n\nGeneral Options");
  general_options.add_options()
    ("output-name,o",         po::value(&opt.url),               "Specify the URL of the platefile.")
    ("transaction-id,t",      po::value(&opt.transaction_id),    "Specify the transaction_id to use for this transaction. If you don't specify one, one will be automatically assigned.\n")
    ("file-type",             po::value(&opt.filetype),          "Output file type")
    ("nodata-value",          po::value(&opt.nodata_value),      "Explicitly set the value to treat as na data (i.e. transparent) in the input file.")
    ("mode,m",                po::value(&opt.mode),              "Output mode [toast, equi, polar]")
    ("tile-size",             po::value(&opt.tile_size),         "Tile size, in pixels")
    ("jpeg-quality",          po::value(&opt.jpeg_quality),      "JPEG quality factor (0.0 to 1.0)")
    ("png-compression",       po::value(&opt.png_compression),   "PNG compression level (0 to 9)")
    ("cache",                 po::value(&opt.cache_size),        "Source data cache size, in megabytes")
    ("terrain",               po::bool_switch(&opt.terrain),     "Tweak a few settings that are best for terrain platefiles. Turns on nearest neighbor sampling in mipmapping and zero out semi-transparent pixels.")
    ("nudge-x",               po::value(&opt.nudge_x),           "Nudge the image, in projected coordinates")
    ("nudge-y",               po::value(&opt.nudge_y),           "Nudge the image, in projected coordinates")
    ("force-lunar-datum",     po::bool_switch(&opt.force_lunar), "Use the lunar spherical datum for the input images' geographic coordinate systems, even if they are not encoded to do so.")
    ("force-mars-datum",      po::bool_switch(&opt.force_mars),  "Use the Mars spherical datum for the input images' geographic coordinate systems, even if they are not encoded to do so.")
    ("force-spherical-datum", po::value(&opt.force_spherical),   "Choose an arbitrary input spherical datum to use for input images', overriding the existing datum.")
    ("force-float",           po::bool_switch(&opt.force_float), "Force the platefile to use a channel type of float.")
    ("north",                 po::value(&opt.north),             "The northernmost latitude in projection units")
    ("south",                 po::value(&opt.south),             "The southernmost latitude in projection units")
    ("east",                  po::value(&opt.east),              "The easternmost longitude in projection units")
    ("west",                  po::value(&opt.west),              "The westernmost longitude in projection units")
    ("global",                po::bool_switch(&opt.global),      "Override image size to global (in lonlat)")
    ("debug",                 po::bool_switch(&opt.debug),       "Display helpful debugging messages.")
    ("help,h",                po::bool_switch(&opt.help),        "Display this help message");

  po::options_description hidden_options("");
  hidden_options.add_options()
    ("input-file", po::value(&opt.image_files));

  po::options_description options("Allowed Options");
  options.add(general_options).add(hidden_options);

  po::positional_options_description p;
  p.add("input-file", -1);

  std::ostringstream usage;
  usage << "Usage: " << argv[0] << " [options] <filename>..." <<std::endl << std::endl;
  usage << general_options << std::endl;

  try {
    po::variables_map vm;
    po::store( po::command_line_parser( argc, argv ).options(options).positional(p).run(), vm );
    po::notify( vm );
    opt.validate();
  } catch (const po::error& e) {
    std::cerr << usage.str() << std::endl
              << "Failed to parse command line arguments:" << std::endl
              << "\t" << e.what() << std::endl;
    return false;
  } catch (const Usage& e) {
    const char* msg = e.what();
    std::cerr << usage.str() << std::endl;
    if (::strlen(msg) > 0)
      std::cerr << std::endl << "Invalid argument:" << std::endl << "\t" << msg << std::endl;
    return false;
  }
  return true;
}

void run(const Options& opt) {
  boost::shared_ptr<PlateFile> platefile;
  {
    std::string filetype = opt.filetype;

    // For new plate files, we adopt the pixel and channel type of the
    // first image.  If the platefile already exists, then this is
    // ignored and the native pixel and channel type of the platefile is
    // used instead.
    boost::scoped_ptr<DiskImageResource> rsrc(DiskImageResource::open(opt.image_files[0]));
    PixelFormatEnum pixel_format = rsrc->pixel_format();
    ChannelTypeEnum channel_type = rsrc->channel_type();

    // Plate files should always have an alpha channel.
    if (pixel_format == VW_PIXEL_GRAY)
      pixel_format = VW_PIXEL_GRAYA;
    if (pixel_format == VW_PIXEL_RGB)
      pixel_format = VW_PIXEL_RGBA;

    // User override for channel type
    if (opt.force_float) {
      vw_out() << "\t--> Processing image using a 32-bit floating point channel type.\n";
      if (boost::to_lower_copy(opt.filetype) != "exr" &&
          boost::to_lower_copy(opt.filetype) != "tif" ) {
        vw_out() << "\t    Overriding channel type for 32-bit float: setting it to OpenEXR.\n";
        filetype = "exr";
      }
      channel_type = VW_CHANNEL_FLOAT32;
    }

    vw_out() << "\nOpening plate file: " << opt.url.get() << " [" << filetype << "]\n";
    platefile.reset( new PlateFile(opt.url.get(), opt.mode, "", opt.tile_size, filetype, pixel_format, channel_type) );
  }

  BOOST_FOREACH(const std::string& filename, opt.image_files) {
    VW_ASSERT(fs::exists(filename), ArgumentErr() << "No such file: " << filename);

    vw_out() << "\t--> Building full-resolution tiles for " << filename << "\n";

    // Load the pixel type, channel type, and nodata value, and
    // georeferencing info for this image.
    boost::scoped_ptr<DiskImageResource> rsrc( DiskImageResource::open(filename) );

    optional<double> nodata_value;
    if (opt.nodata_value) {
      nodata_value = opt.nodata_value.get();
      vw_out() << "\t--> Using user-supplied nodata value: " << nodata_value << ".\n";
    } else if ( rsrc->has_nodata_read() ) {
      nodata_value = rsrc->nodata_read();
      vw_out() << "\t--> Extracted nodata value from file: " << nodata_value << ".\n";
    }

    // Load the georef.  If none is found, assume Plate Caree.
    GeoReference georef;
    {
      DiskImageResourceGDAL diskrsrc( filename );

      bool fail_read_georef = false;
      try {
        fail_read_georef = !read_georeference( georef, diskrsrc );
      } catch ( InputErr const& e ) {
        VW_OUT(ErrorMessage) << "Input " << diskrsrc.filename() << " has malformed georeferencing information.\n";
        fail_read_georef = true;
      }

      if(opt.force_lunar) {
        const double LUNAR_RADIUS = 1737400;
        VW_OUT() << "\t--> Using standard lunar spherical datum: " << LUNAR_RADIUS << "\n";
        cartography::Datum datum("D_MOON", "MOON", "Reference Meridian", LUNAR_RADIUS, LUNAR_RADIUS, 0.0);
        georef.set_datum(datum);
      } else if(opt.force_mars) {
        const double MOLA_PEDR_EQUATORIAL_RADIUS = 3396000.0;
        VW_OUT() << "\t--> Using standard MOLA spherical datum: " << MOLA_PEDR_EQUATORIAL_RADIUS << "\n";
        cartography::Datum datum("D_MARS", "MARS", "Reference Meridian", MOLA_PEDR_EQUATORIAL_RADIUS, MOLA_PEDR_EQUATORIAL_RADIUS, 0.0);
        georef.set_datum(datum);
      } else if(opt.force_spherical) {
        VW_OUT() << "\t--> Using user-supplied spherical datum: " << opt.force_spherical << "\n";
        cartography::Datum datum("USER SUPPLIED DATUM", "SPHERICAL DATUM", "Reference Meridian", opt.force_spherical.get(), opt.force_spherical.get(), 0.0);
        georef.set_datum(datum);
      }

      if( opt.manual ) {
        Matrix3x3 m;
        m(0,0) = double(opt.east.get() - opt.west.get()) / diskrsrc.cols();
        m(0,2) = opt.west.get();
        m(1,1) = double(opt.south.get() - opt.north.get()) / diskrsrc.rows();
        m(1,2) = opt.north.get();
        m(2,2) = 1;
        georef.set_transform( m );
      } else if ( fail_read_georef )
        vw_throw(IOErr() << "Image " << diskrsrc.filename() << " missing input georeference. Please provide --north --south --east and --west.");
    }

    // Apply nudge factors
    if( opt.nudge_x || opt.nudge_y ) {
      Matrix3x3 m = georef.transform();
      m(0,2) += opt.nudge_x;
      m(1,2) += opt.nudge_y;
      georef.set_transform( m );
    }

    PixelFormatEnum rsrc_pixel_frmt = rsrc->pixel_format();
    ChannelTypeEnum rsrc_channel_type = rsrc->channel_type();

    // User override for channel type
    if (opt.force_float) {
      vw_out() << "\t--> Forcing floating point tiles, and disabling image autoscaling.\n";
      rsrc_channel_type = VW_CHANNEL_FLOAT32;

      // Turn off automatic rescaling when reading in images
      DiskImageResource::set_default_rescale(false);
    }

    // Dispatch to the compositer based on the pixel type of this mosaic.
    switch(rsrc_pixel_frmt) {
      case VW_PIXEL_GRAY:
      case VW_PIXEL_GRAYA:
        switch(rsrc_channel_type) {
          case VW_CHANNEL_UINT8:
          case VW_CHANNEL_UINT16:
            if (nodata_value)
              do_mosaic(platefile,
                  mask_to_alpha(create_mask(DiskImageView<PixelGray<uint8> >(filename), boost::numeric_cast<uint8>(nodata_value.get()))),
                  filename, georef, opt);
            else
              do_mosaic(platefile, DiskImageView<PixelGrayA<uint8> >(filename),
                  filename, georef, opt);
            break;
          case VW_CHANNEL_INT16:
            if (nodata_value)
              do_mosaic(platefile,
                  mask_to_alpha(create_mask(DiskImageView<PixelGray<int16> >(filename), boost::numeric_cast<int16>(nodata_value.get()))),
                  filename, georef, opt);
            else
              do_mosaic(platefile, DiskImageView<PixelGrayA<int16> >(filename),
                  filename, georef, opt);
            break;
          case VW_CHANNEL_FLOAT32:
            if (nodata_value) {
              do_mosaic(platefile,
                  mask_to_alpha(create_mask(DiskImageView<PixelGray<float32> >(filename), boost::numeric_cast<float32>(nodata_value.get()))),
                  filename, georef, opt);
            } else
              do_mosaic(platefile, DiskImageView<PixelGrayA<float32> >(filename),
                  filename, georef, opt);
            break;
          default:
            vw_throw(NoImplErr() << "Image contains a channel type not supported by image2plate.");
        }
        break;

      case VW_PIXEL_RGB:
      case VW_PIXEL_RGBA:
        switch(rsrc->channel_type()) {
          case VW_CHANNEL_UINT8:
            if (nodata_value)
              do_mosaic(platefile,
                  pixel_cast<PixelRGBA<uint8> >(mask_to_alpha(create_mask(DiskImageView<PixelGray<uint8> >(filename), boost::numeric_cast<uint8>(nodata_value.get())))),
                  filename, georef, opt);
            else
              do_mosaic(platefile, DiskImageView<PixelRGBA<uint8> >(filename),
                  filename, georef, opt);
            break;
          default:
            vw_throw(NoImplErr() << "Image contains a channel type not supported by image2plate.");
        }
        break;
      default:
        vw_throw(NoImplErr() << "Image contains a channel type not supported by image2plate.");
    }
  }
}

int main(int argc, char **argv) {
  Options opt;
  if (!handle_options(argc, argv, opt))
    return 1;
  try {
    run(opt);
    return 0;
  } catch (const PlatefileErr& e) {
    VW_OUT(ErrorMessage) << "A platefile error occured: " << e.what() << "\n";
  } catch (const Exception& e) {
    VW_OUT(ErrorMessage) << "Runtime Error: " << e.what() << "\n";
  }
  return 1;
}
