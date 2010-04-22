// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file image2plate.cc
///
/// This utility adds a georefernced image to a plate file.
///
/// By default, this tool will add an image to the mosaic, along with
/// mipmapped tiles for that image, to its own layer.  Layers are
/// indexed by transaction IDs.
///

#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Cartography.h>
#include <vw/Plate/PlateFile.h>
#include <vw/Plate/ToastPlateManager.h>
#include <vw/Plate/PlateCarreePlateManager.h>

using namespace vw;
using namespace vw::platefile;
using namespace vw::mosaic;
using namespace vw::cartography;

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>
namespace fs = boost::filesystem;

// Global variables
bool g_debug;

// Erases a file suffix if one exists and returns the base string
static std::string prefix_from_filename(std::string const& filename) {
  std::string result = filename;
  int index = result.rfind(".");
  if (index != -1) 
    result.erase(index, result.size());
  return result;
}

// --------------------------------------------------------------------------
//                                DO_MOSAIC
// --------------------------------------------------------------------------

template <class ViewT>
void do_mosaic(boost::shared_ptr<PlateFile> platefile, 
               ImageViewBase<ViewT> const& view,
               std::string filename, int transaction_id_override, 
               GeoReference const& georef, std::string output_mode) {

  std::ostringstream status_str;
  status_str << "\t    " << filename << " : ";

  if (output_mode == "toast") {

    boost::shared_ptr<ToastPlateManager<typename ViewT::pixel_type> > pm( 
      new ToastPlateManager<typename ViewT::pixel_type> (platefile) );

    pm->insert(view.impl(), filename, transaction_id_override, georef, g_debug,
               TerminalProgressCallback( "plate.tools.image2plate", "\t    Processing") );

  }  else if (output_mode == "equi")  {

    boost::shared_ptr<PlateCarreePlateManager<typename ViewT::pixel_type> > pm(
      new PlateCarreePlateManager<typename ViewT::pixel_type> (platefile) );

    pm->insert(view.impl(), filename, transaction_id_override, georef, g_debug,
               TerminalProgressCallback( "plate.tools.image2plate", "\t    Processing") );
    
  }

}

// --------------------------------------------------------------------------
//                                    MAIN
// --------------------------------------------------------------------------

int main( int argc, char *argv[] ) {
  std::string url;
  std::string tile_filetype;
  std::string output_mode;
  int tile_size;
  int transaction_id_override = -1;
  float jpeg_quality;
  int png_compression;
  unsigned cache_size;
  std::vector<std::string> image_files;

  po::options_description general_options("Turns georeferenced image(s) into a TOAST quadtree.\n\nGeneral Options");
  general_options.add_options()
    ("output-name,o", po::value<std::string>(&url), "Specify the URL of the platefile.")
    ("transaction-id,t", po::value<int>(&transaction_id_override), "Specify the transaction_id to use for this transaction. If you don't specify one, one will be automatically assigned.\n")
    ("file-type", po::value<std::string>(&tile_filetype)->default_value("png"), "Output file type")
    ("mode,m", po::value<std::string>(&output_mode)->default_value("toast"), "Output mode [toast, equi]")
    ("tile-size", po::value<int>(&tile_size)->default_value(256), "Tile size, in pixels")
    ("jpeg-quality", po::value<float>(&jpeg_quality)->default_value(0.95), "JPEG quality factor (0.0 to 1.0)")
    ("png-compression", po::value<int>(&png_compression)->default_value(3), "PNG compression level (0 to 9)")
    ("cache", po::value<unsigned>(&cache_size)->default_value(512), "Source data cache size, in megabytes")
    ("debug", "Display helpful debugging messages.")
    ("help", "Display this help message");

  po::options_description hidden_options("");
  hidden_options.add_options()
    ("input-file", po::value<std::vector<std::string> >(&image_files));

  po::options_description options("Allowed Options");
  options.add(general_options).add(hidden_options);

  po::positional_options_description p;
  p.add("input-file", -1);

  std::ostringstream usage;
  usage << "Usage: " << argv[0] << " [options] <filename>..." <<std::endl << std::endl;
  usage << general_options << std::endl;

  po::variables_map vm;
  try { 
    po::store( po::command_line_parser( argc, argv ).options(options).positional(p).run(), vm );
    po::notify( vm );
  } catch (po::error &e) {
    std::cout << "An error occured while parsing command line arguments.\n\n";
    std::cout << usage.str();
    return 0;    
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

  //------------------------- SET DEFAULT OPTIONS -----------------------------

  if( url == "" )
    url = prefix_from_filename(image_files[0]) + "_" + output_mode + ".plate";

  if( tile_size <= 0 || tile_size != pow(2,int(log(tile_size)/log(2))) ) {
    std::cerr << "Error: The tile size must be a power of two!  (You specified: " << tile_size << ")" << std::endl << std::endl;
    std::cout << usage.str();
    return 1;
  }

  DiskImageResourceJPEG::set_default_quality( jpeg_quality );
  DiskImageResourcePNG::set_default_compression_level( png_compression );
  vw_settings().set_system_cache_size( cache_size*1024*1024 );

  //-------------------------- OPEN THE PLATE FILE ------------------------------

  // For new plate files, we adopt the pixel and channel type of the
  // first image.  If the platefile already exists, then this is
  // ignored and the native pixel and channel type of the platefile is
  // used instead.
  PixelFormatEnum pixel_format;
  ChannelTypeEnum channel_type;
  try {

    DiskImageResource *rsrc = DiskImageResource::open(image_files[0]);
    pixel_format = rsrc->pixel_format();
    channel_type = rsrc->channel_type();

    // Plate files should always have an alpha channel.
    if (pixel_format == VW_PIXEL_GRAY)
      pixel_format = VW_PIXEL_GRAYA;
    if (pixel_format == VW_PIXEL_RGB)
      pixel_format = VW_PIXEL_RGBA;
    delete rsrc;

  } catch (vw::Exception &e) {
    vw_out(ErrorMessage) << "An error occured: " << e.what() << "\n";
    exit(1);
  }

  // Set the debug level
  if (vm.count("debug")) {
    g_debug = true;
  } else {
    g_debug = false;
  }
  
  if (vm.count("transaction-id") && image_files.size() != 1) {
    std::cout << "Error: you cannot override the transaction-id while processing multiple images with " << argv[0] << ".\n";
    exit(1);
  }

  if (vm.count("transaction-id") && transaction_id_override < 1) {
    vw_out() << "Error: you must specify a positive transaction-id.\n";
    exit(1);
  }

  // Create both platefile managers (we only end up using one... this
  // just makes the code a little more readable.)
  if (output_mode != "toast" && output_mode != "equi") {
    vw_out() << "Unknown mode type passed in using --mode: " << output_mode << ".  Exiting.\n";
    exit(1);
  }
  try {
      
    std::cout << "\nOpening plate file: " << url << "\n";
    boost::shared_ptr<PlateFile> platefile = 
      boost::shared_ptr<PlateFile>( new PlateFile(url, output_mode, "",
                                                  tile_size, tile_filetype,
                                                  pixel_format, channel_type));

    // Process each image individually
    for ( unsigned i = 0; i < image_files.size(); ++i ) {

      // Check to see if the image exists.
      if ( !fs::exists(image_files[i]) ) {
        vw_out() << "Error: could not open image file named \"" << image_files[i] << "\"";
        exit(1);
      }

      std::cout << "\t--> Building full-resolution tiles for " << image_files[i] << "\n";

      // Load the pixel type, channel type, and nodata value, and
      // georeferencing info for this image.
      boost::shared_ptr<DiskImageResource> rsrc( DiskImageResource::open(image_files[i]) );
      
      double nodata_value = 0;
      bool has_nodata_value = false;
      if ( rsrc->has_nodata_value() ) {
        has_nodata_value = true;
        nodata_value = rsrc->nodata_value();
        std::cout << "\t--> Extracted nodata value from file: " << nodata_value << ".\n";
      }

      // Load the georef.  If none is found, assume Plate Caree.
      GeoReference georef;
      read_georeference( georef, DiskImageResourceGDAL( image_files[i] ) );
      if( georef.transform() == identity_matrix<3>() ) {
        std::cout << "\t    No georeferencing info found for " << image_files[i]
                  << ".  Assuming global plate carree." << std::endl;
        Matrix3x3 M;
        M(0,0) = 360.0 / rsrc->cols();
        M(0,2) = -180.0;
        M(1,1) = -180.0 / rsrc->rows();
        M(1,2) = 90.0;
        M(2,2) = 1;
        georef.set_transform( M );
      }

      // Dispatch to the compositer based on the pixel type of this mosaic.
      switch(rsrc->pixel_format()) {
      case VW_PIXEL_GRAY:
      case VW_PIXEL_GRAYA:
        switch(rsrc->channel_type()) {
        case VW_CHANNEL_UINT8:  
        case VW_CHANNEL_UINT16:  
          if (has_nodata_value)
            do_mosaic(platefile, 
                      mask_to_alpha(create_mask(DiskImageView<PixelGray<uint8> >(image_files[i]), 
                                                nodata_value)),
                      image_files[i], transaction_id_override, georef, output_mode);
          else
            do_mosaic(platefile, DiskImageView<PixelGrayA<uint8> >(image_files[i]), 
                      image_files[i], transaction_id_override, georef, output_mode);
          break;
        case VW_CHANNEL_INT16:  
          if (has_nodata_value)
            do_mosaic(platefile,
                      mask_to_alpha(create_mask(DiskImageView<PixelGray<int16> >(image_files[i]), 
                                                nodata_value)),
                      image_files[i], transaction_id_override, georef, output_mode);
          else
            do_mosaic(platefile, DiskImageView<PixelGrayA<int16> >(image_files[i]), 
                      image_files[i], transaction_id_override, georef, output_mode);
          break;
        case VW_CHANNEL_FLOAT32:
          if (has_nodata_value) {
            do_mosaic(platefile,
                      mask_to_alpha(create_mask(DiskImageView<PixelGray<float32> >(image_files[i]),
                                                nodata_value)),
                      image_files[i], transaction_id_override, georef, output_mode);
          } else
            do_mosaic(platefile, DiskImageView<PixelGrayA<float32> >(image_files[i]),
                      image_files[i], transaction_id_override, georef, output_mode);
          break;
        default:
          vw_out() << "Image contains a channel type not supported by image2plate.\n";
          exit(1);
        }
        break;

      case VW_PIXEL_RGB:
      case VW_PIXEL_RGBA:
        switch(rsrc->channel_type()) {
        case VW_CHANNEL_UINT8:  
          if (has_nodata_value)
            do_mosaic(platefile, 
                      pixel_cast<PixelRGBA<uint8> >(mask_to_alpha(create_mask(DiskImageView<PixelGray<uint8> >(image_files[i]), 
                                                                              nodata_value))),
                      image_files[i], transaction_id_override, georef, output_mode);
          else
            do_mosaic(platefile, DiskImageView<PixelRGBA<uint8> >(image_files[i]), 
                      image_files[i], transaction_id_override, georef, output_mode);
          break;
        default:
          vw_out() << "Platefile contains a channel type not supported by image2plate.\n";
          exit(1);
        }
        break;
      default:
        vw_out() << "Image contains a pixel type not supported by image2plate.\n";
        exit(1);
      }
    }


  } catch (vw::Exception &e) {
    vw_out() << "A vision workbench error occured: " << e.what() << "\n";
    exit(1);
  }
  
}
