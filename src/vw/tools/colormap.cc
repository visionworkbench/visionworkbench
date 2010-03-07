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

#ifdef NDEBUG
#undef NDEBUG
#endif

#include <cstdlib>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <vw/Core/Functors.h>
#include <vw/Image/Algorithms.h>
#include <vw/Image/ImageMath.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Image/PerPixelViews.h>
#include <vw/Image/PixelMask.h>
#include <vw/Image/MaskViews.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Image/Statistics.h>
#include <vw/FileIO/DiskImageView.h>
#include <vw/Cartography/FileIO.h>
#include <vw/Cartography/GeoReference.h>

using namespace vw;

// Global variables
std::string input_file_name, output_file_name = "", shaded_relief_file_name;
float nodata_value;
float min_val = 0, max_val = 0;

// Useful functions

/// Erases a file suffix if one exists and returns the base string
static std::string prefix_from_filename(std::string const& filename) {
  std::string result = filename;
  int index = result.rfind(".");
  if (index != -1) 
    result.erase(index, result.size());
  return result;
}

// Colormap function
class ColormapFunc : public ReturnFixedType<PixelMask<PixelRGB<float> > > {

public:
  ColormapFunc() {}
  
  template <class PixelT>
  PixelMask<PixelRGB<float> > operator() (PixelT const& pix) const {
    if (is_transparent(pix)) 
      return PixelMask<PixelRGB<float> >();
    
    float val = compound_select_channel<const float&>(pix,0);
    if (val > 1.0) val = 1.0;
    if (val < 0.0) val = 0.0;
    
    Vector2 red_range(3/8.0, 1.0);
    Vector2 green_range(2/8.0, 6/8.0);
    Vector2 blue_range(0.0, 5/8.0);

    float red_span = red_range[1] - red_range[0];
    float blue_span = blue_range[1] - blue_range[0];
    float green_span = green_range[1] - green_range[0];

    float red = 0;
    float green = 0;
    float blue = 0;

    // Red
    if (val >= red_range[0] && val <= red_range[0]+red_span/3) 
      red = (val-red_range[0])/(red_span/3);
    if (val >= red_range[0]+red_span/3 && val <= red_range[0]+2*red_span/3) 
      red = 1.0;
    if (val >= red_range[0]+2*red_span/3 && val <= red_range[1]) 
      red = 1-((val-(red_range[0]+2*red_span/3))/(red_span/3));

    // Blue
    if (val >= blue_range[0] && val <= blue_range[0]+blue_span/3) 
      blue = (val-blue_range[0])/(blue_span/3);
    if (val >= blue_range[0]+blue_span/3 && val <= blue_range[0]+2*blue_span/3) 
      blue = 1.0;
    if (val >= blue_range[0]+2*blue_span/3 && val <= blue_range[1]) 
      blue = 1-((val-(blue_range[0]+2*blue_span/3))/(blue_span/3));

     // Green
    if (val >= green_range[0] && val <= green_range[0]+green_span/3) 
      green = (val-green_range[0])/(green_span/3);
    if (val >= green_range[0]+green_span/3 && val <= green_range[0]+2*green_span/3) 
      green = 1.0;
    if (val >= green_range[0]+2*green_span/3 && val <= green_range[1]) 
      green = 1-((val-(green_range[0]+2*green_span/3))/(green_span/3));

    return PixelRGB<float> ( red, green, blue );
  }
};
  
template <class ViewT> 
UnaryPerPixelView<ViewT, ColormapFunc> colormap(ImageViewBase<ViewT> const& view) {
  return UnaryPerPixelView<ViewT, ColormapFunc>(view.impl(), ColormapFunc());
}

// ---------------------------------------------------------------------------------------------------

template <class PixelT>
void do_colorized_dem(po::variables_map const& vm) {
  vw_out() << "Creating colorized DEM.\n";

  cartography::GeoReference georef;
  cartography::read_georeference(georef, input_file_name);

  // Attempt to extract nodata value
  DiskImageResource *disk_dem_rsrc = DiskImageResource::open(input_file_name);
  if (vm.count("nodata-value")) {
    vw_out() << "\t--> Using user-supplied nodata value: " << nodata_value << ".\n";
  } else if ( disk_dem_rsrc->has_nodata_value() ) {
    nodata_value = disk_dem_rsrc->nodata_value();
    vw_out() << "\t--> Extracted nodata value from file: " << nodata_value << ".\n";
  } 

  // Compute min/max
  DiskImageView<PixelT> disk_dem_file(input_file_name);
  ImageViewRef<PixelGray<float> > input_image = pixel_cast<PixelGray<float> >(select_channel(disk_dem_file,0));
  if (min_val == 0 && max_val == 0) {
    min_max_channel_values( create_mask( input_image, nodata_value), min_val, max_val);
    vw_out() << "\t--> DEM color map range: [" << min_val << "  " << max_val << "]\n";
  } else {
    vw_out() << "\t--> Using user-specified color map range: [" << min_val << "  " << max_val << "]\n";
  }

  ImageViewRef<PixelMask<PixelGray<float> > > dem;
  if ( PixelHasAlpha<PixelT>::value ) {
    dem = alpha_to_mask(channel_cast<float>(disk_dem_file) );
  } else if (vm.count("nodata-value")) {
    dem = channel_cast<float>(create_mask(input_image, nodata_value));
  } else if ( disk_dem_rsrc->has_nodata_value() ) {
    dem = create_mask(input_image, nodata_value);
  } else {
    dem = pixel_cast<PixelMask<PixelGray<float> > >(input_image);
  }
  delete disk_dem_rsrc;

  ImageViewRef<PixelMask<PixelRGB<float> > > colorized_image = colormap(normalize(dem,min_val,max_val,0,1.0));

  if (shaded_relief_file_name != "") {
    vw_out() << "\t--> Incorporating hillshading from: " << shaded_relief_file_name << ".\n";
    DiskImageView<PixelMask<float> > shaded_relief_image(shaded_relief_file_name);
    ImageViewRef<PixelMask<PixelRGB<float> > > shaded_image = copy_mask(colorized_image*apply_mask(shaded_relief_image), shaded_relief_image);
    vw_out() << "Writing image color-mapped image: " << output_file_name << "\n";
    write_georeferenced_image(output_file_name, channel_cast_rescale<uint8>(shaded_image), georef,
                              TerminalProgressCallback( "tools.colormap", "Writing:"));
  } else {
    vw_out() << "Writing image color-mapped image: " << output_file_name << "\n";
    write_georeferenced_image(output_file_name, channel_cast_rescale<uint8>(colorized_image), georef,
                              TerminalProgressCallback( "tools.colormap", "Writing:"));
  }
}

void save_legend() {
  min_val = 0.0;
  max_val = 1.0;
  
  ImageView<PixelGray<float> > img(200, 500);
  for (int j = 0; j < img.rows(); ++j) {
    float val = float(j) / img.rows();
    for (int i = 0; i < img.cols(); ++i) {
      img(i,j) = val;
    }
  }
  
  ImageViewRef<PixelMask<PixelRGB<float> > > colorized_image = colormap(img);
  write_image("legend.png", channel_cast_rescale<uint8>(apply_mask(colorized_image)));
}


int main( int argc, char *argv[] ) {

  set_debug_level(VerboseDebugMessage-1);
  
  po::options_description desc("Description: Produces a colorized image of a DEM \n\nUsage: colormap [options] <input file> \n\nOptions");
  desc.add_options()
    ("help,h", "Display this help message")
    ("input-file", po::value<std::string>(&input_file_name), "Explicitly specify the input file")
    ("shaded-relief-file,s", po::value<std::string>(&shaded_relief_file_name)->default_value(""), "Specify a shaded relief image (grayscale) to apply to the colorized image.")
    ("output-file,o", po::value<std::string>(&output_file_name), "Specify the output file")
    ("nodata-value", po::value<float>(&nodata_value), "Remap the DEM default value to the min altitude value.")
    ("min", po::value<float>(&min_val), "Minimum height of the color map.")
    ("max", po::value<float>(&max_val), "Maximum height of the color map.")
    ("moon", "Set the min and max values to [-8499 10208] meters, which is suitable for covering elevations on the Moon.")
    ("mars", "Set the min and max values to [-8208 21249] meters, which is suitable for covering elevations on Mars.")
    ("legend", "Generate the colormap legend.  This image is saved (without labels) as \'legend.png\'")
    ("verbose", "Verbose output");
  po::positional_options_description p;
  p.add("input-file", 1);

  po::variables_map vm;
  try {
    po::store( po::command_line_parser( argc, argv ).options(desc).positional(p).run(), vm );
    po::notify( vm );
  } catch (po::error &e) {
    std::cout << "An error occured while parsing command line arguments.\n";
    std::cout << "\t" << e.what() << "\n\n";
    std::cout << desc << std::endl;
    return 1;
  }

  if( vm.count("help") ) {
    std::cout << desc << std::endl;
    return 1;
  }

  if( vm.count("legend") ) {
    std::cout << "\t--> Saving legend file to disk as \'legend.png\'\n";
    save_legend();
    exit(0);
  }

  // This is a reasonable range of elevation values to cover global
  // lunar topography.
  if( vm.count("moon") ) {
    min_val = -8499;  
    max_val = 10208;
  }

  // This is a reasonable range of elevation values to cover global
  // mars topography.
  if( vm.count("mars") ) {
    min_val = -8208;
    max_val = 21249;
  }

  if( vm.count("input-file") != 1 ) {
    std::cout << "Error: Must specify exactly one input file!\n" << std::endl;
    std::cout << desc << std::endl;
    return 1;
  }

  if( output_file_name == "" ) {
    output_file_name = prefix_from_filename(input_file_name) + "_CMAP.tif";
  }

  if( vm.count("verbose") ) {
    set_debug_level(VerboseDebugMessage);
  }

  try {

    // Get the right pixel/channel type.
    DiskImageResource *rsrc = DiskImageResource::open(input_file_name);
    ChannelTypeEnum channel_type = rsrc->channel_type();
    PixelFormatEnum pixel_format = rsrc->pixel_format();
    delete rsrc;


    switch(pixel_format) {
    case VW_PIXEL_GRAY:
      switch(channel_type) {
      case VW_CHANNEL_UINT8:  do_colorized_dem<PixelGray<uint8>   >(vm); break;
      case VW_CHANNEL_INT16:  do_colorized_dem<PixelGray<int16>   >(vm); break;
      case VW_CHANNEL_UINT16: do_colorized_dem<PixelGray<uint16>  >(vm); break;
      default:                do_colorized_dem<PixelGray<float32> >(vm); break;
      }
      break;
    case VW_PIXEL_GRAYA:
      switch(channel_type) {
      case VW_CHANNEL_UINT8:  do_colorized_dem<PixelGrayA<uint8>   >(vm); break;
      case VW_CHANNEL_INT16:  do_colorized_dem<PixelGrayA<int16>   >(vm); break;
      case VW_CHANNEL_UINT16: do_colorized_dem<PixelGrayA<uint16>  >(vm); break;
      default:                do_colorized_dem<PixelGrayA<float32> >(vm); break;
      }
      break;
    default:
      std::cout << "Error: Unsupported pixel format.  The DEM image must have only one channel.";
      exit(0);
    }
  }
  catch( Exception& e ) {
    std::cerr << "Error: " << e.what() << std::endl;
  }

  return 0;
}
