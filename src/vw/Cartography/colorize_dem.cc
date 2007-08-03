#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#ifdef NDEBUG
#undef NDEBUG
#endif

#include <stdlib.h>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <vw/vw.h>
#include <vw/Image/PerPixelAccessorViews.h>

using namespace vw;

class ColormapFunc : public ReturnFixedType<PixelRGB<float> > {

public:
  ColormapFunc() {}
  
  template <class PixelT>
  PixelRGB<float> operator() (PixelT const& pix) const {
    const float val = compound_select_channel<const float&>(pix,0);
    
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

template <class PixelT>
class RemapPixelFunc : public ReturnFixedType<PixelT> {

  PixelT m_src_pix, m_dst_pix;

public:
  RemapPixelFunc(PixelT const& src_pix, PixelT const& dst_pix) : 
    m_src_pix(src_pix), m_dst_pix(dst_pix) {}
  
  PixelT operator() (PixelT const& pix) const {
    if (pix == m_src_pix) 
      return m_dst_pix;
    else 
      return pix;
  }
};
  
template <class ViewT> 
UnaryPerPixelView<ViewT, RemapPixelFunc<typename ViewT::pixel_type> > remap_pixel(ImageViewBase<ViewT> const& view, typename ViewT::pixel_type src_value, typename ViewT::pixel_type dst_value) {
  return UnaryPerPixelView<ViewT, RemapPixelFunc<typename ViewT::pixel_type> >(view.impl(), RemapPixelFunc<typename ViewT::pixel_type>(src_value, dst_value));
}


int main( int argc, char *argv[] ) {

  set_debug_level(VerboseDebugMessage-1);

  std::string input_file_name, output_file_name, shaded_relief_file_name;
  float offset;
  float dem_default_value;
  float clamp_low, clamp_high;
  
  po::options_description desc("Options");
  desc.add_options()
    ("help", "Display this help message")
    ("input-file", po::value<std::string>(&input_file_name), "Explicitly specify the input file")
    ("shaded-relief-file", po::value<std::string>(&shaded_relief_file_name)->default_value(""), "Specify a shaded relief image (grayscale) to apply to the colorized image.")
    ("output-file,o", po::value<std::string>(&output_file_name)->default_value("output.tif"), "Specify the output file")
    ("offset", po::value<float>(&offset)->default_value(0), "Add this offset to the height values in the file.")
    ("clamp-low", po::value<float>(&clamp_low), "Clamp the DEM to the specified range.  Must be used in conjunction with clamp-high")
    ("clamp-high", po::value<float>(&clamp_high), "Clamp the DEM to the specified range.  Must be used in conjunction with clamp-low")
    ("dem-default-value", po::value<float>(&dem_default_value), "Remap the DEM default value to the min altitude value.")
    ("verbose", "Verbose output");
  po::positional_options_description p;
  p.add("input-file", 1);

  po::variables_map vm;
  po::store( po::command_line_parser( argc, argv ).options(desc).positional(p).run(), vm );
  po::notify( vm );

  if( vm.count("help") ) {
    std::cout << desc << std::endl;
    return 1;
  }

  if( vm.count("input-file") != 1 ) {
    std::cout << "Error: Must specify exactly one input file!" << std::endl;
    std::cout << desc << std::endl;
    return 1;
  }

  if( vm.count("verbose") ) {
    set_debug_level(VerboseDebugMessage);
  }

  try {
    cartography::GeoReference georef;
    cartography::read_georeference(georef, input_file_name);

    std::cout << georef << "\n";

    DiskImageView<float> disk_dem_file(input_file_name);
    ImageViewRef<float> dem = channel_cast<float>(disk_dem_file) + offset;

    if (vm.count("clamp-low") && vm.count("clamp-high")) {
      std::cout << "Clamping DEM pixel values to range [" << clamp_low << ", " << clamp_high << "].\n";
      dem = clamp(dem, clamp_low, clamp_high);
    }

    if (vm.count("dem-default-value")) {
      std::cout << "Remapping default pixel value " << dem_default_value << " to 0.\n";
      dem = remap_pixel(dem, dem_default_value, 0);
    }

    float min_val, max_val;
    min_max_channel_values(dem, min_val, max_val);
    std::cout << "\n\nDEM Pixel Range: [" << min_val << "  " << max_val << "]\n\n";

    ImageViewRef<PixelRGB<float> > colorized_image = colormap(normalize(dem));

    std::cout << "Writing image.\n";
    if (shaded_relief_file_name != "") {
      DiskImageView<float> shaded_relief_image(shaded_relief_file_name);
      write_georeferenced_image(output_file_name, channel_cast_rescale<uint8>(colorized_image*shaded_relief_image), georef, TerminalProgressCallback());
    } else {
      write_georeferenced_image(output_file_name, channel_cast_rescale<uint8>(colorized_image), georef, TerminalProgressCallback());
    }
  }
  catch( Exception& e ) {
    std::cerr << "Error: " << e.what() << std::endl;
  }

  return 0;
}
