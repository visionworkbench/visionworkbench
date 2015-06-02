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


#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <cstdlib>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/foreach.hpp>
namespace fs = boost::filesystem;
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
#include <vw/Cartography/GeoReference.h>
#include <vw/tools/Common.h>

using namespace vw;

struct Options {
  // Input
  std::string input_file_name;
  std::string shaded_relief_file_name;

  // Settings
  std::string output_file_name, colormap_style;
  float       nodata_value, min_val, max_val;
  bool        draw_legend;

  typedef Vector<uint8,3>                 Vector3u;
  typedef std::pair<std::string,Vector3u> lut_element;
  typedef std::vector<lut_element>        lut_type;

  lut_type lut;
  std::map<float,Vector3u> lut_map;
};

// Colormap function
class ColormapFunc : public ReturnFixedType<PixelMask<PixelRGB<uint8> > > {
  typedef std::map<float,Options::Vector3u> map_type;
  map_type m_colormap;

public:
  ColormapFunc( std::map<float,Options::Vector3u> const& map) : m_colormap(map) {}

  template <class PixelT>
  PixelMask<PixelRGB<uint8> > operator() (PixelT const& pix) const {
    if (is_transparent(pix))
      return PixelMask<PixelRGB<uint8> >(); // Skip transparent pixels

    float val = compound_select_channel<const float&>(pix,0);
    if (val > 1.0) val = 1.0;
    if (val < 0.0) val = 0.0;

    // Get locations on sparse colormap that bound this pixel value
    map_type::const_iterator bot = m_colormap.upper_bound( val ); bot--;
    map_type::const_iterator top = m_colormap.upper_bound( val );

    if ( top == m_colormap.end() ) // If this is above the top colormap value
      return PixelRGB<uint8>(bot->second[0], bot->second[1], bot->second[2]); // Use the max colormap value

    // Otherwise determine a proportional color between the bounding colormap values
    Options::Vector3u output = bot->second + // Lower colormap value
                                ( ((val - bot->first)/(top->first - bot->first)) * // Fraction of distance to next colormap value
                                  (Vector3i(top->second) - Vector3i(bot->second))  ); // Difference between successive colormap values
    return PixelRGB<uint8>(output[0], output[1], output[2]);
  }
};

template <class ViewT>
UnaryPerPixelView<ViewT, ColormapFunc> colormap(ImageViewBase<ViewT> const& view,
                                                std::map<float,Options::Vector3u> const& map) {
  return UnaryPerPixelView<ViewT, ColormapFunc>(view.impl(), ColormapFunc(map));
}

// -------------------------------------------------------------------------------------

template <class PixelT>
void do_colorized_dem(Options& opt) {
  vw_out() << "Creating colorized DEM.\n";

  cartography::GeoReference georef;
  cartography::read_georeference(georef, opt.input_file_name);

  // Attempt to extract nodata value
  SrcImageResource *disk_dem_rsrc =
    DiskImageResource::open(opt.input_file_name);
  if (opt.nodata_value != std::numeric_limits<float>::max()) {
    vw_out() << "\t--> Using user-supplied nodata value: " << opt.nodata_value << ".\n";
  } else if ( disk_dem_rsrc->has_nodata_read() ) {
    opt.nodata_value = disk_dem_rsrc->nodata_read();
    vw_out() << "\t--> Extracted nodata value from file: " << opt.nodata_value << ".\n";
  }

  // Compute min/max of input image values
  DiskImageView<PixelT> disk_dem_file(opt.input_file_name);
  ImageViewRef<PixelGray<float> > input_image =
    pixel_cast<PixelGray<float> >(select_channel(disk_dem_file,0));
  if (opt.min_val == 0 && opt.max_val == 0) {
    min_max_channel_values( create_mask( input_image, opt.nodata_value),
                            opt.min_val, opt.max_val);
    vw_out() << "\t--> DEM color map range: ["
             << opt.min_val << "  " << opt.max_val << "]\n";
  } else {
    vw_out() << "\t--> Using user-specified color map range: ["
             << opt.min_val << "  " << opt.max_val << "]\n";
  }

  // Convert lut to lut_map ( converts altitudes to relative percent )
  opt.lut_map.clear();
  BOOST_FOREACH( Options::lut_element const& pair, opt.lut ) {
    try {
      if ( boost::contains(pair.first,"%") ) {
        float key = boost::lexical_cast<float>(boost::erase_all_copy(pair.first,"%"))/100.0;
        opt.lut_map[ key ] = pair.second;
      } else {
        float key = boost::lexical_cast<float>(pair.first);
        opt.lut_map[ ( key - opt.min_val ) / ( opt.max_val - opt.min_val ) ] =
          pair.second;
      }
    } catch ( const boost::bad_lexical_cast& e ) {
      continue;
    }
  }

  // Mask input
  ImageViewRef<PixelMask<PixelGray<float> > > dem;
  if ( PixelHasAlpha<PixelT>::value )
    dem = alpha_to_mask(channel_cast<float>(disk_dem_file) );
  else if (opt.nodata_value != std::numeric_limits<float>::max())
    dem = channel_cast<float>(create_mask(input_image, opt.nodata_value));
  else if ( disk_dem_rsrc->has_nodata_read() )
    dem = create_mask(input_image, opt.nodata_value);
  else
    dem = pixel_cast<PixelMask<PixelGray<float> > >(input_image);

  delete disk_dem_rsrc;

  // Apply colormap
  ImageViewRef<PixelMask<PixelRGB<uint8> > > colorized_image =
    colormap(normalize(dem, opt.min_val, opt.max_val, 0, 1.0), opt.lut_map);

  if (!opt.shaded_relief_file_name.empty()) { // Using a hillshade file
    vw_out() << "\t--> Incorporating hillshading from: "
             << opt.shaded_relief_file_name << ".\n";
    DiskImageView<PixelGray<float> >
      shaded_relief_image(opt.shaded_relief_file_name); // It's okay to throw  away the
                                                        // second channel if it exists.
    ImageViewRef<PixelMask<PixelRGB<uint8> > > shaded_image =
      copy_mask(channel_cast<uint8>(colorized_image*pixel_cast<float>(shaded_relief_image)),
                colorized_image);
    vw_out() << "Writing color-mapped image: " << opt.output_file_name << "\n";
    boost::scoped_ptr<DiskImageResource> r(DiskImageResource::create(opt.output_file_name,shaded_image.format()));
    if ( r->has_block_write() )
      r->set_block_write_size( Vector2i( vw_settings().default_tile_size(),
                                         vw_settings().default_tile_size() ) );
    write_georeference( *r, georef );
    block_write_image( *r, shaded_image,
                       TerminalProgressCallback( "tools.colormap", "Writing:") );
  } else { // Not using a hillshade file
    vw_out() << "Writing color-mapped image: " << opt.output_file_name << "\n";

    boost::scoped_ptr<DiskImageResource> r(DiskImageResource::create(opt.output_file_name,colorized_image.format()));
    if ( r->has_block_write() )
      r->set_block_write_size( Vector2i( vw_settings().default_tile_size(),
                                         vw_settings().default_tile_size() ) );
    write_georeference( *r, georef );
    block_write_image( *r, colorized_image,
                       TerminalProgressCallback( "tools.colormap", "Writing:") );
  }
}

void save_legend( Options const& opt) {
  ImageView<PixelGray<float> > img(100, 500);
  for (int j = 0; j < img.rows(); ++j) {
    float val = float(j) / img.rows();
    for (int i = 0; i < img.cols(); ++i) {
      img(i,j) = val;
    }
  }

  ImageViewRef<PixelMask<PixelRGB<uint8> > > colorized_image =
    colormap(img, opt.lut_map);
  write_image("legend.png", channel_cast_rescale<uint8>(apply_mask(colorized_image)));
}

void handle_arguments( int argc, char *argv[], Options& opt ) {
  po::options_description general_options("");
  general_options.add_options()
    ("shaded-relief-file,s", po::value(&opt.shaded_relief_file_name),
                      "Specify a shaded relief image (grayscale) to apply to the colorized image.")
    ("output-file,o", po::value(&opt.output_file_name),
                      "Specify the output file.")
    ("colormap-style",po::value(&opt.colormap_style)->default_value("jet"),
                      "Specify the colormap style. Options: jet (default), binary-red-blue, or the name of a file having the colormap, similar to the file used by gdaldem.")
    ("nodata-value",  po::value(&opt.nodata_value)->default_value(std::numeric_limits<float>::max()),
                      "Remap the DEM default value to the min altitude value.")
    ("min",           po::value(&opt.min_val)->default_value(0),
                      "Minimum height of the color map.")
    ("max",           po::value(&opt.max_val)->default_value(0),
                      "Maximum height of the color map.")
    ("moon",          "Set the min and max values to [-8499 10208] meters, which is suitable for covering elevations on the Moon.")
    ("mars",          "Set the min and max values to [-8208 21249] meters, which is suitable for covering elevations on Mars.")
    ("legend",        "Generate the colormap legend.  This image is saved (without labels) as \'legend.png\'.")
    ("help,h",        "Display this help message.");

  po::options_description positional("");
  positional.add_options()
    ("input-file", po::value(&opt.input_file_name), "Input disparity map.");

  po::positional_options_description positional_desc;
  positional_desc.add("input-file", 1);

  po::options_description all_options;
  all_options.add(general_options).add(positional);

  po::variables_map vm;
  try {
    po::store( po::command_line_parser( argc, argv ).options(all_options).positional(positional_desc).run(), vm );
    po::notify( vm );
  } catch (const po::error& e) {
    vw_throw( ArgumentErr() << "Error parsing input:\n\t"
              << e.what() << general_options );
  }

  std::ostringstream usage;
  usage << "Usage: " << argv[0] << " [options] <input dem> \n";

  if ( vm.count("help") )
    vw_throw( ArgumentErr() << usage.str() << general_options );
  if ( opt.input_file_name.empty() )
    vw_throw( ArgumentErr() << "Missing input file!\n"
              << usage.str() << general_options );
  if ( vm.count("moon") ) {
    opt.min_val = -8499;
    opt.max_val = 10208;
  }
  if ( vm.count("mars") ) {
    opt.min_val = -8208;
    opt.max_val = 21249;
  }
  if ( opt.output_file_name.empty() )
    opt.output_file_name =
      fs::path(opt.input_file_name).replace_extension().string()+"_CMAP.tif";
  opt.draw_legend = vm.count("legend");
}

int main( int argc, char *argv[] ) {

  Options opt;
  try {
    handle_arguments( argc, argv, opt );

    // Decide legend
    if ( opt.colormap_style.empty() || opt.colormap_style == "jet" ) { // default
      opt.lut.clear();
      opt.lut.push_back( Options::lut_element("0%",   Options::Vector3u(  0,   0,   0)) ); // Black
      opt.lut.push_back( Options::lut_element("20.8%",Options::Vector3u(  0,   0, 255)) ); // Blue
      opt.lut.push_back( Options::lut_element("25%",  Options::Vector3u(  0,   0, 255)) ); // Blue
      opt.lut.push_back( Options::lut_element("37.5%",Options::Vector3u(  0, 191, 255)) ); // Light blue
      opt.lut.push_back( Options::lut_element("41.7%",Options::Vector3u(  0, 255, 255)) ); // Teal
      opt.lut.push_back( Options::lut_element("58.3%",Options::Vector3u(255, 255,  51)) ); // Yellow
      opt.lut.push_back( Options::lut_element("62.5%",Options::Vector3u(255, 191,   0)) ); // Orange
      opt.lut.push_back( Options::lut_element("75%",  Options::Vector3u(255,   0,   0)) ); // Red
      opt.lut.push_back( Options::lut_element("79.1%",Options::Vector3u(255,   0,   0)) ); // Red
      opt.lut.push_back( Options::lut_element("100%", Options::Vector3u(  0,   0,   0)) ); // Black
    }
    if ( opt.colormap_style == "binary-red-blue" ) {
      opt.lut.push_back( Options::lut_element( "0%",      Options::Vector3u( 59,  76, 192)) );
      opt.lut.push_back( Options::lut_element( "3.13%",   Options::Vector3u( 68,  90, 204)) );
      opt.lut.push_back( Options::lut_element( "6.25%",   Options::Vector3u( 78, 104, 215)) );
      opt.lut.push_back( Options::lut_element( "9.38%",   Options::Vector3u( 88, 117, 225)) );
      opt.lut.push_back( Options::lut_element("12.50%",   Options::Vector3u( 98, 130, 234)) );
      opt.lut.push_back( Options::lut_element("15.63%",   Options::Vector3u(108, 142, 241)) );
      opt.lut.push_back( Options::lut_element("18.75%",   Options::Vector3u(119, 154, 247)) );
      opt.lut.push_back( Options::lut_element("21.88%",   Options::Vector3u(130, 165, 251)) );
      opt.lut.push_back( Options::lut_element("25.00%",   Options::Vector3u(141, 176, 254)) );
      opt.lut.push_back( Options::lut_element("28.13%",   Options::Vector3u(152, 185, 255)) );
      opt.lut.push_back( Options::lut_element("31.25%",   Options::Vector3u(163, 194, 255)) );
      opt.lut.push_back( Options::lut_element("34.38%",   Options::Vector3u(174, 201, 253)) );
      opt.lut.push_back( Options::lut_element("37.50%",   Options::Vector3u(184, 208, 249)) );
      opt.lut.push_back( Options::lut_element("40.63%",   Options::Vector3u(194, 213, 244)) );
      opt.lut.push_back( Options::lut_element("43.75%",   Options::Vector3u(204, 217, 238)) );
      opt.lut.push_back( Options::lut_element("46.88%",   Options::Vector3u(213, 219, 230)) );
      opt.lut.push_back( Options::lut_element("50.00%",   Options::Vector3u(221, 221, 221)) );
      opt.lut.push_back( Options::lut_element("53.13%",   Options::Vector3u(229, 216, 209)) );
      opt.lut.push_back( Options::lut_element("56.25%",   Options::Vector3u(236, 211, 198)) );
      opt.lut.push_back( Options::lut_element("59.38%",   Options::Vector3u(241, 204, 185)) );
      opt.lut.push_back( Options::lut_element("62.50%",   Options::Vector3u(245, 196, 173)) );
      opt.lut.push_back( Options::lut_element("65.63%",   Options::Vector3u(247, 187, 160)) );
      opt.lut.push_back( Options::lut_element("68.75%",   Options::Vector3u(247, 177, 148)) );
      opt.lut.push_back( Options::lut_element("71.88%",   Options::Vector3u(247, 166, 135)) );
      opt.lut.push_back( Options::lut_element("75.00%",   Options::Vector3u(244, 154, 123)) );
      opt.lut.push_back( Options::lut_element("78.13%",   Options::Vector3u(241, 141, 111)) );
      opt.lut.push_back( Options::lut_element("81.25%",   Options::Vector3u(236, 127,  99)) );
      opt.lut.push_back( Options::lut_element("84.38%",   Options::Vector3u(229, 112,  87)) );
      opt.lut.push_back( Options::lut_element("87.50%",   Options::Vector3u(222,  97,  76)) );
      opt.lut.push_back( Options::lut_element("90.63%",   Options::Vector3u(213,  80,  66)) );
      opt.lut.push_back( Options::lut_element("93.75%",   Options::Vector3u(203,  62,  56)) );
      opt.lut.push_back( Options::lut_element("96.88%",   Options::Vector3u(192,  40,  47)) );
      opt.lut.push_back( Options::lut_element("100.00%",  Options::Vector3u(180,   4,  38)) );
    }
    else {
      // Read input LUT
      typedef boost::tokenizer<> tokenizer;
      boost::char_delimiters_separator<char> sep(false,",:");

      std::ifstream lut_file( opt.colormap_style.c_str() );
      if ( !lut_file.is_open() )
        vw_throw( IOErr() << "Unable to open colormap style file: "
                  << opt.colormap_style );
      std::string line;
      std::getline( lut_file, line );
      while ( lut_file.good() ) {

        // Skip lines containing spaces only
        bool only_spaces = true;
        for (unsigned s = 0; s < line.size(); s++) if (!isspace(line[s])) only_spaces = false;
        if (only_spaces){
          std::getline( lut_file, line );
          continue;
        }

        tokenizer tokens(line,sep);
        tokenizer::iterator iter = tokens.begin();

        std::string key;
        Options::Vector3u value;

        try {
          if ( iter == tokens.end()) vw_throw( IOErr() << "Unable to read input LUT" );
          key = *iter; iter++;
          if ( iter == tokens.end()) vw_throw( IOErr() << "Unable to read input LUT" );
          value[0] = boost::numeric_cast<uint8>(boost::lexical_cast<int>(*iter)); iter++;
          if ( iter == tokens.end()) vw_throw( IOErr() << "Unable to read input LUT" );
          value[1] = boost::numeric_cast<uint8>(boost::lexical_cast<int>(*iter)); iter++;
          if ( iter == tokens.end()) vw_throw( IOErr() << "Unable to read input LUT" );
          value[2] = boost::numeric_cast<uint8>(boost::lexical_cast<int>(*iter));
        } catch ( const boost::bad_lexical_cast& e ) {
          std::getline( lut_file, line );
          continue;
        }
        opt.lut.push_back( Options::lut_element(key, value) );
        std::getline( lut_file, line );
      }
      lut_file.close();
    }

    // Get the right pixel/channel type.
    ImageFormat fmt = tools::image_format(opt.input_file_name);

    switch(fmt.pixel_format) {
    case VW_PIXEL_GRAY:
      switch(fmt.channel_type) {
      case VW_CHANNEL_UINT8:  do_colorized_dem<PixelGray<uint8>   >(opt); break;
      case VW_CHANNEL_INT16:  do_colorized_dem<PixelGray<int16>   >(opt); break;
      case VW_CHANNEL_UINT16: do_colorized_dem<PixelGray<uint16>  >(opt); break;
      default:                do_colorized_dem<PixelGray<float32> >(opt); break;
      }
      break;
    case VW_PIXEL_GRAYA:
      switch(fmt.channel_type) {
      case VW_CHANNEL_UINT8:  do_colorized_dem<PixelGrayA<uint8>   >(opt); break;
      case VW_CHANNEL_INT16:  do_colorized_dem<PixelGrayA<int16>   >(opt); break;
      case VW_CHANNEL_UINT16: do_colorized_dem<PixelGrayA<uint16>  >(opt); break;
      default:                do_colorized_dem<PixelGrayA<float32> >(opt); break;
      }
      break;
    default:
      vw_throw( ArgumentErr() << "Unsupported pixel format. The DEM image must have only one channel." );
    }

    // Draw legend
    if ( opt.draw_legend )
      save_legend(opt);

  } catch ( const ArgumentErr& e ) {
    vw_out() << e.what() << std::endl;
    return 1;
  } catch ( const Exception& e ) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }
  return 0;
}
