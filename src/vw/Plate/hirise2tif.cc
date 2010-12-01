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

#include <stdlib.h>
#include <iostream>
#include <fstream>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <vw/Core/Cache.h>
#include <vw/Core/ProgressCallback.h>
#include <vw/Math/Matrix.h>
#include <vw/Image/Transform.h>
#include <vw/Image/Palette.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/FileIO/DiskImageResourceJPEG.h>
#include <vw/FileIO/DiskImageResourcePNG.h>
#include <vw/FileIO/DiskImageResourceGDAL.h>
#include <vw/FileIO/DiskImageView.h>
#include <vw/Cartography/GeoReference.h>
#include <vw/Cartography/GeoTransform.h>
#include <vw/Mosaic/ImageComposite.h>
#include <vw/Mosaic/QuadTreeGenerator.h>
#include <vw/Mosaic/KMLQuadTreeConfig.h>

using namespace vw;
using namespace vw::cartography;
using namespace vw::mosaic;

// --------------------------------------------------------------------------
//                     BIT DEPTH CONVERSION ROUTINES
// --------------------------------------------------------------------------


/// Parses arguments for the min/max values of the color channels in
/// the HiRISE image.  If no min/max values are supplied this class
/// computes suitable min/max values by opening up the images and
/// running the statistics itself. (this takes a while...)
class ImageStats {

  void error(std::string arg, std::string const& params) {
    vw_out(ErrorMessage) << "Error parsing arguments for --" << arg << " : " << params << "\n";
    exit(1);
  }

public:

  // Public variables.  You access these directly.
  uint16 min_gray;
  uint16 max_gray;
  uint16 min_i;
  uint16 max_i;
  uint16 min_r;
  uint16 max_r;
  uint16 min_b;
  uint16 max_b;

  // Constructor
  ImageStats(std::string const& range_string,
             std::string const& gray_filename,
             std::string const& color_filename) {

    typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
    boost::char_separator<char> sep(",;");

    // STEP 1 : PARSE RANGE STRING

    if (range_string.empty()) {
      std::cout << "Error: you must specify pixel ranges using --range...\n";
      exit(1);
    }

    // If the range string is not empty, we attempt to parse the
    // two parameters out from the range string.
    tokenizer tokens(range_string, sep);
    tokenizer::iterator tok_iter = tokens.begin();

    if (tok_iter == tokens.end()) this->error("stats", range_string);
    min_gray = boost::lexical_cast<uint16>(*tok_iter);
    ++tok_iter;

    if (tok_iter == tokens.end()) this->error("stats", range_string);
    max_gray = boost::lexical_cast<uint16>(*tok_iter);
    ++tok_iter;

    if (tok_iter == tokens.end()) this->error("stats", range_string);
    min_i = boost::lexical_cast<uint16>(*tok_iter);
    ++tok_iter;

    if (tok_iter == tokens.end()) this->error("stats", range_string);
    max_i = boost::lexical_cast<uint16>(*tok_iter);
    ++tok_iter;

    if (tok_iter == tokens.end()) this->error("stats", range_string);
    min_r = boost::lexical_cast<uint16>(*tok_iter);
    ++tok_iter;

    if (tok_iter == tokens.end()) this->error("stats", range_string);
    max_r = boost::lexical_cast<uint16>(*tok_iter);
    ++tok_iter;

    if (tok_iter == tokens.end()) this->error("stats", range_string);
    min_b = boost::lexical_cast<uint16>(*tok_iter);
    ++tok_iter;

    if (tok_iter == tokens.end()) this->error("stats", range_string);
    max_b = boost::lexical_cast<uint16>(*tok_iter);
    ++tok_iter;

    if (tok_iter != tokens.end()) this->error("snapshot", range_string);

    if (min_gray == min_i && min_gray == min_r && min_gray == min_b &&
        max_gray == max_i && max_gray == max_r && max_gray == max_b) {
      std::cout << "\t--> The color ranges all match, which is not likely.  Double checking the color ranges...\n";

      min_gray=1024;
      max_gray=0;
      min_i=1024;
      max_i=0;
      min_r=1024;
      max_r=0;
      min_b=1024;
      max_b=0;

      int filter = 64;
      std::cout << "Scanning IRB image for autoscaling..." << std::endl;
      DiskImageView<PixelRGB<uint16> > image( color_filename );
      ImageView<PixelRGB<uint16> > stripe( image.cols(), filter );
      TerminalProgressCallback progress("plate.tools.hirise2tif", "");
      for( int y=0; y+filter-1<image.rows(); y+=filter ) {
        progress.report_progress( double(y) / image.rows() );
        stripe = crop( image, 0, y, image.cols(), filter );
        ImageView<PixelRGB<uint16> >::pixel_accessor col = stripe.origin();
        for( int x=0; x+filter-1<image.cols(); x+=filter ) {
          int i=0,r=0,b=0,n=0;
          ImageView<PixelRGB<uint16> >::pixel_accessor row = col;
          for( int yy=0; yy<filter; ++yy ) {
            ImageView<PixelRGB<uint16> >::pixel_accessor pos = row;
            for( int xx=0; xx<filter; ++xx ) {
              PixelRGB<uint16> &val = *pos;
              pos.next_col();
              if( val.r()==0 || val.g()==0 || val.b()==0 ) continue;
              i += val.r();
              r += val.g();
              b += val.b();
              n ++;
            }
            row.next_row();
          }
          col.advance(filter,0);
          if( n < filter*filter/4 ) continue;
          if( i/n > max_i ) max_i = i/n;
          if( i/n < min_i ) min_i = i/n;
          if( r/n > max_r ) max_r = r/n;
          if( r/n < min_r ) min_r = r/n;
          if( b/n > max_b ) max_b = b/n;
          if( b/n < min_b ) min_b = b/n;
        }
      }
      progress.report_finished();

      {
        std::cout << "Scanning red image for autoscaling..." << std::endl;
        DiskImageView<uint16> image( gray_filename );
        ImageView<uint16> stripe( image.cols(), filter );
        TerminalProgressCallback progress("plate.tools.hirise2tif", "");
        for( int y=0; y+filter-1<image.rows(); y+=filter ) {
          progress.report_progress( double(y) / image.rows() );
          stripe = crop( image, 0, y, image.cols(), filter );
          ImageView<uint16>::pixel_accessor col = stripe.origin();
          for( int x=0; x+filter-1<image.cols(); x+=filter ) {
            int r=0,n=0;
            ImageView<uint16>::pixel_accessor row = col;
            for( int yy=0; yy<filter; ++yy ) {
              ImageView<uint16>::pixel_accessor pos = row;
              for( int xx=0; xx<filter; ++xx ) {
                uint16 &val = *pos;
                pos.next_col();
                if( val == 0 ) continue;
                r += val;
                n ++;
              }
              row.next_row();
            }
            col.advance(filter,0);
            if( n < filter*filter/4 ) continue;
            if( r/n > max_gray ) max_gray = r/n;
            if( r/n < min_gray ) min_gray = r/n;
          }
        }
        progress.report_finished();
      }

    }
  }
};


// --------------------------------------------------------------------------
//                     GEOREFERENCE PARSER/BUILDER
// --------------------------------------------------------------------------


/// Parses georeference related arguments and turs them into useful
/// georefernce objects for the red and color images.  If no arguments
/// are supplied, we attempt to fall back to reading the georeference
/// from the files themselves.
class ImageGeorefs {

  void error(std::string arg, std::string const& params) {
    vw_out(ErrorMessage) << "Error parsing arguments for --" << arg << " : " << params << "\n";
    exit(1);
  }

public:

  // Public variables.  You access these directly.
  cartography::GeoReference gray_georef;
  cartography::GeoReference color_georef;

  // Constructor
  ImageGeorefs(std::string const& gray_wkt,
             std::string const& color_wkt,
             std::string const& gray_ullr,
             std::string const& color_ullr,
             std::string const& gray_filename,
             std::string const& color_filename) {

    typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
    boost::char_separator<char> sep(",;");

    // STEP 1 : PARSE GEOREFERENCES

    // If we don't have a full complement of command line options,
    // then we fall back to reading georeferences from the files.
    if (gray_wkt.empty() && color_wkt.empty()) {
      std::cout << "\t-->Reading georeferences from input images.\n";
      read_georeference( gray_georef, gray_filename );
      read_georeference( color_georef, color_filename );

    } else {
      std::cout << "\t--> Parsing georeferences from the command line.\n";

      if (!gray_wkt.empty() && !gray_ullr.empty()) {

        // Initialize the georef objects using the well known text
        // strings supplied on the command line.
        gray_georef.set_wkt(gray_wkt);

        // Parse through the ullr strings to extract upper left and
        // lower right coordinates.
        Vector2 ul_gray, lr_gray;

        // Gray
        tokenizer tokens_gray(gray_ullr, sep);
        tokenizer::iterator tok_iter = tokens_gray.begin();
        if (tok_iter == tokens_gray.end()) this->error("ullr-gray", gray_ullr);
        ul_gray[0] = boost::lexical_cast<float>(*tok_iter);
        ++tok_iter;

        if (tok_iter == tokens_gray.end()) this->error("ullr-gray", gray_ullr);
        ul_gray[1] = boost::lexical_cast<float>(*tok_iter);
        ++tok_iter;

        if (tok_iter == tokens_gray.end()) this->error("ullr-gray", gray_ullr);
        lr_gray[0] = boost::lexical_cast<float>(*tok_iter);
        ++tok_iter;

        if (tok_iter == tokens_gray.end()) this->error("ullr-gray", gray_ullr);
        lr_gray[1] = boost::lexical_cast<float>(*tok_iter);
        ++tok_iter;
        if (tok_iter != tokens_gray.end()) this->error("ullr-gray", gray_ullr);

        // Now use the ullr bounds to construct proper georeference information.
        Matrix3x3 T_gray;
        T_gray.set_identity();
        boost::shared_ptr<DiskImageResource> gray_rsrc(DiskImageResource::open(gray_filename));

        T_gray(0,0) = (lr_gray[0] - ul_gray[0]) / gray_rsrc->cols();
        T_gray(1,1) = (lr_gray[1] - ul_gray[1]) / gray_rsrc->rows();
        T_gray(0,2) = ul_gray[0];
        T_gray(1,2) = ul_gray[1];

        // Check to make sure that the longitude is within the range
        // [-180, 180].
        if (T_gray(0,2) > 180)
          T_gray(0,2) -= 360;
        if (T_gray(0,2) < -180)
          T_gray(0,2) += 360;

        std::cout << "T_GRAY=" << T_gray << "\n";
        gray_georef.set_transform(T_gray);
      }

      if (!color_wkt.empty() && !color_ullr.empty()) {

        // Initialize the georef objects using the well known text
        // strings supplied on the command line.
        color_georef.set_wkt(color_wkt);

        // Parse through the ullr strings to extract upper left and
        // lower right coordinates.
        Vector2 ul_color, lr_color;

        // Color
        tokenizer tokens_color(color_ullr, sep);
        tokenizer::iterator tok_iter = tokens_color.begin();
        if (tok_iter == tokens_color.end()) this->error("ullr-color", color_ullr);
        ul_color[0] = boost::lexical_cast<float>(*tok_iter);
        ++tok_iter;

        if (tok_iter == tokens_color.end()) this->error("ullr-color", color_ullr);
        ul_color[1] = boost::lexical_cast<float>(*tok_iter);
        ++tok_iter;

        if (tok_iter == tokens_color.end()) this->error("ullr-color", color_ullr);
        lr_color[0] = boost::lexical_cast<float>(*tok_iter);
        ++tok_iter;

        if (tok_iter == tokens_color.end()) this->error("ullr-color", color_ullr);
        lr_color[1] = boost::lexical_cast<float>(*tok_iter);
        ++tok_iter;
        if (tok_iter != tokens_color.end()) this->error("ullr-color", color_ullr);

        Matrix3x3 T_color;
        T_color.set_identity();
        boost::shared_ptr<DiskImageResource> color_rsrc(DiskImageResource::open(color_filename));

        T_color(0,0) = (lr_color[0] - ul_color[0]) / color_rsrc->cols();
        T_color(1,1) = (lr_color[1] - ul_color[1]) / color_rsrc->rows();
        T_color(0,2) = ul_color[0];
        T_color(1,2) = ul_color[1];


        // Check to make sure that the longitude is within the range
        // [-180, 180].
        if (T_color(0,2) > 180)
          T_color(0,2) -= 360;
        if (T_color(0,2) < -180)
          T_color(0,2) += 360;

        std::cout << "T_COLOR=" << T_color << "\n";
        color_georef.set_transform(T_color);
      }

    }
  }
};

// --------------------------------------------------------------------------
//                     BIT DEPTH CONVERSION ROUTINES
// --------------------------------------------------------------------------


class uint16_to_rgba8 : public ReturnFixedType<PixelRGBA<uint8> > {
  uint16 minval, maxval;
  double scale;
public:
  uint16_to_rgba8( uint16 minval, uint16 maxval ) : minval(minval), maxval(maxval), scale(255.0/(maxval-minval)) {}
  PixelRGBA<uint8> operator()( uint16 pix ) const {
    if( pix == 0 ) return PixelRGBA<uint8>();
    uint8 val = (uint8) ( scale * (pix - minval) );
    if( pix < minval ) val = 0;
    else if( pix > maxval ) val = 255;
    return PixelRGBA<uint8>( val, val, val, 255 );
  }
};

class rgb16_to_rgba8 : public ReturnFixedType<PixelRGBA<uint8> > {
  uint16 minval_ir, maxval_ir;
  double scale_ir;
  uint16 minval_red, maxval_red;
  double scale_red;
  uint16 minval_bg, maxval_bg;
  double scale_bg;
  uint16 minval_overall, maxval_overall;
  double scale_overall;
  bool m_rgb, m_uniform_strech;
public:
  rgb16_to_rgba8( uint16 minval_ir, uint16 maxval_ir,
                  uint16 minval_red, uint16 maxval_red,
                  uint16 minval_bg, uint16 maxval_bg,
                  bool rgb, bool uniform_stretch)
    : minval_ir(minval_ir), maxval_ir(maxval_ir), scale_ir(255.0/(maxval_ir-minval_ir)),
      minval_red(minval_red), maxval_red(maxval_red), scale_red(255.0/(maxval_red-minval_red)),
      minval_bg(minval_bg), maxval_bg(maxval_bg), scale_bg(255.0/(maxval_bg-minval_bg)),
      m_rgb(rgb), m_uniform_strech(uniform_stretch)
  {
    if (uniform_stretch)  {
      minval_overall = std::min(minval_ir, std::min(minval_red, minval_bg));
      maxval_overall = std::max(maxval_ir, std::max(maxval_red, maxval_bg));
      scale_overall =255.0/(maxval_overall-minval_overall);
      std::cout << "\t--> Using uniform stretch for color image: [ "
                << minval_overall << " " << maxval_overall << " ]\n";
      minval_ir = minval_overall; maxval_ir = maxval_overall; scale_ir = scale_overall;
      minval_red = minval_overall; maxval_red = maxval_overall; scale_red = scale_overall;
      minval_bg = minval_overall; maxval_bg = maxval_overall; scale_bg = scale_overall;
    }

    // If the user has requested RGB color, we repurpose the IR
    // channel to store the synthetic blue color information.  The
    // synthetic blue channel is created by multiplying the bg channel
    // by two and subtracting one third of the red channel, as
    // perscibed in the HiRISE paper entitled: "Color imaging of Mars
    // by the High Resolution Imaging Science Experiment (HiRISE)"
    if (rgb && !uniform_stretch) {
      if (2*minval_bg < minval_red * 3 / 10)
        // Handle the case where the min red is significantly larger than min bg.
        minval_ir = 1;
      else
        minval_ir = 2 * minval_bg - minval_red * 3 / 10;
      maxval_ir = 2 * maxval_bg - maxval_red * 3 / 10;
      scale_ir = 255.0/(maxval_ir-minval_ir);
      std::cout << "\t--> Computed range for synthetic blue channel: [ "
                << minval_ir << " " << maxval_ir << " ]\n";
    }

  }

  PixelRGBA<uint8> operator()( PixelRGB<uint16> const& pix ) const {
    // Check for NULL pixels
    if( pix.r()==0 || pix.g()==0 || pix.b()==0 ) return PixelRGBA<uint8>();

    PixelRGBA<uint8> result;
    if (m_rgb) {
      // This is a little confusing: scale_ir contains that scale

      int16 synthetic_blue = 2 * pix.b() - pix.g() * 3 / 10;
      result = PixelRGBA<uint8>( (uint8)(scale_red*(pix.g()-minval_red)),      // Red
                                 (uint8)(scale_bg*(pix.b()-minval_bg)),        // Blue-Green
                                 (uint8)(scale_ir*(synthetic_blue-minval_ir)), // Synthetic Blue
                                 255 );
      if( synthetic_blue < int16(minval_ir) ) result.b() = 0;
      else if( synthetic_blue > maxval_ir ) result.b() = 255;
      if( pix.g() < minval_red ) result.r() = 0;
      else if( pix.g() > maxval_red ) result.r() = 255;
      if( pix.b() < minval_bg ) result.g() = 0;
      else if( pix.b() > maxval_bg ) result.g() = 255;

    } else {
      result = PixelRGBA<uint8>( (uint8)(scale_ir*(pix.r()-minval_ir)),    // IR
                                 (uint8)(scale_red*(pix.g()-minval_red)),  // Red
                                 (uint8)(scale_bg*(pix.b()-minval_bg)),    // Blue-Green
                                 255 );
      if( pix.r() < minval_ir ) result.r() = 0;
      else if( pix.r() > maxval_ir ) result.r() = 255;
      if( pix.g() < minval_red ) result.g() = 0;
      else if( pix.g() > maxval_red ) result.g() = 255;
      if( pix.b() < minval_bg ) result.b() = 0;
      else if( pix.b() > maxval_bg ) result.b() = 255;
    }
    return result;
  }
};

PixelRGB<uint8> rgba8_to_rgb8( PixelRGBA<uint8> const& rgba ) {
  if( rgba.a() != 255 ) return PixelRGB<uint8>();
  else return PixelRGB<uint8>( rgba );
}


// --------------------------------------------------------------------------
//                            MAIN FUNCTION
// --------------------------------------------------------------------------


int main( int argc, char *argv[] ) {

  std::vector<std::string> image_files;
  std::string output_file_name;
  unsigned cache_size;
  unsigned gdal_cache_size;
  std::string wkt_gray, wkt_color;
  std::string ullr_gray, ullr_color;
  std::string stats_string;

  po::options_description general_options("General Options");
  general_options.add_options()
    ("output-name,o", po::value<std::string>(&output_file_name)->default_value("output"), "Specify the base output filename")
    ("wkt-gray", po::value<std::string>(&wkt_gray),
     "Specify a well-known text string for the red channel input image.")
    ("wkt-color", po::value<std::string>(&wkt_color),
     "Specify a well-known text string for the color input image.")
    ("ullr-gray", po::value<std::string>(&ullr_gray),
     "Specify a well-known text string for the red channel input image.")
    ("ullr-color", po::value<std::string>(&ullr_color),
     "Specify a well-known text string for the color input image.")
    ("stats", po::value<std::string>(&stats_string),
     "where arg = <gray_min>,<gray_max>;<ir_min>,<ir_max>;<red_min>,<red_max>;<bg_min>,<bg_max>")
    ("rgb", "Compute \"true\" color using the RGB algorithm rather than the straight-forward IRB color scheme.")
    ("uniform-stretch", "Stretch all color channels using the same range.")
    ("cache", po::value<unsigned>(&cache_size)->default_value(512), "Cache size, in megabytes")
    ("gdal-cache", po::value<unsigned>(&gdal_cache_size)->default_value(256), "GDAL internal cache size, in megabytes")
    ("help,h", "Display this help message");

  po::options_description output_options("Output Options");
  output_options.add_options()
    ("alpha", "Output an alpha-masked geotiff");

  po::options_description hidden_options("");
  hidden_options.add_options()
    ("input-file", po::value<std::vector<std::string> >(&image_files));

  po::options_description options("Allowed Options");
  options.add(general_options).add(output_options).add(hidden_options);

  po::positional_options_description p;
  p.add("input-file", -1);

  po::variables_map vm;
  po::store( po::command_line_parser( argc, argv ).options(options).positional(p).run(), vm );
  po::notify( vm );

  std::ostringstream usage;
  usage << "Usage: hirise2tif [options] <filename>..." << std::endl << std::endl;
  usage << general_options << std::endl;
  usage << output_options << std::endl;

  if( vm.count("help") ) {
    std::cout << usage.str();
    return 1;
  }

  if( vm.count("input-file") < 1 ) {
    std::cout << "Error: Must specify at least one input file!" << std::endl << std::endl;
    std::cout << usage.str();
    return 1;
  }

  vw_settings().set_system_cache_size( cache_size*1024*1024 );
  DiskImageResourceGDAL::set_gdal_cache_size( gdal_cache_size*1024*1024 );

  // Compute (or parse from input) image statistics.
  ImageStats stats(stats_string, image_files[0], image_files[1]);
  std::cout << "Pixel ranges:\n";
  std::cout << "\t-->GRAY: " << stats.min_gray << " - " << stats.max_gray << std::endl;
  std::cout << "\t-->IR:   " << stats.min_i << " - " << stats.max_i << std::endl;
  std::cout << "\t-->RED:  " << stats.min_r << " - " << stats.max_r << std::endl;
  std::cout << "\t-->BLUE: " << stats.min_b << " - " << stats.max_b << std::endl << std::endl;

  // Compute georeferences
  ImageGeorefs georefs(wkt_gray, wkt_color, ullr_gray, ullr_color,
                       image_files[0], image_files[1]);
  std::cout << "Georeferences:\n";
  std::cout << "GRAY: " << georefs.gray_georef << std::endl;
  std::cout << "COLOR: " << georefs.color_georef << std::endl;

  // Add the transformed input files to the composite
  ImageComposite<PixelRGBA<uint8> > composite;
  GeoReference master_georef = georefs.gray_georef;

  // Set options for JP2 decoding.
  //
  // TODO: These aren't used yet. We don't have a way of passing them into GDAL...
  // DiskImageResourceGDAL::Options jp2_options;
  // jp2_options["JP2KAK_THREADS"] = "2";
  // jp2_options["GDAL_CACHEMAX"]  = "512";

  // Add the grayscale image
  DiskImageResourceGDAL *gray_rsrc = new DiskImageResourceGDAL( image_files[0] );

  composite.insert( per_pixel_filter( DiskImageView<uint16>( gray_rsrc ),
                                      uint16_to_rgba8(stats.min_gray,stats.max_gray) ), 0, 0 );

  // Add the color image
  if (image_files.size() == 2) {
    DiskImageResourceGDAL *color_rsrc = new DiskImageResourceGDAL( image_files[1] );
    Vector2 position =
      master_georef.lonlat_to_pixel( georefs.color_georef.pixel_to_lonlat( Vector2() ) );
    std::cout << "\t--> Adding color image at " << position << std::endl;
    composite.insert( per_pixel_filter( DiskImageView<PixelRGB<uint16> >( color_rsrc ),
                                        rgb16_to_rgba8(stats.min_i,stats.max_i,
                                                       stats.min_r,stats.max_r,
                                                       stats.min_b,stats.max_b,
                                                       vm.count("rgb"),
                                                       vm.count("uniform-stretch")) ),
                      math::impl::_round(position.x()), math::impl::_round(position.y()) );
  }

  // Update the georeference
  BBox2i pixel_bbox = composite.bbox();
  Vector2 offset = master_georef.pixel_to_point( pixel_bbox.min() ) - master_georef.pixel_to_point( Vector2() );
  Matrix3x3 M = master_georef.transform();
  M(0,2) += offset.x();
  M(1,2) += offset.y();
  if( ! master_georef.is_projected() ) {
    if( M(0,2) > 180 ) M(0,2) -= 360;
  }
  master_georef.set_transform( M );

  // Prepare the composite
  composite.set_draft_mode( true );
  composite.prepare();

  // Heuristic: If the image area is bigger than 100K by 100K, die before we trigger an OOM.
  VW_ASSERT(uint64(composite.cols()) * uint64(composite.rows()) < 10000000000ULL,
    LogicErr() << "Composite is WAY too big (" << composite.cols() << " x " << composite.rows() << "). Something's wrong!");

  std::cout << "\t--> Generating output GeoTIFF.  Size = " << composite.cols() << " x " << composite.rows() << "\n";
  DiskImageResourceGDAL::Options gdal_options;
  gdal_options["COMPRESS"] = "LZW";
  gdal_options["BIGTIFF"]  = "YES";

  std::cout << "NAME: " << output_file_name << std::endl;
  ImageFormat format = composite.format();
  if( vm.count("alpha") ) {
    DiskImageResourceGDAL output_resource( output_file_name, format,
                                           Vector2i(256,256), gdal_options );
    write_georeference( output_resource, master_georef );
    std::cout << "\t--> Retaining alpha channel!" << std::endl;
    write_image( output_resource, composite, TerminalProgressCallback("plate.tools.hirise2tif", "") );
  } else {
    format.pixel_format = VW_PIXEL_RGB;
    DiskImageResourceGDAL output_resource( output_file_name, format,
                                           Vector2i(256,256), gdal_options );
    write_georeference( output_resource, master_georef );
    std::cout << "\t--> Replacing alpha with black!" << std::endl;
    write_image( output_resource, per_pixel_filter(composite,&rgba8_to_rgb8),
                 TerminalProgressCallback("plate.tools.hirise2tif", "") );
  }

  return 0;
}
