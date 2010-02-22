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
#include <vw/Cartography/FileIO.h>
#include <vw/Mosaic/ImageComposite.h>
#include <vw/Mosaic/QuadTreeGenerator.h>
#include <vw/Mosaic/KMLQuadTreeConfig.h>

using namespace vw;
using namespace vw::cartography;
using namespace vw::mosaic;

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
  uint16 minval_r, maxval_r;
  double scale_r;
  uint16 minval_g, maxval_g;
  double scale_g;
  uint16 minval_b, maxval_b;
  double scale_b;
public:
  rgb16_to_rgba8( uint16 minval_r, uint16 maxval_r, uint16 minval_g, uint16 maxval_g, uint16 minval_b, uint16 maxval_b )
    : minval_r(minval_r), maxval_r(maxval_r), scale_r(255.0/(maxval_r-minval_r)),
      minval_g(minval_g), maxval_g(maxval_g), scale_g(255.0/(maxval_g-minval_g)),
      minval_b(minval_b), maxval_b(maxval_b), scale_b(255.0/(maxval_b-minval_b))
  {}
  PixelRGBA<uint8> operator()( PixelRGB<uint16> const& pix ) const {
  if( pix.r()==0 || pix.g()==0 || pix.b()==0 ) return PixelRGBA<uint8>();
  PixelRGBA<uint8> result( (uint8)(scale_r*(pix.r()-minval_r)), (uint8)(scale_g*(pix.g()-minval_g)), (uint8)(scale_b*(pix.b()-minval_b)), 255 );
  if( pix.r() < minval_r ) result.r() = 0;
  else if( pix.r() > maxval_r ) result.r() = 255;
  if( pix.g() < minval_g ) result.g() = 0;
  else if( pix.g() > maxval_g ) result.g() = 255;
  if( pix.b() < minval_b ) result.b() = 0;
  else if( pix.b() > maxval_b ) result.b() = 255;
  return result;
  }
};

PixelRGB<uint8> rgba8_to_rgb8( PixelRGBA<uint8> const& rgba ) {
  if( rgba.a() != 255 ) return PixelRGB<uint8>();
  else return PixelRGB<uint8>( rgba );
}


int main( int argc, char *argv[] ) {

  std::vector<std::string> image_files;
  std::string output_file_name;
  std::string output_file_type;
  float jpeg_quality;
  int png_compression;
  unsigned cache_size;
  int draw_order_offset;
  int max_lod_pixels;
  int levels=0;
  int aspect_ratio=0;

  po::options_description general_options("General Options");
  general_options.add_options()
    ("output-name,o", po::value<std::string>(&output_file_name)->default_value("output"), "Specify the base output filename")
    ("quiet,q", "Quiet output")
    ("verbose,v", "Verbose output")
    ("cache", po::value<unsigned>(&cache_size)->default_value(512), "Cache size, in megabytes")
    ("help", "Display this help message");

  po::options_description output_options("Output Options");
  output_options.add_options()
    ("file-type", po::value<std::string>(&output_file_type)->default_value("auto"), "Output file type")
    ("jpeg-quality", po::value<float>(&jpeg_quality)->default_value(0.95), "JPEG quality factor (0.0 to 1.0)")
    ("png-compression", po::value<int>(&png_compression)->default_value(3), "PNG compression level (0 to 9)")
    ("draw-order-offset", po::value<int>(&draw_order_offset)->default_value(0), "Set an offset for the KML <drawOrder> tag for this overlay")
    ("max-lod-pixels", po::value<int>(&max_lod_pixels)->default_value(1024), "Max LoD in pixels, or -1 for none (kml only)")
    ("levels", po::value<int>(&levels), "Number of levels in the global quadtree")
    ("alpha", "Output an alpha-masked geotiff")
    ("aspect-ratio", po::value<int>(&aspect_ratio)->default_value(1),"Pixel aspect ratio (for polar overlays; should be a power of two)");

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
  usage << "Usage: image2kml [options] <filename>..." << std::endl << std::endl;
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

  TerminalProgressCallback tpc("plate.tools.hirise2tif", "");
  const ProgressCallback *progress = &tpc;
  if( vm.count("verbose") ) {
    set_debug_level(VerboseDebugMessage);
    progress = &ProgressCallback::dummy_instance();
  }
  else if( vm.count("quiet") ) {
    set_debug_level(WarningMessage);
  }

  DiskImageResourceJPEG::set_default_quality( jpeg_quality );
  DiskImageResourcePNG::set_default_compression_level( png_compression );
  vw_settings().set_system_cache_size( cache_size*1024*1024 );

  uint16 min_gray=1024, max_gray=0, min_i=1024, max_i=0, min_r=1024, max_r=0, min_b=1024, max_b=0;

  int filter = 128;
  if( image_files.size() > 1 ) {
    std::cout << "Scanning IRB image for autoscaling..." << std::endl;
    DiskImageView<PixelRGB<uint16> > image( image_files[1] );
    ImageView<PixelRGB<uint16> > stripe( image.cols(), filter );
    TerminalProgressCallback progress("plate.tools.hirise2tif", "");
    for( int y=0; y+filter-1<image.rows(); y+=filter ) {
      progress.report_progress( (double)y / image.rows() );
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
  }
  {
    std::cout << "Scanning red image for autoscaling..." << std::endl;
    DiskImageView<uint16> image( image_files[0] );
    ImageView<uint16> stripe( image.cols(), filter );
    TerminalProgressCallback progress;
    for( int y=0; y+filter-1<image.rows(); y+=filter ) {
      progress.report_progress( (double)y / image.rows() );
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

  std::cout << "Pixel ranges:\n";
  std::cout << "\tGRAY: " << min_gray << " - " << max_gray << std::endl;
  std::cout << "\tIR:   " << min_i << " - " << max_i << std::endl;
  std::cout << "\tRED:  " << min_r << " - " << max_r << std::endl;
  std::cout << "\tBLUE: " << min_b << " - " << max_b << std::endl;

  // Add the transformed input files to the composite
  ImageComposite<PixelRGBA<uint8> > composite;
  GeoReference master_georef;
  for( unsigned i=0; i<image_files.size(); ++i ) {
    DiskImageResourceGDAL *r = new DiskImageResourceGDAL( image_files[i] );
    GeoReference georef;
    read_georeference( georef, *r );
    if( i==0 ) master_georef = georef;
    Vector2 position = master_georef.lonlat_to_pixel( georef.pixel_to_lonlat( Vector2() ) );
    std::cout << position << std::endl;
    if( r->pixel_format() == VW_PIXEL_GRAY ) {
      composite.insert( per_pixel_filter( DiskImageView<uint16>( r ), uint16_to_rgba8(min_gray,max_gray) ),
                        math::impl::_round(position.x()), math::impl::_round(position.y()) );
    }
    else {
      composite.insert( per_pixel_filter( DiskImageView<PixelRGB<uint16> >( r ), rgb16_to_rgba8(min_i,max_i,min_r,max_r,min_b,max_b) ),
                        math::impl::_round(position.x()), math::impl::_round(position.y()) );
    }
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

  std::cout << "Generating output GeoTIFF..." << std::endl;
  DiskImageResourceGDAL::Options gdal_options;
  gdal_options["COMPRESS"] = "NONE";
  gdal_options["BIGTIFF"]  = "YES";

  std::cout << "NAME: " << output_file_name << std::endl;
  ImageFormat format = composite.format();
  if( vm.count("alpha") ) {
    DiskImageResourceGDAL output_resource( output_file_name, format, Vector2i(256,256), gdal_options );
    write_georeference( output_resource, master_georef );
    std::cout << "Retaining alpha channel!" << std::endl;
    write_image( output_resource, composite, TerminalProgressCallback() );
  } else {
    format.pixel_format = VW_PIXEL_RGB;
    DiskImageResourceGDAL output_resource( output_file_name, format, Vector2i(256,256), gdal_options );
    write_georeference( output_resource, master_georef );
    std::cout << "Replacing alpha with black!" << std::endl;
    write_image( output_resource, per_pixel_filter(composite,&rgba8_to_rgb8), TerminalProgressCallback() );
  }

  return 0;
}
