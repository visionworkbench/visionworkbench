// __BEGIN_LICENSE__
// 
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
// 
// Copyright 2006 Carnegie Mellon University. All rights reserved.
// 
// This software is distributed under the NASA Open Source Agreement
// (NOSA), version 1.3.  The NOSA has been approved by the Open Source
// Initiative.  See the file COPYING at the top of the distribution
// directory tree for the complete NOSA document.
// 
// THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
// KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
// LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
// SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
// A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
// THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
// DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
// 
// __END_LICENSE__

#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif


#include <string>
#include <fstream>
#include <list>
#include <sys/types.h>
#include <sys/stat.h>

#include <boost/operators.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/fstream.hpp>
namespace fs = boost::filesystem;

#include <vw/Core.h>
#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Mosaic.h>
#include <vw/Cartography.h>
#include <vw/Math.h>
using namespace vw;
using namespace vw::math;
using namespace vw::cartography;


//  mask_zero_pixels()
//
// Creates an alpha channel based on pixels with a value of zero.
template <class PixelT> struct AlphaTypeFromPixelType { typedef PixelT type; };
template<class ChannelT> struct AlphaTypeFromPixelType<PixelGray<ChannelT> > { typedef PixelGrayA<ChannelT> type; };
template<class ChannelT> struct AlphaTypeFromPixelType<PixelRGB<ChannelT> > { typedef PixelRGBA<ChannelT> type; };

struct MaskZeroPixelFunc: public vw::UnaryReturnTemplateType<AlphaTypeFromPixelType> {
  template <class PixelT>
  typename AlphaTypeFromPixelType<PixelT>::type operator() (PixelT const& pix) const {
    typedef typename AlphaTypeFromPixelType<PixelT>::type result_type;
    if (pix[0] == 0) 
      return result_type();  // Mask pixel
    else
      return result_type(pix);
  }
};

template <class ViewT>
vw::UnaryPerPixelView<ViewT, MaskZeroPixelFunc> 
mask_zero_pixels(vw::ImageViewBase<ViewT> const& view) {
  return vw::per_pixel_filter(view.impl(), MaskZeroPixelFunc());
}

// do_blend()
//
template <class PixelT>
void do_blend( std::vector<std::string> const& image_files, std::string const& mosaic_name, std::string const& output_file_type, bool draft, bool qtree, int patch_size, int patch_overlap, int draw_order_offset, int max_lod_pixels ) {
  vw::mosaic::ImageComposite<PixelT> composite;
    
  if( draft ) composite.set_draft_mode( true );

  double smallest_x_scale = vw::ScalarTypeLimits<float>::highest();
  double smallest_y_scale = vw::ScalarTypeLimits<float>::highest();

  // First pass, read georeferencing information and build an output
  // georef.
  for(unsigned i = 0; i < image_files.size(); ++i) {
    std::cout << "Adding file " << image_files[i] << std::endl;
    DiskImageResourceGDAL file_resource( image_files[i] );
    GeoReference input_georef;
    file_resource.read_metadata( input_georef );
    std::cout << "\t" << input_georef.transform() << "\n";

    Matrix3x3 affine = input_georef.transform();
    if (fabs(affine(0,0)) < smallest_x_scale || fabs(affine(1,1)) < smallest_y_scale) {
      smallest_x_scale = affine(0,0);
      smallest_y_scale = affine(1,1);
    }

    if( input_georef.transform() == identity_matrix<3>() ) {
      vw_out(InfoMessage) << "No georeferencing info found for image: \"" << image_files[i] << "\".  Aborting.\n";
      exit(0);
    }
  }

  // Adopt the scale of the highest resolution being composited.
  Matrix3x3 output_affine = identity_matrix<3>();
  output_affine(0,0) = smallest_x_scale;
  output_affine(1,1) = smallest_y_scale;
  std::cout << "\nScale for mosaic: [" << smallest_x_scale << " " << smallest_y_scale << "] units/pixel\n";

  // Second pass: add files to the image composite.
  for(unsigned i = 0; i < image_files.size(); ++i) {
    DiskImageResourceGDAL file_resource( image_files[i] );
    GeoReference input_georef;
    file_resource.read_metadata( input_georef );
    Matrix3x3 affine = input_georef.transform();

    DiskImageView<PixelT> source( image_files[i] );
    //     HomographyTransform trans(inverse(output_affine)*affine);
    //     BBox2i bbox = trans.forward_bbox( BBox2i(0,0,source.cols(),source.rows()) );
    //    composite.insert( crop( transform( mask_zero_pixels(source), trans ), bbox ), bbox.min().x(), bbox.min().y() );
    composite.insert( mask_zero_pixels(source), affine(0,2), affine(1,2) );
  }

  vw::vw_out(vw::InfoMessage) << "Preparing the composite..." << std::endl;
  composite.prepare(TerminalProgressCallback() );
  if( qtree ) {
    vw::vw_out(vw::InfoMessage) << "Preparing the quadtree..." << std::endl;
    BBox2 ll_bbox( -30,-30,30.0,30.0*(float)(composite.rows())/(float)(composite.cols()) );
    vw_out(0) << "\tOverlay will appear in this bounding box in Google Earth: " << ll_bbox << "\n";
    vw_out(0) << "\tComposite bbox: " << composite.source_data_bbox() << "\n";
    vw::mosaic::KMLQuadTreeGenerator<PixelT > quadtree( mosaic_name, composite, ll_bbox );
    quadtree.set_output_image_file_type( output_file_type );
    quadtree.set_patch_size( patch_size );
    quadtree.set_patch_overlap( patch_overlap );
    quadtree.set_draw_order_offset( draw_order_offset );
    quadtree.set_max_lod_pixels(max_lod_pixels);
    vw::vw_out(vw::InfoMessage) << "Generating..." << std::endl;
    quadtree.generate(TerminalProgressCallback());
    vw::vw_out(vw::InfoMessage) << "Done!" << std::endl;
  } else {
    std::string mosaic_filename = mosaic_name+".blend."+output_file_type;
    vw::vw_out(vw::InfoMessage) << "Blending..." << std::endl;
    //    vw::ImageViewRef<PixelT> im = block_rasterize( composite );
    write_image( mosaic_filename, composite, TerminalProgressCallback() );
    vw::vw_out(vw::InfoMessage) << "Done!" << std::endl;
  }
}

int main( int argc, char *argv[] ) {
  try {
    std::string mosaic_name, input_file_type, output_file_type;
    std::vector<std::string> image_files;
    int patch_size, patch_overlap;
    unsigned cache_size;
    int draw_order_offset;
    int max_lod_pixels;

    po::options_description general_options("Options");
    general_options.add_options()
      ("help", "Display this help message")
      ("mosaic-name,o", po::value<std::string>(&mosaic_name)->default_value("mosaic"), "Explicitly specify the input directory")
      ("input-file-type", po::value<std::string>(&input_file_type)->default_value("tif"), "Output file type")
      ("output-file-type,t", po::value<std::string>(&output_file_type)->default_value("tif"), "Output file type")
      ("size", po::value<int>(&patch_size)->default_value(256), "Patch size, in pixels")
      ("overlap", po::value<int>(&patch_overlap)->default_value(0), "Patch overlap, in pixels (must be even)")
      ("cache", po::value<unsigned>(&cache_size)->default_value(1024), "Cache size, in megabytes")
      ("draft", "Draft mode (no blending)")
      ("qtree", "Output in quadtree format")
      ("grayscale", "Process in grayscale only")
      ("max-lod-pixels", po::value<int>(&max_lod_pixels)->default_value(1024), "Max LoD in pixels, or -1 for none")
      ("draw-order-offset", po::value<int>(&draw_order_offset)->default_value(100), "Set an offset for the KML <drawOrder> tag for this overlay")
      ("verbose", "Verbose output");

    po::options_description hidden_options("");
    hidden_options.add_options()
      ("input-files", po::value<std::vector<std::string> >(&image_files));
    
    po::options_description options("Allowed Options");
    options.add(general_options).add(hidden_options);

    po::positional_options_description p;
    p.add("input-files", -1);

    po::variables_map vm;
    po::store( po::command_line_parser( argc, argv ).options(options).positional(p).run(), vm );
    po::notify( vm );

    std::ostringstream usage;
    usage << "Usage: geoblend [options] <filename>..." << std::endl << std::endl;
    usage << general_options << std::endl;

    if( vm.count("help") ) {
      std::cout << usage << std::endl;
      return 1;
    }
    
    if( vm.count("input-files") < 1 ) {
      std::cout << "Error: Must specify at least one input file!" << std::endl << std::endl;
      std::cout << usage.str();
      return 1;
    }
    
    if( vm.count("verbose") ) {
      vw::set_debug_level(vw::VerboseDebugMessage);
    }

    if( patch_size <= 0 ) {
      std::cerr << "Error: The patch size must be a positive number!  (You specified " << patch_size << ".)" << std::endl;
      std::cout << usage << std::endl;
      return 1;
    }
    
    if( patch_overlap<0 || patch_overlap>=patch_size || patch_overlap%2==1 ) {
      std::cerr << "Error: The patch overlap must be an even number nonnegative number" << std::endl;
      std::cerr << "smaller than the patch size!  (You specified " << patch_overlap << ".)" << std::endl;
      std::cout << usage << std::endl;
      return 1;
    }
    
    vw::Cache::system_cache().resize( cache_size*1024*1024 );

    if( vm.count("grayscale") ) {
      do_blend<vw::PixelGrayA<float> >( image_files, mosaic_name, output_file_type, vm.count("draft"), vm.count("qtree"), patch_size, patch_overlap, draw_order_offset, max_lod_pixels );
    }
    else {
      do_blend<vw::PixelRGBA<float> >( image_files, mosaic_name, output_file_type, vm.count("draft"), vm.count("qtree"), patch_size, patch_overlap, draw_order_offset, max_lod_pixels );
    }

  }
  catch( std::exception &err ) {
    vw::vw_out(vw::ErrorMessage) << "Error: " << err.what() << std::endl << "Aborting!" << std::endl;
    return 1;
  }
  return 0;
}
