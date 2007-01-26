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

#include <vw/Core/Debugging.h>
#include <vw/Image/BlockRasterize.h>
#include <vw/FileIO/DiskImageView.h>
#include <vw/Mosaic/QuadTreeGenerator.h>
#include <vw/Mosaic/ImageComposite.h>

template <class PixelT>
void do_blend( std::string const& mosaic_name, std::string const& file_type, bool draft, bool qtree, int patch_size, int patch_overlap ) {
  vw::mosaic::ImageComposite<PixelT> composite;
    
  if( draft ) composite.set_draft_mode( true );
    
  std::map<std::string,fs::path> image_files;
  std::map<std::string,fs::path> offset_files;
  fs::path source_dir_path( mosaic_name, fs::native );
  fs::directory_iterator pi( source_dir_path ), pend;
  for( ; pi != pend; ++pi ) {
    if( extension(*pi) == ".offset" )
      offset_files[basename(*pi)] = *pi;
    else image_files[basename(*pi)] = *pi;
  }
  std::map<std::string,fs::path>::iterator ofi=offset_files.begin(), ofend=offset_files.end();
  for( ; ofi != ofend; ++ofi ) {
    std::map<std::string,fs::path>::iterator ifi = image_files.find( ofi->first );
    if( ifi != image_files.end() ) {
      fs::ifstream offset( ofi->second );
      int x, y;
      offset >> x >> y;
      std::cout << "Importing image file " << ifi->second.string() << " at offet (" << x << "," << y << ")" << std::endl;
      composite.insert( vw::DiskImageView<PixelT>( ifi->second.string(), false ), x, y );
    }
  }
    
  vw::vw_out(vw::InfoMessage) << "Preparing the composite..." << std::endl;
  composite.prepare();
  if( qtree ) {
    vw::vw_out(vw::InfoMessage) << "Preparing the quadtree..." << std::endl;
    vw::mosaic::ImageQuadTreeGenerator<PixelT > quadtree( mosaic_name, composite );
    quadtree.set_output_image_file_type( file_type );
    quadtree.set_patch_size( patch_size );
    quadtree.set_patch_overlap( patch_overlap );
    vw::vw_out(vw::InfoMessage) << "Generating..." << std::endl;
    quadtree.generate();
    vw::vw_out(vw::InfoMessage) << "Done!" << std::endl;
  }
  else {
    vw::vw_out(vw::InfoMessage) << "Blending..." << std::endl;
    vw::ImageView<PixelT> im = block_rasterize( composite );
    write_image( mosaic_name+".blend."+file_type, im );
    vw::vw_out(vw::InfoMessage) << "Done!" << std::endl;
  }
}

int main( int argc, char *argv[] ) {
  try {
    std::string mosaic_name, file_type;
    int patch_size, patch_overlap;
    unsigned cache_size;
    
    po::options_description desc("Options");
    desc.add_options()
      ("help", "Display this help message")
      ("input-dir", po::value<std::string>(&mosaic_name), "Explicitly specify the input directory")
      ("file-type", po::value<std::string>(&file_type)->default_value("png"), "Output file type")
      ("size", po::value<int>(&patch_size)->default_value(256), "Patch size, in pixels")
      ("overlap", po::value<int>(&patch_overlap)->default_value(0), "Patch overlap, in pixels (must be even)")
      ("cache", po::value<unsigned>(&cache_size)->default_value(1024), "Cache size, in megabytes")
      ("draft", "Draft mode (no blending)")
      ("qtree", "Output in quadtree format")
      ("grayscale", "Process in grayscale only")
      ("verbose", "Verbose output");
    po::positional_options_description p;
    p.add("input-dir", 1);
    
    po::variables_map vm;
    po::store( po::command_line_parser( argc, argv ).options(desc).positional(p).run(), vm );
    po::notify( vm );
    
    if( vm.count("help") ) {
      std::cout << desc << std::endl;
      return 1;
    }
    
    if( vm.count("input-dir") != 1 ) {
      std::cout << "Error: Must specify one (and only one) input directory!" << std::endl;
      std::cout << desc << std::endl;
      return 1;
    }
    
    if( vm.count("verbose") ) {
      vw::set_debug_level(vw::VerboseDebugMessage);
    }

    if( patch_size <= 0 ) {
      std::cerr << "Error: The patch size must be a positive number!  (You specified " << patch_size << ".)" << std::endl;
      std::cout << desc << std::endl;
      return 1;
    }
    
    if( patch_overlap<0 || patch_overlap>=patch_size || patch_overlap%2==1 ) {
      std::cerr << "Error: The patch overlap must be an even number nonnegative number" << std::endl;
      std::cerr << "smaller than the patch size!  (You specified " << patch_overlap << ".)" << std::endl;
      std::cout << desc << std::endl;
      return 1;
    }
    
    vw::Cache::system_cache().resize( cache_size*1024*1024 );

    if( vm.count("grayscale") ) {
      do_blend<vw::PixelGrayA<float> >( mosaic_name, file_type, vm.count("draft"), vm.count("qtree"), patch_size, patch_overlap );
    }
    else {
      do_blend<vw::PixelRGBA<float> >( mosaic_name, file_type, vm.count("draft"), vm.count("qtree"), patch_size, patch_overlap );
    }

  }
  catch( std::exception &err ) {
    vw::vw_out(vw::ErrorMessage) << "Error: " << err.what() << std::endl << "Aborting!" << std::endl;
    return 1;
  }
  return 0;
}
