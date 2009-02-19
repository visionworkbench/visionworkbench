// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <string>
#include <fstream>
#include <vector>
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

using namespace std;
using namespace vw;

string mosaic_name;
string file_type;
unsigned int cache_size;
int tile_size;
bool draft;
bool qtree;

template <class PixelT>
void do_blend() {
  mosaic::ImageComposite<PixelT> composite;
  if( draft ) composite.set_draft_mode( true );

  map<string,fs::path> image_files;
  map<string,fs::path> offset_files;
  fs::path source_dir_path( mosaic_name, fs::native );
  fs::directory_iterator pi( source_dir_path ), pend;
  for( ; pi != pend; ++pi ) {
    if( extension(*pi) == ".offset" )
      offset_files[basename(*pi)] = *pi;
    else image_files[basename(*pi)] = *pi;
  }
  map<string,fs::path>::iterator ofi=offset_files.begin(), ofend=offset_files.end();
  for( ; ofi != ofend; ++ofi ) {
    map<string,fs::path>::iterator ifi = image_files.find( ofi->first );
    if( ifi != image_files.end() ) {
      fs::ifstream offset( ofi->second );
      int x, y;
      offset >> x >> y;
      cout << "Importing image file " << ifi->second.string() << " at offet (" << x << "," << y << ")" << endl;
      composite.insert( DiskImageView<PixelT>( ifi->second.string(), false ), x, y );
    }
  }

  vw_out(InfoMessage) << "Preparing the composite..." << endl;
  composite.prepare();
  if( qtree ) {
    vw_out(InfoMessage) << "Preparing the quadtree..." << endl;
    mosaic::QuadTreeGenerator quadtree( composite, mosaic_name );
    quadtree.set_file_type( file_type );
    quadtree.set_tile_size( tile_size );
    vw_out(InfoMessage) << "Generating..." << endl;
    quadtree.generate();
    vw_out(InfoMessage) << "Done!" << endl;
  }
  else {
    vw_out(InfoMessage) << "Blending..." << endl;
    write_image( mosaic_name+".blend."+file_type, composite );
    vw_out(InfoMessage) << "Done!" << endl;
  }
}

int main( int argc, char *argv[] ) {
  try {
    po::options_description desc("Options");
    desc.add_options()
      ("help", "Display this help message")
      ("input-dir", po::value<string>(&mosaic_name), "Explicitly specify the input directory")
      ("file-type", po::value<string>(&file_type)->default_value("png"), "Output file type")
      ("tile-size", po::value<int>(&tile_size)->default_value(256), "Tile size, in pixels")
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
      cout << desc << endl;
      return 1;
    }

    if( vm.count("input-dir") != 1 ) {
      cout << "Error: Must specify one (and only one) input directory!" << endl;
      cout << desc << endl;
      return 1;
    }

    if( vm.count("verbose") ) {
      set_debug_level(VerboseDebugMessage);
    }

    if( vm.count("draft") ) draft = true; else draft = false;
    if( vm.count("qtree") ) qtree = true; else qtree = false;

    if( tile_size <= 0 ) {
      cerr << "Error: The tile size must be a positive number!  (You specified " << tile_size << ".)" << endl;
      cout << desc << endl;
      return 1;
    }

    vw_system_cache().resize( cache_size*1024*1024 );

    if( vm.count("grayscale") ) {
      do_blend<PixelGrayA<float> >();
    }
    else {
      do_blend<PixelRGBA<float> >();
    }

  }
  catch( exception &err ) {
    vw_out(ErrorMessage) << "Error: " << err.what() << endl << "Aborting!" << endl;
    return 1;
  }
  return 0;
}
