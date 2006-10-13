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
#include <vw/FileIO/DiskImageView.h>
#include <vw/Mosaic/QuadTreeGenerator.h>
#include <vw/Mosaic/ImageComposite.h>

int main( int argc, char *argv[] ) {
  try {
    int patch_size, patch_overlap;
    std::string mosaic_name;
    unsigned cache_size;
    
    po::options_description desc("Options");
    desc.add_options()
      ("help", "Display this help message")
      ("input-dir", po::value<std::string>(&mosaic_name), "Explicitly specify the input directory")
      ("size", po::value<int>(&patch_size)->default_value(256), "Patch size, in pixels")
      ("overlap", po::value<int>(&patch_overlap)->default_value(0), "Patch overlap, in pixels (must be even)")
      ("cache", po::value<unsigned>(&cache_size)->default_value(1024), "Cache size, in megabytes")
      ("draft", "Draft mode (no blending)")
      ("qtree", "Output in quadtree format")
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

    vw::mosaic::ImageComposite composite;
    
    if( vm.count("draft") > 0 ) {
      composite.set_draft_mode( true );
    }
    
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
        composite.insert( vw::DiskImageView<vw::PixelRGBA<float> >( ifi->second.string(), false ), x, y );
      }
    }
    
    vw::vw_out(vw::InfoMessage) << "Preparing the composite..." << std::endl;
    composite.prepare();
    if( vm.count("qtree") ) {
      vw::vw_out(vw::InfoMessage) << "Preparing the quadtree..." << std::endl;
      vw::mosaic::ImageQuadTreeGenerator<vw::PixelRGBA<float> > quadtree( mosaic_name, composite );
      quadtree.set_patch_size( patch_size );
      quadtree.set_patch_overlap( patch_overlap );
      vw::vw_out(vw::InfoMessage) << "Generating..." << std::endl;
      quadtree.generate();
      vw::vw_out(vw::InfoMessage) << "Done!" << std::endl;
    }
    else {
      vw::vw_out(vw::InfoMessage) << "Blending..." << std::endl;
      vw::ImageView<vw::PixelRGBA<float> > im = composite;
      write_image( mosaic_name+".blend.png", im );
      vw::vw_out(vw::InfoMessage) << "Done!" << std::endl;
    }
  }
  catch( std::exception &err ) {
    vw::vw_out(vw::ErrorMessage) << "Error: " << err.what() << std::endl << "Aborting!" << std::endl;
    return 1;
  }
  return 0;
}
