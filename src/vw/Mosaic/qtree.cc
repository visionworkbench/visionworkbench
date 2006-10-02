#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#ifdef NDEBUG
#undef NDEBUG
#endif

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <vw/Core/Debugging.h>
#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Mosaic/QuadTreeGenerator.h>
using namespace vw;
using namespace vw::mosaic;

int main( int argc, char *argv[] ) {

  ImageQuadTreeGenerator<PixelRGBA<uint8> >::Settings settings;
  settings.crop_images = true;
  settings.output_image_type = "png";
  std::string input_file_name;

  po::options_description desc("Options");
  desc.add_options()
    ("help", "Display this help message")
    ("input-file", po::value<std::string>(&input_file_name), "Explicitly specify the input file")
    ("size", po::value<int>(&settings.patch_size)->default_value(256), "Patch size, in pixels")
    ("overlap", po::value<int>(&settings.patch_overlap)->default_value(0), "Patch overlap, in pixels (must be even)")
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
    std::cout << "Error: Must specify one (and only one) input file!" << std::endl;
    std::cout << desc << std::endl;
    return 1;
  }

  if( vm.count("verbose") ) {
    set_debug_level(VerboseDebugMessage);
  }

  if( settings.patch_size <= 0 ) {
    std::cerr << "Error: The patch size must be a positive number!  (You specified " << settings.patch_size << ".)" << std::endl;
    std::cout << desc << std::endl;
    return 1;
  }

  if( settings.patch_overlap<0 || settings.patch_overlap>=settings.patch_size || settings.patch_overlap%2==1 ) {
    std::cerr << "Error: The patch overlap must be an even number nonnegative number" << std::endl;
    std::cerr << "smaller than the patch size!  (You specified " << settings.patch_overlap << ".)" << std::endl;
    std::cout << desc << std::endl;
    return 1;
  }

  settings.tree_name = input_file_name.substr( 0, input_file_name.rfind( '.' ) );

  try {
    DiskImageView<PixelRGBA<uint8> > input( input_file_name );
    ImageViewPatchSource<DiskImageView<PixelRGBA<uint8> > > source( input );
    ImageQuadTreeGenerator<PixelRGBA<uint8> >( settings, source ).generate();
  }
  catch( Exception& e ) {
    std::cerr << "Error: " << e.what() << std::endl;
  }

  return 0;
}
