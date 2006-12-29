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
#include <vw/Camera/FisheyeModel.h>
#include <vw/Camera/CAHVModel.h>
#include <vw/Camera/CameraTransform.h>
using namespace vw;

int main( int argc, char *argv[] ) {

  std::string input_file_name, output_file_name;
  float focal_length;
  float field_of_view;

  po::options_description desc("Options");
  desc.add_options()
    ("help", "Display this help message")
    ("input-file", po::value<std::string>(&input_file_name), "Explicitly specify the input file")
    ("output-file,o", po::value<std::string>(&output_file_name)->default_value("output.png"), "Specify the output file")
    ("focal-length,f", po::value<float>(&focal_length), "Focal length of the lens (in meters)")
    ("fov", po::value<float>(&field_of_view), "Field of view of the fisheye lens (in degrees)")
    ("normalize", "Normalize the output image")
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

  if( vm.count("focal_length") != 1 ) {
    std::cout << "Error: Must specify the focal length!" << std::endl;
    std::cout << desc << std::endl;
    return 1;    
  }

  if( vm.count("fov") != 1 ) {
    std::cout << "Error: Must specify the field of view!" << std::endl;
    std::cout << desc << std::endl;
    return 1;    
  }

  if( vm.count("verbose") ) {
    set_debug_level(VerboseDebugMessage);
  }

  try {
    ImageView<PixelRGB<float> > image;
    read_image( image, input_file_name );

    if( vm.count("normalize") ) {
      image = normalize(image);
    }

    write_image( output_file_name, image );
  }
  catch( Exception& e ) {
    std::cerr << "Error: " << e.what() << std::endl;
  }

  return 0;
}
