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
using namespace vw;

int main( int argc, char *argv[] ) {

  std::string input_file_name, output_file_name;
  float rotate_angle;
  float noise_level;

  po::options_description desc("Options");
  desc.add_options()
    ("help", "Display this help message")
    ("input-file", po::value<std::string>(&input_file_name), "Explicitly specify the input file")
    ("output-file,o", po::value<std::string>(&output_file_name)->default_value("output.png"), "Specify the output file")
    ("rotate", po::value<float>(&rotate_angle), "Rotate by the given angle (degrees)")
    ("noise", po::value<float>(&noise_level), "Additive noise (percentage)")
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

  if( vm.count("verbose") ) {
    set_debug_level(VerboseDebugMessage);
  }

  try {
    ImageView<PixelGray<float> > image;
    read_image( image, input_file_name );

    if( vm.count("rotate") ) {
      Matrix3x3 R;
      R.set_identity();
      float theta = - rotate_angle * (M_PI/180);
      R(0,0) = cos(theta);     R(0,1) = sin(theta);
      R(1,0) = -sin(theta);    R(1,1) = cos(theta);
      Matrix3x3 T;
      T.set_identity();
      T(0,2) = 0.5*image.cols();
      T(1,2) = 0.5*image.rows();
      image = transform( image, HomographyTransform(T*R*inverse(T)), ConstantEdgeExtend() );
    }

    if( vm.count("noise") ) {
      ImageView<PixelGray<float> >::iterator i=image.begin(), end=image.end();
      while( i != end ) {
        *(i++) += PixelGray<float>( noise_level * ((random()%2000000)/1000000.0-1.0) );
      }
    }

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
