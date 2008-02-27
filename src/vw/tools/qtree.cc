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

std::string intToString(int i)
{
  std::ostringstream oss;
  oss << i;
  return oss.str();
}

template <class PixelT>
class FlatQuadTreeGenerator : public vw::mosaic::ImageQuadTreeGenerator<PixelT>
{
  typedef vw::mosaic::ImageQuadTreeGenerator<PixelT> base_type;

public:
  template <class ImageT>
  FlatQuadTreeGenerator( std::string const& tree_name, ImageViewBase<ImageT> const& source ) :
    ImageQuadTreeGenerator<PixelT>( tree_name, source )
  {}

  virtual std::string compute_image_path( std::string const& name ) const {
    Vector2i pos(0,0);
    for ( int i=1; i<(int)name.length(); ++i ) {
      pos *= 2;
      switch (name[i]) {
      case '0': pos += Vector2i(0,1); break;
      case '1': pos += Vector2i(1,1); break;
      case '2': pos += Vector2i(0,0); break;
      case '3': pos += Vector2i(1,0); break;
      default:
	assert(0); // never reach this point
      }
    }

    fs::path path( base_type::m_base_dir, fs::native );
    path /= intToString( name.length() - 1 );
    path /= intToString( pos.x() );
    path /= intToString( pos.y() );
    
    return path.native_file_string();
  }
};

int main( int argc, char *argv[] ) {

  std::string input_file_name;
  std::string output_image_file_type;
  int patch_size;
  int patch_overlap;

  po::options_description desc("Options");
  desc.add_options()
    ("help", "Display this help message")
    ("input-file", po::value<std::string>(&input_file_name), "Explicitly specify the input file")
    ("size", po::value<int>(&patch_size)->default_value(256), "Patch size, in pixels")
    ("overlap", po::value<int>(&patch_overlap)->default_value(0), "Patch overlap, in pixels (must be even)")
    ("file-type", po::value<std::string>(&output_image_file_type)->default_value("png"), "Output image file type")
    ("pad", "Do not crop output images to available data;  pad images to full quadtree size")
    ("tms", "Output TMS-compatible tile paths (but try image2tms for full TMS support)")
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

  std::string tree_name = input_file_name.substr( 0, input_file_name.rfind( '.' ) );

  try {
    DiskImageView<PixelRGBA<uint8> > input( input_file_name );
    ImageQuadTreeGenerator<PixelRGBA<uint8> >* quadtree;
    if ( vm.count("tms") ) {
      // pad image to have dimension 2**i x 2**i for integer i, with the valid
      // data anchored at the lower left corner rather than the upper left
      int max_dim = std::max( input.cols(), input.rows() );
      int pad_dim = 1 << ((int) ceil( log2( max_dim ) ));
      ImageView<PixelRGBA<uint8> > padded( pad_dim, pad_dim );
      crop( padded, 0, pad_dim - input.rows(), input.cols(), input.rows() ) = input;
      // use the padded image and flat tile path naming
      quadtree = new FlatQuadTreeGenerator<PixelRGBA<uint8> >( tree_name, padded );
    } else {
      quadtree = new ImageQuadTreeGenerator<PixelRGBA<uint8> >( tree_name, input );
    }
    quadtree->set_patch_size( patch_size );
    quadtree->set_patch_overlap( patch_overlap );
    quadtree->set_output_image_file_type( output_image_file_type );
    if( vm.count("pad") ) {
      quadtree->set_crop_images( false );
    }
    quadtree->generate();
  }
  catch( Exception& e ) {
    std::cerr << "Error: " << e.what() << std::endl;
  }

  return 0;
}
