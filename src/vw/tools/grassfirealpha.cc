#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <vw/FileIO.h>
#include <vw/Image.h>
#include <vw/Math.h>
#include <vw/Cartography.h>
using namespace vw;

// Function for highlighting spots of data
template<class PixelT>
class NotNoDataFunctor {
  typedef typename CompoundChannelType<PixelT>::type channel_type;
  channel_type m_nodata;
  typedef ChannelRange<channel_type> range_type;
public:
  NotNoDataFunctor( channel_type nodata ) : m_nodata(nodata) {}

  template <class Args> struct result {
    typedef channel_type type;
  };

  inline channel_type operator()( channel_type const& val ) const {
    return val != m_nodata ? range_type::max() : range_type::min();
  }
};

template <class ImageT, class NoDataT>
UnaryPerPixelView<ImageT,UnaryCompoundFunctor<NotNoDataFunctor<typename ImageT::pixel_type>, typename ImageT::pixel_type>  >
inline notnodata( ImageViewBase<ImageT> const& image, NoDataT nodata ) {
  typedef UnaryCompoundFunctor<NotNoDataFunctor<typename ImageT::pixel_type>, typename ImageT::pixel_type> func_type;
  func_type func( nodata );
  return UnaryPerPixelView<ImageT,func_type>( image.impl(), func );
}

// Cosine Transfer Function
template <class PixelT>
class CosineTransFunc: public vw::UnaryReturnSameType {
  typedef typename CompoundChannelType<PixelT>::type channel_type;
  typedef ChannelRange<channel_type> range_type;
public:
  CosineTransFunc() {}

  template <class ArgT>
  inline ArgT operator()( ArgT const& value ) const {
    return range_type::max()*((1.0-cos(float(value)/float(range_type::max())*M_PI))/2.0);;
    //return (1-cos(value*M_PI))/2;
  }
};

// Function to zip 2 images into content and alpha image
template <class Arg1T, class Arg2T>
struct CreateAlphaFunc: public ReturnFixedType<typename PixelWithAlpha<Arg1T>::type > {
  typedef typename PixelWithAlpha<Arg1T>::type result_type;
  inline result_type operator()(Arg1T const& arg1,
                                Arg2T const& arg2 ) const {
    result_type t(arg1);
    t.a() = arg2;
    return t;
  }
};

template <class Image1T, class Image2T>
inline BinaryPerPixelView<Image1T, Image2T, CreateAlphaFunc<typename Image1T::pixel_type, typename Image2T::pixel_type> >
create_alpha( ImageViewBase<Image1T> const& image1,
              ImageViewBase<Image2T> const& image2 ) {
  typedef CreateAlphaFunc<typename Image1T::pixel_type, typename Image2T::pixel_type> func_type;
  return BinaryPerPixelView<Image1T, Image2T, func_type >(image1.impl(), image2.impl(), func_type());
}

// Operation code for data that uses nodata
template <class PixelT>
void grassfire_nodata( std::string input, typename CompoundChannelType<PixelT>::type nodata,
                       std::string output ) {
  DiskImageView<PixelT> input_image(input);
  ImageView<int32> distance = grassfire(notnodata(input_image,nodata));
  int32 max = max_pixel_value( distance );
  vw_out() << "\tMax distance: " << max << "\n";
  typedef typename CompoundChannelType<PixelT>::type inter_type;
  typedef ChannelRange<typename CompoundChannelType<PixelT>::type> range_type;
  ImageViewRef<inter_type> norm_dist = pixel_cast<inter_type>(range_type::max()*pixel_cast<float>(distance)/float(max));
  ImageViewRef<typename PixelWithAlpha<PixelT>::type> result = create_alpha(input_image,per_pixel_filter(norm_dist,CosineTransFunc<inter_type>()));
  write_image(output, result);
}

// Same as above but modified for alpha input
template <class PixelT>
void grassfire_already_alpha( std::string input,
                              std::string output ) {
  DiskImageView<PixelT> input_image(input);
  ImageView<int32> distance = grassfire(apply_mask(invert_mask(alpha_to_mask(input_image)),1));
  int32 max = max_pixel_value(distance);
  vw_out() << "\tMax distance: " << max << "\n";
  typedef typename CompoundChannelType<PixelT>::type inter_type;
  typedef ChannelRange<typename CompoundChannelType<PixelT>::type> range_type;
  ImageViewRef<inter_type> norm_dist = pixel_cast<inter_type>(range_type::max()*pixel_cast<float>(distance)/float(max));
  ImageViewRef<PixelT> result = create_alpha(input_image,per_pixel_filter(norm_dist,CosineTransFunc<inter_type>()));
  write_image(output, result);
}

// Standard interface
int main( int argc, char *argv[] ) {
  std::vector<std::string> input_file_names;
  float nodata_value;

  po::options_description general_options("Options");
  general_options.add_options()
    ("nodata-value", po::value<float>(&nodata_value), "Value that is no data in input image. Not used if input has alpha.")
    ("help,h", "Display this help message");

  po::options_description hidden_options("");
  hidden_options.add_options()
    ("input-files", po::value<std::vector<std::string> >(&input_file_names));

  po::options_description options("Allowed Options");
  options.add(general_options).add(hidden_options);

  po::positional_options_description p;
  p.add("input-files", -1);

  po::variables_map vm;
  po::store( po::command_line_parser( argc, argv ).options(options).positional(p).run(), vm );
  po::notify( vm );

  std::ostringstream usage;
  usage << "Usage: " << argv[0] << "[options] <image-files> ... " << std::endl << std::endl;
  usage << general_options << std::endl;

  if( vm.count("help") ) {
    vw_out() << usage.str() << std::endl;
    return 1;
  } else if ( input_file_names.size() == 0 ) {
    vw_out() << "ERROR! Require an input file.\n";
    vw_out() << "\n" << usage.str() << "\n";
    return 1;
  }

  // iterate through input files
  for ( uint i = 0; i < input_file_names.size(); i++ ) {

    // First determine the format of the input
    DiskImageResource *first_resource = DiskImageResource::open(input_file_names[i]);
    ChannelTypeEnum channel_type = first_resource->channel_type();
    PixelFormatEnum pixel_format = first_resource->pixel_format();
    delete first_resource;

    vw_out() << "Loading: " << input_file_names[i] << "\n";
    size_t pt_idx = input_file_names[i].rfind(".");
    std::string output_file = input_file_names[i].substr(0,pt_idx)+"_grass"+input_file_names[i].substr(pt_idx,input_file_names[i].size()-pt_idx);

    switch (pixel_format) {
    case VW_PIXEL_GRAY:
      switch (channel_type) {
      case VW_CHANNEL_UINT8:
        grassfire_nodata<PixelGray<uint8> >(input_file_names[i],
                                            nodata_value,output_file); break;
      case VW_CHANNEL_INT16:
        grassfire_nodata<PixelGray<int16> >(input_file_names[i],
                                            nodata_value,output_file); break;
      case VW_CHANNEL_UINT16:
        grassfire_nodata<PixelGray<uint16> >(input_file_names[i],
                                             nodata_value,output_file); break;
      default:
        grassfire_nodata<PixelGray<float32> >(input_file_names[i],
                                              nodata_value,output_file); break;
      }
      break;
    case VW_PIXEL_GRAYA:
      switch (channel_type) {
      case VW_CHANNEL_UINT8:
        grassfire_already_alpha<PixelGrayA<uint8> >(input_file_names[i],
                                                    output_file); break;
      case VW_CHANNEL_INT16:
        grassfire_already_alpha<PixelGrayA<int16> >(input_file_names[i],
                                                    output_file); break;
      case VW_CHANNEL_UINT16:
        grassfire_already_alpha<PixelGrayA<uint16> >(input_file_names[i],
                                                     output_file); break;
      default:
        grassfire_already_alpha<PixelGrayA<float32> >(input_file_names[i],
                                                      output_file); break;
      }
      break;
    case VW_PIXEL_RGB:
      switch (channel_type) {
      case VW_CHANNEL_UINT8:
        grassfire_nodata<PixelRGB<uint8> >(input_file_names[i],
                                            nodata_value,output_file); break;
      case VW_CHANNEL_INT16:
        grassfire_nodata<PixelRGB<int16> >(input_file_names[i],
                                            nodata_value,output_file); break;
      case VW_CHANNEL_UINT16:
        grassfire_nodata<PixelRGB<uint16> >(input_file_names[i],
                                             nodata_value,output_file); break;
      default:
        grassfire_nodata<PixelRGB<float32> >(input_file_names[i],
                                              nodata_value,output_file); break;
      }
      break;
    default:
      switch (channel_type) {
      case VW_CHANNEL_UINT8:
        grassfire_already_alpha<PixelRGBA<uint8> >(input_file_names[i],
                                                   output_file); break;
      case VW_CHANNEL_INT16:
        grassfire_already_alpha<PixelRGBA<int16> >(input_file_names[i],
                                                   output_file); break;
      case VW_CHANNEL_UINT16:
        grassfire_already_alpha<PixelRGBA<uint16> >(input_file_names[i],
                                                    output_file); break;
      default:
        grassfire_already_alpha<PixelRGBA<float32> >(input_file_names[i],
                                                     output_file); break;
      }
      break;
    }
  } // end of looping through files
}
