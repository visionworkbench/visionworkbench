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

#include <vw/tools/Common.h>
#include <vw/Mosaic/ImageComposite.h>
#include <vw/Cartography/GeoReference.h>
#include <vw/Cartography/GeoTransform.h>
#include <vw/FileIO/DiskImageView.h>
#include <vw/Image/PerPixelViews.h>
#include <vw/Image/Filter.h>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

using namespace vw;
using namespace vw::math;
using namespace vw::cartography;

// Variables from the command line.
std::vector<std::string> image_files;
std::string mosaic_name;
std::string output_file_type;
std::string channel_type_str;
bool draft;
unsigned int tilesize;
bool tile_output = false;
unsigned int patch_size, patch_overlap;
float nodata_value;
bool has_nodata_value = false;

// Creates an alpha channel based on pixels with a value of zero.  The
// first version covers scalar types.  The remaining versions cover
// compound pixel types.
template <class PixelT> struct AlphaTypeFromPixelType { typedef PixelGrayA<PixelT> type; };
template<class ChannelT> struct AlphaTypeFromPixelType<PixelGray<ChannelT> > { typedef PixelGrayA<ChannelT> type; };
template<class ChannelT> struct AlphaTypeFromPixelType<PixelGrayA<ChannelT> > { typedef PixelGrayA<ChannelT> type; };
template<class ChannelT> struct AlphaTypeFromPixelType<PixelRGB<ChannelT> > { typedef PixelRGBA<ChannelT> type; };
template<class ChannelT> struct AlphaTypeFromPixelType<PixelRGBA<ChannelT> > { typedef PixelRGBA<ChannelT> type; };

template <class PixelT> struct NonAlphaTypeFromPixelType { typedef PixelGrayA<PixelT> type; };
template<class ChannelT> struct NonAlphaTypeFromPixelType<PixelGray<ChannelT> > { typedef PixelGray<ChannelT> type; };
template<class ChannelT> struct NonAlphaTypeFromPixelType<PixelGrayA<ChannelT> > { typedef PixelGray<ChannelT> type; };
template<class ChannelT> struct NonAlphaTypeFromPixelType<PixelRGB<ChannelT> > { typedef PixelRGB<ChannelT> type; };
template<class ChannelT> struct NonAlphaTypeFromPixelType<PixelRGBA<ChannelT> > { typedef PixelRGB<ChannelT> type; };

template <class PixelT>
class NodataToMaskFunctor: public vw::UnaryReturnTemplateType<AlphaTypeFromPixelType> {
  PixelT m_nodata_value;

public:
  NodataToMaskFunctor(typename PixelChannelType<PixelT>::type nodata_value = 0) : m_nodata_value(nodata_value) {}

  typename AlphaTypeFromPixelType<PixelT>::type operator() (PixelT const& pix) const {
    typedef typename AlphaTypeFromPixelType<PixelT>::type result_type;
    if (pix == m_nodata_value)
      return result_type();  // Mask pixel
    else
      return result_type(pix);
  }
};

template <class ViewT>
vw::UnaryPerPixelView<ViewT, NodataToMaskFunctor<typename ViewT::pixel_type> >
nodata_to_mask(vw::ImageViewBase<ViewT> const& view,
               typename PixelChannelType<typename ViewT::pixel_type>::type const& nodata_value = 0 ) {
  return vw::per_pixel_filter(view.impl(), NodataToMaskFunctor<typename ViewT::pixel_type>(nodata_value));
}

template <class PixelT>
class MaskToNodataFunctor: public vw::UnaryReturnTemplateType<NonAlphaTypeFromPixelType> {
  PixelT m_nodata_value;
  typedef typename PixelChannelType<PixelT>::type channel_type;

public:
  MaskToNodataFunctor(float nodata_value = 0) : m_nodata_value((channel_type)nodata_value) {}

  typename NonAlphaTypeFromPixelType<PixelT>::type operator() (PixelT const& pix) const {
    typedef typename NonAlphaTypeFromPixelType<PixelT>::type result_type;
    if (is_transparent(pix))
      return result_type(m_nodata_value);
    else
      return result_type(pix);
  }
};

template <class ViewT>
vw::UnaryPerPixelView<ViewT, MaskToNodataFunctor<typename ViewT::pixel_type> >
mask_to_nodata(vw::ImageViewBase<ViewT> const& view, float nodata_value = 0 ) {
  return vw::per_pixel_filter(view.impl(), MaskToNodataFunctor<typename ViewT::pixel_type>(nodata_value));
}

// do_blend()
//
template <class PixelT>
void do_blend() {

  typedef typename AlphaTypeFromPixelType<PixelT>::type alpha_pixel_type;
  typedef typename PixelChannelCast<alpha_pixel_type,float32>::type float_pixel_type;

  TerminalProgressCallback tpc( "tools.geoblend", "");

  vw::mosaic::ImageComposite<float_pixel_type> composite;
  if( draft ) composite.set_draft_mode( true );

  double smallest_x_scale = vw::ScalarTypeLimits<float>::highest();
  double smallest_y_scale = vw::ScalarTypeLimits<float>::highest();
  double smallest_x_val = vw::ScalarTypeLimits<float>::highest();
  double largest_y_val = vw::ScalarTypeLimits<float>::lowest();

  // First pass, read georeferencing information and build an output
  // georef.
  tpc.set_progress_text( "Status (scanning):   " );
  SubProgressCallback scanning_pc( tpc, 0, 0.05 );
  for(unsigned i = 0; i < image_files.size(); ++i) {
    vw_out(vw::VerboseDebugMessage) << "Adding file " << image_files[i] << std::endl;
    scanning_pc.report_fractional_progress(i,image_files.size());

    GeoReference input_georef;
    read_georeference( input_georef, image_files[i] );
    DiskImageView<PixelT> source_disk_image( image_files[i] );
    vw_out(vw::VerboseDebugMessage) << "\tTransform: " << input_georef.transform()
                                    << "\t\tBBox: " << input_georef.bounding_box(source_disk_image) << std::endl;

    // Check to make sure the image has valid georeferencing
    // information.
    if( input_georef.transform() == identity_matrix<3>() ) {
      vw_out(InfoMessage) << "No georeferencing info found for image: \"" << image_files[i] << "\".  Aborting." << std::endl;
      exit(0);
    }

    Matrix3x3 affine = input_georef.transform();
    if (fabs(affine(0,0)) < smallest_x_scale || fabs(affine(1,1)) < smallest_y_scale) {
      smallest_x_scale = affine(0,0);
      smallest_y_scale = affine(1,1);
    }

    if (affine(0,2) < smallest_x_val)
      smallest_x_val = affine(0,2);

    // Note: since the y coordinates are typically flipped in a DEM,
    // we look here for the _largest_ y value, since it corresponds to
    // the upper left hand pixel.
    if (affine(1,2) > largest_y_val)
      largest_y_val = affine(1,2);

  }
  scanning_pc.report_finished();

  // Convert all of the images so they share the same scale factor.
  // Adopt the scale of the highest resolution being composited.
  Matrix3x3 output_affine = identity_matrix<3>();
  output_affine(0,0) = smallest_x_scale;
  output_affine(1,1) = smallest_y_scale;
  output_affine(0,2) = smallest_x_val;
  output_affine(1,2) = largest_y_val;

  vw_out(VerboseDebugMessage) << "Output affine transform: " << output_affine << std::endl
                              << int(smallest_x_val) << "  " << int(largest_y_val) << std::endl;

  // Take the georef from the first file (this ensures that the
  // projection and datum information is preserved...), but update the
  // affine transform.
  GeoReference output_georef;
  read_georeference( output_georef, image_files[0] );
  output_georef.set_transform(output_affine);

  tpc.set_progress_text( "Status (assembling): " );
  SubProgressCallback assembling_pc( tpc, 0.05, 0.1 );
  // Second pass: add files to the image composite.
  for(unsigned i = 0; i < image_files.size(); ++i) {
    assembling_pc.report_fractional_progress(i, image_files.size() );
    GeoReference input_georef;
    read_georeference(input_georef, image_files[i]);
    DiskImageView<PixelT> source_disk_image( image_files[i] );

    GeoTransform trans(input_georef, output_georef);
    BBox2 output_bbox = trans.forward_bbox( BBox2(0,0,source_disk_image.cols(),source_disk_image.rows()) );
    vw_out(vw::VerboseDebugMessage) << "output_bbox = " << output_bbox << std::endl;

    // I've hardwired this to use nearest pixel interpolation for now
    // until we have a chance to sit down and develop a better
    // strategy for intepolating and filtering in the presence of
    // missing pixels in DEMs. -mbroxton
    if (has_nodata_value) {
      ImageViewRef<alpha_pixel_type> masked_source = crop( transform( nodata_to_mask(source_disk_image, (typename PixelChannelType<PixelT>::type)(nodata_value) ), trans, ZeroEdgeExtension(), NearestPixelInterpolation() ), output_bbox );
      composite.insert( channel_cast_rescale<float32>(masked_source), (int)output_bbox.min().x(), (int)output_bbox.min().y() );
    } else {
     ImageViewRef<alpha_pixel_type> masked_source = crop( transform( pixel_cast<alpha_pixel_type>(source_disk_image), trans, ZeroEdgeExtension(), NearestPixelInterpolation() ), output_bbox );
     composite.insert( channel_cast_rescale<float32>(masked_source), (int)output_bbox.min().x(), (int)output_bbox.min().y() );
    }

  }
  assembling_pc.report_finished();

  tpc.set_progress_text( "Status (preparing):  " );
  vw_out(vw::VerboseDebugMessage) << std::endl;
  composite.prepare( SubProgressCallback( tpc, 0.1, 0.5 ) );
  vw_out(vw::VerboseDebugMessage) << "Composite dimensions: " << composite.cols() << "  " << composite.rows() << std::endl;

  tpc.set_progress_text( "Status (blending):   " );
  SubProgressCallback blending_pc( tpc, 0.5, 1.0 );
  // Output the image in tiles, or one large image.
  if(tile_output) {
    const int dim = patch_size - patch_overlap;
    const int tile_width = composite.cols() / dim;
    const int tile_height = composite.rows() / dim;
    vw_out(vw::VerboseDebugMessage) << "Outputting composite in " << tile_width * tile_height << " tiles." << std::endl;

    for(int i=0; i < composite.rows(); i += dim) {
      for(int j=0; j < composite.cols(); j += dim) {
        BBox2i tile_bbox(j, i, dim, dim);
        if(tile_bbox.max().x() >= composite.cols()) tile_bbox.max().x() = composite.cols();
        if(tile_bbox.max().y() >= composite.rows()) tile_bbox.max().y() = composite.cols();
        ImageView<PixelT> tile_view = crop(channel_cast_rescale<typename PixelChannelType<PixelT>::type>(composite), tile_bbox);
        GeoReference tile_georef = output_georef;

        // Adjust the affine transformation's offset to point to the upper
        // left of this tile.
        Vector2 upper_left = output_georef.pixel_to_point( Vector2(j, i) );
        Matrix3x3 tile_transform = tile_georef.transform();
        tile_transform(0,2) += upper_left[0];
        tile_transform(1,2) += upper_left[1];
        tile_georef.set_transform(tile_transform);

        // Filename for this tile.
        std::stringstream tile_filename;
        tile_filename << mosaic_name;
        tile_filename << '.' << j << '.' << i << '.';
        tile_filename << output_file_type;

        // Finally, write.
        write_georeferenced_image( tile_filename.str(), tile_view, tile_georef, blending_pc );
      }
    }
  } else {
    vw_out(vw::VerboseDebugMessage) << "Output image:" << std::endl
                                    << "\tTransform: " << output_affine << std::endl
                                    << "\t\tBBox: " << output_georef.bounding_box(composite) << " [ W: " << output_georef.bounding_box(composite).width() << " H: " << output_georef.bounding_box(composite).height() << " ]" << std::endl << std::endl;

    std::string mosaic_filename = mosaic_name+".blend."+output_file_type;
    DiskImageResourceGDAL *out_resource;
    ImageViewRef<PixelT> out_image;

    // Specify the output image resource.
    if (has_nodata_value) {
      out_image = pixel_cast<PixelT>( mask_to_nodata( channel_cast_rescale<typename PixelChannelType<PixelT>::type>(composite), nodata_value ) );
    } else {
      out_image = pixel_cast<PixelT>( channel_cast_rescale<typename PixelChannelType<PixelT>::type>(composite) );
    }

    // Set up tiled TIFF output, if it's specified.
    if(tilesize > 0)
      out_resource = new DiskImageResourceGDAL( mosaic_filename, out_image.format(), Vector2i(tilesize, tilesize) );
    else
      out_resource = new DiskImageResourceGDAL( mosaic_filename, out_image.format() );

    // Finally, write.
    write_georeference(*out_resource, output_georef);
    write_image(*out_resource, out_image, blending_pc);
    delete out_resource;
  }
  blending_pc.report_finished();

  tpc.report_finished();
}

int main( int argc, char *argv[] ) {
  try {

    po::options_description general_options("Options");
    general_options.add_options()
      ("mosaic-name,o", po::value<std::string>(&mosaic_name)->default_value("mosaic"), "Specify base output directory")
      ("output-file-type,t", po::value<std::string>(&output_file_type)->default_value("tif"), "Output file type")
      ("tile-output", "Output the leaf tiles of a quadtree, instead of a single blended image.")
      ("tiled-tiff", po::value<unsigned int>(&tilesize)->default_value(0), "Output a tiled TIFF image, with given tile size (0 disables, TIFF only)")
      ("patch-size", po::value<unsigned int>(&patch_size)->default_value(256), "Patch size for tiled output, in pixels")
      ("patch-overlap", po::value<unsigned int>(&patch_overlap)->default_value(0), "Patch overlap for tiled output, in pixels")
      ("draft", "Draft mode (no blending)")
      ("ignore-alpha", "Ignore the alpha channel of the input images, and don't write an alpha channel in output.")
      ("nodata-value", po::value<float>(&nodata_value), "Pixel value to use for nodata in input and output (when there's no alpha channel)")
      ("channel-type", po::value<std::string>(&channel_type_str), "Images' channel type. One of [uint8, uint16, int16, float].")
      ("help,h", "Display this help message");

    po::options_description hidden_options("");
    hidden_options.add_options()
      ("input-files", po::value<std::vector<std::string> >(&image_files));

    po::options_description options("Allowed Options");
    options.add(general_options).add(hidden_options);

    po::positional_options_description p;
    p.add("input-files", -1);

    std::ostringstream usage;
    usage << "Description: merges several DEMs" << std::endl << std::endl;
    usage << "Usage: geoblend [options] <filename1> <filename2> ..." << std::endl << std::endl;
    usage << general_options << std::endl;

    po::variables_map vm;
    try {
      po::store( po::command_line_parser( argc, argv ).options(options).positional(p).run(), vm );
      po::notify( vm );
    } catch(po::error &e) {
      std::cout << "An error occured while parsing command line arguments.\n";
      std::cout << "\t" << e.what() << "\n\n";
      std::cout << usage.str();
      return 1;
    }

    if( vm.count("help") ) {
      std::cerr << usage.str() << std::endl;
      return 1;
    }

    draft = false;
    if( vm.count("draft") ) {
      draft = true;
    }

    if( vm.count("input-files") < 1 ) {
      std::cerr << "Error: Must specify at least one input file!" << std::endl << std::endl;
      std::cerr << usage.str();
      return 1;
    }

    if(vm.count("tile-output")) tile_output = true;

    if( patch_size <= 0 ) {
      std::cerr << "Error: The patch size must be a positive number!  (You specified " << patch_size << ".)" << std::endl;
      std::cerr << usage << std::endl;
      return 1;
    }

    if( patch_overlap>=patch_size || patch_overlap%2==1 ) {
      std::cerr << "Error: The patch overlap must be an even number nonnegative number" << std::endl;
      std::cerr << "smaller than the patch size!  (You specified " << patch_overlap << ".)" << std::endl;
      std::cerr << usage << std::endl;
      return 1;
    }

    if(tilesize > 0 && vm.count("tile-output")) {
        std::cerr << "Error: Cannot output both a unified tiled TIFF (single image output) and individual tiles (multiple image output)." << std::endl;
        std::cerr << "\tThe two options --tiled-tiff (> 0) and --tile-output conflict." << std::endl;
        return 1;
    }

    if(output_file_type != "tif") {
        std::cerr << "Error: Cannot output a unified tiled TIFF if output file type is not set to tif." << std::endl;
        std::cerr << usage << std::endl;
        return 1;
    }

    if(vm.count("nodata-value")) has_nodata_value = true;

    ImageFormat fmt = tools::taste_image(image_files[0]);

    if (vm.count("channel-type")) {
      fmt.channel_type = channel_name_to_enum(channel_type_str);
      switch (fmt.channel_type) {
        case VW_CHANNEL_UINT8:  case VW_CHANNEL_INT16:
        case VW_CHANNEL_UINT16: case VW_CHANNEL_FLOAT32:
          break;
        default:
          std::cerr << "Error: Bad channel type specified." << std::endl;
          std::cerr << usage.str() << std::endl;
          exit(1);
      }
    }

    if (vm.count("ignore-alpha")) {
      if (fmt.pixel_format == VW_PIXEL_RGBA)  fmt.pixel_format = VW_PIXEL_RGB;
      if (fmt.pixel_format == VW_PIXEL_GRAYA) fmt.pixel_format = VW_PIXEL_GRAY;
    }

    vw_out(vw::VerboseDebugMessage) << "Using pixel type " << pixel_format_name(fmt.pixel_format) << ":" << channel_type_name(fmt.channel_type) << std::endl;

    switch (fmt.pixel_format) {
    case VW_PIXEL_GRAY:
      switch (fmt.channel_type) {
      case VW_CHANNEL_UINT8:   do_blend<vw::PixelGray<vw::uint8> >(); break;
      case VW_CHANNEL_INT16:   do_blend<vw::PixelGray<vw::int16> >(); break;
      case VW_CHANNEL_UINT16:  do_blend<vw::PixelGray<vw::uint16> >(); break;
      default:                 do_blend<vw::PixelGray<vw::float32> >(); break;
      }
      break;
    case VW_PIXEL_GRAYA:
      switch (fmt.channel_type) {
      case VW_CHANNEL_UINT8:   do_blend<vw::PixelGrayA<vw::uint8> >(); break;
      case VW_CHANNEL_INT16:   do_blend<vw::PixelGrayA<vw::int16> >(); break;
      case VW_CHANNEL_UINT16:  do_blend<vw::PixelGrayA<vw::uint16> >(); break;
      default:                 do_blend<vw::PixelGrayA<vw::float32> >(); break;
      }
      break;
    case VW_PIXEL_RGB:
      switch (fmt.channel_type) {
      case VW_CHANNEL_UINT8:   do_blend<vw::PixelRGB<vw::uint8> >(); break;
      case VW_CHANNEL_INT16:   do_blend<vw::PixelRGB<vw::int16> >(); break;
      case VW_CHANNEL_UINT16:  do_blend<vw::PixelRGB<vw::uint16> >(); break;
      default:                 do_blend<vw::PixelRGB<vw::float32> >(); break;
      }
      break;
    default:
      switch (fmt.channel_type) {
      case VW_CHANNEL_UINT8:   do_blend<vw::PixelRGBA<vw::uint8> >(); break;
      case VW_CHANNEL_INT16:   do_blend<vw::PixelRGBA<vw::int16> >(); break;
      case VW_CHANNEL_UINT16:  do_blend<vw::PixelRGBA<vw::uint16> >(); break;
      default:                 do_blend<vw::PixelRGBA<vw::float32> >(); break;
      }
      break;
    }

  }
  catch( std::exception &err ) {
    vw::vw_out(vw::ErrorMessage) << "Error: " << err.what() << std::endl << "Aborting!" << std::endl;
    return 1;
  }
  return 0;
}
