// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// deplate.cc
///
/// Converts a plate file to GeoTIFF tiles on disk.

#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Math.h>
#include <vw/Plate/PlateView.h>
#include <vw/Plate/PlateCarreePlateManager.h>
#include <vw/Cartography/GeoReference.h>

using namespace vw;
using namespace vw::platefile;

#include <boost/filesystem/path.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;
namespace fs = boost::filesystem;

struct Options {
  Options() : west(0), east(0), north(0), south(0), tile_size(0),
              tile_size_deg(0), pds_dem_mode(false), pds_imagery_mode(false) {}

  // Input for project file
  std::string plate_file_name;

  // Settings
  int32 west, east, north, south, tile_size, tile_size_deg;
  double tile_ppd;
  bool pds_dem_mode, pds_imagery_mode;
  int level;

  // Output
  std::string output_datum, output_prefix;
};

// ConvertToPDSImagery
//  .. coverts an image into the range of 1-254 where 0 is invalid and 255 is reserved.
template <class PixelT>
class ConvertToPDSImagery : public ReturnFixedType<typename CompoundChannelCast<typename PixelWithoutAlpha<PixelT>::type,uint8>::type> {
  typedef typename CompoundChannelCast<typename PixelWithoutAlpha<PixelT>::type,uint8>::type return_type;
  typedef typename CompoundChannelCast<typename PixelWithoutAlpha<PixelT>::type,double>::type return_double_type;
  typedef typename CompoundChannelType<PixelT>::type input_channel_type;
public:
  return_type operator()( PixelT value ) const {
    if ( is_transparent(value) ) {
      return return_type();
    } else {
      return_double_type result
        = channel_cast<double>(non_alpha_channels(value))*253.0/ChannelRange<input_channel_type>::max();
      for ( uint8 i = 0; i < CompoundNumChannels<return_double_type>::value; i++ ) {
        if (result[i] < 0)
          result[i] = 0;
        if (result[i] > 253)
          result[i] = 253;
      }
      return channel_cast<uint8>(result)+1;
    }
  }
};

template <class ImageT>
UnaryPerPixelView<ImageT, ConvertToPDSImagery<typename ImageT::pixel_type> >
inline convert_to_pds_imagery( ImageViewBase<ImageT> const& image ) {
  typedef ConvertToPDSImagery<typename ImageT::pixel_type> func_type;
  return UnaryPerPixelView<ImageT, func_type>( image.impl(), func_type() );
}

template <class PixelT>
void do_tiles(boost::shared_ptr<PlateFile> platefile, Options& opt) {

  PlateCarreePlateManager<PixelT> pm(platefile);
  cartography::GeoReference output_georef;
  output_georef = pm.georeference(platefile->num_levels()-1);

  cartography::Datum datum;
  datum.set_well_known_datum( opt.output_datum );
  output_georef.set_datum( datum );

  PlateView<PixelT> plate_view(opt.plate_file_name);
  if ( opt.level != -1 )
    plate_view.set_level( opt.level );
  ImageViewRef<PixelT> plate_view_ref = plate_view;
  double scale_change = 1;
  if ( opt.tile_ppd > 0 ) {
    // Finding out our current PPD and attempting to match
    double curr_ppd = norm_2(output_georef.lonlat_to_pixel(Vector2(0,0))-
                             output_georef.lonlat_to_pixel(Vector2(1,0)));
    scale_change = opt.tile_ppd / curr_ppd;
    plate_view_ref = resample( plate_view, scale_change, scale_change,
                               ZeroEdgeExtension());
    Matrix3x3 scale;
    scale.set_identity();
    scale(0,0) /= scale_change;
    scale(1,1) /= scale_change;
    output_georef.set_transform( output_georef.transform()*scale );

    // Double check
    curr_ppd = norm_2(output_georef.lonlat_to_pixel(Vector2(0,0))-
                      output_georef.lonlat_to_pixel(Vector2(1,0)));
  }
  std::cout << "Converting " << opt.plate_file_name << " to " << opt.output_prefix << "\n";
  std::cout << output_georef << "\n";

  // Get the output georeference.
  vw::BBox2i output_bbox;
  output_bbox.grow(output_georef.lonlat_to_pixel(Vector2(opt.west, opt.north)));
  output_bbox.grow(output_georef.lonlat_to_pixel(Vector2(opt.east, opt.north)));
  output_bbox.grow(output_georef.lonlat_to_pixel(Vector2(opt.west, opt.south)));
  output_bbox.grow(output_georef.lonlat_to_pixel(Vector2(opt.east, opt.south)));
  std::cout << "\t--> Output bbox: " << output_bbox << "\n";

  if ( opt.tile_size_deg > 0 ) {
    if ( opt.tile_ppd > 0 ) {
      opt.tile_size =
        boost::numeric_cast<int32>(opt.tile_size_deg * opt.tile_ppd);
    } else {
      // User must have specified out to be sized in degrees
      opt.tile_size =
        boost::numeric_cast<int32>(norm_2(output_georef.lonlat_to_pixel(Vector2(0,0)) -
                                          output_georef.lonlat_to_pixel(Vector2(opt.tile_size_deg,0))));
    }
  }

  // Compute the bounding box for each tile.
  std::vector<BBox2i> crop_bboxes = image_blocks(crop(plate_view_ref, output_bbox),
                                                 opt.tile_size, opt.tile_size);

  for (unsigned i = 0; i < crop_bboxes.size(); ++i) {
    // The crop bboxes start at (0,0), and we want them to start at
    // the upper left corner of the output_bbox.
    crop_bboxes[i].min() += output_bbox.min();
    crop_bboxes[i].max() += output_bbox.min();

    { // Checking to see if this section is transparent
      std::list<TileHeader> theaders =
        plate_view.search_for_tiles( crop_bboxes[i] / scale_change );
      if (theaders.empty())
        continue;
    }

    cartography::GeoReference tile_georef = output_georef;
    Vector2 top_left_ll = output_georef.pixel_to_lonlat(crop_bboxes[i].min());
    Matrix3x3 T = tile_georef.transform();
    T(0,2) = top_left_ll(0);
    T(1,2) = top_left_ll(1);
    tile_georef.set_transform(T);

    std::cout << "\t--> Generating tile " << (i+1) << " / " << crop_bboxes.size()
              << " : " << crop_bboxes[i] << "\n"
              << "\t    with transform  "
              << tile_georef.transform() << "\n";

    std::ostringstream output_filename;
    output_filename << opt.output_prefix << "_"
                    << abs(int32(top_left_ll[0]));
    if ( top_left_ll[0] < 0 )
      output_filename << "W_";
    else
      output_filename << "E_";
    output_filename << abs(int32(top_left_ll[1]));
    if ( top_left_ll[1] >= 0 )
      output_filename << "N.tif";
    else
      output_filename << "S.tif";

    ImageView<PixelT> cropped_view = crop(plate_view_ref, crop_bboxes[i]);
    DiskImageResourceGDAL::Options gdal_options;
    gdal_options["COMPRESS"] = "LZW";

    if ( opt.pds_dem_mode ) {
      ImageViewRef<typename CompoundChannelCast<typename PixelWithoutAlpha<PixelT>::type ,int16>::type > dem_image = apply_mask(alpha_to_mask(channel_cast<int16>(cropped_view)),-32767);
      DiskImageResourceGDAL rsrc(output_filename.str(), dem_image.format(),
                                 Vector2i(256,256), gdal_options);
      rsrc.set_nodata_write( -32767 );
      write_georeference(rsrc, tile_georef);
      write_image(rsrc, dem_image,
                  TerminalProgressCallback( "plate.tools", "\t    Writing: "));
    } else if ( opt.pds_imagery_mode ) {
      ImageViewRef<typename CompoundChannelCast<typename PixelWithoutAlpha<PixelT>::type ,uint8>::type > dem_image = convert_to_pds_imagery( cropped_view );
      DiskImageResourceGDAL rsrc(output_filename.str(), dem_image.format(),
                                 Vector2i(256,256), gdal_options);
      rsrc.set_nodata_write( 0 );
      write_georeference(rsrc, tile_georef);
      write_image(rsrc, dem_image,
                  TerminalProgressCallback( "plate.tools", "\t    Writing: "));
    } else {
      DiskImageResourceGDAL rsrc(output_filename.str(), cropped_view.format(),
                                 Vector2i(256,256), gdal_options);
      write_georeference(rsrc, tile_georef);
      write_image(rsrc, cropped_view,
                  TerminalProgressCallback( "plate.tools", "\t    Writing: "));
    }

  }
}

void handle_arguments( int argc, char *argv[], Options& opt ) {
  po::options_description general_options("Chops platefile into georeferenced squares.\n\nGeneral Options");
  general_options.add_options()
    ("output-prefix,o", po::value(&opt.output_prefix), "Specify the base output directory")
    ("level,l", po::value(&opt.level)->default_value(-1), "Level inside the plate in which to process. -1 will error out at show the number of levels available.")
    ("west,w", po::value(&opt.west)->default_value(-180), "Specify west edge of the region to extract (deg).")
    ("east,e", po::value(&opt.east)->default_value(180), "Specify east edge of the region to extract (deg).")
    ("north,n", po::value(&opt.north)->default_value(90), "Specify north edge of the region to extract (deg).")
    ("south,s", po::value(&opt.south)->default_value(-90), "Specify south edge of the region to extract (deg).")
    ("tile-size-px,t", po::value(&opt.tile_size)->default_value(4096), "Specify the size of each output dem in pixels.")
    ("tile-size-deg,d", po::value(&opt.tile_size_deg), "Specify the size of each output dem in degrees.")
    ("tile-px-per-degree,p", po::value(&opt.tile_ppd), "Specify the output tiles' pixel per degrees.")
    ("output-datum", po::value(&opt.output_datum)->default_value("WGS84"), "Specify the output datum to use, [WGS84, WGS72, D_MOON, D_MARS]")
    ("export-pds-dem", "Export using int16 channel value with a -32767 nodata value")
    ("export-pds-imagery", "Export using uint8 channel value with a 0 nodata value")
    ("help,h", "Display this help message");

  po::options_description positional("");
  positional.add_options()
    ("plate-file", po::value(&opt.plate_file_name));

  po::positional_options_description positional_desc;
  positional_desc.add("plate-file", 1);

  po::options_description all_options;
  all_options.add(general_options).add(positional);

  po::variables_map vm;
  try {
    po::store( po::command_line_parser( argc, argv ).options(all_options).positional(positional_desc).run(), vm );
    po::notify( vm );
  } catch (po::error &e) {
    vw_throw( ArgumentErr() << "Error parsing input:\n\t"
              << e.what() << general_options );
  }

  std::ostringstream usage;
  usage << "Usage: " << argv[0] << " <platefile-url> <optional project settings> \n";

  if ( vm.count("help") )
    vw_throw( ArgumentErr() << usage.str() << general_options );
  if ( opt.plate_file_name.empty() )
    vw_throw( ArgumentErr() << "Missing input platefile!\n"
              << usage.str() << general_options );

  opt.pds_dem_mode = vm.count("export-pds-dem");
  opt.pds_imagery_mode = vm.count("export-pds-imagery");

  if( opt.output_prefix == "" ) {
    opt.output_prefix = fs::path(opt.plate_file_name).stem();
    size_t indx = opt.output_prefix.rfind("/");
    if ( indx != std::string::npos )
      opt.output_prefix = opt.output_prefix.substr(indx+1);
  }
}

int main( int argc, char *argv[] ) {

  Options opt;
  try {
    handle_arguments( argc, argv, opt );

    // XXX: should make plate_file_name a Url
    boost::shared_ptr<PlateFile> platefile(new PlateFile(Url(opt.plate_file_name)));

    std::cout << "Opened " << opt.plate_file_name << ".     Depth: "
              << platefile->num_levels() << " levels.\n";

    PixelFormatEnum pixel_format = platefile->pixel_format();
    ChannelTypeEnum channel_type = platefile->channel_type();

    switch(pixel_format) {
    case VW_PIXEL_GRAYA:
      switch(channel_type) {
      case VW_CHANNEL_UINT8:
        do_tiles<PixelGrayA<uint8> >(platefile, opt);
        break;
      case VW_CHANNEL_INT16:
        do_tiles<PixelGrayA<int16> >(platefile, opt);
        break;
      case VW_CHANNEL_FLOAT32:
        do_tiles<PixelGrayA<float32> >(platefile, opt);
        break;
      default:
        vw_throw(ArgumentErr() << "Platefile contains a channel type not supported by image2plate.\n");
      }
      break;
    case VW_PIXEL_RGB:
    case VW_PIXEL_RGBA:
    default:
      switch(channel_type) {
      case VW_CHANNEL_UINT8:
        do_tiles<PixelRGBA<uint8> >(platefile, opt);
        break;
      default:
        vw_throw(ArgumentErr() << "Platefile contains a channel type not supported by image2plate.\n");
      }
      break;
    }

  } catch ( const ArgumentErr& e ) {
    vw_out() << e.what() << std::endl;
    return 1;
  } catch ( const Exception& e ) {
    std::cerr << "VW Error: " << e.what() << std::endl;
    return 1;
  } catch ( const std::bad_alloc& e ) {
    std::cerr << "Error: Ran out of Memory!" << std::endl;
    return 1;
  } catch ( const std::exception& e ) {
    std::cerr << "Error: " << e.what() <<  std::endl;
    return 1;
  }
}
