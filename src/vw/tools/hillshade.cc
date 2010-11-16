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

#ifdef NDEBUG
#undef NDEBUG
#endif

#include <cstdlib>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <boost/filesystem/path.hpp>
namespace fs = boost::filesystem;

#include <vw/Math/EulerAngles.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/Algorithms.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Image/Filter.h>
#include <vw/Image/PixelMask.h>
#include <vw/Image/MaskViews.h>
#include <vw/Image/PerPixelAccessorViews.h>
#include <vw/FileIO/DiskImageView.h>
#include <vw/FileIO/DiskImageResourceGDAL.h>
#include <vw/Cartography/GeoReference.h>
#include <vw/tools/Common.h>

using namespace vw;


// Global Variables
std::string input_file_name, output_file_name = "";
double azimuth, elevation, scale;
double nodata_value;
double blur_sigma;

// Allows FileIO to correctly read/write these pixel types
namespace vw {
  template<> struct PixelFormatID<Vector3>   { static const PixelFormatEnum value = VW_PIXEL_XYZ; };
}

// ---------------------------------------------------------------------------------------------------

//  compute_normals()
//
// Compute a vector normal to the surface of a DEM for each given
// pixel.  The normal is computed by forming a plane with three points
// in the vicinity of the requested pixel, and then finding the vector
// normal to that plane.  The user must specify the scale in the [u,v]
// directions so that the direction of the vector in physical space
// can be properly ascertained.  This is often contained in the (0,0)
// and (1,1) entry of the georeference transform.
class ComputeNormalsFunc : public ReturnFixedType<PixelMask<Vector3> >
{
  double m_u_scale, m_v_scale;

public:
  ComputeNormalsFunc(double u_scale, double v_scale) :
    m_u_scale(u_scale), m_v_scale(v_scale) {}

  BBox2i work_area() const { return BBox2i(Vector2i(0, 0), Vector2i(1, 1)); }

  template <class PixelAccessorT>
  PixelMask<Vector3> operator() (PixelAccessorT const& accessor_loc) const {
    PixelAccessorT acc = accessor_loc;

    // Pick out the three altitude values.
    if (is_transparent(*acc))
      return PixelMask<Vector3>();
    double alt1 = *acc;

    acc.advance(1,0);
    if (is_transparent(*acc))
      return PixelMask<Vector3>();
    double alt2 = *acc;

    acc.advance(-1,1);
    if (is_transparent(*acc))
      return PixelMask<Vector3>();
    double alt3 = *acc;

    // Form two orthogonal vectors in the plane containing the three
    // altitude points
    Vector3 n1(m_u_scale, 0, alt2-alt1);
    Vector3 n2(0, m_v_scale, alt3-alt1);

    // Return the vector normal to the local plane.
    return normalize(cross_prod(n1,n2));
  }
};

template <class ViewT>
UnaryPerPixelAccessorView<EdgeExtensionView<ViewT,ConstantEdgeExtension>, ComputeNormalsFunc> compute_normals(ImageViewBase<ViewT> const& image,
                                                                                                              double u_scale, double v_scale) {
  return UnaryPerPixelAccessorView<EdgeExtensionView<ViewT,ConstantEdgeExtension>, ComputeNormalsFunc>(edge_extend(image.impl(), ConstantEdgeExtension()),
                                                                                                       ComputeNormalsFunc (u_scale, v_scale));
}

class DotProdFunc : public ReturnFixedType<PixelMask<PixelGray<double> > > {
  Vector3 m_vec;
public:
  DotProdFunc(Vector3 const& vec) : m_vec(vec) {}
  PixelMask<PixelGray<double> > operator() (PixelMask<Vector3> const& pix) const {
    if (is_transparent(pix))
      return PixelMask<PixelGray<double> >();
    else {
//       std::cout << "Vec1 : " << pix.child() << "   " << norm_2(pix.child()) << "\n";
//       std::cout << "Vec2 : " << m_vec << "   " << norm_2(m_vec) << "\n";
//       std::cout << "OVerall: " << dot_prod(pix.child(),m_vec)/(norm_2(pix.child()) * norm_2(m_vec)) << "\n\n";
      return dot_prod(pix.child(),m_vec)/(norm_2(pix.child()) * norm_2(m_vec));
    }
  }
};

template <class ViewT>
UnaryPerPixelView<ViewT, DotProdFunc> dot_prod(ImageViewBase<ViewT> const& view, Vector3 const& vec) {
  return UnaryPerPixelView<ViewT, DotProdFunc>(view.impl(), DotProdFunc(vec));
}

// ---------------------------------------------------------------------------------------------------

template <class PixelT>
void do_hillshade(po::variables_map const& vm) {

  cartography::GeoReference georef;
  cartography::read_georeference(georef, input_file_name);

  // Select the pixel scale.
  double u_scale, v_scale;
  if (scale == 0) {
    if (georef.is_projected()) {
      u_scale = georef.transform()(0,0);
      v_scale = georef.transform()(1,1);
    } else {
      double meters_per_degree = 2*M_PI*georef.datum().semi_major_axis()/360.0;
      u_scale = georef.transform()(0,0) * meters_per_degree;
      v_scale = georef.transform()(1,1) * meters_per_degree;
    }
  } else {
    u_scale = scale;
    v_scale = -scale;
  }
  // For debugging:
  //  std::cout << "\t--> Scale: " << u_scale << "  " << v_scale << "\n";

  // Set the direction of the light source.
  Vector3 light_0(1,0,0);
  Vector3 light = vw::math::euler_to_rotation_matrix(elevation*M_PI/180, azimuth*M_PI/180, 0, "yzx") * light_0;

  // Compute the surface normals
  std::cout << "Loading: " << input_file_name << ".\n";
  DiskImageView<PixelT> disk_dem_file(input_file_name);
  ImageViewRef<PixelGray<double> > input_image = channel_cast<double>(disk_dem_file);

  ImageViewRef<PixelMask<PixelGray<double> > > dem;
  SrcImageResource *disk_dem_rsrc = DiskImageResource::open(input_file_name);
  if (vm.count("nodata-value")) {
    std::cout << "\t--> Masking pixel value: " << nodata_value << ".\n";
    dem = create_mask(input_image, nodata_value);
  } else if ( disk_dem_rsrc->has_nodata_read() ) {
    nodata_value = disk_dem_rsrc->nodata_read();
    std::cout << "\t--> Extracted nodata value from file: " << nodata_value << ".\n";
    dem = create_mask(input_image, nodata_value);
  } else {
    dem = pixel_cast<PixelMask<PixelGray<double> > >(input_image);
  }
  delete disk_dem_rsrc;

  if (vm.count("blur")) {
    std::cout << "\t--> Blurring pixel with gaussian kernal.  Sigma = " << blur_sigma << "\n";
    dem = gaussian_filter(dem, blur_sigma);
  }

  // The final result is the dot product of the light source with the normals
  ImageViewRef<PixelMask<PixelGray<uint8> > > shaded_image = channel_cast_rescale<uint8>(clamp(dot_prod(compute_normals(dem, u_scale, v_scale), light)));

  // Save the result
  std::cout << "Writing shaded relief image: " << output_file_name << "\n";

  DiskImageResourceGDAL rsrc(output_file_name, shaded_image.format());
  rsrc.set_block_write_size(Vector2i(1024,1024));
  write_georeference(rsrc, georef);
  write_image(rsrc, shaded_image,
              TerminalProgressCallback( "tools.hillshade", "Writing:"));
}

int main( int argc, char *argv[] ) {

  po::options_description desc("Description: Outputs image of a DEM lighted as specified\n\nUsage: hillshade [options] <input file> \n\nOptions");
  desc.add_options()
    ("help,h", "Display this help message")
    ("input-file", po::value<std::string>(&input_file_name), "Explicitly specify the input file")
    ("output-file,o", po::value<std::string>(&output_file_name), "Specify the output file")
    ("azimuth,a", po::value<double>(&azimuth)->default_value(0), "Sets the direction tha the light source is coming from (in degrees).  Zero degrees is to the right, with positive degree counter-clockwise.")
    ("elevation,e", po::value<double>(&elevation)->default_value(45), "Set the elevation of the light source (in degrees).")
    ("scale,s", po::value<double>(&scale)->default_value(0), "Set the scale of a pixel (in the same units as the DTM height values.")
    ("nodata-value", po::value<double>(&nodata_value), "Remap the DEM default value to the min altitude value.")
    ("blur", po::value<double>(&blur_sigma), "Pre-blur the DEM with the specified sigma.");
  po::positional_options_description p;
  p.add("input-file", 1);

  po::variables_map vm;
  try {
    po::store( po::command_line_parser( argc, argv ).options(desc).positional(p).run(), vm );
    po::notify( vm );
  } catch (po::error &e) {
    std::cout << "An error occured while parsing command line arguments.\n";
    std::cout << "\t" << e.what() << "\n\n";
    std::cout << desc << std::endl;
    return 1;
  }

  if( vm.count("help") ) {
    std::cout << desc << std::endl;
    return 1;
  }

  if( vm.count("input-file") != 1 ) {
    std::cout << "Error: Must specify exactly one input file!\n" << std::endl;
    std::cout << desc << std::endl;
    return 1;
  }

  if( output_file_name == "" ) {
    output_file_name = fs::path(input_file_name).replace_extension().string() + "_HILLSHADE.tif";
  }

  try {
    // Get the right pixel/channel type.
    ImageFormat fmt = tools::taste_image(input_file_name);

    switch(fmt.pixel_format) {
    case VW_PIXEL_GRAY:
    case VW_PIXEL_GRAYA:
    case VW_PIXEL_RGB:
    case VW_PIXEL_RGBA:
      switch(fmt.channel_type) {
      case VW_CHANNEL_UINT8:  do_hillshade<PixelGray<uint8>   >(vm); break;
      case VW_CHANNEL_INT16:  do_hillshade<PixelGray<int16>   >(vm); break;
      case VW_CHANNEL_UINT16: do_hillshade<PixelGray<uint16>  >(vm); break;
      case VW_CHANNEL_FLOAT64:do_hillshade<PixelGray<float64> >(vm); break;
      default:                do_hillshade<PixelGray<float32> >(vm); break;
      }
      break;
    default:
      std::cout << "Error: Unsupported pixel format.\n";
      exit(0);
    }
  } catch( Exception& e ) {
    std::cout << "Error: " << e.what() << std::endl;
  }

  return 0;
}
