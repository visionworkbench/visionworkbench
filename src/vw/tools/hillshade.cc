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

// ----------------------------------------------------------------------------

//  compute_normals()
//
// Compute a vector normal to the surface of a DEM for each given
// pixel.  The normal is computed by forming a plane with three points
// in the vicinity of the requested pixel, and then finding the vector
// normal to that plane.  The user must specify the scale in the [u,v]
// directions so that the direction of the vector in physical space
// can be properly ascertained.  This is often contained in the (0,0)
// and (1,1) entry of the georeference transform.
class ComputeNormalsFunc : public ReturnFixedType<PixelMask<Vector3f> >
{
  float m_u_scale, m_v_scale;

public:
  ComputeNormalsFunc(float u_scale, float v_scale) :
    m_u_scale(u_scale), m_v_scale(v_scale) {}

  BBox2i work_area() const { return BBox2i(Vector2i(0, 0), Vector2i(1, 1)); }

  template <class PixelAccessorT>
  PixelMask<Vector3f> operator() (PixelAccessorT const& accessor_loc) const {
    PixelAccessorT acc = accessor_loc;

    // Pick out the three altitude values.
    if (is_transparent(*acc))
      return PixelMask<Vector3f>();
    float alt1 = *acc;

    acc.advance(1,0);
    if (is_transparent(*acc))
      return PixelMask<Vector3f>();
    float alt2 = *acc;

    acc.advance(-1,1);
    if (is_transparent(*acc))
      return PixelMask<Vector3f>();
    float alt3 = *acc;

    // Form two orthogonal vectors in the plane containing the three
    // altitude points
    Vector3f n1(m_u_scale, 0, alt2-alt1);
    Vector3f n2(0, m_v_scale, alt3-alt1);

    // Return the vector normal to the local plane.
    return normalize(cross_prod(n1,n2));
  }
};

template <class ViewT>
UnaryPerPixelAccessorView<EdgeExtensionView<ViewT,ConstantEdgeExtension>, ComputeNormalsFunc> compute_normals(ImageViewBase<ViewT> const& image,
                                                                                                              float u_scale, float v_scale) {
  return UnaryPerPixelAccessorView<EdgeExtensionView<ViewT,ConstantEdgeExtension>, ComputeNormalsFunc>(edge_extend(image.impl(), ConstantEdgeExtension()),
                                                                                                       ComputeNormalsFunc (u_scale, v_scale));
}

class DotProdFunc : public ReturnFixedType<PixelMask<PixelGray<float> > > {
  Vector3f m_vec;
public:
  DotProdFunc(Vector3f const& vec) : m_vec(vec) {}
  PixelMask<PixelGray<float> > operator() (PixelMask<Vector3f> const& pix) const {
    if (is_transparent(pix))
      return PixelMask<PixelGray<float> >();
    else
      return dot_prod(pix.child(),m_vec)/(norm_2(pix.child()) * norm_2(m_vec));
  }
};

template <class ViewT>
UnaryPerPixelView<ViewT, DotProdFunc> dot_prod(ImageViewBase<ViewT> const& view, Vector3f const& vec) {
  return UnaryPerPixelView<ViewT, DotProdFunc>(view.impl(), DotProdFunc(vec));
}

// ----------------------------------------------------------------------------

void do_hillshade(po::variables_map const& vm) {

  cartography::GeoReference georef;
  cartography::read_georeference(georef, input_file_name);

  // Select the pixel scale.
  float u_scale, v_scale;
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

  // Set the direction of the light source.
  Vector3f light_0(1,0,0);
  Vector3f light = vw::math::euler_to_rotation_matrix(elevation*M_PI/180, azimuth*M_PI/180, 0, "yzx") * light_0;

  // Compute the surface normals
  std::cout << "Loading: " << input_file_name << ".\n";
  DiskImageView<PixelGray<float> > disk_dem_file(input_file_name);

  ImageViewRef<PixelMask<PixelGray<float> > > dem;
  SrcImageResource *disk_dem_rsrc = DiskImageResource::open(input_file_name);
  if (vm.count("nodata-value")) {
    vw_out() << "\t--> Masking pixel value: " << nodata_value << ".\n";
    dem = create_mask(disk_dem_file, nodata_value);
  } else if ( disk_dem_rsrc->has_nodata_read() ) {
    nodata_value = disk_dem_rsrc->nodata_read();
    vw_out() << "\t--> Extracted nodata value from file: "
             << nodata_value << ".\n";
    dem = create_mask(disk_dem_file, nodata_value);
  } else {
    dem = pixel_cast<PixelMask<PixelGray<float> > >(disk_dem_file);
  }
  delete disk_dem_rsrc;

  if (vm.count("blur")) {
    vw_out() << "\t--> Blurring pixel with gaussian kernal.  Sigma = "
             << blur_sigma << "\n";
    dem = gaussian_filter(dem, blur_sigma);
  }

  // The final result is the dot product of the light source with the normals
  ImageViewRef<PixelMask<PixelGray<uint8> > > shaded_image =
    channel_cast_rescale<uint8>(clamp(dot_prod(compute_normals(dem, u_scale, v_scale), light)));

  // Save the result
  vw_out() << "Writing shaded relief image: " << output_file_name << "\n";

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
    ("input-file", po::value(&input_file_name), "Explicitly specify the input file")
    ("output-file,o", po::value(&output_file_name), "Specify the output file")
    ("azimuth,a", po::value(&azimuth)->default_value(0), "Sets the direction tha the light source is coming from (in degrees).  Zero degrees is to the right, with positive degree counter-clockwise.")
    ("elevation,e", po::value(&elevation)->default_value(45), "Set the elevation of the light source (in degrees).")
    ("scale,s", po::value(&scale)->default_value(0), "Set the scale of a pixel (in the same units as the DTM height values.")
    ("nodata-value", po::value(&nodata_value), "Remap the DEM default value to the min altitude value.")
    ("blur", po::value(&blur_sigma), "Pre-blur the DEM with the specified sigma.");
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
    do_hillshade(vm);
  } catch( Exception& e ) {
    std::cout << "Error: " << e.what() << std::endl;
  }

  return 0;
}
