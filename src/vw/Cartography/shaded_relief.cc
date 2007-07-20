// __BEGIN_LICENSE__
// 
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
// 
// Copyright 2006 Carnegie Mellon University. All rights reserved.
// 
// This software is distributed under the NASA Open Source Agreement
// (NOSA), version 1.3.  The NOSA has been approved by the Open Source
// Initiative.  See the file COPYING at the top of the distribution
// directory tree for the complete NOSA document.
// 
// THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
// KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
// LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
// SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
// A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
// THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
// DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
// 
// __END_LICENSE__
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

#include <vw/Image/ImageView.h>
#include <vw/Image/Algorithms.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Image/PerPixelAccessorViews.h>
#include <vw/Math/EulerAngles.h>
#include <vw/FileIO/DiskImageView.h>
#include <vw/Cartography/GeoReference.h>
#include <vw/Cartography/FileIO.h>

using namespace vw;

//  compute_normals()
//
// Compute a vector normal to the surface of a DEM for each given
// pixel.  The normal is computed by forming a plane with three points
// in the vicinity of the requested pixel, and then finding the vector
// normal to that plane.  The user must specify the scale in the [u,v]
// directions so that the direction of the vector in physical space
// can be properly ascertained.  This is often contained in the (0,0)
// and (1,1) entry of the georeference transform.
class ComputeNormalsFunc : public ReturnFixedType<Vector3> 
{
  double m_u_scale, m_v_scale;

public:
  ComputeNormalsFunc(double u_scale, double v_scale) : 
    m_u_scale(u_scale), m_v_scale(v_scale) {}

  BBox2i work_area() const { return BBox2i(Vector2i(0, 0), Vector2i(1, 1)); }
    
  template <class PixelAccessorT>
  Vector3 operator() (PixelAccessorT const& accessor_loc) const {
    PixelAccessorT acc = accessor_loc;
    // Pick out the three altitude values.
    double alt1 = *acc;
    acc.advance(1,0);
    double alt2 = *acc;
    acc.advance(-1,1);
    double alt3 = *acc;
    
    // Form two orthogonal vectors in the plane containing the three
    // altitude points
    Vector3 n1(m_u_scale, 0, alt2-alt1);
    Vector3 n2(0, m_v_scale, alt3-alt1);

    // Return the vector normal to the local plane.
    return cross_prod(n1,n2);
  }
};

template <class ViewT>
UnaryPerPixelAccessorView<EdgeExtensionView<ViewT,ConstantEdgeExtension>, ComputeNormalsFunc> compute_normals(ImageViewBase<ViewT> const& image,
                                                                                                              double u_scale, double v_scale) {
  return UnaryPerPixelAccessorView<EdgeExtensionView<ViewT,ConstantEdgeExtension>, ComputeNormalsFunc>(edge_extend(image.impl(), ConstantEdgeExtension()),
                                                                                                       ComputeNormalsFunc (u_scale, v_scale)); 
}

int main( int argc, char *argv[] ) {

  set_debug_level(InfoMessage);

  std::string input_file_name, output_file_name;
  float azimuth, elevation, scale, clamp_range;

  po::options_description desc("Options");
  desc.add_options()
    ("help", "Display this help message")
    ("input-file", po::value<std::string>(&input_file_name), "Explicitly specify the input file")
    ("output-file,o", po::value<std::string>(&output_file_name)->default_value("shaded-relief.tif"), "Specify the output file")
    ("azimuth,a", po::value<float>(&azimuth)->default_value(0), "Sets the direction tha the light source is coming from.  Zero degrees is to the right, with positive degree counter-clockwise.")
    ("elevation,e", po::value<float>(&elevation)->default_value(0), "Set the elevation of the light source.")
    ("scale,s", po::value<float>(&scale)->default_value(0), "Set the scale of a pixel (in the same units as the DTM height values.")
    ("clamp-range", po::value<float>(&clamp_range)->default_value(4), "Set the range of floating point values to clamp to prior to normalizing.  You can normally leave this setting untouched.")
    ("no-normalize", "Don't normalize the result -- save the original values as a floating point file (if possible).  This is most often used for debugging.");
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

  try {
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
      v_scale = scale;
    }

    // Set the direction of the light source.
    Vector3 light_0(1,0,0);
    Vector3 light = math::euler_to_rotation_matrix(azimuth*M_PI/180, elevation*M_PI/180, 0, "zyx") * light_0;  

    // Compute the surface normals
    std::cout << "Computing normals..." << std::flush;
    DiskImageView<float> disk_dem_file(input_file_name);
    ImageView<Vector3> normals = compute_normals(disk_dem_file, u_scale, v_scale);

    // The final result is the dot product of the light source with the normals
    ImageView<float> result(normals.cols(), normals.rows());
    for (int j = 0; j < result.rows(); ++j) 
      for (int i = 0; i < result.cols(); ++i) 
        result(i,j) = dot_prod(light, normals(i,j));

    // Save the result
    std::cout << " done.\nWriting shaded relief image.\n";
    if (vm.count("no-normalize"))
      write_image(output_file_name, result, TerminalProgressCallback());
    else
      write_image(output_file_name, channel_cast<uint8>(normalize(clamp(result,-clamp_range,clamp_range),0,255)), TerminalProgressCallback());

  } catch( Exception& e ) {
    std::cerr << "Error: " << e.what() << std::endl;
  }

  return 0;
}
