// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NASA Vision Workbench is licensed under the Apache License,
//  Version 2.0 (the "License"); you may not use this file except in
//  compliance with the License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__


#ifndef __VW_TOOLS_HILLSHADE_H__
#define __VW_TOOLS_HILLSHADE_H__

#include <vw/Core/System.h>
#include <vw/Core/Log.h>
#include <vw/Core/FundamentalTypes.h>
#include <vw/Math/EulerAngles.h>
#include <vw/Image/Algorithms.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Image/Filter.h>
#include <vw/Image/PixelMask.h>
#include <vw/Image/MaskViews.h>
#include <vw/Image/PerPixelAccessorViews.h>
#include <vw/FileIO/DiskImageView.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/Cartography/GeoReference.h>
#include <vw/tools/Common.h>

namespace vw{

  template <class PixelT>
  void do_hillshade(std::string const& input_file_name,
                    std::string const& output_file_name,
                    double azimuth, double elevation, double scale,
                    double nodata_value, double blur_sigma) {

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
    vw_out() << "Loading: " << input_file_name << ".\n";
    DiskImageView<PixelT > disk_dem_file(input_file_name);

    ImageViewRef<PixelMask<PixelT > > dem;
    boost::scoped_ptr<SrcImageResource> disk_dem_rsrc(DiskImageResource::open(input_file_name));
    if ( !std::isnan(nodata_value) ) {
      vw_out() << "\t--> Masking pixel value: " << nodata_value << ".\n";
      dem = create_mask(disk_dem_file, nodata_value);
    } else if ( disk_dem_rsrc->has_nodata_read() ) {
      nodata_value = disk_dem_rsrc->nodata_read();
      vw_out() << "\t--> Extracted nodata value from file: "
               << nodata_value << ".\n";
      dem = create_mask(disk_dem_file, nodata_value);
    } else {
      dem = pixel_cast<PixelMask<PixelT > >(disk_dem_file);
    }

    if ( !std::isnan(blur_sigma) ) {
      vw_out() << "\t--> Blurring pixel with gaussian kernal.  Sigma = "
               << blur_sigma << "\n";
      dem = gaussian_filter(dem, blur_sigma);
    }

    // The final result is the dot product of the light source with the normals
    ImageViewRef<PixelMask<PixelGray<uint8> > > shaded_image =
      channel_cast_rescale<uint8>(clamp(dot_prod(compute_normals(dem, u_scale, v_scale), light)));

    // Save the result
    vw_out() << "Writing shaded relief image: " << output_file_name << "\n";

    boost::scoped_ptr<DiskImageResource> r(DiskImageResource::create(output_file_name,
                                                                     shaded_image.format()));
    if ( r->has_block_write() )
      r->set_block_write_size( Vector2i( vw_settings().default_tile_size(),
                                         vw_settings().default_tile_size() ) );
    write_georeference( *r, georef );
    block_write_image( *r, shaded_image,
                       TerminalProgressCallback( "tools.hillshade", "Writing:") );
  }

  void do_multitype_hillshade(std::string const& input_file,
                              std::string const& output_file,
                              double azimuth, double elevation, double scale,
                              double nodata_value, double blur_sigma) {

    ImageFormat fmt = vw::image_format(input_file);

    switch(fmt.pixel_format) {
    case VW_PIXEL_GRAY:
    case VW_PIXEL_GRAYA:
      switch (fmt.channel_type) {

      case VW_CHANNEL_UINT8:
        do_hillshade<PixelGray<uint8>  >(input_file, output_file,
                                         azimuth, elevation, scale,
                                         nodata_value, blur_sigma);
        break;
      case VW_CHANNEL_INT16:
        do_hillshade<PixelGray<int16>  >(input_file, output_file,
                                         azimuth, elevation, scale,
                                         nodata_value, blur_sigma);
        break;
      case VW_CHANNEL_UINT16:
        do_hillshade<PixelGray<uint16> >(input_file, output_file,
                                         azimuth, elevation, scale,
                                         nodata_value, blur_sigma);
        break;
      default:
        do_hillshade<PixelGray<float>  >(input_file, output_file,
                                         azimuth, elevation, scale,
                                         nodata_value, blur_sigma);
        break;
      }
      break;
    default:
      vw_throw( ArgumentErr() << "Unsupported pixel format. The DEM image must have only one channel." );
    }
  }

}

#endif
