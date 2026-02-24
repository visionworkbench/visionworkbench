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


#include <vw/Image/ImageView.h>
#include <vw/Math/BresenhamLine.h>
#include <vw/Cartography/GeoReference.h>
#include <vw/Cartography/GeoReferenceUtils.h>
#include <vw/FileIO/GdalWriteOptions.h>
using namespace vw;

void grow_bbox( math::BresenhamLine l,
                cartography::GeoReference const& georef,
                BBox2& box ) {
  while ( l.is_good() ) {
    try {
      box.grow( georef.pixel_to_lonlat( *l ) );
    } catch ( const cartography::ProjectionErr& e ) {}
    ++l;
  }
}

void render_degree( ImageView<PixelGray<uint8> >& image,
                    cartography::GeoReference const& georef ) {

  // Find our lat lon boundary. Notice how all our pixels
  // have 0<= x < image.cols() and 0 <= y < image.rows().
  BBox2 ll_box;
  grow_bbox( math::BresenhamLine( Vector2(),
                                  Vector2(image.cols(),image.rows()) ),
             georef, ll_box );
  grow_bbox( math::BresenhamLine( Vector2(image.cols()-1,0),
                                  Vector2(0,image.rows()-1) ),
             georef, ll_box );
  grow_bbox( math::BresenhamLine( Vector2(),
                                  Vector2(0,image.rows()) ),
             georef, ll_box );
  grow_bbox( math::BresenhamLine( Vector2(),
                                  Vector2(image.cols(),0) ),
             georef, ll_box );
  grow_bbox( math::BresenhamLine( Vector2(image.cols()-1,0),
                                  Vector2(image.cols()-1,image.rows()-1) ),
             georef, ll_box );
  grow_bbox( math::BresenhamLine( Vector2(0,image.rows()-1),
                                  Vector2(image.cols()-1,image.rows()-1) ),
             georef, ll_box );

  BBox2i lli_box;
  lli_box.min() = floor( ll_box.min() );
  lli_box.max() = ceil( ll_box.max() );
  std::cout << "LL Box: " << ll_box << "\n";
  std::cout << "LLi Box: " << lli_box << "\n";

  BBox2i image_box = bounding_box( image );

  const int32 inc = 2;

  // Iterate through longitude lines
  for ( int32 lon = lli_box.min()[0];
        lon <= lli_box.max()[0]; lon += inc ) {
    for ( int32 lat = lli_box.min()[1];
          lat < lli_box.max()[1]; lat += inc ) {
      Vector2i start = georef.lonlat_to_pixel( Vector2(lon,lat) );
      Vector2i end = georef.lonlat_to_pixel( Vector2(lon,lat+inc) );
      math::BresenhamLine l( start, end );
      while ( l.is_good() ) {
        if ( !image_box.contains( *l ) ) {
          ++l;
          continue;
        }
        image( (*l).x(), (*l).y() ) = PixelGray<uint8>(255);
        ++l;
      }
    }
  }

  // Iterate through latitude lines
  for ( int32 lat = lli_box.min()[1];
        lat <= lli_box.max()[1]; lat += inc ) {
    for ( int32 lon = lli_box.min()[0];
          lon < lli_box.max()[0]; lon += inc ) {
      Vector2i start = georef.lonlat_to_pixel( Vector2(lon,lat) );
      Vector2i end = georef.lonlat_to_pixel( Vector2(lon+inc,lat) );
      math::BresenhamLine l( start, end );
      while ( l.is_good() ) {
        if ( !image_box.contains( *l ) ) {
          ++l;
          continue;
        }
        image( (*l).x(), (*l).y() ) = PixelGray<uint8>(255);
        ++l;
      }
    }
  }
}

int main( int argc, char *argv[] ) {

  { vw_out() << "Writing: North Stereographic\n";
    ImageView<PixelGray<uint8> > input =
      constant_view( PixelGray<uint8>(64), 1024, 1024 );
    cartography::GeoReference georef;
    georef.set_stereographic(90,0,1);
    georef.set_transform( Matrix3x3(1000,0,-412000,0,-1000,512000,0,0,1) );
    render_degree( input, georef, GdalWriteOptions());
    write_gdal_image( "n_stereo.tif", input, georef, GdalWriteOptions());
  }

  { vw_out() << "Writing: South Stereographic\n";
    ImageView<PixelGray<uint8> > input =
      constant_view( PixelGray<uint8>(64), 1024, 1024 );
    cartography::GeoReference georef;
    georef.set_stereographic(-90,0,1);
    georef.set_transform( Matrix3x3(1000,0,-512000,0,-1000,612000,0,0,1) );
    render_degree( input, georef, GdalWriteOptions());
    write_gdal_image( "s_stereo.tif", input, georef, GdalWriteOptions());
  }

  { vw_out() << "Writing: Center -90 Trans Equirectangular\n";
    ImageView<PixelGray<uint8> > input =
      constant_view( PixelGray<uint8>(64), 1024, 1024 );
    cartography::GeoReference georef;
    georef.set_transform( Matrix3x3(0.0292968,0,-105,0,-0.0292968,15,0,0,1) );
    render_degree( input, georef, GdalWriteOptions());
    write_gdal_image( "n90_equi.tif", input, georef, GdalWriteOptions());
  }

  { vw_out() << "Writing: Center 1 Trans Equirectangular\n";
    ImageView<PixelGray<uint8> > input =
      constant_view( PixelGray<uint8>(64), 1024, 1024 );
    cartography::GeoReference georef;
    georef.set_transform( Matrix3x3(0.0292968,0,-14,0,-0.0292968,15,0,0,1) );
    render_degree( input, georef, GdalWriteOptions());
    write_gdal_image( "p1_equi.tif", input, georef, GdalWriteOptions());
  }

  { vw_out() << "Writing: Center 90 Trans Equirectangular\n";
    ImageView<PixelGray<uint8> > input =
      constant_view( PixelGray<uint8>(64), 1024, 1024 );
    cartography::GeoReference georef;
    georef.set_transform( Matrix3x3(0.0292968,0,75,0,-0.0292968,7.5,0,0,1) );
    render_degree( input, georef, GdalWriteOptions());
    write_gdal_image( "p90_equi.tif", input, georef, GdalWriteOptions());
  }

  { vw_out() << "Writing: Center 180 Trans Equirectangular\n";
    ImageView<PixelGray<uint8> > input =
      constant_view( PixelGray<uint8>(64), 1024, 1024 );
    cartography::GeoReference georef;
    georef.set_transform( Matrix3x3(0.0292968,0,165,0,-0.0292968,7.5,0,0,1) );
    render_degree( input, georef, GdalWriteOptions());
    write_gdal_image( "p180_equi.tif", input, georef, GdalWriteOptions());
  }

  { vw_out() << "Writing: Center 271 Trans Equirectangular\n";
    ImageView<PixelGray<uint8> > input =
      constant_view( PixelGray<uint8>(64), 1024, 1024 );
    cartography::GeoReference georef;
    georef.set_transform( Matrix3x3(0.0292968,0,256,0,-0.0292968,15,0,0,1) );
    render_degree( input, georef, GdalWriteOptions());
    write_gdal_image( "p270_equi.tif", input, georef, GdalWriteOptions());
  }

  { vw_out() << "Writing: Center -90 Proj Equirectangular\n";
    ImageView<PixelGray<uint8> > input =
      constant_view( PixelGray<uint8>(64), 1024, 1024 );
    cartography::GeoReference georef;
    georef.set_equirectangular(0,-90);
    georef.set_transform( Matrix3x3(5000,0,-2560000,0,-5000,2560000,0,0,1) );
    render_degree( input, georef, GdalWriteOptions());
    write_gdal_image( "n90_proj_equi.tif", input, georef, GdalWriteOptions());
  }

  { vw_out() << "Writing: Center 90 Proj Equirectangular\n";
    ImageView<PixelGray<uint8> > input =
      constant_view( PixelGray<uint8>(64), 1024, 1024 );
    cartography::GeoReference georef;
    georef.set_equirectangular(0,90);
    georef.set_transform( Matrix3x3(5000,0,-2560000,0,-5000,2560000,0,0,1) );
    render_degree( input, georef, GdalWriteOptions());
    write_gdal_image( "p90_proj_equi.tif", input, georef, GdalWriteOptions());
  }

  { vw_out() << "Writing: Center 270 Proj Equirectangular\n";
    ImageView<PixelGray<uint8> > input =
      constant_view( PixelGray<uint8>(64), 1024, 1024 );
    cartography::GeoReference georef;
    georef.set_equirectangular(0,270);
    georef.set_transform( Matrix3x3(5000,0,-2560000,0,-5000,2560000,0,0,1) );
    render_degree( input, georef, GdalWriteOptions());
    write_gdal_image( "p270_proj_equi.tif", input, georef, GdalWriteOptions());
  }

  { vw_out() << "Writing: Global [-180,180] Equirectangular\n";
    ImageView<PixelGray<uint8> > input =
      constant_view( PixelGray<uint8>(64), 2048, 1024 );
    cartography::GeoReference georef;
    georef.set_transform( Matrix3x3(.17578125,0,-180,0,-.17578125,90,0,0,1) );
    render_degree( input, georef, GdalWriteOptions());
    write_gdal_image( "global_180_180_equi.tif", input, georef, GdalWriteOptions());
  }

  { vw_out() << "Writing: Global [0,360] Equirectangular\n";
    ImageView<PixelGray<uint8> > input =
      constant_view( PixelGray<uint8>(64), 2048, 1024 );
    cartography::GeoReference georef;
    georef.set_transform( Matrix3x3(.17578125,0,0,0,-.17578125,90,0,0,1) );
    render_degree( input, georef, GdalWriteOptions());
    write_gdal_image( "global_0_360_equi.tif", input, georef, GdalWriteOptions());
  }

  { vw_out() << "Writing: Center -90 Orthoprojection\n";
    ImageView<PixelGray<uint8> > input =
      constant_view( PixelGray<uint8>(64), 1024, 1024 );
    cartography::GeoReference georef;
    georef.set_orthographic(-75,-90);
    georef.set_transform( Matrix3x3(5000,0,-2560000,0,-5000,2560000,0,0,1) );
    render_degree( input, georef, GdalWriteOptions());
    write_gdal_image( "n90_ortho.tif", input, georef, GdalWriteOptions());
  }

  { vw_out() << "Writing: Center 271 Orthoprojection\n";
    ImageView<PixelGray<uint8> > input =
      constant_view( PixelGray<uint8>(64), 1024, 1024 );
    cartography::GeoReference georef;
    georef.set_orthographic(15,271);
    georef.set_transform( Matrix3x3(5000,0,-2560000,0,-5000,2560000,0,0,1) );
    render_degree( input, georef, GdalWriteOptions());
    write_gdal_image( "p271_ortho.tif", input, georef, GdalWriteOptions());
  }

  { vw_out() << "Writing: Center -90 Lambert Azimuthal\n";
    ImageView<PixelGray<uint8> > input =
      constant_view( PixelGray<uint8>(64), 1024, 1024 );
    cartography::GeoReference georef;
    georef.set_lambert_azimuthal(15,-90);
    georef.set_transform( Matrix3x3(5000,0,-2560000,0,-5000,2560000,0,0,1) );
    render_degree( input, georef, GdalWriteOptions());
    write_gdal_image( "n90_lambert_azi.tif", input, georef, GdalWriteOptions());
  }

  { vw_out() << "Writing: Center 270 Lambert Azimuthal\n";
    ImageView<PixelGray<uint8> > input =
      constant_view( PixelGray<uint8>(64), 1024, 1024 );
    cartography::GeoReference georef;
    georef.set_lambert_azimuthal(-35,270);
    georef.set_transform( Matrix3x3(5000,0,-2560000,0,-5000,2560000,0,0,1) );
    render_degree( input, georef, GdalWriteOptions());
    write_gdal_image( "p270_lambert_azi.tif", input, georef, GdalWriteOptions());
  }

  { vw_out() << "Writing: Center -90 Lambert Conformal\n";
    ImageView<PixelGray<uint8> > input =
      constant_view( PixelGray<uint8>(64), 1024, 1024 );
    cartography::GeoReference georef;
    georef.set_lambert_conformal(0,30,15,-90);
    georef.set_transform( Matrix3x3(5000,0,-2560000,0,-5000,2560000,0,0,1) );
    render_degree( input, georef, GdalWriteOptions());
    write_gdal_image( "n90_lambert_con.tif", input, georef, GdalWriteOptions());
  }

  { vw_out() << "Writing: Center 270 Lambert Conformal\n";
    ImageView<PixelGray<uint8> > input =
      constant_view( PixelGray<uint8>(64), 1024, 1024 );
    cartography::GeoReference georef;
    georef.set_lambert_conformal(-70,-40,-55,270);
    georef.set_transform( Matrix3x3(5000,0,-2560000,0,-5000,2560000,0,0,1) );
    render_degree( input, georef, GdalWriteOptions());
    write_gdal_image( "p270_lambert_con.tif", input, georef, GdalWriteOptions());
  }

  { vw_out() << "Writing: Center -90 Transverse Mercator\n";
    ImageView<PixelGray<uint8> > input =
      constant_view( PixelGray<uint8>(64), 1024, 1024 );
    cartography::GeoReference georef;
    georef.set_transverse_mercator(-50,-90,1);
    georef.set_transform( Matrix3x3(5000,0,-2560000,0,-5000,2560000,0,0,1) );
    render_degree( input, georef, GdalWriteOptions());
    write_gdal_image( "n90_trans_merc.tif", input, georef, GdalWriteOptions());
  }

  { vw_out() << "Writing: Center 270 Transverse Mercator\n";
    ImageView<PixelGray<uint8> > input =
      constant_view( PixelGray<uint8>(64), 1024, 1024 );
    cartography::GeoReference georef;
    georef.set_transverse_mercator(10,270,1);
    georef.set_transform( Matrix3x3(5000,0,-2560000,0,-5000,2560000,0,0,1) );
    render_degree( input, georef, GdalWriteOptions());
    write_gdal_image( "p270_trans_merc.tif", input, georef, GdalWriteOptions());
  }

  { vw_out() << "Writing: UTM 10\n";
    ImageView<PixelGray<uint8> > input =
      constant_view( PixelGray<uint8>(64), 1024, 1024 );
    cartography::GeoReference georef;
    georef.set_UTM( 10, 42 );
    georef.set_transform( Matrix3x3(5000,0,-1000000,0,-5000,2560000,0,0,1) );
    render_degree( input, georef, GdalWriteOptions());
    write_gdal_image( "utm10.tif", input, georef, GdalWriteOptions());
  }

  { vw_out() << "Writing: UTM 60S\n";
    ImageView<PixelGray<uint8> > input =
      constant_view( PixelGray<uint8>(64), 1024, 1024 );
    cartography::GeoReference georef;
    georef.set_UTM( 60, 40 );
    georef.set_transform( Matrix3x3(1000,0,-300000,0,-1000,5e6,0,0,1) );
    render_degree( input, georef, GdalWriteOptions());
    write_gdal_image( "utm60s.tif", input, georef, GdalWriteOptions());
  }

  { vw_out() << "Writing: Center -90 Sinusoidal\n";
    ImageView<PixelGray<uint8> > input =
      constant_view( PixelGray<uint8>(64), 1024, 1024 );
    cartography::GeoReference georef;
    georef.set_sinusoidal(15,-90);
    georef.set_transform( Matrix3x3(5000,0,-2560000,0,-5000,560000,0,0,1) );
    render_degree( input, georef, GdalWriteOptions());
    write_gdal_image( "n90_sin.tif", input, georef, GdalWriteOptions());
  }

  { vw_out() << "Writing: Center 270 Sinusoidal\n";
    ImageView<PixelGray<uint8> > input =
      constant_view( PixelGray<uint8>(64), 1024, 1024 );
    cartography::GeoReference georef;
    georef.set_sinusoidal(-60,270);
    georef.set_transform( Matrix3x3(5000,0,-2560000,0,-5000,560000,0,0,1) );
    render_degree( input, georef, GdalWriteOptions());
    write_gdal_image( "p270_sin.tif", input, georef, GdalWriteOptions());
  }

  return 0;
}
