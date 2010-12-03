// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/Plate/PolarStereoPlateManager.h>
#include <vw/Cartography/GeoReference.h>
#include <vw/Cartography/GeoTransform.h>
#include <vw/Image/Filter.h>

using namespace vw::platefile;
using namespace vw;

template <class PixelT>
cartography::GeoReference
PolarStereoPlateManager<PixelT>::georeference( int level, bool north_pole,
                                               cartography::Datum const& datum) const {
  cartography::GeoReference output_georef( datum );
  output_georef.set_stereographic( north_pole ? 90 : -90, 0, 1.0, 0, 0 );
  Matrix3x3 transform = math::identity_matrix<3>();
  double pixels_per_meters =
    256*pow(2,level)/(2*datum.semi_major_axis());
  transform(0,0) = 1/pixels_per_meters;
  transform(1,1) = -1/pixels_per_meters;
  transform(0,2) = -datum.semi_major_axis();
  transform(1,2) = datum.semi_major_axis();
  output_georef.set_transform( transform );
  return output_georef;
}

template <class PixelT>
cartography::GeoReference
PolarStereoPlateManager<PixelT>::georeference( int level ) const {
  vw_out(WarningMessage, "plate") << "Return PolarStereoGraphic georeference that is north pole regardless of data!";
  return this->georeference( level, true, cartography::Datum("WGS84") );
}

template <class PixelT>
void PolarStereoPlateManager<PixelT>::transform_image(
                          cartography::GeoReference const& georef,
                          ImageViewRef<PixelT>& image,
                          TransformRef& txref, int& level ) const {
  // Determine if input is North or South Pole from points
  Vector2 test_points[5];
  test_points[0] = Vector2( image.cols()/2, image.rows()/2 );
  test_points[1] = Vector2( image.cols()*3/4, image.rows()/2 );
  test_points[2] = Vector2( image.cols()*1/4, image.rows()/2 );
  test_points[3] = Vector2( image.cols()/2, image.rows()*3/4 );
  test_points[4] = Vector2( image.cols()/2, image.rows()*1/4 );
  uint8 north_count = 0;
  for ( uint8 i = 0; i < 5; i++ )
    if ( georef.pixel_to_lonlat( test_points[i] )[1] > 0 )
      north_count++;
  bool is_north = north_count > 2;

  // Work out output resolution from 5 points
  cartography::GeoReference output_georef( georef.datum() );
  output_georef.set_stereographic( is_north ? 90.0 : -90.0, 0, 1.0, 0, 0 );
  {
    Matrix3x3 transform = math::identity_matrix<3>();
    transform(1,1) = -1;
    transform(0,2) = -georef.datum().semi_major_axis();
    transform(1,2) = georef.datum().semi_major_axis();
    output_georef.set_transform( transform );
  }
  cartography::GeoTransform geotx( georef, output_georef );
  // We are seeding pixel_per_meters with the lowest resolution possible
  double requested_pixels_per_meters=256.0/(2*georef.datum().semi_major_axis());
  for ( uint i = 0; i < 5; i++ ) {
    Vector2 i_pos = geotx.forward( test_points[i] );
    Vector2 x_res = geotx.forward( test_points[i]+Vector2(1,0) ) - i_pos;
    Vector2 y_res = geotx.forward( test_points[i]+Vector2(0,1) ) - i_pos;
    double i_resolution = 1.0 / std::min( norm_2(x_res), norm_2(y_res) );
    if ( i_resolution > requested_pixels_per_meters )
      requested_pixels_per_meters = i_resolution;
  }

  // Fit requested_pixels_per_meters to the nearest (256*2^n) / (2*major)
  level =
    boost::numeric_cast<int>(ceil(log(requested_pixels_per_meters*2*georef.datum().semi_major_axis()/256)/log(2)));
  output_georef = this->georeference(level,is_north,georef.datum());

  geotx = cartography::GeoTransform( georef, output_georef );
  BBox2i output_bbox = geotx.forward_bbox( bounding_box(image) );
  vw_out() << "\t    Placing image at level " << level
           << " with bbox " << output_bbox << "\n"
           << "\t    (Total Stereographic resolution at this level =  "
           << requested_pixels_per_meters*2*georef.datum().semi_major_axis() << " pixels.)\n";
  if ( is_north )
    vw_out() << "\t    This is a North Pole image.\n";
  else
    vw_out() << "\t    This is a South Pole image.\n";

  // Perform transform and rewrite to input
  ImageViewRef<PixelT> holding =
    transform( image, geotx, ZeroEdgeExtension(),
               BicubicInterpolation() );
  image = holding;
  txref = TransformRef( geotx );
}

// Explicit template instantiation
namespace vw {
namespace platefile {

#define VW_INSTANTIATE_POLAR_STEREO_TYPES(PIXELT)                            \
  template void                                                              \
  PolarStereoPlateManager<PIXELT >::transform_image(                         \
                                    cartography::GeoReference const& georef, \
                                    ImageViewRef<PIXELT >& image,            \
                                    TransformRef& txref, int& level ) const; \
  template cartography::GeoReference                                         \
  PolarStereoPlateManager<PIXELT >::georeference(int level,                  \
                    bool north_pole, cartography::Datum const& datum) const; \
  template cartography::GeoReference                                         \
  PolarStereoPlateManager<PIXELT >::georeference(int level) const;

  VW_INSTANTIATE_POLAR_STEREO_TYPES(PixelGrayA<uint8>)
  VW_INSTANTIATE_POLAR_STEREO_TYPES(PixelGrayA<int16>)
  VW_INSTANTIATE_POLAR_STEREO_TYPES(PixelGrayA<float32>)
  VW_INSTANTIATE_POLAR_STEREO_TYPES(PixelRGBA<uint8>)

}}

