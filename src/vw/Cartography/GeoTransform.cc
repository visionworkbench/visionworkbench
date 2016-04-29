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


#include <vw/Math/BresenhamLine.h>
#include <vw/Cartography/GeoTransform.h>

// Vision Workbench
#include <vw/Image/ImageView.h>

// Proj.4
#include <proj_api.h>

// Macro for checking Proj.4 output, something we do a lot of.
#define CHECK_PROJ_ERROR(ctx_input) if(ctx_input.error_no()) vw_throw(ProjectionErr() << "Proj.4 error: " << pj_strerrno(ctx_input.error_no()))

namespace vw {
namespace cartography {

  using vw::math::BresenhamLine;

  // Constructor
  GeoTransform::GeoTransform(GeoReference const& src_georef, GeoReference const& dst_georef,
                             BBox2i const& src_bbox, BBox2i const& dst_bbox) :
    m_src_georef(src_georef), m_dst_georef(dst_georef), m_offset(Vector2(0, 0)) {

    // Deal with the fact that longitudes could differ by 360 degrees
    // between src and dst.
    if (src_bbox.empty() || dst_bbox.empty()) {
      // If we don't know the image areas, simply use the 0, 0 corner
      try {
	Vector2 src_origin = src_georef.pixel_to_lonlat(Vector2(0, 0));
	Vector2 dst_origin = dst_georef.pixel_to_lonlat(Vector2(0, 0));
	m_offset = Vector2( 360.0*round( (dst_origin[0] - src_origin[0])/360.0 ), 0.0 );
      }catch ( const std::exception & e ) {
	m_offset = Vector2(0, 0);
      }
    } else{
	
      // Try to offset by 360 degrees until the lon-lat boxes are most compatible,
      // measured by how close the centers of the boxes are.
	
      BBox2 src_lonlat_box = src_georef.pixel_to_lonlat_bbox(src_bbox);
      BBox2 dst_lonlat_box = dst_georef.pixel_to_lonlat_bbox(dst_bbox);
      std::vector<double> shift, dist_vec;
      for (int val = -360; val <= 360; val+= 360) {
	shift.push_back(val);
	BBox2 shifted_src = src_lonlat_box + Vector2(val, 0);
	  
	// Distance between centers of boxes. The smaller it is, the better.
	double dist =
	  norm_2((shifted_src.min() + shifted_src.max())/2.0
		 - (dst_lonlat_box.min() + dst_lonlat_box.max())/2.0);
	dist_vec.push_back(dist);
      }
      int min_index = std::distance(dist_vec.begin(),
				    min_element(dist_vec.begin(), dist_vec.end()));
      m_offset = Vector2(shift[min_index], 0);
    }
    
    const std::string src_datum = m_src_georef.datum().proj4_str();
    const std::string dst_datum = m_dst_georef.datum().proj4_str();

    // This optimizes in the common case where the two images are
    // already in the same map projection, and we need only apply
    // the affine transform.  This will break, of course, as soon as
    // we have mare than one type of GeoReference object, but it
    // makes life faster for now. -mbroxton
    if (m_src_georef.proj4_str() == m_dst_georef.proj4_str())
      m_skip_map_projection = true;
    else
      m_skip_map_projection = false;

    // Bugfix: If one of the projections is lon-lat, and the offset is 360 degrees,
    // don't skip map-projection, as have to correct for the offset.
    bool has_lonat = (m_src_georef.proj4_str().find("+proj=longlat") != std::string::npos ||
                      m_dst_georef.proj4_str().find("+proj=longlat") != std::string::npos );
    if (has_lonat && m_offset != Vector2()) {
      m_skip_map_projection = false;
    }

    // This optimizes the case where the two datums are the same,
    // and thus we don't need to call proj to convert between them
    // as we transform.
    if(src_datum == dst_datum) {
      m_skip_datum_conversion = true;
    } else {
      // Set up the various variables for proj.
      m_skip_datum_conversion = false;

      // The source proj4 context.
      std::stringstream ss_src;
      // We convert lat/long to lat/long regardless of what the
      // source or destination georef uses.
      ss_src << "+proj=latlong " << src_datum;
      m_src_proj = ProjContext( ss_src.str() );

      // The destination proj4 context.
      std::stringstream ss_dst;
      ss_dst << "+proj=latlong " << dst_datum;
      m_dst_proj = ProjContext( ss_dst.str() );
    }
    // Because GeoTransform is typically very slow, we default to a tolerance
    // of 0.1 pixels to allow ourselves to be approximated.
    set_tolerance( 0.1 );
  }

  void GeoTransform::set_offset(Vector2 const& offset){
    m_offset = offset;
  }

  // Performs a forward or reverse datum conversion.
  Vector2 GeoTransform::datum_convert(Vector2 const& v, bool forward) const {
    double x = v[0];
    double y = v[1];
    double z = 0;

    if(forward)
      pj_transform(m_src_proj.proj_ptr(), m_dst_proj.proj_ptr(), 1, 0, &x, &y, &z);
    else
      pj_transform(m_dst_proj.proj_ptr(), m_src_proj.proj_ptr(), 1, 0, &x, &y, &z);
    CHECK_PROJ_ERROR( m_src_proj );
    CHECK_PROJ_ERROR( m_dst_proj );

    return Vector2(x, y);
  }

  BBox2i GeoTransform::forward_bbox( BBox2i const& bbox ) const {
    BBox2 r = TransformHelper<GeoTransform,ContinuousFunction,ContinuousFunction>::forward_bbox(bbox);
    BresenhamLine l1( bbox.min(), bbox.max() );
    while ( l1.is_good() ) {
      try {
        r.grow( this->forward( *l1 ) );
      } catch ( const cartography::ProjectionErr& e ) {}
      ++l1;
    }
    BresenhamLine l2( bbox.min() + Vector2i(bbox.width(),0),
        bbox.max() + Vector2i(-bbox.width(),0) );
    while ( l2.is_good() ) {
      try {
        r.grow( this->forward( *l2 ) );
      } catch ( const cartography::ProjectionErr& e ) {}
      ++l2;
    }

    return grow_bbox_to_int(r);
  }

  BBox2i GeoTransform::reverse_bbox( BBox2i const& bbox ) const {
    BBox2 r = TransformHelper<GeoTransform,ContinuousFunction,ContinuousFunction>::reverse_bbox(bbox);
    BresenhamLine l1( bbox.min(), bbox.max() );
    while ( l1.is_good() ) {
      try {
        r.grow( this->reverse( *l1 ) );
      } catch ( const cartography::ProjectionErr& e ) {}
      ++l1;
    }
    BresenhamLine l2( bbox.min() + Vector2i(bbox.width(),0),
        bbox.max() + Vector2i(-bbox.width(),0) );
    while ( l2.is_good() ) {
      try {
        r.grow( this->reverse( *l2 ) );
      } catch ( const cartography::ProjectionErr& e ) {}
      ++l2;
    }

    return grow_bbox_to_int(r);
  }

  // Convert a pixel in respect to src_georef to a point (hence in projected coordinates)
  // in respect to dst_georef.
  Vector2 GeoTransform::forward_pixel_to_point( Vector2 const& v ) const {
    if (m_skip_map_projection)
      return m_src_georef.pixel_to_point(v);
    Vector2 dst_lonlat = m_src_georef.pixel_to_lonlat(v) + m_offset;
    if(m_skip_datum_conversion)
      return m_dst_georef.lonlat_to_point(dst_lonlat);
    dst_lonlat = datum_convert(dst_lonlat, true);
    return m_dst_georef.lonlat_to_point(dst_lonlat);
  }

  // Convert a pixel box in respect to src_georef to a point box
  // in respect to dst_georef.
  BBox2 GeoTransform::forward_pixel_to_point_bbox( BBox2i const& pixel_bbox ) const {

    BBox2 point_bbox;

    // Go along the perimeter of the pixel bbox.
    for ( int32 x=pixel_bbox.min().x(); x<pixel_bbox.max().x(); ++x ) {
      try {
        point_bbox.grow(this->forward_pixel_to_point( Vector2(x,pixel_bbox.min().y())   ));
        point_bbox.grow(this->forward_pixel_to_point( Vector2(x,pixel_bbox.max().y()-1) ));
      } catch ( const cartography::ProjectionErr& e ) {}
    }
    for ( int32 y=pixel_bbox.min().y()+1; y<pixel_bbox.max().y()-1; ++y ) {
      try {
        point_bbox.grow(this->forward_pixel_to_point( Vector2(pixel_bbox.min().x(),y)   ));
        point_bbox.grow(this->forward_pixel_to_point( Vector2(pixel_bbox.max().x()-1,y) ));
      } catch ( const cartography::ProjectionErr& e ) {}
    }

    // Draw an X inside the bbox. This covers the poles. It will
    // produce a lonlat boundary that is within at least one pixel of
    // the pole. This will also help catch terminator boundaries from
    // orthographic projections.
    BresenhamLine l1( pixel_bbox.min(), pixel_bbox.max() );
    while ( l1.is_good() ) {
      try {
        point_bbox.grow( this->forward_pixel_to_point( *l1 ) );
      } catch ( const cartography::ProjectionErr& e ) {}
      ++l1;
    }
    BresenhamLine l2( pixel_bbox.min() + Vector2i(pixel_bbox.width(),0),
                      pixel_bbox.max() - Vector2i(pixel_bbox.width(),0) );
    while ( l2.is_good() ) {
      try {
        point_bbox.grow( this->forward_pixel_to_point( *l2 ) );
      } catch ( const cartography::ProjectionErr& e ) {}
      ++l2;
    }

    return point_bbox;
  }

  void reproject_point_image(ImageView<Vector3> const& point_image,
                             GeoReference const& src_georef,
                             GeoReference const& dst_georef) {

    GeoTransform gtx(src_georef, dst_georef);

    // Iterate over the image, transforming the first two coordinates
    // in the Vector one at a time.  The third coordinate is taken to
    // be the altitude value, and this value is not touched.
    for (int32 j=0; j < point_image.rows(); ++j) {
      for (int32 i=0; i < point_image.cols(); ++i) {
        if (point_image(i,j) != Vector3()) {
          Vector2 in(point_image(i,j)[0], point_image(i,j)[1]);
          Vector2 out = gtx.forward(in);
          point_image(i,j).x() = out[0];
          point_image(i,j).y() = out[1];
        }
      }
    }
  }
}} // namespace vw::cartography

#undef CHECK_PROJ_ERROR
