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


#ifndef __VW_CARTOGRAPHY_ORTHOIMAGEVIEW_H__
#define __VW_CARTOGRAPHY_ORTHOIMAGEVIEW_H__

#include <vw/Math/BresenhamLine.h>
#include <vw/Image/ImageViewBase.h>
#include <vw/Image/EdgeExtension.h>
#include <vw/Image/Interpolation.h>
#include <vw/Cartography/GeoReference.h>
#include <vw/Camera/CameraModel.h>

#include <boost/shared_ptr.hpp>


/// \file OrthoImageView.h Image view that projects a camera image onto a given digital elevation model.

// TODO: We don't use this class for anything, why not??

namespace vw {
namespace cartography {

  /// This image view projects a camera image onto a given digital
  /// elevation model.
  ///
  /// This image view assumes the dimensions and georeferencing of the
  /// Terrain image (i.e. the DTM), but it assumes the pixel type of
  /// the camera image.
  template <class TerrainImageT, class CameraImageT, class InterpT, class EdgeT, bool markNoProcessedData>
  class OrthoImageView : public ImageViewBase<OrthoImageView<TerrainImageT, CameraImageT, InterpT, EdgeT, markNoProcessedData> > {

    typedef typename boost::mpl::if_<IsFloatingPointIndexable<TerrainImageT>, double, int32>::type offset_type;

    TerrainImageT        m_terrain;
    GeoReference         m_georef;
    camera::CameraModel* m_camera_model; // There is a big assumption
                                         // here that the user's
                                         // camera model is thread
                                         // safe.
    InterpolationView<EdgeExtensionView<CameraImageT, EdgeT>, InterpT> m_camera_image;
    CameraImageT m_camera_image_ref;
    InterpT      m_interp_func;
    EdgeT        m_edge_func;

    // Provide safe interaction with DEMs that are scalar or compound
    template <class PixelT>
    typename boost::enable_if< IsScalar<PixelT>, double >::type
    inline Helper( offset_type x, offset_type y ) const {
      return m_terrain(x,y);
    }

    template <class PixelT>
    typename boost::enable_if< IsCompound<PixelT>, double>::type
    inline Helper( offset_type x, offset_type y ) const {
      return m_terrain(x,y)[0];
    }

    inline Vector2 find_camera_coordinates( offset_type i, offset_type j, double height ) const {
      Vector2 lon_lat( m_georef.pixel_to_lonlat(Vector2(i,j)) );
      return m_camera_model->point_to_pixel( m_georef.datum().geodetic_to_cartesian( Vector3( lon_lat.x(), lon_lat.y(), height ) ) );
    }

    inline void apply_bresen( math::BresenhamLine line, double min, double max, BBox2i& camera_bbox ) const {
      while ( line.is_good() ) {
        Vector2i pt( *line );
        try { camera_bbox.grow( find_camera_coordinates( pt.x(), pt.y(), min ) ); }
        catch (...) {} // PointToPixelErr, MathErr, etc.
        try { camera_bbox.grow( find_camera_coordinates( pt.x(), pt.y(), max ) ); }
        catch (...) {} // PointToPixelErr, MathErr, etc.
        ++line;
      }
    }

  public:
    typedef typename CameraImageT::pixel_type pixel_type;
    typedef const pixel_type result_type;
    typedef ProceduralPixelAccessor<OrthoImageView> pixel_accessor;

    OrthoImageView(TerrainImageT const& terrain, GeoReference const& georef,
                   CameraImageT const& camera_image, camera::CameraModel* camera_model,
                   InterpT const& interp_func, EdgeT const& edge_func) :
      m_terrain(terrain), m_georef(georef), m_camera_model(camera_model),
      m_camera_image(interpolate(camera_image, m_interp_func, m_edge_func)),
      m_camera_image_ref(camera_image), m_interp_func(interp_func),
      m_edge_func(edge_func) {}

    inline int32 cols  () const { return m_terrain.cols();        }
    inline int32 rows  () const { return m_terrain.rows();        }
    inline int32 planes() const { return m_camera_image.planes(); }

    inline pixel_accessor origin() const { return pixel_accessor(*this); }

    /// Querying the point in the OrthoImageView is straight-forward.
    /// First, the georefencing and altitude information is used to
    /// compute a 3D point corresponding to this location in the DTM.
    /// This point is then "imaged" by the camera model and the
    /// resulting pixel location is returned from the camera image.

    // We need to convert the georefernced positions into a
    // cartesian coordinate system so that they can be imaged by the
    // camera model.  Doing so require we proceed through 3 steps:
    //
    // 1. Convert from the projection used for the terrain into
    //    lon,lat,altitude.
    // 2. Add in the offset from the datum that was used, which
    //    converts from altitude to planetary radius.
    // 3. Convert to cartesian (xyz) coordinates.

    inline result_type operator()(offset_type i, offset_type j, int32 p=0) const {

      if (!markNoProcessedData){ // This 'if' will be evaluated at compile time

        // Default behavior, don't mark no-processed-data separately from no-data

        // Check for missing DEM pixels.
        if ( is_transparent(m_terrain(i,j)) ) {
          return result_type();
        }

        Vector2 pix;
        try{
          pix = find_camera_coordinates( i, j, Helper<typename TerrainImageT::pixel_type>(i,j) );
        }catch (...) { // PointToPixelErr, MathErr, etc.
          return result_type();
        }

        return m_camera_image(pix[0], pix[1], p);

      }else{

        // Do mark no-processed-data separately from no-data

        // Check for missing DEM pixels.
        double height;
        if (is_transparent(m_terrain(i,j)) ){
          height = 0.0;
        }else{
          height = Helper<typename TerrainImageT::pixel_type>(i,j);
        }

        Vector2 pix;
        try{
          pix = find_camera_coordinates( i, j, height );
        }catch (...) { // PointToPixelErr, MathErr, etc.
          // No data, return a transparent pixel
          return result_type();
        }

        result_type ans = m_camera_image(pix[0], pix[1], p);

        if ( is_transparent(m_terrain(i,j)) ){
          if (0 <= pix[0] && pix[0] < m_camera_image.cols() &&
              0 <= pix[1] && pix[1] < m_camera_image.rows()
              ){
            // No processed data, return a black pixel
            return result_type(0.0);
          }else{
            // No data, return a transparent pixel
            return result_type();
          }
        }
        return ans;
      }

      return result_type();
    }

    /// \cond INTERNAL
    typedef OrthoImageView<typename TerrainImageT::prerasterize_type,
                           typename CameraImageT::prerasterize_type, InterpT, EdgeT, markNoProcessedData> prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i const& bbox ) const {
      // Prerasterize the terrain that we'll be using
      typename TerrainImageT::prerasterize_type terrain_preraster = m_terrain.prerasterize(bbox);

      double terrain_min = std::numeric_limits<double>::max(),
             terrain_max = std::numeric_limits<double>::min();

      // Determine min max
      for ( int32 j = bbox.min().y(); j < bbox.max().y(); j++ ) {
        for ( int32 i = bbox.min().x(); i < bbox.max().x(); i++ ) {
          if ( !is_transparent( terrain_preraster( i,j ) ) ) {
            double val = terrain_preraster(i,j);
            terrain_min = std::min( val, terrain_min );
            terrain_max = std::max( val, terrain_max );
          }
        }
      }

      // Work out what the active area of the camera image that we'll
      // be using. Unfortunately this is no linear. We need to project
      // all the pixels really. For speed I'm only processing an X
      BBox2i camera_bbox;
      apply_bresen( math::BresenhamLine( bbox.min().x(), bbox.min().y(), bbox.max().x(), bbox.min().y() ),
                    terrain_min, terrain_max, camera_bbox );
      apply_bresen( math::BresenhamLine( bbox.min().x(), bbox.max().y(), bbox.max().x(), bbox.max().y() ),
                    terrain_min, terrain_max, camera_bbox );
      apply_bresen( math::BresenhamLine( bbox.min().x(), bbox.min().y(), bbox.min().x(), bbox.max().y() ),
                    terrain_min, terrain_max, camera_bbox );
      apply_bresen( math::BresenhamLine( bbox.max().x(), bbox.min().y(), bbox.max().x(), bbox.max().y() ),
                    terrain_min, terrain_max, camera_bbox );
      apply_bresen( math::BresenhamLine( bbox.min().x(), bbox.min().y(), bbox.max().x(), bbox.max().y() ),
                    terrain_min, terrain_max, camera_bbox );
      apply_bresen( math::BresenhamLine( bbox.min().x(), bbox.max().y(), bbox.max().x(), bbox.min().y() ),
                    terrain_min, terrain_max, camera_bbox );
      camera_bbox.max() += Vector2i(1,1);        // Because grow is
                                                 // inclusive and we
                                                 // need exclusive
      camera_bbox.expand(InterpT::pixel_buffer); // Fudge factor
      camera_bbox.crop( bounding_box( m_camera_image_ref ) );
      if ( camera_bbox.width() * camera_bbox.height() == 0 )
        camera_bbox = BBox2i(0,0,0,0);

      // Prerasterize the parts of the camera image we'll be
      // using. This is important as otherwise we'll just be waiting
      // on cache repeatedly.
      return prerasterize_type( m_terrain.prerasterize(bbox),
                                m_georef, m_camera_image_ref.prerasterize(camera_bbox),
                                m_camera_model, m_interp_func, m_edge_func);
    }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const { vw::rasterize( prerasterize(bbox), dest, bbox ); }
    /// \endcond
  };

  // --------------------------------------------------------------------------
  // Functional API
  // --------------------------------------------------------------------------
  template <class TerrainImageT, class CameraImageT, class InterpT, class EdgeT>
  OrthoImageView<TerrainImageT, CameraImageT, InterpT, EdgeT, false>
  orthoproject(ImageViewBase<TerrainImageT> const& terrain_image,
               GeoReference const& georef,
               ImageViewBase<CameraImageT> const& camera_image,
               camera::CameraModel* camera_model,
               InterpT const& interp_func,
               EdgeT const& edge_extend_func) {
    return OrthoImageView<TerrainImageT, CameraImageT, InterpT, EdgeT, false>( terrain_image.impl(), georef, camera_image.impl(), camera_model, interp_func, edge_extend_func);
  }

  // A special version of orthoproject which will distinguish no-processed-data
  // from no-data.
  template <class TerrainImageT, class CameraImageT, class InterpT, class EdgeT>
  OrthoImageView<TerrainImageT, CameraImageT, InterpT, EdgeT, true>
  orthoproject_markNoProcessedData(ImageViewBase<TerrainImageT> const& terrain_image,
                                   GeoReference const& georef,
                                   ImageViewBase<CameraImageT> const& camera_image,
                                   camera::CameraModel* camera_model,
                                   InterpT const& interp_func,
                                   EdgeT const& edge_extend_func) {
    return OrthoImageView<TerrainImageT, CameraImageT, InterpT, EdgeT, true>( terrain_image.impl(), georef, camera_image.impl(), camera_model, interp_func, edge_extend_func);
  }

} // namespace cartography

  /// \cond INTERNAL
  // Type traits
  template <class TerrainImageT, class CameraImageT, class InterpT, class EdgeT, bool markNoProcessedData>
  struct IsFloatingPointIndexable< cartography::OrthoImageView<TerrainImageT, CameraImageT, InterpT, EdgeT, markNoProcessedData> > : public IsFloatingPointIndexable<TerrainImageT> {};
  /// \endcond

} // namespace vw

#endif // __VW_CARTOGRAPHY_ORTHOIMAGEVIEW_H__
