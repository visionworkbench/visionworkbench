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

#ifndef __VW_CARTOGRAPHY_ORTHOIMAGEVIEW_H__
#define __VW_CARTOGRAPHY_ORTHOIMAGEVIEW_H__

#include <vw/Image/ImageViewBase.h>
#include <vw/Image/EdgeExtension.h>
#include <vw/Cartography/GeoReference.h>
#include <vw/Cartography/Datum.h>
#include <vw/Cartography/PointImageManipulation.h>
#include <vw/Camera/CameraModel.h>
#include <vw/Core/CompoundTypes.h>

#include <boost/shared_ptr.hpp>

namespace vw {
namespace cartography {

  /// This image view projects a camera image onto a given digital
  /// elevation model.
  /// 
  /// This image view assumes the dimensions and georeferencing of the
  /// Terrain image (i.e. the DTM), but it assumes the pixel type of
  /// the camera image. 
  template <class TerrainImageT, class CameraImageT, class InterpT, class EdgeT>
  class OrthoImageView : public ImageViewBase<OrthoImageView<TerrainImageT, CameraImageT, InterpT, EdgeT> > {

    TerrainImageT m_terrain;
    CameraImageT m_camera_image_ref;
    GeoReference m_georef;
    boost::shared_ptr<vw::camera::CameraModel> m_camera_model;
    InterpolationView<EdgeExtensionView<CameraImageT, EdgeT>, InterpT> m_camera_image;
    InterpT m_interp_func;
    EdgeT m_edge_func;

  public:
    typedef typename CameraImageT::pixel_type pixel_type;
    typedef const pixel_type result_type;
    typedef ProceduralPixelAccessor<OrthoImageView> pixel_accessor;
    
    OrthoImageView(TerrainImageT const& terrain, GeoReference const& georef,
                   CameraImageT const& camera_image, boost::shared_ptr<vw::camera::CameraModel> camera_model,
                   InterpT const& interp_func, EdgeT const& edge_func) :
      m_terrain(terrain), 
      m_camera_image_ref(camera_image),
      m_georef(georef), 
      m_camera_model(camera_model),
      m_camera_image(interpolate(camera_image, m_interp_func, m_edge_func)),
      m_interp_func(interp_func), 
      m_edge_func(edge_func) {}

    inline int32 cols() const { return m_terrain.cols(); }
    inline int32 rows() const { return m_terrain.rows(); }
    inline int32 planes() const { return m_camera_image.planes(); }
    
    inline pixel_accessor origin() const { return pixel_accessor(*this); }

    /// Querying the point in the OrthoImageView is straight-forward.
    /// First, the georefencing and altitude information is used to
    /// compute a 3D point corresponding to this location in the DTM.
    /// This point is then "imaged" by the camera model and the
    /// resulting pixel location is returned from the camera image.
    inline result_type operator()( int32 i, int32 j, int p=0 ) const { 

      // We need to convert the georefernced positions into a
      // cartesian coordinate system so that they can be imaged by the
      // camera model.  Doing so require we proceed through 3 steps:
      //
      // 1. Convert from the projection used for the terrain into 
      //    lon,lat,altitude.
      // 2. Add in the offset from the datum that was used, which
      //    converts from altitude to planetary radius.
      // 3. Convert to cartesian (xyz) coordinates.
      Vector2 lon_lat( m_georef.pixel_to_lonlat(Vector2(i,j)) );
      Vector3 xyz = m_georef.datum().geodetic_to_cartesian( Vector3( lon_lat.x(), lon_lat.y(), m_terrain(i,j) ) );
                                                           
      // Now we can image the point using the camera model and return
      // the resulting pixel from the camera image.
      Vector2 pix = m_camera_model->point_to_pixel(xyz);
      return m_camera_image(pix[0], pix[1], p);
    }

    /// \cond INTERNAL
    typedef OrthoImageView<typename TerrainImageT::prerasterize_type,
                           typename CameraImageT::prerasterize_type,
                           InterpT, EdgeT> prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i const& bbox ) const {
      return prerasterize_type( m_terrain.prerasterize(bbox), m_georef, 
                                m_camera_image_ref.prerasterize(BBox2i(0,0,m_camera_image_ref.cols(),m_camera_image_ref.rows())), m_camera_model,
                                m_interp_func, m_edge_func);
    }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const { vw::rasterize( prerasterize(bbox), dest, bbox ); }
    /// \endcond

  };


  // -------------------------------------------------------------------------------
  // Functional API
  // -------------------------------------------------------------------------------

  template <class TerrainImageT, class CameraImageT, class InterpT, class EdgeT>
  OrthoImageView<TerrainImageT, CameraImageT, InterpT, EdgeT> orthoproject( ImageViewBase<TerrainImageT> const& terrain_image, 
                                                                            GeoReference const& georef,
                                                                            ImageViewBase<CameraImageT> const& camera_image,
                                                                            boost::shared_ptr<vw::camera::CameraModel> camera_model,
                                                                            InterpT const& interp_func,
                                                                            EdgeT const& edge_extend_func) {
    return OrthoImageView<TerrainImageT, CameraImageT, InterpT, EdgeT>( terrain_image.impl(), georef, camera_image.impl(), camera_model, interp_func, edge_extend_func);
  }
  
} // namespace cartography
  
  /// \cond INTERNAL
  // Type traits 
  template <class TerrainImageT, class CameraImageT, class InterpT, class EdgeT>
  struct IsFloatingPointIndexable< cartography::OrthoImageView<TerrainImageT, CameraImageT, InterpT, EdgeT> > : public true_type {};
  /// \endcond

} // namespace vw

#endif // __VW_CARTOGRAPHY_ORTHOIMAGEVIEW_H__
