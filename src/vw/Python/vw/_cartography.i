// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


%module cartography

%{
#define SWIG_FILE_WITH_INIT
#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Cartography.h>
%}

%include "std_string.i"
%include "carrays.i"
%include "numpy.i"

%import "_image.i"
%import "_fileio.i"
%import "_transform.i"

%init %{
  import_array();
%}

%array_class(double, doublea);

%pythoncode {
  from numpy import array

  def array_to_doublea(M):
    m = doublea(9)
    for i in xrange(3):
      for j in xrange(3):
        m[3*i+j] = M[i,j]
    return m

  def doublea_to_array(M):
    return array(((M[0],M[1],M[2]),(M[3],M[4],M[5]),(M[6],M[7],M[8])))
}

%inline %{
namespace vw {
namespace cartography {

  enum PixelInterpretation {
    PixelAsArea = GeoReference::PixelAsArea,
    PixelAsPoint = GeoReference::PixelAsPoint
  };

}
}
%}

namespace vw {
namespace cartography {

  class Datum {
  public:
    Datum();
    Datum( std::string const& name );
    Datum( std::string const& name,
           std::string const& spheroid_name,
           std::string const& meridian_name,
           double semi_major_axis,
           double semi_minor_axis,
           double meridian_offset );

    void set_semi_major_axis(double val);
    double semi_major_axis() const;

    %extend {
      void _geodetic_to_cartesian( double arg[3], double result[3] ) {
        (*(vw::Vector3*)result) = self->geodetic_to_cartesian( (*(vw::Vector3*)arg) );
      }
      void _cartesian_to_geodetic( double arg[3], double result[3] ) {
        (*(vw::Vector3*)result) = self->cartesian_to_geodetic( (*(vw::Vector3*)arg) );
      }
    }
    %pythoncode {
      semi_major_axis = property(semi_major_axis,set_semi_major_axis)

      def geodetic_to_cartesian(self, geodetic):
        g = doublea(3)
        c = doublea(3)
        for i in range(0,3):
          g[i] = geodetic[i]
        self._geodetic_to_cartesian(g, c)
        return [c[i] for i in range(0,3)]

      def cartesian_to_geodetic(self, cartesian):
        c = doublea(3)
        g = doublea(3)
        for i in range(0,3):
          c[i] = cartesian[i]
        self._cartesian_to_geodetic(c, g)
        return [g[i] for i in range(0,3)]
    }
  };

  class GeoReference {
  public:
    GeoReference();
    GeoReference( Datum const& datum );

    //void set_spatial_ref(void* spatial_ref_ptr);
    //void set_proj4_str( std::string const& proj4_str );
    //void set_wkt_str( std::string const& wkt_str );

    //std::string proj4_str() const;
    //std::string wkt_str() const;

    //const void*         spatial_ref_ptr() const;
    //GeoDatum datum() const;

    //std::string projection_name() const;

    //Matrix<double,3,3> transform() const { return m_transform; }

    //bool is_projected() const;

    //BBox2 bounding_box(int width, int height) const;

    void set_well_known_geogcs( std::string const& name );

    void set_sinusoidal( double center_longitude, double false_easting=0, double false_northing=0 );
    void set_mercator( double center_latitude, double center_longitude, double scale, double false_easting=0, double false_northing=0 );
    void set_orthographic( double center_latitude, double center_longitude, double false_easting=0, double false_northing=0 );
    void set_stereographic( double center_latitude, double center_longitude, double scale, double false_easting=0, double false_northing=0 );

    void set_UTM( int zone, bool north=true );

    %extend {
      void _set_transform( double arg[9] ) {
        self->set_transform( *(vw::Matrix3x3*)arg );
      }
      void _get_transform( double arg[9] ) {
        *(vw::Matrix3x3*)arg = self->transform();
      }
      void _point_to_pixel( double arg[2], double result[2] ) {
        (*(vw::Vector2*)result) = self->point_to_pixel( (*(vw::Vector2*)arg) );
      }
      void _pixel_to_point( double arg[2], double result[2] ) {
        (*(vw::Vector2*)result) = self->pixel_to_point( (*(vw::Vector2*)arg) );
      }
      void _point_to_lonlat( double arg[2], double result[2] ) {
        (*(vw::Vector2*)result) = self->point_to_lonlat( (*(vw::Vector2*)arg) );
      }
      void _lonlat_to_point( double arg[2], double result[2] ) {
        (*(vw::Vector2*)result) = self->lonlat_to_point( (*(vw::Vector2*)arg) );
      }
      void _pixel_to_lonlat( double arg[2], double result[2] ) {
        (*(vw::Vector2*)result) = self->pixel_to_lonlat( (*(vw::Vector2*)arg) );
      }
      void _lonlat_to_pixel( double arg[2], double result[2] ) {
        (*(vw::Vector2*)result) = self->lonlat_to_pixel( (*(vw::Vector2*)arg) );
      }
      void _set_pixel_interpretation( PixelInterpretation interp ) {
        self->set_pixel_interpretation( (vw::cartography::GeoReference::PixelInterpretation) interp );
      }
      PixelInterpretation _get_pixel_interpretation() {
        return (vw::cartography::PixelInterpretation) self->pixel_interpretation();
      }
    }
    %pythoncode {
      def set_transform(self,M):
        self._set_transform(array_to_doublea(M))

      def get_transform(self):
        m = doublea(9)
        self._get_transform(m)
        M = doublea_to_array(m)
        return M

      transform = property(get_transform,set_transform)

      def _wrap_converter(func):
        def converter(self,arg):
          a = doublea(2)
          r = doublea(2)
          for i in range(0,2):
            a[i] = arg[i]
          func(self,a,r)
          return (r[0],r[1])
        return converter
      point_to_pixel = _wrap_converter(_point_to_pixel)
      pixel_to_point = _wrap_converter(_pixel_to_point)
      point_to_lonlat = _wrap_converter(_point_to_lonlat)
      lonlat_to_point = _wrap_converter(_lonlat_to_point)
      pixel_to_lonlat = _wrap_converter(_pixel_to_lonlat)
      lonlat_to_pixel = _wrap_converter(_lonlat_to_pixel)

      pixel_interpretation = property(_get_pixel_interpretation, _set_pixel_interpretation)

      # This nonsense is to override SWIG to enable setting and getting properties
      _old_getattr = __getattr__
      def __getattr__(self, name):
        if name in self.__class__.__dict__ and self.__class__.__dict__[name].__class__ is property and self.__class__.__dict__[name].fget is not None:
          return self.__class__.__dict__[name].fget(self)
        else:
          return self._old_getattr(name)
      _old_setattr = __setattr__
      def __setattr__(self, name, value):
        if name in self.__class__.__dict__ and self.__class__.__dict__[name].__class__ is property and self.__class__.__dict__[name].fset is not None:
          self.__class__.__dict__[name].fset(self,value)
        else:
          self._old_setattr(name,value)
    }

    %extend {
      std::string __repr__() { std::ostringstream oss; oss << *self; return oss.str(); }
    }
  };

} // namespace cartography
} // namespace vw

%inline %{
  template <class T>
  vw::cartography::GeoReference read_georeference( T const& filename ) {
    vw::cartography::GeoReference georef;
    read_georeference( georef, filename );
    return georef;
  }

  vw::BBox2i _forward_bbox( vw::BBox2i const& bbox, vw::cartography::GeoReference const& src, vw::cartography::GeoReference const& dest ) {
    return vw::cartography::GeoTransform(src,dest).forward_bbox(bbox);
  }

  template <class PixelT, class InterpT>
  vw::ImageViewRef<PixelT> _geotransform( vw::ImageViewRef<PixelT> const& image, vw::cartography::GeoReference const& src, vw::cartography::GeoReference const& dest, InterpT const& interp, vw::int32 cols, vw::int32 rows ) {
    vw::InterpolationView<vw::ImageViewRef<PixelT>, InterpT> interpolated( image, interp );
    vw::TransformRef transform( vw::cartography::GeoTransform(src,dest) );
    vw::BBox2i bbox = transform.forward_bbox( vw::BBox2i(0,0,image.cols(),image.rows()) );
    return vw::TransformView<vw::InterpolationView<vw::ImageViewRef<PixelT>, InterpT>, vw::TransformRef>( interpolated, transform, cols?cols:bbox.max().x(), rows?rows:bbox.max().y() );
  }
%}

%template(read_georeference) read_georeference<std::string>;
%template(read_georeference) read_georeference<vw::DiskImageResource>;

%define %instantiate_geotransform(cname,ctype,pname,ptype,in,...)
  %template(_geotransform) _geotransform<ptype, in>;
%enddef

%instantiate_for_pixel_types_and_interpolations(instantiate_geotransform)

%pythoncode {
  from vwmath import BBox2i
  from image import edge_extend, ZeroEdgeExtension
  from transform import BicubicInterpolation

  def geotransform(image,src,dest,edge=None,interp=None,cols=0,rows=0):
    if image.__class__ is BBox2i:
      return _forward_bbox(image,src,dest)
    if edge is None: edge = ZeroEdgeExtension()
    if interp is None: interp = BicubicInterpolation()
    return _geotransform(edge_extend(image,edge=edge),src,dest,interp,cols,rows)
}
