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

namespace vw {
namespace cartography {

  class Datum {
  public:
    Datum();
    Datum( std::string const& name );

    %extend {
      void _geodetic_to_cartesian( double arg[3], double result[3] ) {
        (*(vw::Vector3*)result) = self->geodetic_to_cartesian( (*(vw::Vector3*)arg) );
      }
      void _cartesian_to_geodetic( double arg[3], double result[3] ) {
        (*(vw::Vector3*)result) = self->cartesian_to_geodetic( (*(vw::Vector3*)arg) );
      }
    }
    %pythoncode {
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
    
    //void set_spatial_ref(void* spatial_ref_ptr);
    //void set_transform(Matrix<double,3,3> transform) { m_transform = transform; }
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
    }
    %pythoncode {
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

  template <class PixelT, class InterpT>
  vw::ImageViewRef<PixelT> _geotransform( vw::ImageViewRef<PixelT> const& image, vw::cartography::GeoReference const& src, vw::cartography::GeoReference const& dest, InterpT const& interp ) {
    vw::InterpolationView<vw::ImageViewRef<PixelT>, InterpT> interpolated( image, interp );
    vw::TransformRef transform( vw::cartography::GeoTransform(src,dest) );
    return vw::TransformView<vw::InterpolationView<vw::ImageViewRef<PixelT>, InterpT>, vw::TransformRef>( interpolated, transform );
  }
%}

%template(read_georeference) read_georeference<std::string>;
%template(read_georeference) read_georeference<vw::DiskImageResource>;

%define %instantiate_geotransform(cname,ctype,pname,ptype,in,...)
  %template(_geotransform) _geotransform<ptype, in>;
%enddef
   
%instantiate_for_pixel_types_and_interpolations(instantiate_geotransform)

%pythoncode {
  from image import edge_extend, ZeroEdgeExtension
  from transform import BilinearInterpolation

  def geotransform(image,src,dest,edge=None,interp=None):
    if edge is None: edge = ZeroEdgeExtension()
    if interp is None: interp = BilinearInterpolation()
    return _geotransform(edge_extend(image,edge=edge),src,dest,interp)
}
