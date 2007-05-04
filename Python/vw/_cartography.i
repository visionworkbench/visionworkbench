%module cartography

%{
#define SWIG_FILE_WITH_INIT
#include <vw/Image.h>
#include <vw/Cartography.h>
%}

%include "std_string.i"
%include "numpy.i"

%import "_image.i"
%import "_transform.i"

%init %{
  import_array();
%}

namespace vw {
namespace cartography {

  class GeoReference {
  public:
    GeoReference();
    
    //void set_spatial_ref(void* spatial_ref_ptr);
    //void set_transform(Matrix<double,3,3> transform) { m_transform = transform; }
    void set_proj4_str( std::string const& proj4_str );
    void set_wkt_str( std::string const& wkt_str );
    
    std::string proj4_str() const;
    std::string wkt_str() const;

    //const void*         spatial_ref_ptr() const;
    //GeoDatum datum() const;

    std::string projection_name() const;

    //Matrix<double,3,3> transform() const { return m_transform; }

    bool is_projected() const;

    //BBox2 bounding_box(int width, int height) const;    

    void set_well_known_geogcs( std::string const& name );

    void set_sinusoidal( double center_longitude, double false_easting=0, double false_northing=0 );
    void set_mercator( double center_latitude, double center_longitude, double scale, double false_easting=0, double false_northing=0 );
    void set_orthographic( double center_latitude, double center_longitude, double false_easting=0, double false_northing=0 );
    void set_stereographic( double center_latitude, double center_longitude, double scale, double false_easting=0, double false_northing=0 );

    void set_UTM( int zone, bool north=true );

    %extend {
      std::string __repr__() { std::ostringstream oss; oss << *self; return oss.str(); }
    }
  };

} // namespace cartography
} // namespace vw

%inline %{
  template <class PixelT, class InterpT>
  vw::ImageViewRef<PixelT> _geotransform( vw::ImageViewRef<PixelT> const& image, vw::cartography::GeoReference const& src, vw::cartography::GeoReference const& dest, InterpT const& interp ) {
    vw::InterpolationView<vw::ImageViewRef<PixelT>, InterpT> interpolated( image, interp );
    vw::TransformRef transform( vw::cartography::GeoTransform(src,dest) );
    return vw::TransformView<vw::InterpolationView<vw::ImageViewRef<PixelT>, InterpT>, vw::TransformRef>( interpolated, transform );
  }
%}

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
