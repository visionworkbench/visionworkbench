#include <iostream>
#include <vw/Image.h>
#include <vw/FileIO.h>

#include <vw/Cartography/GeoReference.h>
#include <vw/Cartography/GeoTransform.h>
#include <vw/Cartography/FileIO.h>
#include <vw/Cartography/Datum.h>
#include <vw/Image/ImageMath.h>
#include <vw/Image/Algorithms.h>
#include <vw/Math/Matrix.h>
#include <vw/Math/Vector.h>
#include <vw/Math/LinearAlgebra.h>
#include <boost/program_options.hpp>
#include <vw/FileIO/DiskImageView.h>

/*
Implements modified versions of finite difference algorithms + fitting a plane to 9 points of a 3x3 window
*/

namespace po = boost::program_options;

using namespace vw;
using namespace vw::math;
using namespace vw::cartography;
using namespace std;

//Global variables
double pi=4*atan(1.0);
std::string input_file_name, output_prefix, algorithm_string = "";

bool output_gradient=true;
bool output_aspect=true;
bool output_pretty=false; //probably more for debugging purposes than for anything else

enum Algorithm { HORN, SA, FH, PLANEFIT };
Algorithm algorithm;

bool spherically_defined=false;

// Erases a file suffix if one exists and returns the base string
static std::string prefix_from_filename(std::string const& filename) {
  std::string result=filename;
  int index = result.rfind(".");
  if (index != -1)
    result.erase(index, result.size());
  return result;
}

//basic utilities

Vector3 pixel_to_cart (Vector2 pos, double alt, GeoReference GR) {
  Vector2 loc_longlat2=GR.point_to_lonlat(GR.pixel_to_point(pos));
  Vector3 loc_longlat3(loc_longlat2(0),loc_longlat2(1),alt);
  Vector3 loc_cartesian=GR.datum().geodetic_to_cartesian(loc_longlat3); 
  return loc_cartesian;
}

template <class ImageT>
Vector3 pixel_to_cart (Vector2 pos, DiskImageView<ImageT> img, GeoReference GR) {
  return pixel_to_cart(pos,img((int)pos[0],(int)pos[1]),GR);
}

double dist_from_2pi(double n) {
  if( abs(n+2*pi) < abs(n) ) return n+2*pi;
  if( abs(n-2*pi) < abs(n) ) return n-2*pi;
  return n;
}

//utilities...

Vector2 gradient_aspect_from_normals(Vector3 center_normal, Vector3 plane_normal) {
  //gradient angle is absval
  double dotprod=dot_prod(plane_normal,center_normal);
  double gradient_angle=acos(dotprod);
  Vector3 surface_normal_on_sphere_tangent_plane=normalize(plane_normal-dot_prod(plane_normal,center_normal)*center_normal);

  //get projection of (0,0,1) onto sphere tangent plane. 
  Vector3 north(0,0,1);
  Vector3 north_projected=normalize(north-dot_prod(north,center_normal)*center_normal);

  //find the angle between those two
  double dotprod2=dot_prod(surface_normal_on_sphere_tangent_plane, north_projected);

  //figure out which angle it is. 
  double aspect=acos(dotprod2);
  if(dotprod2>1) { dotprod2=1; aspect=0;}
  if(dotprod2<-1) { dotprod2=-1; aspect=0;}

  if( dot_prod((cross_prod(north_projected,surface_normal_on_sphere_tangent_plane)),center_normal) < 0)
    aspect=pi+aspect;
  else
    aspect=pi-aspect;
  
  if(aspect>=2*pi)
    aspect=aspect-2*pi;
  return Vector2(aspect,gradient_angle);
}

Vector2 gradient_aspect_from_dx_dy(double& dx, double& dy) {
  //total gradient:
  double gradient=norm_2(Vector2(dx,dy));    
  double gradient_angle=atan(gradient);
  //aspect: dy/dx=tan(angle)
  double aspect=atan(dy/dx);
  if(dx<0) aspect=aspect+pi;
  aspect=2*pi-aspect+pi/2;
  if(aspect<0) aspect=2*pi+aspect;
  if(aspect>=2*pi) aspect=aspect-2*pi*(int)(aspect/2/pi);  
  return Vector2(aspect,gradient_angle);
}

Vector2 gradient_aspect_from_dtheta_dphi(double rho, double theta, double phi, double& dtheta, double& dphi, Vector3 center) {
  double r_comp=1;  
  double theta_comp=dtheta/rho;
  double phi_comp=1/(rho*sin(theta))*dphi;

  double n_x=r_comp*sin(theta)*cos(phi)-phi_comp*sin(phi)+theta_comp*cos(theta)*cos(phi);
  double n_y=r_comp*sin(theta)*sin(phi)+phi_comp*cos(phi)+theta_comp*cos(theta)*sin(phi);
  double n_z=r_comp*cos(theta)-theta_comp*sin(theta);

  Vector3 plane_normal=normalize(Vector3(n_x,n_y,n_z));
  center=normalize(center);
  return gradient_aspect_from_normals(center, plane_normal);
}

template <class ImageT>
Vector2 uneven_grid (int x, int y, DiskImageView<ImageT> img, GeoReference GR) {

  Vector3 center=pixel_to_cart(Vector2(x,y),img,GR);
  Vector3 center_normal=normalize(center);
  Vector3 center_below=pixel_to_cart(Vector2(x,y),0,GR);

  Vector3 north=Vector3(0,0,1);
        Vector3 north_projected=normalize(north-dot_prod(north,center_normal)*center_normal);
  
  //define temporary axis
  Vector3 up_normal=normalize(north_projected);        //also theta hat
  Vector3 third=normalize(cross_prod(up_normal,center_normal));     //also phi hat

  //spherical
  Vector2 lonlat=GR.pixel_to_lonlat(Vector2(x,y));
  double lat=lonlat(1);
  double lon=lonlat(0);

  //rise= run*slope
  Matrix<double> rises;  //rise over vector distance
  Matrix<double> runs;  //components in either direction
  
  //different ways of weighting neighbors. 
  if(algorithm==HORN) { rises=Matrix<double>(12,1); runs=Matrix<double>(12,2); }
  if(algorithm==SA)   { rises=Matrix<double>(8,1); runs=Matrix<double>(8,2); }
  if(algorithm==FH)   { rises=Matrix<double>(4,1); runs=Matrix<double>(4,2); }
    
  int ct=0;
  for(int i=-1;i<=1;i++) {
    for(int j=-1;j<=1;j++) {
      if(i==0 && j==0) continue;
      int repeat=1;
      if(algorithm==HORN) {
        if(i==0 && j!=0 || i!=0 && j==0) //for horn, weight direct neighbors twice
          repeat=2;
      }
      if(algorithm==FH) {
        if(i!=0 && j!=0) //ignore diagonals
          continue;
      }
      
      for(int k=0;k<repeat;k++) {
        if(!spherically_defined) {
          Vector3 neighbor=pixel_to_cart(Vector2(x+i,y+j),img,GR);
          Vector3 neighbor_below_rescale=normalize(neighbor)*(norm_2(center_below)/dot_prod(normalize(neighbor),center_normal));  
          Vector2 lonlat=GR.pixel_to_lonlat(Vector2(x,y));
          Vector3 v=neighbor_below_rescale-center_below;
          //or, project both onto tangent plane...
          rises(ct,0)=(img(x+i,y+j)-img(x,y))/norm_2(v);
          v=normalize(v);
          runs(ct,0)=dot_prod(v,up_normal);
          runs(ct,1)=dot_prod(v,third);        
        } else {
          rises(ct,0)=(img(x+i,y+j)-img(x,y));    
          Vector2 neighbor_lonlat=GR.pixel_to_lonlat(Vector2(x+i,y+j));
          runs(ct,0)=(neighbor_lonlat[1]-lonlat[1])*pi/180.0;
          runs(ct,1)=dist_from_2pi(neighbor_lonlat[0]-lonlat[0])*pi/180.0;      
        }  
        ct++;          
      }
    }
  }
  //solve using least squares
  //gradient components=(VtV)-1Vt * rises
  Matrix<double> VtV(2,2);
  VtV=transpose(runs)*runs;
  Matrix<double> VtVinv(2,2);
  VtVinv(0,0)=VtV(1,1);
  VtVinv(1,1)=VtV(0,0);
  VtVinv(0,1)=-VtV(0,1);
  VtVinv(1,0)=-VtV(1,0);
  VtVinv=VtVinv/(VtV(0,0)*VtV(1,1)-VtV(0,1)*VtV(1,0));

  Matrix<double> ans=VtVinv*transpose(runs)*rises;

  if(spherically_defined) {
    double rho=norm_2(center);
    double phi=lon/180.0*pi;
    double theta=(-lat+90.0)/180.0*pi;
    double dtheta=ans(0,0);
    double dphi=-ans(1,0);
    return gradient_aspect_from_dtheta_dphi(rho, theta, phi, dtheta, dphi, center);
  } 
  //otherwise,
  return gradient_aspect_from_dx_dy(ans(0,1), ans(0,0));
}

template <class ImageT>
Vector2 interpolate_plane (int x, int y, DiskImageView<ImageT> img, GeoReference GR) { 
  Matrix<double> A(9,4);
  int i=0;
  int j=0;  
  int ct=0;
  Vector3 center_normal=pixel_to_cart(Vector2(x,y),img,GR);
  
  for(i=-1;i<=1;i++) {
    for(j=-1;j<=1;j++) {
      Vector3 tmp=pixel_to_cart(Vector2(x+i,y+j),img,GR);
      A(ct,0)=tmp(0);
      A(ct,1)=tmp(1);
      A(ct,2)=tmp(2);
      A(ct,3)=1;
      ct++;
    }
  }

  Matrix<double> U;
  Matrix<double> VT;
  Vector<double> s;
  
  svd(A, U, s, VT);
  Vector<double> plane_normal(3);
  plane_normal(0)=VT(3,0);
  plane_normal(1)=VT(3,1);
  plane_normal(2)=VT(3,2);
  //normalize sphere normal
  center_normal=normalize(center_normal);
  plane_normal=normalize(plane_normal);
  if(dot_prod(plane_normal,center_normal) <0) plane_normal=plane_normal*-1;
  return gradient_aspect_from_normals(center_normal, plane_normal);  
}

template <class imageT>
void do_slopemap (po::variables_map const& vm) { //not sure what the arguments are

  GeoReference GR;
  read_georeference( GR, input_file_name );

  DiskImageView<imageT> img(input_file_name);

  int x;
  int y;

  ImageView<double> gradient_angle;
  ImageView<double> aspect;
  ImageView<PixelHSV<double> > pretty;

  if(output_gradient) 
    gradient_angle.set_size(img.cols(),img.rows());
  if(output_aspect) 
    aspect.set_size(img.cols(),img.rows());
  if(output_pretty) pretty.set_size(img.cols(),img.rows());

  for(x=1;x<img.cols()-1;x++) {
    for(y=1;y<img.rows()-1;y++) {  
      Vector2 res;
      //these are pretty similar...
      if(algorithm==PLANEFIT) res=interpolate_plane(x,y,img,GR);
      else res=uneven_grid(x,y,img,GR);

      if(output_aspect)   aspect(x,y) = res(0);
      if(output_gradient) gradient_angle(x,y) = res(1);
      if(output_pretty)   pretty(x,y) = PixelHSV<double>(res(0),res(1),(res(1))+0.2*abs(pi-res(0)));//(res(1)/pi*2)*abs(pi-res(0))); 
     } 
  }  
  ImageView<PixelRGB<uint8> > pretty2;

  if(output_pretty) {
    select_channel(pretty,0)=normalize(select_channel(pretty,0),0,2*pi,0,1);
    select_channel(pretty,1)=normalize(select_channel(pretty,1),0,pi/2,0.1,1);
    select_channel(pretty,2)=normalize(select_channel(pretty,2),0.3,0.6);
  
    pretty2=pixel_cast_rescale<PixelRGB<uint8> >( copy(pretty) );
    pretty2=PixelRGB<uint8>(255,255,255)-pretty2;
  }
  //save everything to file
  if(output_gradient) write_georeferenced_image( output_prefix + "_gradient.tif" , gradient_angle, GR);
  if(output_aspect)   write_georeferenced_image( output_prefix + "_aspect.tif"   , aspect, GR);
  if(output_pretty)   write_image( output_prefix + "_pretty.tif"   , pretty2);   
}


int main( int argc, char *argv[] ) {

  set_debug_level(InfoMessage);

  po::options_description desc("Description: Outputs gradient and/or aspect at each point of an input DEM with altitude values\n\nUsage: slopemap [options] <input file> \n\nOptions");
  desc.add_options()
    ("help", "Display this help messsage")
    ("input-file", po::value<std::string>(&input_file_name), "Explicitly specify the input file")
    ("output-prefix,o", po::value<std::string>(&output_prefix), "Specify the output prefix") //should add more description...
    ("no-aspect", "Do not output aspect")
    ("no-gradient", "Do not output gradient")
    ("pretty", "Output colored image.") 
    ("algorithm", po::value<std::string>(&algorithm_string)->default_value("horn"), "Choose an algorithm to calculate slope/aspect from [ horn, fh, sa, planefit ]. Horn: Horn's algorithm; FH: Fleming & Hoffer's (rook's case); SA: Sharpnack & Akin's (queen's case)")
    ("spherical", po::value<bool>(&spherically_defined)->default_value(true), "Spherical/elliptical datum (recommended); otherwise, a flat grid");

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
    std::cout << "Error: Must specify exactly one input file!\n" << std::endl;
    std::cout << desc << std::endl;
    return 1;
  }

  if( output_prefix == "" ) { output_prefix=prefix_from_filename(input_file_name); }

  if( vm.count("verbose") ) {
    set_debug_level(VerboseDebugMessage);  
  }

  //checking strings

  boost::to_lower(algorithm_string);
  if( !(  algorithm_string == "horn" ||
    algorithm_string == "fh" ||
    algorithm_string == "sa" ||
    algorithm_string == "planefit" ||

    algorithm_string == "" ) ) { //it's okay if it isn't set?
    vw_out(0) << "Unknown algorithm: " << algorithm_string << ". Options are : [ horn, fh, sa, planefit ]\n";
    exit(0); 
  }
  else {
    if(algorithm_string=="horn")
      algorithm=HORN;
    else if(algorithm_string=="fh")
      algorithm=FH;
    else if(algorithm_string=="sa")
      algorithm=SA;
    else if(algorithm_string=="planefit")
      algorithm=PLANEFIT;
  }

  if(vm.count("no-aspect"))   output_aspect=false;
  if(vm.count("no-gradient"))   output_gradient=false;
  if(vm.count("pretty")) output_pretty=true;

  if(!output_aspect && !output_gradient && !output_pretty) {
    vw_out(0) << "No output specified. Select at least one of [ gradient, output, pretty ].\n" << endl;
  }

  try {
    // Get the right pixel/channel type.
    DiskImageResource *rsrc = DiskImageResource::open(input_file_name);
    ChannelTypeEnum channel_type = rsrc->channel_type();
    PixelFormatEnum pixel_format = rsrc->pixel_format();
    delete rsrc;
    
    switch(pixel_format) {
    case VW_PIXEL_GRAY:
    case VW_PIXEL_GRAYA:
    case VW_PIXEL_RGB:
    case VW_PIXEL_RGBA:
      switch(channel_type) {
      case VW_CHANNEL_UINT8:  do_slopemap<PixelGray<uint8>   >(vm); break;
      case VW_CHANNEL_INT16:  do_slopemap<PixelGray<int16>   >(vm); break;
      case VW_CHANNEL_UINT16: do_slopemap<PixelGray<uint16>  >(vm); break;
      case VW_CHANNEL_FLOAT32:do_slopemap<PixelGray<float32> >(vm); break;
      case VW_CHANNEL_FLOAT64:do_slopemap<PixelGray<float64> >(vm); break;
      default:                do_slopemap<PixelGray<float32> >(vm); break;
      }
      break;
    default:
      std::cout << "Error: Unsupported pixel format.\n";
      exit(0);
    }
  } catch( Exception& e ) {
    std::cout << "Error: " << e.what() << std::endl;
  }
  return 0;

}
