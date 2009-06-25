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

namespace po = boost::program_options;

using namespace vw;
using namespace vw::math;
using namespace vw::cartography;
using namespace std;

//Global variables
double pi=4*atan(1.0);
std::string input_file_name, output_prefix = "", algorithm;
double radius=1737400.0;
bool output_gradient=true;
bool output_aspect=true;
bool output_pretty=true;
bool use_horn=true; //arbitrary default

// Erases a file suffix if one exists and returns the base string
static std::string prefix_from_filename(std::string const& filename) {
  std::string result=filename;
  int index = result.rfind(".");
  if (index != -1)
    result.erase(index, result.size());
  return result;
}

template <class ImageT>
Vector3 pixel_to_cart (Vector2 pos, ImageView<ImageT> img, GeoReference GR) {

	Vector2 loc_longlat2=GR.point_to_lonlat(GR.pixel_to_point(pos));
	Vector3 loc_longlat3(loc_longlat2(0),loc_longlat2(1),img((int)pos(0),(int)pos(1))+radius);
	Vector3 loc_cartesian=GR.datum().geodetic_to_cartesian(loc_longlat3); 
	return loc_cartesian;
}

template <class ImageT>
double pixel_cartesian_dist(Vector2 pt1, Vector2 pt2, ImageView<ImageT> img, GeoReference GR) {

	Vector3 loc1=pixel_to_cart(pt1, img, GR); 
	Vector3 loc2=pixel_to_cart(pt2, img, GR);
	unsigned int i;
	double res=0;
	//cartesian distance. sqrt(dx^2+dy^2+dz^2).
	for(i=0;i<3;i++)
		res+= pow(loc1[i]-loc2[i],2);
	res=pow(res,0.5);	
	return res;
}

template <class ImageT>
Vector3 naive_horn (int x, int y, ImageView<ImageT> img, GeoReference GR) { 

	double pi=(4.0*atan(1.0));//ok, this should be somewhere else

	Vector2 Nc(x,y-1);
	Vector2 Wc(x-1,y);
	Vector2 Sc(x,y+1);
	Vector2 Ec(x+1,y);

	double N=img(x,y-1);
	double NW=img(x-1,y-1);
	double W=img(x-1,y);
	double SW=img(x-1,y+1);
	double S=img(x,y+1);
	double SE=img(x+1,y+1);
	double E=img(x+1,y);
	double NE=img(x+1,y-1);

	//calculate an east-west gradient and a north-south gradient
	//p: west to east
	double p=( (NE+2*E+SE)-(NW+2*W+SW) )/(4*pixel_cartesian_dist(Wc,Ec, img, GR) );			
	
	//q: south to north
	double q=( (NW+2*N+NE)-(SW+2*S+SE) )/(4*pixel_cartesian_dist(Nc,Sc, img, GR) );
	
	//total gradient:
	double gradient=pow(p*p+q*q, 0.5);		
	
	//aspect: q/p=tan(angle)
	double aspect=atan(q/p);

	if(p<0) aspect=aspect+pi;

	aspect=2*pi-aspect+pi/2;
	if(aspect<0) aspect=2*pi+aspect;
	if(aspect>=2*pi) aspect=aspect-2*pi*(int)(aspect/2/pi);

	if(p==0) return Vector3(0,0,0);		
			
		
	return Vector3(aspect,gradient,1);
}

template <class ImageT>
Vector3 interpolate_plane (int x, int y, ImageView<ImageT> img, GeoReference GR) { 

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
	
	//normalize
	center_normal=normalize(center_normal);
	plane_normal=normalize(plane_normal);

	if(dot_prod(plane_normal,center_normal) <0) plane_normal=plane_normal*-1;
		
	double dotprod=dot_prod(plane_normal,center_normal);

	//gradient angle is absval
	double gradient_angle=acos(dotprod);
	if(gradient_angle>pi/2)
		gradient_angle=pi-gradient_angle;
	
	//calculating aspect...
	Vector3 surface_normal_on_sphere_tangent_plane=normalize(plane_normal-dot_prod(plane_normal,center_normal)*center_normal);
	//get projection of (0,0,1) onto sphere tangent plane. 
	Vector3 north(0,0,1);
	Vector3 north_projected=normalize(north-dot_prod(north,center_normal)*center_normal);
	//find the angle between those two
	double dotprod2=dot_prod(surface_normal_on_sphere_tangent_plane, north_projected);

	//figure out which angle it is. 
	double aspect=acos(dotprod2);

	if( dot_prod((cross_prod(north_projected,surface_normal_on_sphere_tangent_plane)),center_normal) < 0)
		aspect=pi+aspect;
	else
		aspect=pi-aspect;

	return Vector3(aspect,abs(gradient_angle),1);
}

template <class imageT>
void do_slopemap (po::variables_map const& vm) { //not sure what the arguments are

	GeoReference GR;
	read_georeference( GR, input_file_name );

	ImageView<imageT> img;
	read_image( img, input_file_name );

	int x;
	int y;

	ImageView<double> gradient_angle;
	ImageView<double> aspect;
	ImageView<PixelHSV<double> > pretty;

	if(output_gradient) gradient_angle.set_size(img.cols(),img.rows());
	if(output_aspect) aspect.set_size(img.cols(),img.rows());
	if(output_pretty) pretty.set_size(img.cols(),img.rows());
	
	for(x=1;x<img.cols()-1;x++) {

		for(y=1;y<img.rows()-1;y++) {	
			//does the calculations for both gradient and aspect but outputs one or the other or both
			Vector3 res;
			if(use_horn) res=naive_horn(x,y,img,GR);
			else res=interpolate_plane(x,y,img,GR);

			if(output_aspect)   aspect(x,y) = res(0);
			if(output_gradient) gradient_angle(x,y) = res(1); 
			if(output_pretty)   pretty(x,y) = PixelHSV<double>(res(0),res(1)+0.1*abs(pi-res(0))/pi,res(1)+0.2*abs(pi-res(0))); 
			//pretty things arbitrary
		}	
	}	

	ImageView<PixelRGB<uint8> > pretty2;

	if(output_pretty) {
		select_channel(pretty,0)=normalize(select_channel(pretty,0));
		select_channel(pretty,1)=normalize(select_channel(pretty,1),0,0.25,0.3,1);
		select_channel(pretty,2)=normalize(select_channel(pretty,2),0.5,1);
	
		pretty2=pixel_cast_rescale<PixelRGB<uint8> >( copy(pretty) );
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
    ("radius", po::value<double>(&radius), "Set radius in meters as specified. [default: 1737400 meters (moon radius)]") 
    ("no-aspect", "Do not output aspect")
    ("no-gradient", "Do not output gradient")
    ("algorithm", po::value<std::string>(&algorithm)->default_value("horn"), "Choose an algorithm to calculate slope/aspect from [horn, planefit]") //could use better names
    ("no-pretty", "Do not output colored image."); 
  
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
  
  if( output_prefix == "" ) { output_prefix=prefix_from_filename(input_file_name);
  }
  
  if( vm.count("verbose") ) {
    set_debug_level(VerboseDebugMessage);
  }


//checking strings
boost::to_lower(algorithm);
if( !(  algorithm == "horn" ||
	algorithm == "planefit" ||
	algorithm == "" ) ) { //it's okay if it isn't set?
	vw_out(0) << "Unknown algorithm: " << algorithm << ". Options are : [ horn, planefit ]\n";
	exit(0);
}
else {
if(algorithm=="horn")
	use_horn=true;
else
	use_horn=false;
}

if(vm.count("no-aspect")) output_aspect=false;
if(vm.count("no-gradient")) output_gradient=false;
if(!output_aspect && !output_gradient && !output_pretty) {
	vw_out(0) << "No output specified. Select at least one of [gradient, output, pretty].\n" << endl;
}

if(vm.count("no-pretty")) output_pretty=false;


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
