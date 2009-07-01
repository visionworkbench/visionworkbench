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

/*

Implements a modified Horn's algorithm for calculating slope, where df/dx and df/dy on the lon/lat map correspond to df/dphi and df/dtheta. 

*/

namespace po = boost::program_options;

using namespace vw;
using namespace vw::math;
using namespace vw::cartography;
using namespace std;

//Global variables
double pi=4*atan(1.0);
std::string input_file_name, output_prefix = "";

double total_lat=0;
double total_lon=0;

bool output_gradient=true;
bool output_aspect=true;
bool output_pretty=true;

bool use_horn=true; //which is default? is this the best way?


// Erases a file suffix if one exists and returns the base string
static std::string prefix_from_filename(std::string const& filename) {
  std::string result=filename;
  int index = result.rfind(".");
  if (index != -1)
    result.erase(index, result.size());
  return result;
}

Vector3 pixel_to_cart (Vector2 pos, double alt, GeoReference GR) {

	Vector2 loc_longlat2=GR.point_to_lonlat(GR.pixel_to_point(pos));
	Vector3 loc_longlat3(loc_longlat2(0),loc_longlat2(1),alt);
	Vector3 loc_cartesian=GR.datum().geodetic_to_cartesian(loc_longlat3); 
	return loc_cartesian;
}

template <class ImageT>
Vector3 pixel_to_cart (Vector2 pos, ImageView<ImageT> img, GeoReference GR) {
	return pixel_to_cart(pos,img((int)pos[0],(int)pos[1]),GR);
}

double pixel_cart_dist(Vector3 pos1, Vector3 pos2) { //was used with original horn's
	Vector3 diff=pos1-pos2;
	return pow(diff[0]*diff[0]+diff[1]*diff[1]+diff[2]*diff[2],0.5);
}

template <class ImageT>
Vector3 modified_horn (int x, int y, ImageView<ImageT> img, GeoReference GR) { 

	Vector3 center=pixel_to_cart(Vector2(x,y),img,GR);
	Vector2 Nc(x,y-1);
	Vector2 NWc(x-1,y-1);
	Vector2 Wc(x-1,y);
	Vector2 SWc(x-1,y+1);
	Vector2 Sc(x,y+1);
	Vector2 SEc(x+1,y+1);
	Vector2 Ec(x+1,y);
	Vector2 NEc(x+1,y-1);

	double N=img(x,y-1);
	double NW=img(x-1,y-1);
	double W=img(x-1,y);
	double SW=img(x-1,y+1);
	double S=img(x,y+1);
	double SE=img(x+1,y+1);
	double E=img(x+1,y);
	double NE=img(x+1,y-1);

	//p: west to east
	double p=( (NE+2*E+SE)-(NW+2*W+SW) )/(4);			
	
	//q: south to north
	double q=( (NW+2*N+NE)-(SW+2*S+SE) )/(4);

	Vector2 lonlat=GR.pixel_to_lonlat(Vector2(x,y));

	double lat=lonlat(1);
	double lon=lonlat(0);

	double phi=lon/180.0*pi;
	double theta=(-lat+90.0)/180.0*pi;

	double rho=pow(center[0]*center[0]+center[1]*center[1]+center[2]*center[2],0.5);
	double df_dtheta=-q*img.cols()*180.0/2.0/total_lat/pi;
	double df_dphi=p*img.rows()*180.0/2.0/total_lon/pi;

	double r_comp=1;	
	double theta_comp=df_dtheta/rho;
	double phi_comp=1/(rho*sin(theta))*df_dphi;

	double n_x=r_comp*sin(theta)*cos(phi)-phi_comp*sin(phi)+theta_comp*cos(theta)*cos(phi);
	double n_y=r_comp*sin(theta)*sin(phi)+phi_comp*cos(phi)+theta_comp*cos(theta)*sin(phi);
	double n_z=r_comp*cos(theta)-theta_comp*sin(theta);

	Vector3 center_normal=normalize(center);
	Vector3 plane_normal=normalize(Vector3(n_x,n_y,n_z));

	if(dot_prod(plane_normal,center_normal) <0) plane_normal=plane_normal*-1;

        double dotprod=dot_prod(plane_normal,center_normal);
	double gradient_angle=pi/2-acos(dotprod);

	Vector3 surface_normal_on_sphere_tangent_plane=normalize(plane_normal-dot_prod(plane_normal,center_normal)*center_normal);
	Vector3 north(0,0,1);
	Vector3 north_projected=normalize(north-dot_prod(north,center_normal)*center_normal);
	double dotprod2=dot_prod(surface_normal_on_sphere_tangent_plane, north_projected);

	double aspect=acos(dotprod2);

	if(dotprod2>1)	{aspect=0;}
	if(dotprod2<-1) {aspect=0;}

        if( dot_prod((cross_prod(north_projected,surface_normal_on_sphere_tangent_plane)),center_normal) < 0)
	        aspect=pi+aspect;
	else
        	aspect=pi-aspect;
	if(aspect>=2*pi)
		aspect=aspect-2*pi;
        return Vector3(aspect,gradient_angle,1);
}

template <class imageT>
void do_slopemap (po::variables_map const& vm) { //not sure what the arguments are

	GeoReference GR;
	read_georeference( GR, input_file_name );

	ImageView<imageT> img;
	read_image( img, input_file_name );

	int x;
	int y;

	//setting up global variables for modified_horn
	Vector2 lonlatdiff=GR.pixel_to_lonlat(Vector2(img.cols(),img.rows()))-GR.pixel_to_lonlat(Vector2(0,0));
	total_lat=abs(lonlatdiff[1]);
	total_lon=abs(lonlatdiff[0]);

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
			Vector3 res=modified_horn(x,y,img,GR);
			if(output_aspect)   aspect(x,y) = res(0);
			if(output_gradient) gradient_angle(x,y) = res(1);
			if(output_pretty)   pretty(x,y) = PixelHSV<double>(res(0),res(1),res(1)+0.1*abs(pi-res(0))); 
		 } 
	}	
	ImageView<PixelRGB<uint8> > pretty2;

	if(output_pretty) {
		select_channel(pretty,0)=normalize(select_channel(pretty,0),0,2*pi,0,1);
		select_channel(pretty,1)=normalize(select_channel(pretty,1),0,pi/2,0,1);
		select_channel(pretty,2)=normalize(select_channel(pretty,2),0.5,1);
	
		pretty2=pixel_cast_rescale<PixelRGB<uint8> >( copy(pretty) );
	}
	//save everything to file
	if(output_gradient) write_georeferenced_image( output_prefix + "_gradient.tif" , gradient_angle, GR);
	if(output_aspect)  write_georeferenced_image( output_prefix + "_aspect.tif"   , aspect, GR);
	if(output_pretty)  write_image( output_prefix + "_pretty.tif"   , pretty2); 	
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
