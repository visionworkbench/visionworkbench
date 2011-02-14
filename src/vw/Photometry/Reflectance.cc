// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <string>
#include <fstream>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>

#include <vw/Core.h>
#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Cartography.h>
#include <vw/Math.h>
using namespace vw;
using namespace vw::cartography;

#include <vw/Photometry/Reflectance.h>
#include <vw/Photometry/ShapeFromShading.h>
#include <vw/Photometry/Shape.h>
#include <vw/Photometry/Albedo.h>
#include <vw/Photometry/Exposure.h>
#include <vw/Photometry/Reconstruct.h>
using namespace vw::photometry;

//generates the normal of a point p1 from the 3D coordinates of p1, p2, p3
//this function is currently buggy. Needs to be fixed and this will replace
//computeNormalFrom3DPoints
Vector3
vw::photometry::computeNormalFrom3DPointsGeneral(Vector3 p1,
                                                 Vector3 p2,
                                                 Vector3 p3) {
  return -normalize(cross_prod(p2-p1,p3-p1));
}

//generates the normal of a point p1 from the 3D coordinates of p1, p2, p3
//assumes unit dstance between adjacent pixels in x, y coordinates
//this function will be grandfathered
//author Ara Nefian
Vector3
vw::photometry::computeNormalFrom3DPoints(Vector3 p1, Vector3 p2,
                                          Vector3 p3) {
  Vector3 normal;

  normal[0] = p1[2]-p2[2];
  normal[1] = p1[2]-p3[2];
  normal[2] = -1;

  return normalize(normal);
}

/*
//generates the 3D coordinates of a point from longitude and latitude on the Moon
//Vector2 lon_lat is a 2D vector. First element is the longitude and the second the latitude
//author Ara Nefian
Vector3 get_normal(Vector2 lon_lat){

  Vector3 xyz;
  float x, y, z;
  float rho;
  float longitude, latitude;

  longitude = lon_lat[0]*3.14/180;
  latitude = lon_lat[1]*3.14/180;
  rho = 1738; //kilometers
  x = rho*cos(longitude)*cos(latitude);
  y = rho*sin(longitude)*cos(latitude);
  z = rho*sin(latitude);

  xyz[0] = x;
  xyz[1] = y;
  xyz[2] = z;

  return xyz;
}
*/

//reads the 3D position of the Sun in Moon centric cordinates
//char *filename: the filename with the sun positions (sunposition.txt)
//int numEntres: number of entries in the sun position file
std::vector<Vector3>
vw::photometry::ReadSunPosition( std::string const& filename,
                                 int const& numEntries ) {
  std::vector<Vector3> sunPositions(numEntries);

  std::ifstream infile( filename.c_str() );
  if ( infile.is_open() ) {
    for ( int i = 0; i < numEntries; i++ ) {
      if ( !infile.good() )
        vw_throw( IOErr() << "Sun Position file does not have enough entries for the "
                  << numEntries << " requested." );
      std::string imgname;
      infile >> imgname >> sunPositions[i][0]
             >> sunPositions[i][1] >> sunPositions[i][2];
    }
  } else {
    vw_throw( IOErr() << "Unable to open Sun Position file: " << filename );
  }
  infile.close();

  return sunPositions;
}

//reads the 3D position of the Spacecraft in Moon centric cordinates
//char *filename: the filename with the sun positions (sunposition.txt)
//int numEntres: number of entries in the sun position file
std::vector<Vector3>
vw::photometry::ReadSpacecraftPosition(std::string const& filename,
                                       int const& numEntries) {
  std::vector<Vector3> spacecraftPositions(numEntries);

  std::ifstream infile( filename.c_str() );
  if ( infile.is_open() ) {
    for ( int i = 0; i < numEntries; i++ ) {
      if ( !infile.good() )
        vw_throw( IOErr() << "Spacecraft Position file does not have enough entries for the "
                  << numEntries << " requested." );
      std::string imgname;
      infile >> imgname >> spacecraftPositions[i][0]
             >> spacecraftPositions[i][1] >> spacecraftPositions[i][2];
    }
  } else {
    vw_throw( IOErr() << "Unable to open Spacecraft Position file: " << filename );
  }
  infile.close();

  return spacecraftPositions;
}

/*
//computes the cosine of the light direction and the normal to the Moon
//Vector3  sunpos: the 3D coordinates of the Sun relative to the center of the Moon
//Vector2 lon_lat is a 2D vector. First element is the longitude and the second the latitude
//author Ara Nefian
float computeReflectance(Vector3 sunPos, Vector2 lon_lat)
{
  float reflectance;

  Vector3 xyz = get_normal(lon_lat);

  Vector3 sunDirection; // sun coordinates relative to the xyz point on the Moon surface
  sunDirection[0] = sunPos[0]-xyz[0];
  sunDirection[1] = sunPos[1]-xyz[1];
  sunDirection[2] = sunPos[2]-xyz[2];

  float sunDirectionNorm;
  sunDirectionNorm = sqrt(sunDirection[0]*sunDirection[0] + sunDirection[1]*sunDirection[1] + sunDirection[2]*sunDirection[2]);

    //printf("sun_direction norm: %f, 0: %f, 1: %f, 2: %f\n", sun_direction_norm, sun_direction[0], sun_direction[1], sun_direction[2]);
    sunDirection[0]=sunDirection[0]/sunDirectionNorm;
    sunDirection[1]=sunDirection[1]/sunDirectionNorm;
    sunDirection[2]=sunDirection[2]/sunDirectionNorm;
    //printf("sun_direction norm: %f, 0: %f, 1: %f, 2: %f\n", sun_direction_norm, sun_direction[0], sun_direction[1], sun_direction[2]);

    float xyzNorm;
    xyzNorm = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[2]*xyz[2]);

    xyz[0] = xyz[0]/xyzNorm;
    xyz[1] = xyz[1]/xyzNorm;
    xyz[2] = xyz[2]/xyzNorm;
    //printf("xyz 0: %f, 1: %f, 2: %f\n", xyz[0], xyz[1], xyz[2]);

    reflectance = sunDirection[0]*xyz[0] + sunDirection[1]*xyz[1] + sunDirection[2]*xyz[2];

    return reflectance;
}
*/

//computes the Lambertian reflectance model (cosine of the light direction and the normal to the Moon)
//Vector3  sunpos: the 3D coordinates of the Sun relative to the center of the Moon
//Vector2 lon_lat is a 2D vector. First element is the longitude and the second the latitude
//author Ara Nefian
float
vw::photometry::computeLambertianReflectanceFromNormal(Vector3 sunPos,
                                                       Vector3 xyz,
                                                       Vector3 normal) {
  float reflectance;

  //Vector3 xyz = get_normal(lon_lat);
  //printf("xyz[0] = %f, xyz[1] = %f, xyz[2] = %f\n", xyz[0], xyz[1], xyz[2]);
        // sun coordinates relative to the xyz point on the Moon surface
        Vector3 sunDirection = normalize(sunPos-xyz);

  reflectance = sunDirection[0]*normal[0] + sunDirection[1]*normal[1] + sunDirection[2]*normal[2];

  return reflectance;
}

//computes a simple naive reflectance model
//this function will be grandfathered and replaced by computeLambertianReflectanceFromNormal
//Vector3  sunpos: the 3D coordinates of the Sun relative to the center of the Moon
//Vector2 lon_lat is a 2D vector. First element is the longitude and the second the latitude
//author Ara Nefian
float
vw::photometry::computeReflectanceFromNormal(Vector3 /*sunPos*/,
                                             Vector3 /*xyz*/,
                                             Vector3 /*normal*/) {
  float reflectance;

  reflectance = -1;
  return reflectance;
}

//computes the Lunar-Lambertian reflectance mode
//Vector3  sunpos: the 3D coordinates of the Sun relative to the center of the Moon
//Vector2 lon_lat is a 2D vector. First element is the longitude and the second the latitude
//author Ara Nefian
float
vw::photometry::computeLunarLambertianReflectanceFromNormal(Vector3 sunPos, Vector3 viewPos, Vector3 xyz,  Vector3 normal, float B_0, float L) {
  float reflectance;

  //compute /mu_0 = cosine of the angle between the light direction and the surface normal.
  // sun coordinates relative to the xyz point on the Moon surface
  Vector3 sunDirection = normalize(sunPos-xyz);
  float mu_0 = dot_prod(sunDirection,normal);

  //compute  /mu = cosine of the angle between the viewer direction and the surface normal.
  // viewer coordinates relative to the xyz point on the Moon surface
  Vector3 viewDirection = normalize(viewPos-xyz);
  float mu = dot_prod(viewDirection,normal);

  reflectance = B_0*((2*L*mu_0/(mu_0+mu)) + (1-L)*mu_0 );
  return reflectance;
}

float
vw::photometry::computeLunarLambertianReflectanceFromNormal(Vector3 sunPos, Vector3 viewPos, Vector3 xyz, Vector3 normal) {
  float reflectance;
  float L;

  //compute /mu_0 = cosine of the angle between the light direction and the surface normal.
  //sun coordinates relative to the xyz point on the Moon surface
  //Vector3 sunDirection = -normalize(sunPos-xyz);
  Vector3 sunDirection = normalize(sunPos-xyz);
  float mu_0 = dot_prod(sunDirection,normal);

  //compute  /mu = cosine of the angle between the viewer direction and the surface normal.
  //viewer coordinates relative to the xyz point on the Moon surface
  Vector3 viewDirection = normalize(viewPos-xyz);
  float mu = dot_prod(viewDirection,normal);

  //compute the phase angle /alpha between the viewing direction and the light source direction
  float rad_alpha, deg_alpha;
  float cos_alpha;

  cos_alpha = dot_prod(sunDirection,viewDirection);
  if ((cos_alpha > 1)||(cos_alpha< -1)){
    printf("cos_alpha error\n");
  }

  rad_alpha = acos(cos_alpha);
  deg_alpha = rad_alpha*180/3.141592;

  //printf("deg_alpha = %f\n", deg_alpha);

  //Bob Gaskell's model
  //L = exp(-deg_alpha/60.0);

  if (deg_alpha > 90){
    //printf("Error!!: rad_alpha = %f, deg_alpha = %f\n", rad_alpha, deg_alpha);
    return(0.0);
  }
  if (deg_alpha < -90){
    //printf("Error!!: rad_alpha = %f, deg_alpha = %f\n", rad_alpha, deg_alpha);
    return(0.0);
  }

  //Alfred McEwen's model
  float A = -0.019;
  float B =  0.000242;//0.242*1e-3;
  float C = -0.00000146;//-1.46*1e-6;

  L = 1.0 + A*deg_alpha + B*deg_alpha*deg_alpha + C*deg_alpha*deg_alpha*deg_alpha;

  //        std::cout << " sun direction " << sunDirection << " view direction " << viewDirection << " normal " << normal;
  //        std::cout << " cos_alpha " << cos_alpha << " incident " << mu_0 << " emission " << mu;
  //printf(" deg_alpha = %f, L = %f\n", deg_alpha, L);

  //if (mu_0 < 0.15){ //incidence angle is close to 90 deg
  if (mu_0 < 0.0){
    //mu_0 = 0.15;
    return (0.0);
  }

  if (mu < 0.0){ //emission angle is > 90
    mu = 0.0;
    //return (0.0);
  }

  if (mu_0 + mu == 0){
    //printf("negative reflectance\n");
    reflectance = 0.0;
  }
  else{
    reflectance = 2*L*mu_0/(mu_0+mu) + (1-L)*mu_0;
  }
  if (reflectance < 0){
    //printf("negative reflectance\n");
    reflectance = 0;
  }
  return reflectance;
}

float
vw::photometry::ComputeReflectance(Vector3 normal, Vector3 xyz,
                                   ModelParams input_img_params,
                                   GlobalParams globalParams) {
  float input_img_reflectance;

  switch ( globalParams.reflectanceType )
    {
    case LUNAR_LAMBERT:
      //printf("Lunar Lambert\n");
      input_img_reflectance = computeLunarLambertianReflectanceFromNormal(input_img_params.sunPosition,
                                                                          input_img_params.spacecraftPosition,
                                                                          xyz,  normal);
      break;
    case LAMBERT:
      //printf("Lambert\n");
      input_img_reflectance = computeLambertianReflectanceFromNormal(input_img_params.sunPosition,
                                                                     xyz,  normal);
      break;

    default:
      //printf("No reflectance model\n");
      input_img_reflectance = 1;
    }

  return input_img_reflectance;

}


//computes a reflectance image
//author: Ara Nefian
float vw::photometry::computeImageReflectance(ModelParams input_img_params,
                                              GlobalParams globalParams) {
  int l, k;
  int count = 0;
  float avg_reflectance = 0.0;

  std::string input_img_file = input_img_params.inputFilename;
  std::string DEM_file = input_img_params.meanDEMFilename;
  std::string shadow_file = input_img_params.shadowFilename;
  std::string output_img_file = input_img_params.reliefFilename;;

  DiskImageView<PixelMask<PixelGray<uint8> > >  input_img(input_img_file);
  GeoReference input_img_geo;
  read_georeference(input_img_geo, input_img_file);

  DiskImageView<PixelGray<float> >  input_dem_image(DEM_file);
  GeoReference input_dem_geo;
  read_georeference(input_dem_geo, DEM_file);

  DiskImageView<PixelMask<PixelGray<uint8> > >  shadowImage(shadow_file);

  ImageView<PixelMask<PixelGray<float> > > output_img (input_img.cols(), input_img.rows());


  Vector3 xyz;
  Vector3 xyz_prior;
  //int x, y;
  float x, y;

  ImageViewRef<PixelGray<float> >  interp_dem_image = interpolate(edge_extend(input_dem_image.impl(),
                                                                              ConstantEdgeExtension()),
                                                                  BilinearInterpolation());
  //interpolate(image, BilinearInterpolation());

  //initialize the nominator and denomitor images
  for (k = 0 ; k < (int)input_img.rows(); ++k) {
    for (l = 0; l < (int)input_img.cols(); ++l) {

      Vector2 input_image_pix(l,k);

      if ( is_valid(input_img(l,k)) ) {

        //get the corresponding DEM value

        Vector2 lon_lat = input_img_geo.pixel_to_lonlat(input_image_pix);
        Vector2 input_dem_pix = input_dem_geo.lonlat_to_pixel(input_img_geo.pixel_to_lonlat(input_image_pix));

        //x = (int)input_dem_pix[0];
        //y = (int)input_dem_pix[1];
	x = (float)input_dem_pix[0];
        y = (float)input_dem_pix[1];

        //check for valid DEM coordinates to be within the DEM boundaries
        if ((x>=0) && (x < input_dem_image.cols()) && (y>=0) && (y< input_dem_image.rows())){

          Vector3 longlat3(lon_lat(0),lon_lat(1),(interp_dem_image)(x, y));
          Vector3 xyz = input_img_geo.datum().geodetic_to_cartesian(longlat3);

          Vector2 input_img_left_pix;
          input_img_left_pix(0) = l-1;
          input_img_left_pix(1) = k;

          Vector2 input_img_top_pix;
          input_img_top_pix(0) = l;
          input_img_top_pix(1) = k-1;

          //check for valid DEM pixel value and valid left and top coordinates
          if ((input_img_left_pix(0) >= 0) && (input_img_top_pix(1) >= 0) && (input_dem_image(x,y) != globalParams.noDEMDataValue)){

            //determine the 3D coordinates of the pixel left of the current pixel
            Vector2 input_dem_left_pix = input_dem_geo.lonlat_to_pixel(input_img_geo.pixel_to_lonlat(input_img_left_pix));
            Vector2 lon_lat_left = input_img_geo.pixel_to_lonlat(input_img_left_pix);
            Vector3 longlat3_left(lon_lat_left(0),lon_lat_left(1),(interp_dem_image)(input_dem_left_pix(0), input_dem_left_pix(1)));
            Vector3 xyz_left = input_img_geo.datum().geodetic_to_cartesian(longlat3_left);

            //determine the 3D coordinates of the pixel top of the current pixel
            Vector2 input_dem_top_pix = input_dem_geo.lonlat_to_pixel(input_img_geo.pixel_to_lonlat(input_img_top_pix));
            Vector2 lon_lat_top = input_img_geo.pixel_to_lonlat(input_img_top_pix);
            Vector3 longlat3_top(lon_lat_top(0),lon_lat_top(1),(interp_dem_image)(input_dem_top_pix(0), input_dem_top_pix(1)));
            Vector3 xyz_top = input_img_geo.datum().geodetic_to_cartesian(longlat3_top);

            Vector3 normal = computeNormalFrom3DPointsGeneral(xyz, xyz_left, xyz_top);
            //This part is the only image depedent part - START

            float input_img_reflectance;
            input_img_reflectance = ComputeReflectance(normal, xyz, input_img_params, globalParams);

            output_img(l, k) = input_img_reflectance;

            if ((input_img_reflectance != 0.0) && (shadowImage(l, k) == 0)){ //valid non zero reflectance
              avg_reflectance = avg_reflectance + output_img(l,k);
              count++;
            }

          }
        }
      }
    }
  }

  avg_reflectance = avg_reflectance/count;
  printf("avg_reflectance = %f\n", avg_reflectance);
  write_georeferenced_image(output_img_file,
                            output_img,
                            input_img_geo, TerminalProgressCallback("{Core}","Processing:"));

  return avg_reflectance;
}


//computes a reflectance image
//author: Ara Nefian
float vw::photometry::computeImageReflectance(ModelParams input_img_params,
                                              ModelParams overlap_img_params,
                                              GlobalParams globalParams) {
  int l, k;
  int count = 0;
  float avg_reflectance = 0.0;
  float overlap_avg_reflectance = 0.0;

  std::string input_img_file = input_img_params.inputFilename;
  std::string overlap_img_file = overlap_img_params.inputFilename;
  std::string DEM_file = input_img_params.meanDEMFilename;
  std::string shadow_file = input_img_params.shadowFilename;
  std::string overlap_shadow_file = overlap_img_params.shadowFilename;
  std::string output_img_file = input_img_params.reliefFilename;

  DiskImageView<PixelMask<PixelGray<uint8> > >  input_img(input_img_file);
  GeoReference input_img_geo;
  read_georeference(input_img_geo, input_img_file);

  DiskImageView<PixelGray<float> >  input_dem_image(DEM_file);
  GeoReference input_dem_geo;
  read_georeference(input_dem_geo, DEM_file);

  DiskImageView<PixelMask<PixelGray<uint8> > >  shadowImage(shadow_file);

  ImageView<PixelMask<PixelGray<float> > > output_img (input_img.cols(), input_img.rows());

  Vector3 xyz;
  Vector3 xyz_prior;
  int x, y;

  ImageViewRef<PixelGray<float> >  interp_dem_image = interpolate(edge_extend(input_dem_image.impl(),
                                                                              ConstantEdgeExtension()),
                                                                  BilinearInterpolation());



  // start
  DiskImageView<PixelMask<PixelGray<uint8> > >  overlap_img(overlap_img_file);
  GeoReference overlap_geo;
  read_georeference(overlap_geo, overlap_img_file);

  DiskImageView<PixelMask<PixelGray<uint8> > >  overlapShadowImage(overlap_shadow_file);

  ImageViewRef<PixelMask<PixelGray<uint8> > >  interp_overlap_img = interpolate(edge_extend(overlap_img.impl(),
                                                                                            ConstantEdgeExtension()),
                                                                                BilinearInterpolation());

  ImageViewRef<PixelMask<PixelGray<uint8> > >  interpOverlapShadowImage = interpolate(edge_extend(overlapShadowImage.impl(),
                                                                                                  ConstantEdgeExtension()),
                                                                                      BilinearInterpolation());
  //end

  //initialize the nominator and denomitor images
  for (k = 0 ; k < (int)input_img.rows(); ++k) {
    for (l = 0; l < (int)input_img.cols(); ++l) {

      Vector2 input_image_pix(l,k);

      if ( is_valid(input_img(l,k)) ) {

        //get the corresponding DEM value

        Vector2 lon_lat = input_img_geo.pixel_to_lonlat(input_image_pix);
        Vector2 input_dem_pix = input_dem_geo.lonlat_to_pixel(input_img_geo.pixel_to_lonlat(input_image_pix));

        x = (int)input_dem_pix[0];
        y = (int)input_dem_pix[1];

        //check for valid DEM coordinates to be within the DEM boundaries
        if ((x>=0) && (x < input_dem_image.cols()) && (y>=0) && (y< input_dem_image.rows())){

          Vector3 longlat3(lon_lat(0),lon_lat(1),(interp_dem_image)(x, y));
          Vector3 xyz = input_img_geo.datum().geodetic_to_cartesian(longlat3);

          Vector2 input_img_left_pix;
          input_img_left_pix(0) = l-1;
          input_img_left_pix(1) = k;

          Vector2 input_img_top_pix;
          input_img_top_pix(0) = l;
          input_img_top_pix(1) = k-1;

          //check for valid DEM pixel value and valid left and top coordinates
          if ((input_img_left_pix(0) >= 0) && (input_img_top_pix(1) >= 0) && (input_dem_image(x,y) != globalParams.noDEMDataValue/*-10000*/)){

            //determine the 3D coordinates of the pixel left of the current pixel
            Vector2 input_dem_left_pix = input_dem_geo.lonlat_to_pixel(input_img_geo.pixel_to_lonlat(input_img_left_pix));
            Vector2 lon_lat_left = input_img_geo.pixel_to_lonlat(input_img_left_pix);
            Vector3 longlat3_left(lon_lat_left(0),lon_lat_left(1),(interp_dem_image)(input_dem_left_pix(0), input_dem_left_pix(1)));
            Vector3 xyz_left = input_img_geo.datum().geodetic_to_cartesian(longlat3_left);

            //determine the 3D coordinates of the pixel top of the current pixel
            Vector2 input_dem_top_pix = input_dem_geo.lonlat_to_pixel(input_img_geo.pixel_to_lonlat(input_img_top_pix));
            Vector2 lon_lat_top = input_img_geo.pixel_to_lonlat(input_img_top_pix);
            Vector3 longlat3_top(lon_lat_top(0),lon_lat_top(1),(interp_dem_image)(input_dem_top_pix(0), input_dem_top_pix(1)));
            Vector3 xyz_top = input_img_geo.datum().geodetic_to_cartesian(longlat3_top);

            Vector3 normal = computeNormalFrom3DPointsGeneral(xyz, xyz_left, xyz_top);
            //This part is the only image depedent part - START

            float input_img_reflectance;
            input_img_reflectance = ComputeReflectance(normal, xyz, input_img_params, globalParams);

            output_img(l, k) = input_img_reflectance;

            //--------------------------------------------------------------
            //check for overlap between the output image and the input DEM image
            Vector2 input_img_pix(l,k);
            Vector2 overlap_pix = overlap_geo.lonlat_to_pixel(input_img_geo.pixel_to_lonlat(input_img_pix));
            x = (int)overlap_pix[0];
            y = (int)overlap_pix[1];

            //image dependent part of the code  - START
            PixelMask<PixelGray<uint8> > overlap_img_pixel = interp_overlap_img(x, y);

            //check for valid overlap_img coordinates
            //remove shadow pixels in the overlap_img.
            if ((x>=0) && (x < overlap_img.cols()) && (y>=0) && (y< overlap_img.rows()) && (interpOverlapShadowImage(x, y) == 0) && (shadowImage(l,k) == 0)) {

              if ( is_valid(overlap_img_pixel) ) { //common area between input_img and overlap_img

                float overlap_img_reflectance;
                overlap_img_reflectance = ComputeReflectance(normal, xyz, overlap_img_params, globalParams);
                if ((overlap_img_reflectance > 0) && (input_img_reflectance > 0)){
                  avg_reflectance = avg_reflectance + output_img(l,k);
                  overlap_avg_reflectance = overlap_avg_reflectance + overlap_img_reflectance;
                  count++;
                }
              }//if
            }//if
            //----------------------------
          }
        }
      }
    }
  }

  avg_reflectance = avg_reflectance/count;
  overlap_avg_reflectance = overlap_avg_reflectance/count;
  float reflectance_ratio = avg_reflectance/overlap_avg_reflectance;
  printf("avg_reflectance = %f, overlap_avg_reflectance = %f\n", avg_reflectance, overlap_avg_reflectance);
  write_georeferenced_image(output_img_file,
                            output_img,
                            input_img_geo, TerminalProgressCallback("{Core}","Processing:"));

  return reflectance_ratio;
}



