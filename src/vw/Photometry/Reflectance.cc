// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
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
#include <vw/Photometry/Misc.h>
#include <vw/Photometry/Weights.h>
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

  double len = dot_prod(normal, normal);
  if (abs(len - 1.0) > 1.0e-4){
    std::cerr << "Error: Expecting unit normal in the reflectance computation, in "
              << __FILE__ << " at line " << __LINE__ << std::endl;
    exit(1);
  }

  //compute /mu_0 = cosine of the angle between the light direction and the surface normal.
  //sun coordinates relative to the xyz point on the Moon surface
  //Vector3 sunDirection = -normalize(sunPos-xyz);
  Vector3 sunDirection = normalize(sunPos-xyz);
  float mu_0 = dot_prod(sunDirection,normal);

  double tol = 0.3;
  if (mu_0 < tol){
    // Sun is too low, reflectance is too close to 0, the albedo will be inaccurate
    return 0;
  }
  
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
  deg_alpha = rad_alpha*180.0/M_PI;

  //printf("deg_alpha = %f\n", deg_alpha);

  //Bob Gaskell's model
  //L = exp(-deg_alpha/60.0);

#if 0 // trey
  // perfectly valid for alpha to be greater than 90?
  if (deg_alpha > 90){
    //printf("Error!!: rad_alpha = %f, deg_alpha = %f\n", rad_alpha, deg_alpha);
    return(0.0);
  }
  if (deg_alpha < -90){
    //printf("Error!!: rad_alpha = %f, deg_alpha = %f\n", rad_alpha, deg_alpha);
    return(0.0);
  }
#endif

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

void vw::photometry::computeXYZandSurfaceNormal(ImageView<PixelGray<float> > const& DEMTile,
                                                cartography::GeoReference const& DEMGeo,
                                                GlobalParams globalParams,
                                                ImageView<Vector3> & dem_xyz,
                                                ImageView<Vector3> & surface_normal){

  int nodata = globalParams.noDEMDataValue;
 
  // convert dem altitude to xyz cartesian pixels
  dem_xyz = geodetic_to_cartesian( dem_to_geodetic( DEMTile, DEMGeo ),
                                   DEMGeo.datum() );

  // transfer nodata values
  for (int y=0; y < (int)dem_xyz.rows(); y++) {
    for (int x=0; x < (int)dem_xyz.cols(); x++) {
      if (DEMTile(x, y) == nodata) {
        dem_xyz(x, y) = Vector3();
      }
    }
  }  

  // Init the surface normal to zero
  surface_normal.set_size(dem_xyz.cols(), dem_xyz.rows());
  for (int y=0; y < (int)surface_normal.rows(); y++) {
    for (int x=0; x < (int)surface_normal.cols(); x++) {
      surface_normal(x, y) = Vector3();
    }
  }
  
  // convert xyz pixels to surface normals
  for (int y=1; y < (int)surface_normal.rows(); y++) {
    for (int x=1; x < (int)surface_normal.cols(); x++) {
      
      Vector3 base = dem_xyz(x, y);
      if (base == Vector3()) {
        continue;
      }
      Vector3 x1 = dem_xyz(x-1, y);
      if (x1 == Vector3()) {
        continue;
      }
      Vector3 y1 = dem_xyz(x, y-1);
      if (y1 == Vector3()) {
        continue;
      }

      Vector3 dx = base - x1;
      Vector3 dy = base - y1;
      surface_normal(x, y) = -normalize(cross_prod(dx, dy));
    }
  }

  return;
}

void vw::photometry::computeReflectanceAux(ImageView<Vector3> const& dem_xyz,
                                           ImageView<Vector3> const& surface_normal,
                                           ModelParams input_img_params,
                                           GlobalParams globalParams,
                                           ImageView<PixelMask<PixelGray<float> > >& outputReflectance) {

  //std::cout << "---sun        " << input_img_params.sunPosition << std::endl;
  //std::cout << "---spacecraft " << input_img_params.spacecraftPosition << std::endl;
  outputReflectance.set_size(surface_normal.cols(), surface_normal.rows());
  for (int y=0; y < (int)outputReflectance.rows(); y++) {
    for (int x=0; x < (int)outputReflectance.cols(); x++) {
      Vector3 normal = surface_normal(x, y);
      if (normal == Vector3()) {
        outputReflectance(x, y) = 0;
        outputReflectance(x, y).invalidate();
      }else{
        Vector3 xyz = dem_xyz(x, y);
        outputReflectance(x, y) = ComputeReflectance(normal, xyz, input_img_params, globalParams);
      }
    }
  }
  
  return;
}

float vw::photometry::computeImageReflectanceNoWrite(ModelParams input_img_params,
                                                     GlobalParams globalParams,
                                                     ImageView<PixelMask<PixelGray<float> > >& output_img) {
  
  std::string input_img_file = input_img_params.inputFilename;
  std::string DEM_file = input_img_params.meanDEMFilename;
  std::string shadow_file = input_img_params.shadowFilename;

  DiskImageView<PixelMask<PixelGray<uint8> > >  input_img(input_img_file);
  GeoReference input_img_geo;
  read_georeference(input_img_geo, input_img_file);

  // warp dem to drg georef
  ImageView<PixelGray<double> > dem_in_drg_georef_0(input_img.cols(), input_img.rows());
  DiskImageView<PixelGray<float> >  input_dem_image(DEM_file);
  GeoReference input_dem_geo;
  read_georeference(input_dem_geo, DEM_file);
  dem_in_drg_georef_0 =
    crop
    (geo_transform
     (input_dem_image,
      input_dem_geo,
      input_img_geo,
      ConstantEdgeExtension(),
      BilinearInterpolation()),
     bounding_box(input_img)
     );
  
#if 0
  std::cout << std::endl;
  std::string tmp  = "./" + prefix_from_filename(sufix_from_filename(input_img_file)) + "_dem.tif";
  std::cout << "Writing " << tmp << std::endl;
  write_georeferenced_image(tmp,
                            dem_in_drg_georef_0,
                            input_img_geo, TerminalProgressCallback("{Core}","Processing:"));
#endif
    
  // grow no-data areas to avoid bogus values created by interpolating
  // data with no-data
  ImageView<PixelGray<double> > dem_in_drg_georef = copy(dem_in_drg_georef_0);
  int nodata = globalParams.noDEMDataValue;
  for (int y=1; y < (int)dem_in_drg_georef.rows()-1; y++) {
    for (int x=1; x < (int)dem_in_drg_georef.cols()-1; x++) {

#define CHECK_PIX(dx,dy) \
  if (dem_in_drg_georef_0(x+dx, y+dy) == nodata) { \
    dem_in_drg_georef(x, y) = nodata; \
    continue; \
  }

      CHECK_PIX(-1, -1);
      CHECK_PIX(-1,  0);
      CHECK_PIX(-1,  1);
      CHECK_PIX( 0, -1);
      CHECK_PIX( 0,  1);
      CHECK_PIX( 1, -1);
      CHECK_PIX( 1,  0);
      CHECK_PIX( 1,  1);
    }
  }

  ImageView<Vector3> dem_xyz;
  ImageView<Vector3> surface_normal;
  computeXYZandSurfaceNormal(dem_in_drg_georef, input_img_geo, globalParams, // Inputs 
                             dem_xyz, surface_normal        // Outputs
                             );
  computeReflectanceAux(dem_xyz, surface_normal,  
                        input_img_params,
                        globalParams,  
                        output_img
                        );
  
  // compute average reflectance
  int count = 0;
  float reflectance_sum = 0.0;
  DiskImageView<PixelMask<PixelGray<uint8> > >  shadowImage(shadow_file);
  for (int y=0; y < (int)output_img.rows(); y++) {
    for (int x=0; x < (int)output_img.cols(); x++) {
      float refl = output_img(x, y);
      if ((refl != 0.0) && (shadowImage(x, y) == 0)){ //valid non zero reflectance
        reflectance_sum += refl;
        count++;
      }
    }
  }
  float avg_reflectance = reflectance_sum / count;

  printf("avg_reflectance = %f\n", avg_reflectance);
  return avg_reflectance;
}

float vw::photometry::computeAvgReflectanceOverTilesOrUpdateExposure(bool compAvgRefl,
                                                                     bool useReflectance,
                                                                     int pixelPadding, double tileSize,
                                                                     std::vector<ImageRecord> & DEMTiles,
                                                                     std::vector<ImageRecord> & albedoTiles,
                                                                     std::vector<int> & overlap,
                                                                     ModelParams input_img_params,
                                                                     GlobalParams globalParams){

  // Compute the average reflectance of the current image, or update
  // the exposure of that image. Both of these operations use very similar
  // logic. Namely, we will iterate over all tiles overlapping with
  // the image, and accumulate the results.
  
  // We will keep in mind that the tiles partially overlap, so when we
  // accumulate the results we will ignore the padded region of the
  // tile (of size pixelPadding).

  // Use double precision for the intermediate computations, even though
  // the input and output are float, as the float type is not precise
  // enough when accumulating a lot of numbers.
  
  std::string input_img_file = input_img_params.inputFilename;
  DiskImageView<PixelMask<PixelGray<uint8> > >  input_img(input_img_file);
  GeoReference input_img_geo;
  read_georeference(input_img_geo, input_img_file);

  //std::string shadow_file = input_img_params.shadowFilename;
  //DiskImageView<PixelMask<PixelGray<uint8> > >  shadowImage(shadow_file);

  bool useWeights = (globalParams.useWeights != 0);
  if (useWeights && (!compAvgRefl)){
    // Read the weights only if we update the exposure
    bool useTiles = true;
    ReadWeightsParamsFromFile(useTiles, &input_img_params);
  }
  
  // Quantities needed for computing reflectance
  int count = 0;
  double reflectance_sum = 0.0;

  // Quantities needed for updating the exposure
  double numerator = 0.0, denominator = 0.0;
  
  // Iterate over the DEM tiles
  for (int i = 0; i < (int)overlap.size(); i++){

    std::string DEMTileFile = DEMTiles[overlap[i]].path;
    std::cout << "Reading file: "<< DEMTileFile << std::endl;
    DiskImageView<PixelGray<float> > DEMTile(DEMTileFile);
    GeoReference DEMGeo;
    read_georeference(DEMGeo, DEMTileFile);
    ImageView<PixelMask<PixelGray<float> > > Reflectance;
    if (useReflectance){
      // Declare dem_xyz and surface_normal in a block to de-allocate
      // fast the quantities dem_xyz and surface_normal which are
      // needed only temporarily. They take a lot of memory being
      // images of vectors.
      ImageView<Vector3> dem_xyz, surface_normal;
      vw::photometry::computeXYZandSurfaceNormal(DEMTile, DEMGeo, globalParams,
                                                 dem_xyz, surface_normal
                                                 );
      computeReflectanceAux(dem_xyz, surface_normal,
                            input_img_params, globalParams,  
                            Reflectance // output
                            );
      //system("echo reflectance_xyz top is $(top -u $(whoami) -b -n 1|grep lt-reconstruct)");
    }else{
      // The reflectance is set to 1.
      Reflectance.set_size(DEMTile.cols(), DEMTile.rows());
      for (int row = 0; row < (int)DEMTile.rows(); row++){
        for (int col = 0; col < (int)DEMTile.cols(); col++){
          Reflectance(col, row) = 1.0;
          Reflectance(col, row).validate();
        }
      }
    }
    
    InterpolationView<EdgeExtensionView<ImageView< PixelMask<PixelGray<float> > >,
                                        ConstantEdgeExtension>, BilinearInterpolation>
      interp_reflectance = interpolate(Reflectance,
                                       BilinearInterpolation(), ConstantEdgeExtension());
    
    
    ImageView<PixelMask<PixelGray<float> > > albedoTile;
    if (!compAvgRefl){
      // Update the exposure. Need the current albedo.
      std::string albedoTileFile = albedoTiles[overlap[i]].path;
      std::cout << "Reading " << albedoTileFile << std::endl;
      albedoTile = copy(DiskImageView<PixelMask<PixelGray<uint8> > >(albedoTileFile));
    }
    InterpolationView<EdgeExtensionView<ImageView< PixelMask<PixelGray<float> > >, ConstantEdgeExtension>, BilinearInterpolation>
      interp_albedo = interpolate(albedoTile,
                                  BilinearInterpolation(), ConstantEdgeExtension());


    // We need to keep in mind that the tile is padded with pixelPadding pixels on each side.
    // To avoid double-counting pixels, skip them if they are part of the padding and not
    // of the tile proper.
    double min_tile_x, max_tile_x, min_tile_y, max_tile_y;
    getTileCornersWithoutPadding(// Inputs
                                 DEMTile.cols(),  DEMTile.rows(), DEMGeo,  
                                 tileSize, pixelPadding,  
                                 // Outputs
                                 min_tile_x, max_tile_x,  
                                 min_tile_y, max_tile_y
                                 );
    
    // Iterate only over the portion of the image intersecting the DEM tile
    Vector2 beg_pixel = input_img_geo.lonlat_to_pixel(Vector2(min_tile_x, max_tile_y));
    Vector2 end_pixel = input_img_geo.lonlat_to_pixel(Vector2(max_tile_x, min_tile_y));
    int beg_row = std::max(0,                (int)floor(beg_pixel(1)));
    int end_row = std::min(input_img.rows(), (int)ceil(end_pixel(1)));
    int beg_col = std::max(0,                (int)floor(beg_pixel(0)));
    int end_col = std::min(input_img.cols(), (int)ceil(end_pixel(0)));
    for (int k = beg_row; k < end_row; ++k) {
      for (int l = beg_col; l < end_col; ++l) {

        Vector2 input_img_pix(l,k);

        Vector2 lonlat = input_img_geo.pixel_to_lonlat(input_img_pix);
        double lx = lonlat(0), ly = lonlat(1);
        bool isInTile = (min_tile_x <= lx && lx < max_tile_x && min_tile_y <= ly && ly < max_tile_y);
        if (!isInTile) continue;

        //if (k%2 != 0 && l%2 != 0) continue; // Skip some points when computing the exposure to save time
          
        //check for valid DEM coordinates
        Vector2 reflectance_pix = DEMGeo.lonlat_to_pixel(lonlat);
        double x = reflectance_pix[0];
        double y = reflectance_pix[1];
        double t = globalParams.shadowThresh; // Check if the image is above the shadow threshold
        if ( is_valid(input_img(l, k)) && ( (double)input_img(l, k) >= t) ) {
          if ( (x>=0) && (x <= Reflectance.cols()-1) && (y>=0) && (y<= Reflectance.rows()-1)){
            // Check that all four grid points used for interpolation are valid
            if ( is_valid(Reflectance( floor(x), floor(y))) && (double)Reflectance( floor(x), floor(y) ) != 0 &&
                 is_valid(Reflectance( floor(x), ceil (y))) && (double)Reflectance( floor(x), ceil(y)  ) != 0 &&
                 is_valid(Reflectance( ceil (x), floor(y))) && (double)Reflectance( ceil (x), floor(y) ) != 0 &&
                 is_valid(Reflectance( ceil (x), ceil (y))) && (double)Reflectance( ceil (x), ceil(y)  ) != 0
                 ){
              double R = interp_reflectance(x, y);

              if (compAvgRefl){
                reflectance_sum += R;
                count++;
              }else{
                // Update the exposure
                // We assume that the DEM, reflectance, and albedo are on the same grid
                double weight = 1.0;
                if (useWeights) weight = ComputeLineWeightsHV(input_img_pix, input_img_params);
                double A  = interp_albedo(x, y); // we assume albedo, reflectance, and DEM are on the same grid
                double RA = R*A;
                numerator   += ((double)input_img(l, k) - input_img_params.exposureTime*RA)*RA*weight;
                denominator += RA*RA*weight;
              }
            }
          }
        }
      }
    }
    //system("echo reflectance top is $(top -u $(whoami) -b -n 1|grep lt-reconstruct)");
  }

  if (compAvgRefl){
    // Compute the average reflectance
    double avg_reflectance = reflectance_sum/count;
    printf("avg_reflectance = %f\n", avg_reflectance);
    return avg_reflectance;
  }else{
    // Update the exposure
    double exposure = input_img_params.exposureTime + numerator/denominator;
    printf("updated exposure is = %f\n", exposure);
    return exposure;
  }

  return 0;
}

//computes a reflectance image
//author: Ara Nefian
float vw::photometry::computeImageReflectance(ModelParams input_img_params,
                                              GlobalParams globalParams) {
  std::string input_img_file = input_img_params.inputFilename;
  std::string output_img_file = input_img_params.reliefFilename;;

  DiskImageView<PixelMask<PixelGray<uint8> > >  input_img(input_img_file);
  GeoReference input_img_geo;
  read_georeference(input_img_geo, input_img_file);

  ImageView<PixelMask<PixelGray<float> > > output_img;
  float avg_reflectance =
    computeImageReflectanceNoWrite(input_img_params,
                                   globalParams,
                                   output_img);


  std::cout << "Writing " << output_img_file << std::endl;
  write_georeferenced_image
    (output_img_file,
     output_img,
     input_img_geo,
     TerminalProgressCallback("{Core}","Processing:"));
  
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

  InterpolationView<EdgeExtensionView<EdgeExtensionView<DiskImageView<PixelGray<float> >, ConstantEdgeExtension>, ConstantEdgeExtension>, BilinearInterpolation>
      interp_dem_image = interpolate(edge_extend(input_dem_image.impl(),
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
