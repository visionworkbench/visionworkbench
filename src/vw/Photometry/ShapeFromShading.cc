// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <iostream>
#include <fstream>
#include <vw/Image.h>
#include <vw/Image/PixelMath.h>
#include <vw/Image/PixelMask.h>
#include <vw/Image/MaskViews.h>
#include <vw/FileIO.h>
#include <math.h>
#include <time.h>

#include <vw/Core.h>
#include <vw/Cartography.h>
#include <vw/Math.h>


#include <vw/Math/Matrix.h>

#include <vw/Math/LinearAlgebra.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Image/Interpolation.h>
#include <vw/Image/ImageMath.h>
#include <vw/Photometry/Reflectance.h>
#include <vw/Photometry/Weights.h>
#include <vw/Photometry/Misc.h>

#include <cmath>
//#include <ctime>

using namespace std;
using namespace vw;
using namespace vw::cartography;

#include <math.h>
#include <vw/Photometry/Reconstruct.h>
#include <vw/Photometry/ReconstructError.h>
#include <vw/Photometry/ShapeFromShading.h>
using namespace vw::photometry;

#define horBlockSize 8 //8 //4
#define verBlockSize 8 //8 //4
//#define HOR_BLOCK_SIZE 16 //8 //4
//#define VER_BLOCK_SIZE 16 //8 //4
#define numJacobianRows (horBlockSize+1)*(verBlockSize+1)
#define numJacobianCols (horBlockSize*verBlockSize)


//compute the elements of the normal derivative
Vector3 ComputeNormalDerivative(int flag,  Vector3 xyz, Vector3 xyzTOP, Vector3 xyzLEFT)
{
  // xo = (xo0,xo1,xo2), xt = (xt0,xt1,xt2), xl = (xl0,xl1,xl2)
  // normal = cross_prod(xt-xo, xl-xo);
  // normal = [(xt1-xo1)*(xl2-xo2)-(xt2-xo2)*(xl1-xo1),
  //           (xt2-xo2)*(xl0-xo0)-(xt0-xo0)*(xl2-xo2),
  //           (xt0-xo0)*(xl1-xo1)-(xt1-xo1)*(xl0-xo0)]
  // d_n/xo2 = (xl1-xt1,xt0-xl0,0)
  // d_n/xl2 = (xt1-xo1,xo0-xt0,0)
  // d_n/xt2 = (xo1-xl1,xl0-xo0,0)

  Vector3 normalDerivative;
  if (flag == 0){ //wrt z_{i,j}
    normalDerivative(0) = xyzLEFT(1) - xyzTOP(1); //dn_x/dz_{ij}
    normalDerivative(1) = xyzTOP(0) - xyzLEFT(0); //dn_y/dz_{ij}
    normalDerivative(2) = 0; //dnz_dz_{ij}
  }
  if (flag == 1){ //wrt z_{i,j-1} //LEFT
    normalDerivative(0) = xyzTOP(1) - xyz(1); //dn_x/dz_{i, j-1}
    normalDerivative(1) = xyz(0) - xyzTOP(0); //dn_y/dz_{i, j-1}
    normalDerivative(2) = 0; //dnz_dz_{i, j-1}
  }
  if (flag == 2){ //wrt z_{i-1,j}//TOP
    normalDerivative(0) = xyz(1) - xyzLEFT(1); //dn_x/dz_{i-1, j}
    normalDerivative(1) = xyzLEFT(0) - xyz(0); //dn_y/dz_{i-1, j}
    normalDerivative(2) = 0; //dnz_dz_{i-1, j}
  }

  return normalDerivative;
}

//computes the cosine derivative wrt to the height map of the angle between two vectors stored in variables normal and direction
//direction vector is normalized
float ComputeCosDerivative(Vector3 normal, Vector3 direction, Vector3 normalDerivative)
{
  float normalNorm = normal(0)*normal(0) + normal(1)*normal(1) + normal(2)*normal(2);
  float denominator = normalNorm*sqrt(normalNorm);
  float nominator = (normalDerivative(0)*direction(0)+normalDerivative(1)*direction(1))*normalNorm - (normalDerivative(0) + normalDerivative(1))*(normal(0)*direction(0) + normal(1)*direction(1) + normal(2)*direction(2));
  float cosEDeriv = nominator/denominator;
  return cosEDeriv;
}

float ComputeReliefDerivative(Vector3 xyz,Vector3 xyzLEFT,Vector3 xyzTOP, Vector3 normal, ModelParams inputImgParams, int flag)
{
  float reliefDeriv;
  Vector3 sunPos = inputImgParams.sunPosition;
  Vector3 viewPos = inputImgParams.spacecraftPosition;

  Vector3 normalDerivative = ComputeNormalDerivative(flag, xyz, xyzTOP, xyzLEFT);

  //compute /mu_0 = cosine of the angle between the light direction and the surface normal.
  // sun coordinates relative to the xyz point on the Moon surface
  Vector3 sunDirection = normalize(sunPos-xyz);
  float mu_0 = dot_prod(sunDirection,normal);

  //compute  /mu = cosine of the angle between the viewer direction and the surface normal.
  // viewer coordinates relative to the xyz point on the Moon surface
  Vector3 viewDirection = normalize(viewPos-xyz);
  float mu = dot_prod(viewDirection,normal);

  float cosEDeriv = ComputeCosDerivative(normal, viewDirection, normalDerivative);
  float cosIDeriv = ComputeCosDerivative(normal, sunDirection, normalDerivative);

  //Alfred McEwen's model
  float A = -0.019*180/3.141592;
  float B =  0.000242*180/3.141592*180/3.141592;//0.242*1e-3;
  float C = -0.00000146*180/3.141592*180/3.141592*180/3.141592;//-1.46*1e-6;

  float cos_alpha = dot_prod(sunDirection,viewDirection);

  if ((cos_alpha > 1)||(cos_alpha< -1)){
    printf("cos_alpha error\n");
  }

  float rad_alpha = acos(cos_alpha);
  float L = 1.0 + A*rad_alpha + B*rad_alpha*rad_alpha + C*rad_alpha*rad_alpha*rad_alpha;

  //reliefDeriv = (mu+mu_0) ? (1-L)*cosIDeriv + L*(cosIDeriv*(mu+mu_0)+(cosEDeriv+cosIDeriv)*mu)/((mu+mu_0)*(mu+mu_0)) : 0;
  reliefDeriv = (mu+mu_0) ? (1-L)*cosIDeriv + 2*L*(cosIDeriv*(mu+mu_0)+(cosEDeriv+cosIDeriv)*mu)/((mu+mu_0)*(mu+mu_0)) : 0;
  /*
  if (mu_0<0 || mu<0) {
    std::cout << " sun direction " << sunDirection << " view direction " << viewDirection << " normal " << normal << std::endl;
    std::cout << " cos_alpha " << cos_alpha << " incident " << mu_0 << " emission " << mu << " L = " << L << " reliefDeriv = " << reliefDeriv << std::endl;
  } else {

  }
  */
  return reliefDeriv;
}

void printJacobian( Matrix<float, numJacobianRows, numJacobianCols>  &jacobian, string filename)
{
  FILE *fp = fopen (filename.c_str(), "w");

  for (int r = 0; r < numJacobianRows; r++){
    for (int c = 0; c < numJacobianCols; c++){
      if (jacobian(r, c)!=0){
        fprintf(fp, "(%d %d) : %f\n", r, c, jacobian(r,c));
      }
    }
  }

  fclose(fp);
}

bool 
ComputeBlockGeometryTiles(ImageView<PixelGray<float> > const& dem, GeoReference const &demGeo,
                          float noDEMDataValue,
                          int extraPixel, int kb, int lb, 
                          vector<Vector3> &xyzArray, vector<Vector3> &xyzLEFTArray, vector<Vector3> &xyzTOPArray,
                          vector<Vector3> &normalArray){

  int eVerBlockSize = verBlockSize+1;
  int eHorBlockSize = horBlockSize+1;

  for (int k = 0 ; k < eVerBlockSize; ++k) {
    for (int l = 0; l < eHorBlockSize; ++l) {
      
      int ii = kb*(verBlockSize)+k+extraPixel; //row index for the entire image
      int jj = lb*(horBlockSize)+l+extraPixel; //col index for the entire image
      //printf("ii = %d, jj = %d, width = %d, height = %d\n", ii, jj, dem.impl().cols(), dem.impl().rows());

      // In order to be able to compute the normal, we need to be able
      // to access the current pixel in the DEM, as well as the
      // neighboring left and top pixels. And those pixels must have
      // valid values.
      if ( ii <= 0 || ii >= dem.rows() || jj <= 0 || jj >= dem.cols() ) return false;
      if ( dem(jj,   ii  ) == noDEMDataValue ||
           dem(jj-1, ii  ) == noDEMDataValue ||
           dem(jj,   ii-1) == noDEMDataValue
           ) return false;

      //local index in the vector that describes the block image; assumes row-wise concatenation.
      int l_index = k*eHorBlockSize+l;

      Vector2 lonlat = demGeo.pixel_to_lonlat(Vector2(jj, ii));
      Vector3 lonlat3(lonlat(0),lonlat(1),(dem.impl())(jj, ii));
      xyzArray[l_index] = demGeo.datum().geodetic_to_cartesian(lonlat3);//3D coordinates in the img coordinates

      //determine the 3D coordinates of the pixel left of the current pixel
      Vector2 lonlat_left = demGeo.pixel_to_lonlat(Vector2(jj-1, ii));
      Vector3 longlat3_left(lonlat_left(0),lonlat_left(1),(dem.impl())(jj-1, ii));
      Vector3 xyz_left = demGeo.datum().geodetic_to_cartesian(longlat3_left);
      xyzLEFTArray[l_index] = xyz_left;
      
      //determine the 3D coordinates of the pixel top of the current pixel
      Vector2 lonlat_top = demGeo.pixel_to_lonlat(Vector2(jj, ii-1));
      Vector3 longlat3_top(lonlat_top(0),lonlat_top(1),(dem.impl())(jj, ii-1));
      Vector3 xyz_top = demGeo.datum().geodetic_to_cartesian(longlat3_top);
      xyzTOPArray[l_index] = xyz_top;

      normalArray[l_index] = cross_prod(xyz_top-xyzArray[l_index], xyz_left-xyzArray[l_index]);
      //std::cout << "xyz_top: " << xyz_top << "xyz_left: " << xyz_left << "xyzArray[l_index] " << xyzArray[l_index] << "normalArray: " << normalArray[l_index] << std::endl;
      //std::cout<<"k = "<< k <<",l = "<< l <<", normalArray: " << normalArray[l_index] << std::endl;
    } //for l
  } //for k

  return true;
}


template <class ViewT1>
bool
ComputeBlockJacobianTiles(ImageViewBase<ViewT1> const& overlapImage, GeoReference const &overlapImageGeo,
                          GeoReference const &DEMGeo, ImageViewBase<ViewT1> const& albedoTile,
                          int extraPixel, int kb, int lb, ModelParams overlapImgParams, GlobalParams globalParams,
                          vector<Vector3> const &xyzArray, vector<Vector3> const &xyzLEFTArray,
                          vector<Vector3> const &xyzTOPArray, vector<Vector3> const &normalArray,
                          Matrix<float>  &jacobianArray,
                          Vector<float>  &errorVectorArray,
                          Matrix<float>  &weightsArray){
  
  // Must cast the image to float from uint8, to be able to do interpolation with floats rather
  // than with integers.
  ImageViewRef<PixelMask<PixelGray<float> > > overlapImage_float = channel_cast<float>(overlapImage.impl());
  InterpolationView<EdgeExtensionView<ImageViewRef<PixelMask<PixelGray<float> > >, ConstantEdgeExtension>, BilinearInterpolation>
    interpOverlapImage = interpolate(overlapImage_float,
                                     BilinearInterpolation(), ConstantEdgeExtension());
  
  for (int r = 0; r < numJacobianRows-1; r++){ //last row is always zero
    
    int k = r/(horBlockSize+1);     //row index in the extended block
    int l = r - k*(horBlockSize+1); //col index in the extended block

    int ii = kb*verBlockSize+k+extraPixel; //row index for the entire image
    int jj = lb*horBlockSize+l+extraPixel; //col index for the entire image

    if ( ii >= albedoTile.impl().rows() || jj >= albedoTile.impl().cols() ) return false;
    if ( !is_valid(albedoTile.impl()(jj,ii)) ) return false;
               
    float albedoVal = (float)albedoTile.impl()(jj,ii);
    float expRefl   = overlapImgParams.exposureTime;

    // We will interpolate into the overlap image. Go from the albedo
    // (DEM) pixel to the image pixel. If the interpolation uses
    // invalid image pixels or pixels in the shadow, then the current
    // block cannot be used in the Jacobian computation.
    Vector2 overlapImgPix = overlapImageGeo.lonlat_to_pixel(DEMGeo.pixel_to_lonlat(Vector2(jj, ii)));
    float x = overlapImgPix(0), y = overlapImgPix(1);
    int x0 = (int)floor(x), x1 = (int)ceil(x);
    int y0 = (int)floor(y), y1 = (int)ceil(y);
    float t = globalParams.shadowThresh; // Check if the image is above the shadow threshold
    bool isGood = ((x >= 0) && (x <= overlapImage.impl().cols()-1) &&
                   (y >= 0) && (y <= overlapImage.impl().rows()-1) && 
                   is_valid(overlapImage.impl()(x0, y0)) && (float)overlapImage.impl()(x0, y0) >= t &&
                   is_valid(overlapImage.impl()(x0, y1)) && (float)overlapImage.impl()(x0, y1) >= t &&
                   is_valid(overlapImage.impl()(x1, y0)) && (float)overlapImage.impl()(x1, y0) >= t &&
                   is_valid(overlapImage.impl()(x1, y1)) && (float)overlapImage.impl()(x1, y1) >= t
                   );
    if (!isGood) return false;

    float imgVal = interpOverlapImage(x, y);
    
    int c = k*horBlockSize + l; //current point
    //not computed for the last row and last column of the extended block
    if ((k < verBlockSize) && (l < horBlockSize)){

      float recDer = ComputeReliefDerivative(xyzArray[r], xyzLEFTArray[r],
                                             xyzTOPArray[r], normalArray[r],
                                             overlapImgParams, 0)*albedoVal*expRefl;
      jacobianArray(r, c) = recDer;
    }

    c = k*horBlockSize + l-1;//left point
    //not computed for the first column and last row
    if ((c >= 0) && (l > 0) && ( k < verBlockSize)){
      float recDerLEFT = ComputeReliefDerivative(xyzArray[r], xyzLEFTArray[r],
                                                 xyzTOPArray[r], normalArray[r],
                                                 overlapImgParams, 1)*albedoVal*expRefl;
      jacobianArray(r, c) = recDerLEFT;
    }
        
    c = (k-1)*horBlockSize + l;//top point
    //not computed for the first row and last column of the extended block
    if ((c >= 0) && (k > 0) && (l < horBlockSize)){
      float recDerTOP = ComputeReliefDerivative(xyzArray[r], xyzLEFTArray[r],
                                                xyzTOPArray[r], normalArray[r], overlapImgParams, 2)*albedoVal*expRefl;
      jacobianArray(r, c) = recDerTOP;
    }

    //compute the weights matrix
    if (globalParams.useWeights == 1){
      weightsArray(r,r) = ComputeLineWeightsHV(overlapImgPix, overlapImgParams);
    }
    else{
      weightsArray(r,r) = 1;
    }

    //compute the reconstruction errors
    if ((k == verBlockSize) && (l == horBlockSize)){
      errorVectorArray(r) = 0;
    }else{
      float relief = ComputeReflectance(normalize(normalArray[r]), xyzArray[r], overlapImgParams, globalParams);
      float recErr = ComputeError(imgVal, expRefl, albedoVal, relief);
      errorVectorArray(r) = recErr;
    }
    //printf("%d %f\n", r, errorVectorArray(r));
  }

  return true; 
}

//call function for the update of the height map. main call function for shape from shading - from multiple images
void vw::photometry::UpdateHeightMapTiles(std::string DEMTileFile,
                                          std::string albedoTileFile,
                                          std::string sfsTileFile,
                                          std::vector<ModelParams> & overlapImgParams,
                                          GlobalParams globalParams){

  float noDEMDataValue;
  if ( !readNoDEMDataVal(DEMTileFile, noDEMDataValue)){
    std::cerr << "ERROR: Could not read the NoData Value from " << DEMTileFile << std::endl;
    exit(1);
  }
  
  int numOverlapImages = overlapImgParams.size();

  // Load the weights if needed
  if (globalParams.useWeights != 0){
    for (int m = 0; m < numOverlapImages; m++){
      bool useTiles = true;
      ReadWeightsParamsFromFile(useTiles, &overlapImgParams[m]);
    }
  }
  
  DiskImageView<PixelGray<float> >  DEMTile(DEMTileFile);
  GeoReference DEM_geo;
  read_georeference(DEM_geo, DEMTileFile);

  //copy DEMTile 2 sfsDEM
  ImageView<PixelGray<float> > sfsDEM = copy(DEMTile);

  DiskImageView<PixelMask<PixelGray<uint8> > >  albedoTile(albedoTileFile);
  GeoReference albedo_geo;
  read_georeference(albedo_geo, albedoTileFile);
  
  if (DEMTile.rows() != albedoTile.rows() || DEMTile.cols() != albedoTile.cols()){
    std::cout << "ERROR: We expect the DEM and albedo tiles to have the same number "
              << "of rows/columns." << std::endl;
    exit(1);
  }

  Vector<float> lhs; lhs.set_size(numJacobianCols);
  Matrix<float> rhs; rhs.set_size(numJacobianCols, numJacobianCols);

  Matrix<float> jacobianArray;
  jacobianArray.set_size(numJacobianRows, numJacobianCols);

  Matrix<float> weightsArray;
  weightsArray.set_size(numJacobianRows, numJacobianRows);
  
  Vector<float> errorVectorArray;
  errorVectorArray.set_size(numJacobianRows);
  
  // Each block must see a 1-pixel padding on each side
  int extraPixel = 1;
  
  int numHorBlocks = (DEMTile.cols()-2*extraPixel)/horBlockSize;
  int numVerBlocks = (DEMTile.rows()-2*extraPixel)/verBlockSize;
  printf("numVerBlocks = %d, numHorBlocks = %d\n", numVerBlocks, numHorBlocks);

  vector<Vector3> normalArray;
  vector<Vector3> xyzArray;
  vector<Vector3> xyzTOPArray;
  vector<Vector3> xyzLEFTArray;
  vector<float>   reliefArray;

  normalArray.resize((verBlockSize+1)*(horBlockSize+1));
  xyzArray.resize((verBlockSize+1)*(horBlockSize+1));
  xyzTOPArray.resize((verBlockSize+1)*(horBlockSize+1));
  xyzLEFTArray.resize((verBlockSize+1)*(horBlockSize+1));

  //create image for error in terms of height and initialize to zero
  ImageView<PixelMask<PixelGray<float> > > errorHeight(DEMTile.cols(), DEMTile.rows());
  for (int k = 0; k < DEMTile.rows(); ++k){
    for (int l = 0; l < DEMTile.cols(); ++l){
      errorHeight(l, k) = 0;
    }
  }

  for (int kb = 0 ; kb < numVerBlocks; ++kb) {
    for (int lb = 0; lb < numHorBlocks; ++lb) {

      printf("kb = %d, lb=%d, numVerBlocks = %d, numHorBlocks = %d\n", kb, lb, numVerBlocks, numHorBlocks);
 
      bool success_b = ComputeBlockGeometryTiles(DEMTile, DEM_geo, noDEMDataValue,
                                                 extraPixel, kb, lb,
                                                 xyzArray, xyzLEFTArray,
                                                 xyzTOPArray, normalArray);
      
      if (!success_b) {
        printf("kb = %d, lb=%d is skipped\n", kb, lb);
        continue;
      }
      
      //initialization of the error vector
      for (int ii = 0; ii < numJacobianRows; ii++) {
        errorVectorArray(ii) = 0.0;
      }
      
      //initialization of the Jacobian matrix
      for (int ii = 0; ii < numJacobianRows; ii++){
        for (int jj = 0; jj < numJacobianCols; jj++){
          jacobianArray(ii, jj) = 0.0;
        }
      }

      // initialization of the weights matrix
      for (int ii = 0; ii < numJacobianRows; ii++) {
        for (int jj = 0; jj < numJacobianRows; jj++){
          weightsArray(ii, jj) = 0.0;
          if (ii == jj){
            weightsArray(ii, jj) = 1.0;
          }
        }
      }
      
      //reset the right hand side - square matrix of size BlockArea x BlockArea: rhs = J^T x J
      for (int ii = 0; ii < numJacobianCols; ii++) {
        for (int jj = 0; jj < numJacobianCols; jj++){
          rhs(ii, jj) = 0.0;
        }
      }
      //reset the left hand side
      for (int ii = 0; ii < numJacobianCols; ii++){
        lhs(ii) = 0.0;
      }
      
      int numValidOverlaps = 0;
      for (int m = 0; m < numOverlapImages; m++){

        printf("overlap_img = %s\n", overlapImgParams[m].inputFilename.c_str());

        DiskImageView<PixelMask<PixelGray<uint8> > >  overlapImg(overlapImgParams[m].inputFilename);
        GeoReference overlapImg_geo;
        read_georeference(overlapImg_geo, overlapImgParams[m].inputFilename);

        //determine invalid blocks in the overlap image - END
        bool success_j = ComputeBlockJacobianTiles(overlapImg, overlapImg_geo, DEM_geo, albedoTile,
                                                   extraPixel, kb, lb, overlapImgParams[m], globalParams,
                                                   xyzArray, xyzLEFTArray, xyzTOPArray, normalArray,
                                                   jacobianArray, errorVectorArray, weightsArray);
        if (!success_j) continue;
        
        //printJacobian(jacobianArray, "jacobian.txt");
        
        //compute lhs and rhs
        rhs = rhs + transpose(jacobianArray)*weightsArray*jacobianArray;
        lhs = lhs + transpose(jacobianArray)*weightsArray*errorVectorArray;
        numValidOverlaps++;
        
      }// end iterations over overlap images

      // If there are too few overlaps, we don't have enough information to update the DEM.
      if (numValidOverlaps <= 2) continue;
        
      //solves lhs = rhs*x and stores results in lhs
      try {
        solve_symmetric_nocopy(rhs,lhs);
        for (int k = 0 ; k < verBlockSize; ++k) {
          for (int l = 0; l < horBlockSize; ++l) {

            int ii = kb*verBlockSize+k+extraPixel; //row index for the entire image
            int jj = lb*horBlockSize+l+extraPixel; //col index for the entire image

            if ((ii < DEMTile.rows()) && (jj < DEMTile.cols())){
              //local index in the vector that describes the block image; assumes row-wise concatenation.
              int l_index = k*horBlockSize+l;
              
              std::cout << "sfs_before( " << jj << "," << ii << ")=" << (float)sfsDEM(jj,ii) <<  std::endl;
              sfsDEM(jj, ii) = sfsDEM(jj, ii) + lhs(l_index);
              errorHeight(jj, ii) = lhs(l_index);
              std::cout << "sfs_after( " << jj << "," << ii << ")=" << (float)sfsDEM(jj,ii) << ", lhs after= " << lhs(l_index) << std::endl;

            }
          }
        }

        printf("Go, kb = %d, lb = %d\n", kb, lb);
        
      } catch (const ArgumentErr& e) {

        std::cout << "Error @ (kb,lb) = (" << kb << "," << lb << ")\n";
        std::cout << "Exception caught: " << e.what() << "\n";
        //std::cout << "PRERHS: " << pre_rhs << "\n";
        //std::cout << "PRELHS: " << pre_lhs << "\n\n";
        //std::cout << "RHS: " << rhs << "\n";
        //std::cout << "LHS: " << lhs << "\n\n";
        //std::cout << "DEBUG: " << rhs(0,1) << "   " << rhs(1,0) << "\n\n";
        printf("Error\n");
      }

      //solve_symmetric_nocopy(rhs, lhs);

      //copy lhs to back DEMTile

    }//lb
  }//kb

  //write in the updated DEM
  std::cout << "Writing: " << sfsTileFile << std::endl;
  write_georeferenced_image(sfsTileFile, sfsDEM,
      DEM_geo, TerminalProgressCallback("photometry","Processing:"));
}

// Old code, using images not tiles

template <class ViewT, class ViewT1>
void
ComputeBlockGeometry(ImageViewBase<ViewT> const& dem, GeoReference const &demGeo,
    ImageViewBase<ViewT1> const& drg, GeoReference const &drgGeo,
    int kb, int lb, ModelParams modelParams, GlobalParams globalParams,
    vector<Vector3> &xyzArray, vector<Vector3> &xyzLEFTArray, vector<Vector3> &xyzTOPArray,
    vector<Vector3> &normalArray)
{

  //GeoTransform trans(demGeo, drgGeo);
  printf("bug fix in error vec init\n");

  //GeoTransform trans(drgGeo, demGeo);
  //transform(dem, trans);

  //cout<<"Compute Block Geometry" <<endl;

  int eVerBlockSize = verBlockSize+1;
  int eHorBlockSize = horBlockSize+1;

  for (int k = 0 ; k < eVerBlockSize; ++k) {
    for (int l = 0; l < eHorBlockSize; ++l) {

      //These statement is suspicious!!
      //int ii = kb*(eVerBlockSize)+k; //row index for the entire image
      //int jj = lb*(eHorBlockSize)+l; //col index for the entire image
      int ii = kb*(verBlockSize)+k; //row index for the entire image
      int jj = lb*(horBlockSize)+l; //col index for the entire image
      //printf("ii = %d, jj = %d, width = %d, height = %d\n", ii, jj, dem.impl().cols(), dem.impl().rows());

      if ((ii < drg.impl().rows()) && (jj < drg.impl().cols())){

        //local index in the vector that describes the block image; assumes row-wise concatenation.
        int l_index = k*eHorBlockSize+l;

        if ( is_valid(drg.impl()(jj,ii)) ) {

          Vector2 input_img_pix (jj,ii);
          Vector2 lon_lat = drgGeo.pixel_to_lonlat(input_img_pix);

          Vector2 dem_pix = demGeo.lonlat_to_pixel(drgGeo.pixel_to_lonlat(input_img_pix));
          Vector3 lonlat3(lon_lat(0),lon_lat(1),(dem.impl())(dem_pix(0), dem_pix(1)));

          //Vector3 lonlat3(lon_lat(0),lon_lat(1),(dem.impl())(jj, ii));
          //xyzArray[l_index] = drgGeo.datum().geodetic_to_cartesian(lonlat3);//3D coordinates in the img coordinates
          xyzArray[l_index] = demGeo.datum().geodetic_to_cartesian(lonlat3);//3D coordinates in the img coordinates
          /*
          //OLD and BAD!
          Vector2 input_img_left_pix;
          input_img_left_pix(0) = ii-1;
          input_img_left_pix(1) = jj;

          Vector2 input_img_top_pix;
          input_img_top_pix(0) = ii;
          input_img_top_pix(1) = jj-1;
           */
          Vector2 input_img_left_pix;
          input_img_left_pix(0) = jj-1;
          input_img_left_pix(1) = ii;

          Vector2 input_img_top_pix;
          input_img_top_pix(0) = jj;
          input_img_top_pix(1) = ii-1;

          //check for valid DEM pixel value and valid left and top coordinates
          if ((input_img_left_pix(0) >= 0) && (input_img_top_pix(1) >= 0) && (dem.impl()(/*jj,ii*/dem_pix(0), dem_pix(1)) != -10000)){

            //determine the 3D coordinates of the pixel left of the current pixel
            Vector2 lon_lat_left = drgGeo.pixel_to_lonlat(input_img_left_pix);
            Vector2 dem_left_pix = demGeo.lonlat_to_pixel(drgGeo.pixel_to_lonlat(input_img_left_pix));
            Vector3 longlat3_left(lon_lat_left(0),lon_lat_left(1),(dem.impl())(dem_left_pix(0), dem_left_pix(1)/*dem_left_pix*/));
            //Vector3 longlat3_left(lon_lat_left(0),lon_lat_left(1),(dem.impl())(input_img_left_pix(0), input_img_left_pix(1)));
            //Vector3 xyz_left = drgGeo.datum().geodetic_to_cartesian(longlat3_left);
            Vector3 xyz_left = demGeo.datum().geodetic_to_cartesian(longlat3_left);
            xyzLEFTArray[l_index] = xyz_left;

            //determine the 3D coordinates of the pixel top of the current pixel
            Vector2 lon_lat_top = drgGeo.pixel_to_lonlat(input_img_top_pix);
            Vector2 dem_top_pix = demGeo.lonlat_to_pixel(drgGeo.pixel_to_lonlat(input_img_top_pix));
            Vector3 longlat3_top(lon_lat_top(0),lon_lat_top(1),(dem.impl())(dem_top_pix(0), dem_top_pix(1)));
            //Vector3 longlat3_top(lon_lat_top(0),lon_lat_top(1),(dem.impl())(input_img_top_pix(0), input_img_top_pix(1)));
            //Vector3 xyz_top = drgGeo.datum().geodetic_to_cartesian(longlat3_top);
            Vector3 xyz_top = demGeo.datum().geodetic_to_cartesian(longlat3_top);
            xyzTOPArray[l_index] = xyz_top;

            normalArray[l_index] = cross_prod(xyz_top-xyzArray[l_index], xyz_left-xyzArray[l_index]);
            //std::cout << "xyz_top: " << xyz_top << "xyz_left: " << xyz_left << "xyzArray[l_index] " << xyzArray[l_index] << "normalArray: " << normalArray[l_index] << std::endl;
            //std::cout<<"k = "<< k <<",l = "<< l <<", normalArray: " << normalArray[l_index] << std::endl;

          }
        }
      }
    } //for l
  } //for k
}

template <class ViewT1, class ViewT2>
void
ComputeBlockJacobian(ImageViewBase<ViewT1> const& inputImage, GeoReference const &inputImageGeo,
    ImageViewBase<ViewT2> const& shadowImage, ImageViewBase<ViewT1> const& albedoImage,
    int kb, int lb, ModelParams inputImgParams, GlobalParams globalParams,
    vector<Vector3> const &xyzArray, vector<Vector3> const &xyzLEFTArray,
    vector<Vector3> const &xyzTOPArray, vector<Vector3> const &normalArray,
    Matrix<float>  &jacobianArray,
    Vector<float>  &errorVectorArray,
    Matrix<float>  &weightsArray)
{

  int r, c;

  //cout<<"Compute Block Jacobian" <<endl;
  for (r = 0; r < numJacobianRows-1; r++){//last row is always zero

    int k = r/(horBlockSize+1); //row index in the extended block
    int l = r - k*(horBlockSize+1); //col index in the extended block

    int ii = kb*verBlockSize+k; //row index for the entire image
    int jj = lb*horBlockSize+l; //col index for the entire image

    if ((ii < inputImage.impl().rows()) && (jj < inputImage.impl().cols())){

      //this is the small DRG and the large DEM. needs a fix.
      Vector2 input_img_pix(jj,ii);

      //update from the main image
      if (is_valid(inputImage.impl()(jj,ii)) && (shadowImage.impl()(jj, ii) == 0)){

        c = k*horBlockSize + l;//same point
        //not computed for the last row and last column of the extended block
        if ((k < verBlockSize) && (l < horBlockSize)){

          float recDer = ComputeReliefDerivative(xyzArray[r], xyzLEFTArray[r],
              xyzTOPArray[r], normalArray[r],
              inputImgParams, 0)
            *(float)albedoImage.impl()(jj,ii)*inputImgParams.exposureTime;
          jacobianArray(r, c) = recDer;
        }

        c = k*horBlockSize + l-1;//left point
        //not computed for the first column and last row
        if ((c >= 0) && (l > 0) && ( k < verBlockSize)){
          float recDerLEFT = ComputeReliefDerivative(xyzArray[r], xyzLEFTArray[r],
              xyzTOPArray[r], normalArray[r],
              inputImgParams, 1)
            *(float)albedoImage.impl()(jj,ii)*inputImgParams.exposureTime;
          jacobianArray(r, c) = recDerLEFT;
        }

        c = (k-1)*horBlockSize + l;//top point
        //not computed for the first row and last column of the extended block
        if ((c >= 0) && (k > 0) && (l < horBlockSize)){
          float recDerTOP = ComputeReliefDerivative(xyzArray[r], xyzLEFTArray[r],
              xyzTOPArray[r], normalArray[r], inputImgParams, 2)
            *(float)albedoImage.impl()(jj,ii)*inputImgParams.exposureTime;
          jacobianArray(r, c) = recDerTOP;
        }


        //------------------------------------------------------------------------------------------------------------------------------------------
        //compute the weights matrix
        if (globalParams.useWeights == 1){
          weightsArray(r,r) = ComputeLineWeightsHV(input_img_pix, inputImgParams);
        }
        else{
          weightsArray(r,r) = 1;
        }

        //-------------------------------------------------------------------------------------------------------------------------------------------
        //compute the reconstruction errors
        if ((k == verBlockSize) && (l == horBlockSize)){
          errorVectorArray(r) = 0;
        }
        else{

          float relief = ComputeReflectance(normalize(normalArray[r]), xyzArray[r], inputImgParams, globalParams);
          float recErr = ComputeError((float)inputImage.impl()(jj, ii), inputImgParams.exposureTime, (float)albedoImage.impl()(jj, ii), relief);
          //float recErr = ComputeReconstructError((float)inputImage.impl()(jj, ii), inputImgParams.exposureTime, (float)albedoImage.impl()(jj, ii), relief);
          errorVectorArray(r) = recErr;

        }
        //printf("%d %f\n", r, errorVectorArray(r));

      }
    }
  }
}


template <class ViewT1, class ViewT2>
void
ComputeBlockJacobianOverlap(ImageViewBase<ViewT1> const& inputImage, GeoReference const &inputImageGeo,
    ImageViewBase<ViewT1> const& overlapImage, GeoReference const &overlapImageGeo,
    ImageViewBase<ViewT2> const& shadowImage, ImageViewBase<ViewT2> const& overlapShadowImage,
    ImageViewBase<ViewT1> const& albedoImage, int kb, int lb,
    ModelParams inputImgParams,  ModelParams overlapImgParams, GlobalParams globalParams,
    vector<Vector3> const &xyzArray, vector<Vector3> const &xyzLEFTArray,
    vector<Vector3> const &xyzTOPArray, vector<Vector3> const &normalArray,
    Matrix<float> &jacobianArray,
    Vector<float> &errorVectorArray,
    Matrix<float> &weightsArray)
{

  int r, c;

  ImageViewRef<PixelMask<PixelGray<uint8> > >  interpOverlapImage = interpolate(edge_extend(overlapImage.impl(),ConstantEdgeExtension()),
      BilinearInterpolation());

  ImageViewRef<PixelMask<PixelGray<uint8> > >  interpOverlapShadowImage = interpolate(edge_extend(overlapShadowImage.impl(),
      ConstantEdgeExtension()), BilinearInterpolation());
  for (r = 0; r < numJacobianRows-1; r++){

    int k = r/(horBlockSize+1); //row index in the extended block
    int l = r - k*(horBlockSize+1); //col index in the extended block

    int ii = kb*verBlockSize+k; //row index for the entire image
    int jj = lb*horBlockSize+l; //col index for the entire image

    if ((ii < inputImage.impl().rows()) && (jj < inputImage.impl().cols())){ //is valid pixel

      Vector2 input_img_pix(jj,ii);

      //determine the corresponding pixel in the overlaping image
      Vector2 overlap_pix = overlapImageGeo.lonlat_to_pixel(inputImageGeo.pixel_to_lonlat(input_img_pix));

      float x = overlap_pix[0];
      float y = overlap_pix[1];

      //compute and update Jacobian for non shadow pixels
      if ((x>=0) && (x < overlapImage.impl().cols()) && (y>=0) && (y< interpOverlapImage.impl().rows()) && (interpOverlapShadowImage.impl()(x, y) == 0)){

        c = k*horBlockSize + l;//same point
        //not computed for the last row and last column of the extended block
        if ((k < verBlockSize) && (l < horBlockSize)){
          float recDer = ComputeReliefDerivative(xyzArray[r], xyzLEFTArray[r],
              xyzTOPArray[r], normalArray[r],
              overlapImgParams, 0)
            *(float)albedoImage.impl()(jj,ii)*overlapImgParams.exposureTime;
          jacobianArray(r,c) = recDer;
        }

        c = k*horBlockSize + l-1;//left point
        //not computed for the first column and last row
        if ((c >= 0) && (l > 0) && ( k < verBlockSize)){
          float recDerLEFT = ComputeReliefDerivative(xyzArray[r], xyzLEFTArray[r],
              xyzTOPArray[r], normalArray[r],
              overlapImgParams, 1)
            *(float)albedoImage.impl()(jj,ii)*overlapImgParams.exposureTime;
          jacobianArray(r,c) = recDerLEFT;
        }

        c = (k-1)*horBlockSize + l;//top point
        //not computed for the first row and last column of the extended block
        if ((c >= 0) && (k > 0) && (l < horBlockSize)){
          float recDerTOP = ComputeReliefDerivative(xyzArray[r], xyzLEFTArray[r],
              xyzTOPArray[r], normalArray[r], overlapImgParams, 2)
            *(float)albedoImage.impl()(jj,ii)*overlapImgParams.exposureTime;
          jacobianArray(r,c) = recDerTOP;

        }

        //-----------------------------------------------------------------------------------------------------------------
        //compute the weights matrix
        if (globalParams.useWeights == 1){
          weightsArray(r,r) = ComputeLineWeightsHV(overlap_pix, overlapImgParams);
        }
        else{
          weightsArray(r,r) = 1;
        }
        //------------------------------------------------------------------------------------------------------------------
        //compute the error vector
        //suspicious statement
        //if ((c == verBlockSize) && (l == horBlockSize)){

        //error vector is computed only for points within the block and not the extended block
        if ((k == verBlockSize) && (l == horBlockSize)){ //last column and last row of the extended block
          errorVectorArray(r) = 0;
        }
        else{
          float relief = ComputeReflectance(normalize(normalArray[r]), xyzArray[r], overlapImgParams, globalParams);
          float recErr = ComputeError((float)interpOverlapImage.impl()(x, y), overlapImgParams.exposureTime, (float)albedoImage.impl()(jj, ii), relief);
          //float recErr = ComputeReconstructError((float)interpOverlapImage.impl()(x, y), overlapImgParams.exposureTime, (float)albedoImage.impl()(jj, ii), relief);
          errorVectorArray(r) = recErr;
        }
      }
    }
  }
}




//call function for the update of the height map. main call function for shape from shading - from multiple images
void vw::photometry::UpdateHeightMap(ModelParams inputImgParams, std::vector<ModelParams> overlapImgParams, GlobalParams globalParams)
{

  std::string inputImgFilename = inputImgParams.inputFilename;//the original DRG
  std::string shadowFilename = inputImgParams.shadowFilename; //shadow map
  std::string albedoFilename = inputImgParams.outputFilename; //albedo map
  std::string meanDEMFilename = inputImgParams.meanDEMFilename; //original mean DEM
  std::string sfsDEMFilename = inputImgParams.sfsDEMFilename; //DEM from sfs
  std::string errorHeightFilename = inputImgParams.errorHeightFilename; //error in terms of height

  int numOverlapImages = overlapImgParams.size();
  //numOverlapImages = 0;//1;

  DiskImageView<PixelGray<float> >  meanDEM(meanDEMFilename);
  GeoReference DEM_geo;
  read_georeference(DEM_geo, meanDEMFilename);

  /*
  //upsample the meanDEM;
  int upsampleFactor = 2;
  upsample_image(sfsDEMFilename, meanDEMFilename, upsampleFactor);
  DiskImageView<PixelGray<float> >  meanDEM(upDEMFilename);
  GeoReference DEM_geo;
  read_georeference(DEM_geo, meanDEMFilename);
  */

  ImageView<PixelGray<float> > sfsDEM(meanDEM.cols(), meanDEM.rows());

  //copy meanDEM 2 sfsDEM
  sfsDEM = copy(meanDEM);

  DiskImageView<PixelMask<PixelGray<uint8> > >  shadowImage(shadowFilename);

  DiskImageView<PixelMask<PixelGray<uint8> > >  inputImage(inputImgFilename);
  GeoReference inputImg_geo;
  read_georeference(inputImg_geo, inputImgFilename);

  DiskImageView<PixelMask<PixelGray<uint8> > >  albedoImage(albedoFilename);

  ImageViewRef<PixelGray<float> >  interp_dem_image = interpolate(edge_extend(meanDEM.impl(),
        ConstantEdgeExtension()),
      BilinearInterpolation());

  //Vector<float, numJacobianCols> lhs;
  //Matrix<float, numJacobianCols, numJacobianCols> rhs;
  Vector<float> lhs; lhs.set_size(numJacobianCols);
  Matrix<float> rhs; rhs.set_size(numJacobianCols, numJacobianCols);

  // vector<Matrix<float, numJacobianRows, numJacobianCols> >jacobianArray; // old static allocation
  vector<Matrix<float> >jacobianArray;
  jacobianArray.resize(numOverlapImages+1);
  for (int m = 0; m < (int)jacobianArray.size(); m++){
    jacobianArray[m].set_size(numJacobianRows, numJacobianCols);
  }
    
  // vector<Matrix<float, numJacobianRows, numJacobianRows> >weightsArray; // old static allocation
  vector<Matrix<float> >weightsArray;
  weightsArray.resize(numOverlapImages+1);
  for (int m = 0; m < (int)weightsArray.size(); m++){
    weightsArray[m].set_size(numJacobianRows, numJacobianRows);
  }

  //  vector<Vector<float, numJacobianRows> > errorVectorArray; // old static allocation
  vector<Vector<float> > errorVectorArray;
  errorVectorArray.resize(numOverlapImages+1);
  for (int m = 0; m < (int)errorVectorArray.size(); m++){
    errorVectorArray[m].set_size(numJacobianRows);
  }
  
  float recDer, recDerLEFT, recDerTOP, recErr;

  int numHorBlocks = meanDEM.cols()/horBlockSize + 1;
  int numVerBlocks = meanDEM.rows()/verBlockSize + 1;
  printf("numVerBlocks = %d, numHorBlocks = %d\n", numVerBlocks, numHorBlocks);

  vector<Vector3> normalArray;
  vector<Vector3> xyzArray;
  vector<Vector3> xyzTOPArray;
  vector<Vector3> xyzLEFTArray;
  vector<float> reliefArray;

  normalArray.resize((verBlockSize+1)*(horBlockSize+1));
  xyzArray.resize((verBlockSize+1)*(horBlockSize+1));
  xyzTOPArray.resize((verBlockSize+1)*(horBlockSize+1));
  xyzLEFTArray.resize((verBlockSize+1)*(horBlockSize+1));

  //josh - create image for error in terms of height and initialize to zero
  ImageView<PixelMask<PixelGray<float> > > errorHeight(meanDEM.cols(), meanDEM.rows());
  for (int k = 0; k < meanDEM.rows(); ++k){
    for (int l = 0; l < meanDEM.cols(); ++l){
      errorHeight(l, k) = 0;
    }
  }

  for (int kb = 0 ; kb < numVerBlocks; ++kb) {
    for (int lb = 0; lb < numHorBlocks; ++lb) {

      printf("kb = %d, lb=%d, numVerBlocks = %d, numHorBlocks = %d\n", kb, lb, numVerBlocks, numHorBlocks);

      //josh - shouldn't we check the extended image for valid points?
      //determine invalid blocks in the input image - START
      int n = 0;
      for (int k = 0 ; k < verBlockSize; ++k){
        for (int l = 0; l < horBlockSize; ++l) {
          int ii = kb*verBlockSize+k; //row index for the entire image
          int jj = lb*horBlockSize+l; //col index for the entire image
          if ((ii < inputImage.rows()) && (jj < inputImage.cols())){
            if ( is_valid(inputImage(jj,ii)) ){
              n++;
            }
          }
        }
      }

      if ( n < numJacobianCols ) {
        printf("kb = %d, lb=%d is skipped\n", kb, lb);
        continue;
      }
      //determine invalid blocks in the input image - END

      //initialization of the error vector
      for (int m = 0; m < numOverlapImages + 1; m++){
        for (int ii = 0; ii < numJacobianRows; ii++) {
          errorVectorArray[m](ii) = 0.0;
        }
      }

      //initialization of the Jacobian matrix
      for (int m = 0; m < numOverlapImages + 1; m++){
        for (int ii = 0; ii < numJacobianRows; ii++){
          for (int jj = 0; jj < numJacobianCols; jj++){
            jacobianArray[m](ii, jj) = 0.0;
          }
        }
      }

      //initialization of the weights matrix
      for (int m = 0; m < numOverlapImages + 1; m++){
        for (int ii = 0; ii < numJacobianRows; ii++) {
          for (int jj = 0; jj < numJacobianRows; jj++){
            weightsArray[m](ii, jj) = 0.0;
            if (ii == jj){
              weightsArray[m](ii, jj) = 1.0;
            }
          }
        }
      }

     //TO DO: reset xyzArray, xyzLEFT and xyzTOP
      ComputeBlockGeometry(interp_dem_image, DEM_geo,
          inputImage, inputImg_geo, kb, lb,
          inputImgParams, globalParams,
          xyzArray, xyzLEFTArray,
          xyzTOPArray, normalArray);

      ComputeBlockJacobian(inputImage, inputImg_geo, shadowImage, albedoImage,
          kb, lb, inputImgParams, globalParams,
          xyzArray, xyzLEFTArray, xyzTOPArray,normalArray,
          jacobianArray[0], errorVectorArray[0], weightsArray[0]);

      //printJacobian(jacobianArray[0], "mujapenas.txt");

      for (int m = 0; m < numOverlapImages; m++){

        printf("overlap_img = %s\n", overlapImgParams[m].inputFilename.c_str());

        DiskImageView<PixelMask<PixelGray<uint8> > >  overlapImg(overlapImgParams[m].inputFilename);
        GeoReference overlapImg_geo;
        read_georeference(overlapImg_geo, overlapImgParams[m].inputFilename);
        DiskImageView<PixelMask<PixelGray<uint8> > >  overlapShadowImage(overlapImgParams[m].shadowFilename);

        //GeoTransform trans(overlapImg_geo, inputImg_geo);
        //transform(overlapImg, trans);

        //determine invalid blocks in the overlap image - START
        int n = 0;

        for (int k = 0 ; k < verBlockSize; ++k){
          for (int l = 0; l < horBlockSize; ++l) {
            int ii = kb*verBlockSize+k; //row index for the entire image
            int jj = lb*horBlockSize+l; //col index for the entire image
            Vector2 input_img_pix(jj, ii);
            Vector2 overlap_pix = overlapImg_geo.lonlat_to_pixel(inputImg_geo.pixel_to_lonlat(input_img_pix));
            float x = overlap_pix[0];
            float y = overlap_pix[1];

            if ((x>=0) && (x < overlapImg.cols()) && (y>=0) && (y< overlapImg.rows()) && (overlapShadowImage(x, y) == 0)){
              if ( is_valid(overlapImg(x,y)) ){
                n++;
              }
            }
          }
        }

        if ( n < numJacobianCols ) {
          printf("OVERLAP kb = %d, lb=%d is skipped\n", kb, lb);
          continue;
        }

        //determine invalid blocks in the overlap image - END

        ComputeBlockJacobianOverlap(inputImage, inputImg_geo,
            overlapImg, overlapImg_geo,
            shadowImage, overlapShadowImage,
            albedoImage, kb, lb,
            inputImgParams, overlapImgParams[m], globalParams,
            xyzArray, xyzLEFTArray, xyzTOPArray,normalArray,
            jacobianArray[m+1], errorVectorArray[m+1], weightsArray[m+1]);
      }// end iterations over overlap images
      
      //reset the right hand side - square matrix of size BlockArea x BlockArea: rhs = J^T x J
      for (int ii = 0; ii < numJacobianCols; ii++) {
        for (int jj = 0; jj < numJacobianCols; jj++){
          rhs(ii, jj) = 0.0;
        }
      }
      //reset the left hand side
      for (int ii = 0; ii < numJacobianCols; ii++){
        lhs(ii) = 0.0;
      }
      //compute lhs and rhs
      for (int m = 0; m < numOverlapImages+1; m++){
        rhs = rhs + transpose(jacobianArray[m])*weightsArray[m]*jacobianArray[m];
        lhs = lhs + transpose(jacobianArray[m])*weightsArray[m]*errorVectorArray[m];
      }
      //solves lhs = rhs*x and stores results in lhs
      try {
        solve_symmetric_nocopy(rhs,lhs);
        for (int k = 0 ; k < verBlockSize; ++k) {
          for (int l = 0; l < horBlockSize; ++l) {

            int ii = kb*verBlockSize+k; //row index for the entire image
            int jj = lb*horBlockSize+l; //col index for the entire image

            if ((ii < meanDEM.rows()) && (jj < meanDEM.cols())){
              //local index in the vector that describes the block image; assumes row-wise concatenation.
              int l_index = k*horBlockSize+l;
              
              std::cout << "sfs_before( " << jj << "," << ii << ")=" << (float)sfsDEM(jj,ii) <<  std::endl;
              sfsDEM(jj, ii) = sfsDEM(jj, ii) + lhs(l_index);
              errorHeight(jj, ii) = lhs(l_index);
              std::cout << "sfs_after( " << jj << "," << ii << ")=" << (float)sfsDEM(jj,ii) << ", lhs after= " << lhs(l_index) << std::endl;

            }
          }
        }

        printf("Go, kb = %d, lb = %d\n", kb, lb);

      } catch (const ArgumentErr& e) {

        std::cout << "Error @ (kb,lb) = (" << kb << "," << lb << ")\n";
        std::cout << "Exception caught: " << e.what() << "\n";
        //std::cout << "PRERHS: " << pre_rhs << "\n";
        //std::cout << "PRELHS: " << pre_lhs << "\n\n";
        //std::cout << "RHS: " << rhs << "\n";
        //std::cout << "LHS: " << lhs << "\n\n";
        //std::cout << "DEBUG: " << rhs(0,1) << "   " << rhs(1,0) << "\n\n";
        printf("Error\n");
      }

      //solve_symmetric_nocopy(rhs, lhs);

      //copy lhs to back meanDEM

    }//lb
  }//kb

  //write in the updated DEM
  std::cout << "Writing: " << sfsDEMFilename << std::endl;
  write_georeferenced_image(sfsDEMFilename, sfsDEM,
      DEM_geo, TerminalProgressCallback("photometry","Processing:"));
  //write in the error in terms of height
  //write_georeferenced_image(errorHeightFilename, errorHeight,
  //    DEM_geo, TerminalProgressCallback("photometry","Processing:"));
}
