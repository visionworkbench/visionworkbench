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

bool 
ComputeBlockGeometry(ImageView<PixelGray<float> > const& dem, GeoReference const &demGeo,
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

template <class ViewT1>
bool
ComputeBlockJacobian(ImageViewBase<ViewT1> const& overlapImage, GeoReference const &overlapImageGeo,
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
void vw::photometry::UpdateHeightMap(std::string DEMTileFile,
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
 
      bool success_b = ComputeBlockGeometry(DEMTile, DEM_geo, noDEMDataValue,
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
        bool success_j = ComputeBlockJacobian(overlapImg, overlapImg_geo, DEM_geo, albedoTile,
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
