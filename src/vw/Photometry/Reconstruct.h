// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file Reconstruct.h

#ifndef __VW_PHOTOMETRY_RECONSTRUCT_H__
#define __VW_PHOTOMETRY_RECONSTRUCT_H__

#define NO_REFL 0
#define LAMBERT 1
#define LUNAR_LAMBERT 2

#include <iostream>
#include <string>
#include <vw/Math/Vector.h>

struct ImageRecord {
  static const double defaultCoord = -999.0;
  bool useImage;
  std::string path;
  double north, west, south, east;
  
  ImageRecord(void) {}
  ImageRecord(const std::string& path_) :
    useImage(true),
    path(path_),
    north(defaultCoord)
  {}
};
 

namespace vw {
namespace photometry {

  typedef struct GlobalParams{

    char drgDir[1000];
    char demDir[1000];
    char sunPosFile[1000];
    char spacecraftPosFile[1000];

    int initialSetup;
    float tileSize;        // in degrees
    int useTiles;          // 1 or 0 
    int pixelPadding;      // in pixels
    Vector4 simulationBox; // lonMin, lonMax, latMin, latMax (If not present the entire albedo will be simulated)
    
    int reflectanceType;
    int saveReflectance;
    int slopeType;
    int exposureInitType;
    int albedoInitType;
    int DEMInitType;
    int shadowInitType;

    float shadowThresh;

    std::string exposureInfoFilename;

    //float exposureInitRefValue;//this will be removed
    //int exposureInitRefIndex;//this will be removed
    float TRConst;
    int updateAlbedo, updateExposure, updateHeight;
    int useWeights;
    int saveWeights;
    int maxNumIter;
    int computeErrors;
    int noDEMDataValue;
    
  };

  typedef struct ModelParams {
    float   exposureTime;

    Vector2 cameraParams; //currently not used
    Vector3 sunPosition; //relative to the center of the Moon
    Vector3 spacecraftPosition;//relative to the center of the planet

    int *hCenterLine;
    int *hMaxDistArray;
    int *hCenterLineDEM;
    int *hMaxDistArrayDEM;

    int *vCenterLine;
    int *vMaxDistArray;
    int *vCenterLineDEM;
    int *vMaxDistArrayDEM;
    
    Vector4 corners; // cached bounds to quickly calculate overlap
    
    /*
    vector<int> hCenterLine;
    vector<int> hMaxDistArray;
    vector<int> hCenterLineDEM;
    vector<int> hMaxDistArrayDEM;
    vector<int> vCenterLine;
    vector<int> vMaxDistArray;
    vector<int> vCenterLineDEM;
    vector<int> vMaxDistArrayDEM;
    */
    std::string infoFilename, DEMFilename, meanDEMFilename,
      var2DEMFilename, reliefFilename, shadowFilename,
      errorFilename, inputFilename, outputFilename, 
      sfsDEMFilename, errorHeightFilename, weightFilename, exposureFilename;

  };

  // Generic Ostream options for Debugging
  std::ostream& operator<<( std::ostream& os, GlobalParams const& global );
  std::ostream& operator<<( std::ostream& os, ModelParams const& model );

  // Written by Taemin Kim - START
  #define DYNAMIC_RANGE  256
  // Written by Taemin Kim - END
}} // end vw::photometry

#endif//__VW_PHOTOMETRY_RECONSTRUCT_H__
