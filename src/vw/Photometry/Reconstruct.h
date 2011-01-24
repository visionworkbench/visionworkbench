// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
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

namespace vw {
namespace photometry {

  typedef struct GlobalParams{
    int reflectanceType;
    int slopeType;
    int exposureInitType;
    int albedoInitType;
    int DEMInitType;
    int shadowInitType;

    float shadowThresh;

    std::string exposureInfoFilename,
      spacecraftPosFilename, sunPosFilename;

    //float exposureInitRefValue;//this will be removed
    //int exposureInitRefIndex;//this will be removed
    float TRConst;
    int updateAlbedo, updateExposure, updateHeight;
    int useWeights;
    int maxNumIter;
    int computeErrors;
    int noDEMDataValue;
  };

  typedef struct ModelParams {
    float   exposureTime;

    Vector2 cameraParams; //currently not used
    Vector3 sunPosition; //relative to the center of the Moon
    Vector3 spacecraftPosition;//relative to the center of the planet
 
    int *centerLine;
    int *maxDistArray;
    int *centerLineDEM;
    int *maxDistArrayDEM;
    int *horCenterLine;
    int *maxVerDistArray;
    int *horCenterLineDEM;
    int *maxVerDistArrayDEM;

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
