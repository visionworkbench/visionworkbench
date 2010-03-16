// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Math/Vector.h>
#include <string>

int* ComputeImageCenterLine(std::string input_img_file, int **r_maxDistArray);
int* ComputeDEMCenterLine(std::string input_DEM_file, int **r_maxDistArray);
int* ComputeImageHorCenterLine(std::string input_img_file, int **r_maxDistArray);
int* ComputeDEMHorCenterLine(std::string input_DEM_file, int **r_maxDistArray);
float ComputeLineWeights(Vector2 pix, int *centerLine, int *maxDistArray);
float ComputeLineWeightsV(Vector2 pix, int *centerLine, int *maxDistArray);
