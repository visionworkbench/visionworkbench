// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file ControlNetwork.cc
///

#include <vw/Photometry/Reconstruct.h>

// Generic Ostream options for Debugging
std::ostream& operator<<( std::ostream& os, GlobalParams const& global ) {
  os << "-- Global Params --\n";
  os << " ReflectanceType: " << global.reflectanceType
     << " SlopeType: " << global.slopeType << "\n";
  os << " Shadow threshold: " << global.shadowThresh << "\n";

  os << " ExposureInfoFilename : " << global.exposureInfoFilename << "\n";
  os << " SpacecraftPosFilename : " << global.spacecraftPosFilename << "\n";
  os << " SunPosFilename : " << global.sunPosFilename << "\n";

  return os;
}

std::ostream& operator<<( std::ostream& os, ModelParams const& model ) {
  os << "-- Model Params --\n";
  os << " Exposure Time: " << model.exposureTime << "\n";
  os << " Sun Position : " << model.sunPosition << "\n";
  os << " Spacecraft Position : " << model.spacecraftPosition << "\n";
  os << " Info File : " << model.infoFilename << "\n";
  os << " DEMFilename : " << model.DEMFilename << "\n";
  os << " meanDEMFilename : " << model.meanDEMFilename << "\n";
  os << " var2DEMFilename : " << model.var2DEMFilename << "\n";
  os << " reliefFilename  : " << model.reliefFilename << "\n";
  os << " shadowFilename  : " << model.shadowFilename << "\n";
  os << " errorFilename   : " << model.errorFilename << "\n";
  os << " inputFilename   : " << model.inputFilename << "\n";
  os << " outputFilename  : " << model.outputFilename << "\n";

  return os;
}
