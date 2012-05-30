// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file Albedo.h

#ifndef __VW_PHOTOMETRY_ALBEDO_H__
#define __VW_PHOTOMETRY_ALBEDO_H__

#include <string>
#include <vector>
#include <vw/Photometry/Reconstruct.h>

namespace vw {
namespace photometry {

  //image mosaic functions
  void InitImageMosaic(ModelParams input_img_params,
                       std::vector<ModelParams> overlap_img_params,
                       GlobalParams globalParams);

  void InitImageMosaicByBlocks(ModelParams input_img_params,
                               std::vector<ModelParams> overlap_img_params,
                               GlobalParams globalParams);

  void UpdateImageMosaic(ModelParams input_img_params,
                         std::vector<ModelParams> overlap_img_params,
                         GlobalParams globalParams);

  struct phaseCoeffsData{
    double phaseCoeffA1_num, phaseCoeffA1_den;
    double phaseCoeffA2_num, phaseCoeffA2_den;
    phaseCoeffsData(){
      phaseCoeffA1_num = phaseCoeffA1_den = phaseCoeffA2_num = phaseCoeffA2_den = 0.0;
    }
    void writeToFile(std::string fileName){
      std::cout << "Writing to: " << fileName << std::endl;
      std::ofstream fh(fileName.c_str());
      fh.precision(16);
      fh << phaseCoeffA1_num << ' ' << phaseCoeffA1_den << ' '
         << phaseCoeffA2_num << ' ' << phaseCoeffA2_den << std::endl;
      fh.close();
    }
    void readFromFile(std::string fileName){
      std::ifstream fh(fileName.c_str());
      if (!fh || !( fh >> phaseCoeffA1_num >> phaseCoeffA1_den >>
                    phaseCoeffA2_num >> phaseCoeffA2_den) ){
        std::cerr << "Could not read phase data from file: " << fileName << std::endl;
        exit(1);
      }
      fh.close();
    }
  };
  
  double actOnTile(bool isLastIter, bool computeErrors,
                   std::string blankTileFile,
                   std::string DEMTileFile,
                   std::string albedoTileFile,
                   std::string errorTileFile, std::string weightsSumFile, 
                   std::vector<ModelParams> & overlap_img_params,
                   GlobalParams globalParams,
                   phaseCoeffsData & PCD
                   );
  
  void AppendCostFunToFile(double costFunVal, std::string fileName);
  
  //albedo mosaic functions
  void InitAlbedoMosaic(ModelParams input_img_params,
                        std::vector<ModelParams> overlap_img_params,
                        GlobalParams globalParams);

  //albedo mosaic functions
  void UpdateAlbedoMosaic(ModelParams input_img_params,
                          std::vector<ModelParams> overlap_img_params,
                          GlobalParams globalParams);

}} // end vw::photometry

#endif//__VW_PHOTOMETRY_ALBEDO_H__
