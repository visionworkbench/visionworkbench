#include <iostream>
#include <fstream>
#include <vw/Image.h>
#include <vw/Image/PixelMath.h>
#include <vw/Image/PixelMask.h>
#include <vw/Image/MaskViews.h>
#include <vw/FileIO.h>
#include <math.h>
#include <time.h>

using namespace std;
using namespace vw;

#include <vw/Photometry/Reconstruct.h>

//image mosaic functions
void InitImageMosaic(std::string input_img_file,
                     modelParams input_img_params,
                     std::string shadow_file,
                     std::string output_img_file,
                     std::vector<std::string> overlap_img_files,
                     std::vector<modelParams> overlap_img_params,
                     GlobalParams globalParams);

void InitImageMosaicByBlocks(modelParams input_img_params,
                             std::vector<modelParams> overlap_img_params,
                             GlobalParams globalParams);

void UpdateImageMosaic(modelParams input_img_params, std::vector<modelParams> overlap_img_params,
                       GlobalParams globalParams);

//albedo mosaic functions

void InitAlbedoMosaic(modelParams input_img_params,
                      std::vector<modelParams> overlap_img_params,
                      GlobalParams globalParams);

void UpdateAlbedoMosaic(modelParams input_img_params,
                         std::vector<modelParams> overlap_img_params,
                         GlobalParams globalParams);

void ComputeReconstructionErrorMap(modelParams input_img_params,
                                   std::vector<modelParams> overlap_img_params,
                                   GlobalParams globalParams,
                                   float *avgError, int *totalNumSamples);







