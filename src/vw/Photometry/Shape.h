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



void InitDEM( std::string input_DEM_file, std::string mean_DEM_file, std::string var2_DEM_file, modelParams input_img_params,
              std::vector<std::string> overlap_DEM_files,  std::vector<modelParams> overlap_img_params, GlobalParams globalParams);

void ComputeSaveDEM(std::string curr_input_file, std::string prev_input_file,
                    std::string curr_output_file,  std::string DEM_file,
                    modelParams currModelParams, modelParams prevModelParams);

