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

// Does conjugate gradient descent on the DEM, keeping all else fixed.
void optimize_conjugate_gradient(ImageView<PixelGray<double> > *image_predicted,  ImageView<PixelGray<double> > *image, ImageView<PixelGray<double> > *dem, ImageView<PixelGray<double> > *init_dem, ImageView<PixelGray<double> > *albedo, Vector<double, 3> *light_direction);
