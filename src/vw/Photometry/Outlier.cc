// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <vw/Image.h>
#include <vw/Image/PixelMath.h>
#include <vw/Image/PixelMask.h>
#include <vw/Image/MaskViews.h>
#include <vw/FileIO.h>
#include <math.h>
#include <time.h>

//using namespace std;
//using namespace vw;
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/convenience.hpp>
namespace fs=boost::filesystem;

#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions/normal.hpp> // for normal_distribution
#include <boost/math/distributions/chi_squared.hpp>

#include <vw/Core.h>
#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Cartography.h>
#include <vw/Math.h>
using namespace std;
using namespace vw;
using namespace vw::math;
using namespace vw::cartography;
#include <vw/Photometry/Reconstruct.h>
#include <vw/Photometry/Outlier.h>
#include <vw/Photometry/Exposure.h>
#include <vw/Photometry/Misc.h>

