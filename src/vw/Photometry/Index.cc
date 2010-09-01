// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
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
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/normal.hpp> // for normal_distribution
#include <boost/math/distributions/gamma.hpp>

#include <vw/Core.h>
#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Cartography.h>
#include <vw/Math.h>
using namespace vw;
using namespace vw::math;
using namespace vw::cartography;
#include <vw/Photometry/Reconstruct.h>
#include <vw/Photometry/Index.h>

template <class ViewT>
std::vector<float> build_histogram(ImageViewBase<ViewT> const& view) {
        typedef typename ViewT::pixel_accessor pixel_accessor;
        std::vector<float> result(DYNAMIC_RANGE);

        int num_valid = 0;
        pixel_accessor row_acc = view.impl().origin();
        for( int32 row=0; row < view.impl().rows(); ++row ) {
                pixel_accessor col_acc = row_acc;
                for( int32 col=0; col < view.impl().cols(); ++col ) {
                        if ( is_valid(*col_acc) ) {
                                result[ (*col_acc)[0] ] += 1.0;
                                ++num_valid;
                        }
                        col_acc.next_col();
                }
                row_acc.next_row();
        }

        // Renormalize the histogram values
        for (size_t i=0; i< result.size(); ++i)
                result[i] /= num_valid;

        return result;
}

