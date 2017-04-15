// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NASA Vision Workbench is licensed under the Apache License,
//  Version 2.0 (the "License"); you may not use this file except in
//  compliance with the License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__

#ifndef __VW_MODIS_WATER_DETECTION_H__
#define __VW_MODIS_WATER_DETECTION_H__

#include <vw/tools/modis_utilities.h>

/**
  Higher level MODIS water detection functions.

  Development of MODIS code was stopped early and this file is incomplete.
  All MODIS algorithms were originally developed using Google Earth Engine as part
  of the Crisis Mapping project.
*/

namespace vw {
namespace modis {

typedef ImageView<PixelMask<uint8> > WaterDetection;

/// Overload this function for this particular case
template <class FuncT>
void for_each_pixel(ModisImage const& modis, ModisProductImage const& modis_products,
                    WaterDetection &result, FuncT const& functor) {
  // Apply the functor to each pixel       
  result.set_size(modis.cols(), modis.rows());
  for (int r=0; r<modis.rows(); ++r) {
    for (int c=0; c<modis.cols(); ++c) {
      ModisPixelType        pixel    = modis(c,r);
      ModisProductPixelType products = modis_products(c,r);
      // TODO: Mask handling!
      bool value = functor(pixel, products);
      if (value)
        result(c,r) = 255;
      else
        result(c,r) = 0;
    }
  }
}


//==============================================================

/// Simple EVI based classifier
struct detect_water_evi_functor {
  bool operator()( ModisPixelType const& pixel, ModisProductPixelType const& products) const {
    bool criteria1 = (products[EVI] <= 0.3 ) && (products[LSWI]-products[EVI] >= 0.05);
    bool criteria2 = (products[EVI] <= 0.05) && (products[LSWI] <= 0.0);
    return (criteria1 || criteria2);
  }
};


/// Method from paper: Xiao, Boles, Frolking, et. al. Mapping paddy rice agriculture in South and Southeast Asia using
///                    multi-temporal MODIS images, Remote Sensing of Environment, 2006.
/// - This method implements a very simple decision tree from several standard MODIS data products.
/// - The default constants were tuned for (wet) rice paddy detection.
struct detect_water_xiao_functor {
  bool operator()( ModisPixelType const& pixel, ModisProductPixelType const& products) const {
    return (products[LSWI] - products[NDVI] >=0.05) || (products[LSWI] - products[EVI] >= 0.05);
  }
};


/// Compute (b2-b1) < threshold, a simple water detection index.
/// - This method may be all that is needed in cases where the threshold can be hand tuned.
struct detect_water_diff_functor {
  detect_water_diff_functor(double threshold) : m_threshold(threshold) {}
  float m_threshold;
  bool operator()( ModisPixelType const& pixel, ModisProductPixelType const& products) const {
    return (pixel[B2] - pixel[B1]) <= m_threshold;
  }
};


/// A flood detection method from the Dartmouth Flood Observatory.
/// - This method is a refinement of the simple b2-b1 detection method.
struct detect_water_dartmouth_functor {
  detect_water_dartmouth_functor(double threshold) : m_threshold(threshold) {}
  float m_threshold;
  bool operator()( ModisPixelType const& pixel, ModisProductPixelType const& products) const {
    float A = 500;
    float B = 2500;
    return (pixel[B2]+A)/(pixel[B1]+B) <= m_threshold;
  }
};


struct detect_water_mod_ndwi_functor {
  detect_water_mod_ndwi_functor(double threshold) : m_threshold(threshold) {}
  float m_threshold;
  bool operator()( ModisPixelType const& pixel, ModisProductPixelType const& products) const {
    if (pixel[B4] + pixel[B6] == 0)
      return false;
    return (pixel[B6] - pixel[B4]) / (pixel[B4] + pixel[B6]) <= m_threshold;
  }
};


/// Floating Algae Index. Method from paper: Feng, Hu, Chen, Cai, Tian, Gan,
/// Assessment of inundation changes of Poyang Lake using MODIS observations
/// between 2000 and 2010. Remote Sensing of Environment, 2012.
struct detect_water_fai_functor {
  detect_water_fai_functor(double threshold) : m_threshold(threshold) {}
  float m_threshold;
  bool operator()( ModisPixelType const& pixel, ModisProductPixelType const& products) const {
    float constant = (859.0 - 645) / (1240 - 645);
    return (pixel[B2] - (pixel[B1] + constant*(pixel[B5] - pixel[B1]))) <= m_threshold;
  }
};


}} // end namespace vw::modis

#endif

