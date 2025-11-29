// __BEGIN_LICENSE__
//  Copyright (c) 2006-2024, United States Government as represented by the
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

/// \file InterestData.cc
///
/// Basic classes and structures for storing image interest points.
///
#include <vw/InterestPoint/InterestData.h>
#include <fstream>
namespace vw {
namespace ip {

#if defined(VW_HAVE_PKG_OPENCV) && VW_HAVE_PKG_OPENCV == 1

    /// Copy IP information from an OpenCV KeyPoint object.
    void InterestPoint::setFromCvKeypoint(cv::KeyPoint const& cvKey) {
      x  = cvKey.pt.x;
      y  = cvKey.pt.y;
      ix = round(x);
      iy = round(y);
      interest    = cvKey.response;
      octave      = cvKey.octave;
      scale_lvl   = cvKey.octave;
      scale       = cvKey.size;
      orientation = cvKey.angle;
      polarity    = false;
    }

    /// Create an OpenCV KeyPoint object from this IP.
    cv::KeyPoint InterestPoint::makeOpenCvKeypoint() const {
      cv::KeyPoint cvKey;
      cvKey.pt.x     = x;
      cvKey.pt.y     = y;
      cvKey.response = interest;
      cvKey.octave   = octave;
      cvKey.size     = scale;
      cvKey.angle    = orientation;
      return cvKey;
    }

#endif

bool InterestPointLessThan (InterestPoint P1, InterestPoint P2){
  if (P1.x < P2.x) 
    return true; 
  if (P1.x > P2.x)
    return false;
  if (P1.y < P2.y)
    return true; 
  if (P1.y > P2.y)
    return false;
  if (P1.scale < P2.scale) 
    return true; 
  if (P1.scale > P2.scale)
    return false;
  if (P1.orientation < P2.orientation)
    return true; 
  if (P1.orientation > P2.orientation)
    return false;
  if (P1.interest < P2.interest)
    return true; 
  if (P1.interest > P2.interest)
    return false;
  if (P1.polarity < P2.polarity)
    return true; 
  if (P1.polarity > P2.polarity)
    return false;
  if (P1.octave < P2.octave)
    return true; 
  if (P1.octave > P2.octave)
    return false;
  if (P1.scale_lvl < P2.scale_lvl)
    return true; 
  if (P1.scale_lvl > P2.scale_lvl)
    return false;
    
  return false;
}

std::vector<Vector3> iplist_to_vectorlist(std::vector<InterestPoint> const& iplist) {
  std::vector<Vector3> result(iplist.size());
  for (size_t i=0; i < iplist.size(); ++i) {
    result[i][0] = iplist[i].x;
    result[i][1] = iplist[i].y;
    result[i][2] = 1; // homogeneous vector
  }
  return result;
}

/// Helpful functors
void remove_descriptor(InterestPoint & ip) {
  ip.descriptor = ip::InterestPoint::descriptor_type(); // this should free up the memory
}

}} // namespace vw::ip
