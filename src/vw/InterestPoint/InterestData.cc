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

  std::vector<Vector3> iplist_to_vectorlist(std::vector<InterestPoint> const& iplist) {
    std::vector<Vector3> result(iplist.size());
    for (size_t i=0; i < iplist.size(); ++i) {
      result[i][0] = iplist[i].x;
      result[i][1] = iplist[i].y;
      result[i][2] = 1; // homogeneous vector
    }
    return result;
  }
/*
  // You are highly discouraged in using this as all descriptor information is lost.
  std::vector<InterestPoint> vectorlist_to_iplist(std::vector<Vector3> const& veclist) {
    std::vector<InterestPoint> result(veclist.size());
    for (size_t i=0; i < veclist.size(); ++i) {
      result[i].x = veclist[i][0];
      result[i].y = veclist[i][1];
    }
    return result;
  }
*/
  /// Helpful functors
  void remove_descriptor(InterestPoint & ip) {
    ip.descriptor = ip::InterestPoint::descriptor_type(); // this should free up the memory
  }

}} // namespace vw::ip
