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


#include <vector>
#include <gtest/gtest_VW.h>
#include <vw/Math/FLANNTree.h>

using std::vector;
using namespace vw;
using namespace vw::math;



// Make sure that vector bounds are respected in small data sets
TEST(FLANNTree, boundsLimiting) {

  // Generate an arbitrary set of 2D points
  const int numPts = 15;
  Matrix<float> locations(numPts, 2);
  for (int i=0; i<numPts; ++i) {
    locations(i, 0) = i*i;
    locations(i, 1) = 1000 - i*i;
  }

  math::FLANNTree<float> tree;
  tree.load_match_data(locations, FLANN_DistType_L2);
  printf("---------\n");
  Vector<int>    indices;
  Vector<double> distance;
  tree.knn_search(select_row(locations, 9),
                  indices, distance, 11); // Request 11 points
  EXPECT_EQ(indices.size(), 11); // Returned the requested number of points
  for (size_t i=0; i<indices.size(); ++i) {
    EXPECT_LT(indices[i], numPts);
  }

  tree.knn_search(select_row(locations, 3),
                  indices, distance, 20); // Request 20 points
  EXPECT_EQ(indices.size(), numPts); // Size limited by number of points
  for (size_t i=0; i<indices.size(); ++i) {
    EXPECT_LT(indices[i], numPts);
  }

}
