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

#ifndef VW_GEOMETRY_CUTPOLYUTILS_H
#define VW_GEOMETRY_CUTPOLYUTILS_H
#include <vector>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <string>
#include <vw/Geometry/geomUtils.h>

namespace vw { namespace geometry {

  struct valIndex{
    double val;
    int    index;
    bool   isOutward;
    int    nextIndexInward; // Useful only when isOutward is true
  };

  void cutPolyLine(// inputs -- the polygonal line
                   int numVerts,
                   const double * xv, const double * yv,
                   // inputs -- the cutting window
                   double xll, double yll, double xur, double yur,
                   // outputs -- the cut polygons
                   std::vector< double> & cutX,
                   std::vector< double> & cutY,
                   std::vector< int>    & cutNumPolys);
  
  void cutPoly(// inputs -- the polygons
               int numPolys, const int * numVerts,
               const double * xv, const double * yv,
               // inputs -- the cutting window
               double xll, double yll, double xur, double yur,
               // outputs -- the cut polygons
               std::vector< double> & cutX,
               std::vector< double> & cutY,
               std::vector< int>    & cutNumPolys);

  inline bool lessThan (valIndex A, valIndex B){ return A.val < B.val; }
  
  void processPointsOnCutline(std::vector<valIndex> & ptsOnCutline);

  void cutToHalfSpace(// inputs 
                      double nx, double ny, double dotH,
                      int numV, 
                      const double * xv, const double * yv,
                      // outputs -- the cut polygons
                      std::vector<double> & cutX,
                      std::vector<double> & cutY,
                      std::vector<int>    & cutNumPolys);

}}

#endif // VW_GEOMETRY_CUTPOLYUTILS_H
  
