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

#include <cassert>
#include <algorithm>
#include <fstream>
#include <cfloat>
#include <vw/Geometry/cutPoly.h>
#include <vw/Geometry/baseUtils.h>
#include <vw/Geometry/edgeUtils.h>

namespace vw { namespace geometry {

#define DEBUG_CUT_POLY 0 // Must be 0 in production code

void cutPolyLine(// inputs -- the polygonal line
		 int numVerts,
		 const double * xv, const double * yv,
		 // inputs -- the cutting window
		 double xll, double yll, double xur, double yur,
		 // outputs -- the cut polygons
		 std::vector<double> & cutX,
		 std::vector<double> & cutY,
		 std::vector<int>    & cutNumPolys){
  
  // Cut a polygonal line. First make it into a polygon by traveling
  // forward and then backward on the polygonal line, then cut the
  // obtained polygon, then remove the backward portion from each
  // obtained polygon.

  std::vector<double> lXv, lYv, lCutX, lCutY;
  std::vector<int> lCutNumPolys;
  
  lXv.clear(); lYv.clear();
  for (int s = 0; s < numVerts; s++){
    lXv.push_back(xv[s]);
    lYv.push_back(yv[s]);
  }
  for (int s = numVerts - 1; s >= 0; s--){
    lXv.push_back(xv[s]);
    lYv.push_back(yv[s]);
  }
  int lNumVerts = lXv.size();
  
  cutPoly(// inputs -- the polygons
          1, &lNumVerts,  
          vecPtr(lXv), vecPtr(lYv),  
          // inputs -- the cutting window
          xll, yll, xur, yur,  
          // outputs -- the cut polygons
          lCutX, lCutY, lCutNumPolys);

  cutX.clear(); cutY.clear(); cutNumPolys.clear();
  
  int start = 0;
  for (int pIter = 0; pIter < (int)lCutNumPolys.size(); pIter++){

    if (pIter > 0) start += lCutNumPolys[pIter - 1];

    // Keep only half of the points of the cut polygon
    int half = lCutNumPolys[pIter]/2;
    cutNumPolys.push_back(half);
    for (int vIter = 0; vIter < half; vIter++){
      cutX.push_back(lCutX[start + vIter]);
      cutY.push_back(lCutY[start + vIter]);
    }
    
  }

  return;
}

void cutPoly(// inputs -- the polygons
                    int numPolys, const int * numVerts,
                    const double * xv, const double * yv,
                    // inputs -- the cutting window
                    double xll, double yll, double xur, double yur,
                    // outputs -- the cut polygons
                    std::vector< double> & cutX,
                    std::vector< double> & cutY,
                    std::vector< int>    & cutNumPolys){
  
  // Cut a given polygon with a box.
  
  // Intersect the polygon with each of the the half-planes
  // nx*x + ny*y <= (nx + ny)*H.
  // There are four values for the triplet (nx, ny, H):
  double cutParams[] = {
   -1,  0, xll, //  -- left cut
    1,  0, xur, //  -- right cut
    0, -1, yll, //  -- bottom cut
    0,  1, yur  //  -- top cut
  };

  int totalNumVerts = 0;
  for (int s = 0; s < numPolys; s++) totalNumVerts += numVerts[s];

  std::vector<double> Xin(xv, xv + totalNumVerts);        // A copy of xv as vector
  std::vector<double> Yin(yv, yv + totalNumVerts);        // A copy of yv as vector
  std::vector<int>    Pin(numVerts, numVerts + numPolys); // A copy of numVerts 
  
  std::vector<double> cutHalfX, cutHalfY, Xout, Yout;
  std::vector<int>    cutHalfP, Pout;

  for (int c = 0; c < 4; c++){

    Pout.clear(); Xout.clear(); Yout.clear();
    
    double nx   = cutParams[3*c + 0];
    double ny   = cutParams[3*c + 1];
    double H    = cutParams[3*c + 2];
    double dotH = (nx + ny)*H; // This formula works only for nx*ny == 0.

    int start = 0;
    for (int pIter = 0; pIter < (int)Pin.size(); pIter++){
      
      if (pIter > 0) start += Pin[pIter - 1];
      
      int numV = Pin[pIter];
      if (numV == 0) continue;
      
      cutToHalfSpace(nx, ny, dotH,
                     numV, vecPtr(Xin) + start, vecPtr(Yin) + start,
                     cutHalfX, cutHalfY, cutHalfP);
      
      for (int pIterCut = 0; pIterCut < (int)cutHalfP.size(); pIterCut++){
        if (cutHalfP[pIterCut] > 0){
          // Append only non-empty polygons
          Pout.push_back( cutHalfP[pIterCut] );
        }
      }
      
      for (int vIter = 0; vIter < (int)cutHalfX.size(); vIter++){
        Xout.push_back( cutHalfX[vIter] );
        Yout.push_back( cutHalfY[vIter] );
      }
      
    }
    
    Pin = Pout; Xin = Xout; Yin = Yout;
    
  } // End iterating over cutting lines

  cutNumPolys = Pout; cutX = Xout; cutY = Yout;
  
  return;
}

void cutToHalfSpace(// inputs 
                           double nx, double ny, double dotH,
                           int numV, 
                           const double * xv, const double * yv,
                           // outputs -- the cut polygons
                           std::vector<double> & cutX,
                           std::vector<double> & cutY,
                           std::vector<int>    & cutNumPolys){


  std::vector<valIndex> ptsOnCutline; ptsOnCutline.clear();
  valIndex C;
  
  cutX.clear(); cutY.clear(); cutNumPolys.clear();

  int cutPtsIndex = 0;
  
  for (int v = 0; v < numV; v++){

    int vnext = (v + 1)%numV;
    
    double xcurr = xv[v],     ycurr = yv[v];
    double xnext = xv[vnext], ynext = yv[vnext];

    double dotCurr = nx*xcurr + ny*ycurr;
    double dotNext = nx*xnext + ny*ynext;
    double cutx = 0.0, cuty = 0.0;
    
    if (dotCurr < dotH){

      // The current point is inside the half-plane
      
      cutX.push_back(xcurr);
      cutY.push_back(ycurr);
      cutPtsIndex++; 

      if (dotNext <= dotH) continue;
      
      cutEdge(xcurr, ycurr, xnext, ynext, nx, ny, dotH, cutx, cuty);

      cutX.push_back(cutx);
      cutY.push_back(cuty);

      C.val       = nx*cuty - ny*cutx;
      C.index     = cutPtsIndex;
      C.isOutward = true;
      ptsOnCutline.push_back(C);

      cutPtsIndex++; 
      
    }else if (dotCurr > dotH){

      // The current point is outside the half-plane
      
      if (dotNext >= dotH) continue;

      cutEdge(xcurr, ycurr, xnext, ynext, nx, ny, dotH, cutx, cuty);
      
      cutX.push_back(cutx);
      cutY.push_back(cuty);
      
      C.val       = nx*cuty - ny*cutx;
      C.index     = cutPtsIndex;
      C.isOutward = false;
      ptsOnCutline.push_back(C);

      cutPtsIndex++; 

    }else if (dotCurr == dotH){

      // The current point is at the edge of the half-plane

      int    vprev   = (v == 0) ? (numV - 1) : (v - 1);
      double xprev   = xv[vprev], yprev = yv[vprev];
      double dotPrev = nx*xprev + ny*yprev;
      
      if (dotPrev >= dotH && dotNext >= dotH) continue;
      
      cutX.push_back(xcurr);
      cutY.push_back(ycurr);

      if (dotPrev >= dotH || dotNext >= dotH){

        C.val       = nx*ycurr - ny*xcurr;
        C.index     = cutPtsIndex;
        C.isOutward = (dotPrev < dotH);
        ptsOnCutline.push_back(C);
        
      }

      cutPtsIndex++; 
      
    }
    
  }
  
  
  int numPtsOnCutline = ptsOnCutline.size();
  if (numPtsOnCutline == 0){
    cutNumPolys.push_back( cutX.size() );
    return;
  }

  processPointsOnCutline(ptsOnCutline);
  
  // Find the connected components in the cut polygons
  // To do: Move this to its own function.
  
  std::vector<double> X, Y;
  std::vector<int> P;
  X.clear(); Y.clear(); P.clear();

#if DEBUG_CUT_POLY
  static int c = -1;
  c++;
  char file[100];
  snprintf(file, sizeof(file), "beforeCleanup%d.xg", c);
  cout << "\nWriting to " << file << endl;
  ofstream before(file);
  for (int s = 0; s < (int)cutX.size(); s++){
    before << cutX[s] << ' ' << cutY[s] << endl;
    before << "anno " << cutX[s] << ' ' << cutY[s]  << ' ' << s << endl;
  }
  before.close();

  for (int s = 0; s < (int)ptsOnCutline.size(); s++){
    cout.precision(20);
    cout << "point on cutline is (index outward val) "
         << ptsOnCutline[s].index     << ' '
         << ptsOnCutline[s].isOutward << ' '
         << ptsOnCutline[s].val       << endl; 
  }
#endif

  std::vector<int> wasVisited;
  int numCutPts = cutX.size();
  wasVisited.assign(numCutPts, 0);

  int ptIter = 0;
  while(1){

    // Stop when all points are visited
    bool success = false;
    for (int v = 0; v < numCutPts; v++){
      if (!wasVisited[v]){
        ptIter  = v;
        success = true;
        break;
      }
    }
    if (!success) break; 

    int numPtsInComp = 0;
    
    // Visit a given connected component
    while(1){
      
      if (wasVisited[ptIter]){
        P.push_back(numPtsInComp);
        break; // Arrived back to the starting point of the given component
      }
      
      X.push_back(cutX[ptIter]);
      Y.push_back(cutY[ptIter]);
      wasVisited[ptIter] = 1;
#if DEBUG_CUT_POLY
      cout << "ptIter = " << ptIter << endl;
#endif
      numPtsInComp++;

      // Decide which point we will visit next
      
      if (nx*cutX[ptIter] + ny*cutY[ptIter] != dotH){
        // The point is not at the cutline
        ptIter = (ptIter + 1)%numCutPts;
        continue;
      }

      // The point is at the cutline. Find where exactly it is in the
      // sorted cutline points. If it is not among those sorted
      // cutline points, it means that the polygon only touches the
      // cutline at that point rather than crossing over to the other
      // side.
      // To do: Use here some faster lookup, such as a map.
      bool success = false;
      int cutlineIter = 0;
      for (cutlineIter = 0; cutlineIter < numPtsOnCutline; cutlineIter++){
        if (ptsOnCutline[cutlineIter].index == ptIter){
          success = true;
          break;
        }
      }
      if (!success){
        ptIter = (ptIter + 1)%numCutPts;
        continue;
      }

      if (ptsOnCutline[cutlineIter].isOutward){

        // We are getting out of the polygon. Find the point at
        // which we come back.
        ptIter = ptsOnCutline[cutlineIter].nextIndexInward;
        continue;
        
      }else{
        
        // The point ptIter is at the cutline on the way in.
        // The next point will be inside the current half-plane.
        ptIter = (ptIter + 1)%numCutPts;
        continue;
        
      }
      
    } // End iterating over all connected components

  } // End iterating over all points
  
  cutX        = X;
  cutY        = Y;
  cutNumPolys = P;

#if DEBUG_CUT_POLY
  snprintf(file, sizeof(file), "afterCleanup%d.xg", c);
  cout << "Writing to " << file << endl;
  ofstream after(file);
  for (int s = 0; s < (int)cutX.size(); s++){
    after << cutX[s] << ' ' << cutY[s] << endl;
    after << "anno " << cutX[s] << ' ' << cutY[s]  << ' ' << s << endl;
  }
  after.close();
#endif
  
}

void processPointsOnCutline(std::vector<valIndex> & ptsOnCutline){

  
  // Sort the cutline points along the cutline (the sort direction
  // does not matter).
  sort( ptsOnCutline.begin(), ptsOnCutline.end(), lessThan );

  // Find the position of each outward and each inward point on the cutline
  std::vector<int> outwardPositions, inwardPositions;
  outwardPositions.clear(); inwardPositions.clear();
  int numPtsOnCutline = ptsOnCutline.size();
  for (int s = 0; s < numPtsOnCutline; s++){
    const valIndex & C = ptsOnCutline[s];
    if (C.isOutward) outwardPositions.push_back(s);
    else             inwardPositions.push_back (s);
    //cout << "val index isOutward "
    //     << C.val << ' ' << C.index << ' ' << C.isOutward << endl; 
  }

  // There must be an even number of points on a cutline. That holds
  // true for any closed curve, with or without self-intersections.
  // Each time we cross the cutline to the other side at some point
  // there must be a corresponding point at which we come back. Match
  // the i-th outward point to the corresponding i-th inward point.
  int numOut = outwardPositions.size();
#ifndef NDEBUG
  int numIn  = inwardPositions.size();
  assert(numIn == numOut);
#endif
  for (int i = 0; i < numOut; i++){

    valIndex & C = ptsOnCutline[outwardPositions[i]]; //alias
    
    if (C.isOutward){
      C.nextIndexInward = ptsOnCutline[inwardPositions[i]].index;
#if DEBUG_CUT_POLY
      cout << "Going from " << C.index << " to " << C.nextIndexInward << endl;
#endif
    }else{
      C.nextIndexInward = -1; // To not leave it uninitialized
    }
    
  }
  
  return;
}

}}

