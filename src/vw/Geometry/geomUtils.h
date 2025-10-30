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

#ifndef VW_GEOMETRY_GEOMUTILS_H
#define VW_GEOMETRY_GEOMUTILS_H

#include <vw/Math/BBox.h>

#include <sstream>
#include <vector>
#include <fstream>
#include <cmath>
#include <set>
#include <cassert>

namespace vw { namespace geometry {

struct dPoint{
  double x, y;
  dPoint(): x(0), y(0){}
  dPoint(double x_in, double y_in): x(x_in), y(y_in){}
};


inline bool operator< (dPoint P, dPoint Q){
  return ( P.x < Q.x ) || (P.x == Q.x && P.y < Q.y);
}

inline bool greaterThan (dPoint P, dPoint Q){
  return ( P.x > Q.x ) || (P.x == Q.x && P.y > Q.y);
}

struct anno {

  double x;
  double y;
  std::string label;

  void appendTo(std::ofstream & outfile) const{
    outfile << "anno " << x << ' ' << y << ' ' << label << std::endl;
  }

};

std::ostream& operator<<(std::ostream& os, const anno& A);

  void snapPolyLineTo45DegAngles(bool isClosedPolyLine,
                                 int numVerts, double * xv, double * yv);
  void snapOneEdgeTo45(int numAngles, double* xs, double* ys,
                       bool snap2ndClosest,
                       double & x0, double & y0,
                       double & x1, double & y1);

  void minDistFromPtToSeg(//inputs
                          double xin, double yin,
                          double x0, double y0,
                          double x1, double y1,
                          // outputs
                          double & minX, double & minY,
                          double & minDist
                          );

  void searchForLayer(std::string   lineStr, // input
                      std::string & layer    // output
                      );

  double signedPolyArea(int numV, const double* xv, const double* yv, bool counter_cw);

  void searchForColor(std::string lineStr,  // input, not a reference on purpose
                      std::string & color); // output

  bool searchForAnnotation(std::string lineStr, anno & annotation);

  BBox2 expandBoxToRatio(vw::BBox2 const& box, double aspect);

  struct dRect{
    dRect(double xl_in = 0.0, double yl_in = 0.0,
          double xh_in = 0.0, double yh_in = 0.0):
      xl(xl_in), yl(yl_in), xh(xh_in), yh(yh_in) {}
    double xl, yl, xh, yh;

    bool isInSide(double x, double y) const {
      return (x >= xl && x <= xh && y >= yl && y <= yh);
    }
    
  };

  struct dRectWithId: public dRect {
    int id;
    dRectWithId(double xl_in = 0.0, double yl_in = 0.0,
                double xh_in = 0.0, double yh_in = 0.0,
                int id_in = 0):
      dRect(xl_in, yl_in, xh_in, yh_in), id(id_in){}
  };

  struct seg{
    double begx, begy, endx, endy;
    seg(double begx_in = 0.0, double begy_in = 0.0,
        double endx_in = 0.0, double endy_in = 0.0):
      begx(begx_in), begy(begy_in), endx(endx_in), endy(endy_in){}
  };


  struct segDist: public seg{
    double dist;
    segDist(double begx_in, double begy_in, double endx_in,
            double endy_in, double dist_in):
      seg(begx_in, begy_in, endx_in, endy_in), dist(dist_in){}
  };

  inline bool segDistGreaterThan(segDist s, segDist t){
    if (s.dist > t.dist) return true;
    if (s.dist < t.dist) return false;
    if (s.begx > t.begx) return true;
    if (s.begx < t.begx) return false;
    if (s.begy > t.begy) return true;
    if (s.begy < t.begy) return false;
    return false;
  }

  inline bool operator==(segDist s, segDist t){
    return ( s.dist == t.dist ) && ( s.begx == t.begx ) && ( s.begy == t.begy );
  }

  inline std::ostream& operator<<(std::ostream & output, const segDist & S) {
    output << S.begx << ' ' << S.begy << ' ' << S.endx << ' ' << S.endy << ' '
           << S.dist;
    return output;  // for multiple << operators
  }

  bool boxesIntersect(double xl1, double yl1, double xh1, double yh1,
                      double xl2, double yl2, double xh2, double yh2
                      );


  bool mergePolys(int an,
                  const double * ax_in, const double * ay_in,
                  int bn,
                  const double * bx_in, const double * by_in,
                  std::vector<double> & mergedX,
                  std::vector<double> & mergedY
                  );
  bool isPointInPolyOrOnEdges(double x, double y,
                              int n, const double* xv, const double*  yv);

  struct linTrans{
    // Linear transform
    double a11, a12, a21, a22, sx, sy;
    linTrans(){
      a11 = a22 = 1.0;
      a21 = a12 = sx = sy = 0.0;
    }
    void reset(){
      *this = linTrans();
    }
    void print(){
      std::cout << "transform " << a11 << ' ' << a12 << ' '
                << a21 << ' ' << a22 << ' ' << sx << ' ' << sy << std::endl;
    }
  };

  linTrans composeTransforms(linTrans P, linTrans Q);

  struct matrix2{
    // A 2x2 matrix
    double a11, a12, a21, a22;
    matrix2(){
      a11 = a22 = 1.0;
      a21 = a12 = 0.0;
    }
    void reset(){
      *this = matrix2();
    }
    void print(){
      std::cout << "matrix " << a11 << ' ' << a12 << ' ' << a21 << ' ' << a22 << std::endl;
    }
  };

  linTrans transAroundPt(const matrix2 & M, dPoint P);

}} // end namespace vw::geometry
                                                                   
#endif // VW_GEOMETRY_GEOMUTILS_H
