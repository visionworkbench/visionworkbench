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

// For MS Windows
#define _USE_MATH_DEFINES 

#include <vw/Geometry/cutPoly.h>
#include <vw/Geometry/dPoly.h>
#include <vw/Core/Log.h>

#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include <cassert>
#include <cfloat>
#include <cassert>
#include <cstring>
#include <string>
#include <map>

namespace vw { namespace geometry {

using namespace std;

// A double precision polygon class

void dPoly::reset() {
  m_isPointCloud  = false;
  m_numPolys      = 0;
  m_totalNumVerts = 0;
  m_numVerts.clear();
  m_xv.clear();
  m_yv.clear();
  m_isPolyClosed.clear();
  m_colors.clear();
  m_layers.clear();
  m_annotations.clear();
  m_vertIndexAnno.clear();
  m_polyIndexAnno.clear();
  m_layerAnno.clear();
}

void dPoly::bdBox(double & xll, double & yll, double & xur, double & yur) const {

  // The bounding box of all polygons

  if (m_totalNumVerts <= 0) {
    xll = DBL_MAX/4.0, xur = -DBL_MAX/4.0; // Use 1/4.0 to avoid overflow when ...
    yll = DBL_MAX/4.0, yur = -DBL_MAX/4.0; // ... finding width and height
    return;
  }

  xll = *min_element( vecPtr(m_xv), vecPtr(m_xv) + m_totalNumVerts );
  yll = *min_element( vecPtr(m_yv), vecPtr(m_yv) + m_totalNumVerts );
  xur = *max_element( vecPtr(m_xv), vecPtr(m_xv) + m_totalNumVerts );
  yur = *max_element( vecPtr(m_yv), vecPtr(m_yv) + m_totalNumVerts );

  return;
};

void dPoly::bdBoxes(std::vector<double> & xll, std::vector<double> & yll,
                    std::vector<double> & xur, std::vector<double> & yur) const {

  // Bounding boxes of individual polygons

  xll.clear(); yll.clear(); xur.clear(); yur.clear();

  int start = 0;
  for (int pIter = 0; pIter < m_numPolys; pIter++) {

    if (pIter > 0) start += m_numVerts[pIter - 1];

    int numV = m_numVerts[pIter];

    double x0, y0, x1, y1;
    if (numV <= 0) {
      x0 = DBL_MAX/4.0, x1 = -DBL_MAX/4.0; // Use 1/4.0 to avoid overflow when ...
      y0 = DBL_MAX/4.0, y1 = -DBL_MAX/4.0; // ... finding width and height
    }else{
      const double * px = vecPtr(m_xv) + start;
      const double * py = vecPtr(m_yv) + start;
      x0 = *min_element( px, px + numV ); x1 = *max_element( px, px + numV );
      y0 = *min_element( py, py + numV ); y1 = *max_element( py, py + numV );
    }
    xll.push_back(x0); xur.push_back(x1);
    yll.push_back(y0); yur.push_back(y1);

  }

  return;
};

void dPoly::bdBoxCenter(double & mx, double & my) const {

  double xll, yll, xur, yur;
  bdBox(xll, yll, xur, yur);

  mx = (xll + xur)/2.0;
  my = (yll + yur)/2.0;

  return;
}

void dPoly::setPolygon(int numVerts,
                       const double * xv,
                       const double * yv,
                       bool isPolyClosed,
                       const std::string & color,
                       const std::string & layer) {
  reset();
  appendPolygon(numVerts, xv, yv, isPolyClosed, color, layer);
  return;
}

void dPoly::appendPolygon(int numVerts,
                          const double * xv,
                          const double * yv,
                          bool isPolyClosed,
                          const std::string & color,
                          const std::string & layer) {

  if (numVerts <= 0) return;

  m_numPolys      += 1;
  m_totalNumVerts += numVerts;

  m_numVerts.push_back(numVerts);
  m_isPolyClosed.push_back(isPolyClosed);
  m_colors.push_back(color);
  m_layers.push_back(layer);
  for (int s = 0; s < numVerts; s++) {
    m_xv.push_back(xv[s]);
    m_yv.push_back(yv[s]);
  }

  return;
}

void dPoly::appendRectangle(double xl, double yl, double xh, double yh,
                            bool isPolyClosed,
                            const std::string & color, const std::string & layer) {
  double xv[4], yv[4];
  xv[0] = xl; xv[1] = xh; xv[2] = xh; xv[3] = xl;
  yv[0] = yl; yv[1] = yl; yv[2] = yh; yv[3] = yh;
  appendPolygon(4, xv, yv, isPolyClosed, color, layer);
  return;
}

void dPoly::setRectangle(double xl, double yl, double xh, double yh,
                         bool isPolyClosed,
                         const std::string & color, const std::string & layer) {
  reset();
  appendRectangle(xl, yl, xh, yh, isPolyClosed, color, layer);
  return;
}

bool dPoly::isXYRect() {

  // Check if the current polygon set is a (perhaps degenerate)
  // rectangle with sides parallel to the x and y axes.

  if (m_numPolys != 1 || m_totalNumVerts != 4) return false;

  double xll, yll, xur, yur;
  bdBox(xll, yll, xur, yur);
  double tol = 1e-15*(std::abs(xll) + std::abs(yll) +
                      std::abs(xur) + std::abs(yur));

  // Midpoint test (will catch things the direction test won't)
  if (std::abs( m_xv[0] + m_xv[2] - m_xv[1] - m_xv[3] ) > tol ) return false;
  if (std::abs( m_yv[0] + m_yv[2] - m_yv[1] - m_yv[3] ) > tol ) return false;

  // Direction test
  for (int i = 0; i < m_totalNumVerts; i++) {
    int i1 = ( i + m_totalNumVerts + 1) % m_totalNumVerts;
    if (std::abs(m_xv[i1] - m_xv[i]) > tol &&
        std::abs(m_yv[i1] - m_yv[i]) > tol )
      return false;
  }

  return true;
}

void dPoly::get_annoByType(std::vector<anno> & annotations, AnnoType annoType) {

  if (annoType == fileAnno) {
    get_annotations(annotations);
  }else if (annoType == vertAnno) {
    get_vertIndexAnno(annotations);
  }else if (annoType == polyAnno) {
    get_polyIndexAnno(annotations);
  }else if (annoType == layerAnno) {
    get_layerAnno(annotations);
  }else{
    vw::vw_out() << "Unknown annotation type.\n";
  }

  return;
}

void dPoly::set_annoByType(const std::vector<anno> & annotations, AnnoType annoType) {

  if (annoType == fileAnno) {
    set_annotations(annotations);
  }else if (annoType == vertAnno) {
    set_vertIndexAnno(annotations);
  }else if (annoType == polyAnno) {
    set_polyIndexAnno(annotations);
  }else if (annoType == layerAnno) {
    set_layerAnno(annotations);
  }else{
    vw::vw_out() << "Unknown annotation type.\n";
  }

  return;
}

void dPoly::clipPoly(// inputs
                     double clip_xll, double clip_yll,
                     double clip_xur, double clip_yur,
                     dPoly & clippedPoly) { // output

  assert(this != &clippedPoly); // source and destination must be different

  clippedPoly.reset();
  clippedPoly.set_isPointCloud(m_isPointCloud);

  const double * xv               = get_xv();
  const double * yv               = get_yv();
  const int    * numVerts         = get_numVerts();
  int numPolys                    = get_numPolys();
  const vector<char> isPolyClosed = get_isPolyClosed();
  const vector<string> colors     = get_colors();
  const vector<string> layers     = get_layers();

  int start = 0;
  for (int pIter = 0; pIter < numPolys; pIter++) {

    if (pIter > 0) start += numVerts[pIter - 1];

    int  isClosed = isPolyClosed [pIter];
    string color  = colors       [pIter];
    string layer  = layers       [pIter];

    vector<double> cutXv, cutYv;
    vector<int> cutNumVerts;
    cutXv.clear(); cutYv.clear(); cutNumVerts.clear();

    if (m_isPointCloud) {

      // To cut a point cloud to a box all is needed is to select
      // which points are in the box
      for (int vIter = 0; vIter < numVerts[pIter]; vIter++) {
        double x = xv[start + vIter];
        double y = yv[start + vIter];
        if (x >= clip_xll && x <= clip_xur &&
            y >= clip_yll && y <= clip_yur) {
          cutXv.push_back(x);
          cutYv.push_back(y);
        }
      }
      cutNumVerts.push_back( cutXv.size() );

    }else if (isClosed) {

      cutPoly(1, numVerts + pIter, xv + start, yv + start,
              clip_xll, clip_yll, clip_xur, clip_yur,
              cutXv, cutYv, cutNumVerts); // outputs

    }else{

      cutPolyLine(numVerts[pIter], xv + start, yv + start,
                  clip_xll, clip_yll, clip_xur, clip_yur,
                  cutXv, cutYv, cutNumVerts); // outputs
    }

    int cstart = 0;
    for (int cIter = 0; cIter < (int)cutNumVerts.size(); cIter++) {

      if (cIter > 0) cstart += cutNumVerts[cIter - 1];
      int cSize = cutNumVerts[cIter];
      clippedPoly.appendPolygon(cSize,
                                vecPtr(cutXv) + cstart,
                                vecPtr(cutYv) + cstart,
                                isClosed, color, layer);
    }

  }

  // Cutting inherits the annotations at the vertices of the uncut
  // polygons which are in the cutting box.
  vector<anno> annotations, annoInBox;

  for (int annoType = fileAnno; annoType < lastAnno; annoType++) {

    get_annoByType(annotations, (AnnoType)annoType);

    annoInBox.clear();
    for (int s = 0; s < (int)annotations.size(); s++) {
      const anno & A = annotations[s];
      if (clip_xll <= A.x && A.x <= clip_xur && clip_yll <= A.y && A.y <= clip_yur) {
        annoInBox.push_back(A);
      }
    }

    clippedPoly.set_annoByType(annoInBox, (AnnoType)annoType);
  }

  return;
}

void dPoly::shift(double shift_x, double shift_y) {

  // To do: Need to integrate the several very similar transform functions

  for (int i = 0; i < (int)m_xv.size(); i++) {
    m_xv[i] += shift_x;
    m_yv[i] += shift_y;
  }

  std::vector<anno> annotations;
  for (int annoType = fileAnno; annoType < lastAnno; annoType++) {
    get_annoByType(annotations, (AnnoType)annoType);
    for (int i = 0; i < (int)annotations.size(); i++) {
      anno & A = annotations[i]; // alias
      A.x += shift_x;
      A.y += shift_y;
    }
    set_annoByType(annotations, (AnnoType)annoType);
  }

  return;
}

void dPoly::rotate(double angle) { // The angle is given in degrees

  // To do: Need to integrate the several very similar transform functions

  double a = angle*M_PI/180.0, c = cos(a), s= sin(a);

  if (angle == round(angle) && int(angle)%90 == 0 ) {
    // The special case of angle multiple of 90 degrees
    c = round(c), s = round(s);
  }

  for (int i = 0; i < (int)m_xv.size(); i++) {
    double tmpx = c*m_xv[i] - s*m_yv[i];
    double tmpy = s*m_xv[i] + c*m_yv[i];
    m_xv[i] = tmpx;
    m_yv[i] = tmpy;
  }

  vector<anno> annotations;
  for (int annoType = fileAnno; annoType < lastAnno; annoType++) {
    get_annoByType(annotations, (AnnoType)annoType);
    for (int i = 0; i < (int)annotations.size(); i++) {
      anno & A = annotations[i]; // alias
      double tmpx = c*A.x - s*A.y;
      double tmpy = s*A.x + c*A.y;
      A.x = tmpx;
      A.y = tmpy;
    }
    set_annoByType(annotations, (AnnoType)annoType);
  }

  return;
}

void dPoly::scale(double scale) {

  // To do: Need to integrate the several very similar transform functions

  for (int i = 0; i < (int)m_xv.size(); i++) {
    m_xv[i] *= scale;
    m_yv[i] *= scale;
  }

  vector<anno> annotations;
  for (int annoType = fileAnno; annoType < lastAnno; annoType++) {

    get_annoByType(annotations, (AnnoType)annoType);
    for (int i = 0; i < (int)annotations.size(); i++) {
      anno & A = annotations[i]; // alias
      A.x *= scale;
      A.y *= scale;
    }

    set_annoByType(annotations, (AnnoType)annoType);
  }

  return;
}

void dPoly::transformMarkedPolys(std::map<int, int> const& mark, const linTrans & T) {

  int start = 0;
  for (int pIter = 0; pIter < m_numPolys; pIter++) {
    if (pIter > 0) start += m_numVerts[pIter - 1];
    if (mark.find(pIter) == mark.end()) continue;

    for (int vIter = 0; vIter < m_numVerts[pIter]; vIter++) {
      int i = start + vIter;
      double x = T.a11*m_xv[i] + T.a12*m_yv[i] + T.sx;
      double y = T.a21*m_xv[i] + T.a22*m_yv[i] + T.sy;
      m_xv[i] = x;
      m_yv[i] = y;
    }

  }

  m_vertIndexAnno.clear();

  return;
}

void dPoly::transformMarkedAnnos(std::map<int, int> const& mark, const linTrans & T) {
  for (size_t it = 0; it < m_annotations.size(); it++) {

    if (mark.find(it) == mark.end()) continue;

    anno & A = m_annotations[it]; // alias
    double x = T.a11*A.x + T.a12*A.y + T.sx;
    double y = T.a21*A.x + T.a22*A.y + T.sy;
    A.x = x;
    A.y = y;
  }
}

void dPoly::transformMarkedPolysAroundPt(std::map<int, int> const& mark, const matrix2 & M,
                                         dPoint P) {
  linTrans T = transAroundPt(M, P);
  transformMarkedPolys(mark, T);
  return;
}

void dPoly::transformMarkedAnnosAroundPt(std::map<int, int> const& mark, const matrix2 & M,
                                         dPoint P) {
  linTrans T = transAroundPt(M, P);
  transformMarkedAnnos(mark, T);
  return;
}

void dPoly::applyTransform(double a11, double a12, double a21, double a22,
                           double sx, double sy,
                           linTrans & T) { // save the transform here

  // To do: Need to integrate the several very similar transform functions

  // Save the transform before applying it
  T.a11 = a11; T.a12 = a12; T.a21 = a21; T.a22 = a22; T.sx = sx; T.sy = sy;

  for (int i = 0; i < (int)m_xv.size(); i++) {
    double x = a11*m_xv[i] + a12*m_yv[i] + sx;
    double y = a21*m_xv[i] + a22*m_yv[i] + sy;
    m_xv[i] = x;
    m_yv[i] = y;
  }

  vector<anno> annotations;
  for (int annoType = fileAnno; annoType < lastAnno; annoType++) {
    get_annoByType(annotations, (AnnoType)annoType);
    for (int i = 0; i < (int)annotations.size(); i++) {
      anno & A = annotations[i]; // alias
      double x = a11*A.x + a12*A.y + sx;
      double y = a21*A.x + a22*A.y + sy;
      A.x = x;
      A.y = y;
    }
    set_annoByType(annotations, (AnnoType)annoType);
  }

  return;
}


void dPoly::applyTransformAroundBdBoxCenter(double a11, double a12,
                                            double a21, double a22,
                                            linTrans & T) {

  if (m_totalNumVerts == 0) return;

  double mx, my;
  bdBoxCenter(mx, my);
  applyTransform(a11, a12, a21, a22, mx - a11*mx - a12*my, my - a21*mx - a22*my, T);

  return;
}

void dPoly::appendPolygons(const dPoly & poly) {

  const double * xv         = poly.get_xv();
  const double * yv         = poly.get_yv();
  const int    * numVerts   = poly.get_numVerts();
  int numPolys              = poly.get_numPolys();
  vector<char> isPolyClosed = poly.get_isPolyClosed();
  vector<string> colors     = poly.get_colors();
  vector<string> layers     = poly.get_layers();
  vector<anno> annotations;  poly.get_annotations(annotations);

  int start = 0;
  for (int pIter = 0; pIter < numPolys; pIter++) {

    if (pIter > 0) start += numVerts[pIter - 1];

    bool isClosed = isPolyClosed [pIter];
    string color  = colors       [pIter];
    string layer  = layers       [pIter];
    int pSize     = numVerts     [pIter];

    appendPolygon(pSize, xv + start, yv + start, isClosed, color, layer);

  }

  for (int s = 0; s < (int)annotations.size(); s++) {
    addAnno(annotations[s]);
  }

  return;
}

void dPoly::get_annotations(std::vector<anno> & annotations) const {
  annotations = m_annotations;
}

void dPoly::set_annotations(const std::vector<anno> & A) {
  m_annotations = A;
}

void dPoly::set_vertIndexAnno(const std::vector<anno> & annotations) {
  m_vertIndexAnno = annotations;
}

void dPoly::get_vertIndexAnno(std::vector<anno> & annotations) const {
  annotations = m_vertIndexAnno;
}

void dPoly::set_polyIndexAnno(const std::vector<anno> & annotations) {
  m_polyIndexAnno = annotations;
}

void dPoly::get_polyIndexAnno(std::vector<anno> & annotations) const {
  annotations = m_polyIndexAnno;
}

void dPoly::set_layerAnno(const std::vector<anno> & annotations) {
  m_layerAnno = annotations;
}

void dPoly::get_layerAnno(std::vector<anno> & annotations) const {
  annotations = m_layerAnno;
}

void dPoly::set_color(std::string color) {

  m_colors.resize(m_numPolys);
  for (int s = 0; s < (int)m_colors.size(); s++) {
    m_colors[s] = color;
  }

  return;
}

void dPoly::set_isPolyClosed(bool isPolyClosed) {

  m_isPolyClosed.resize(m_numPolys);
  for (int s = 0; s < (int)m_isPolyClosed.size(); s++) {
    m_isPolyClosed[s] = isPolyClosed;
  }

  return;
}

void dPoly::compVertFullIndexAnno() {

  m_vertIndexAnno.clear();

  const double * xv = get_xv();
  const double * yv = get_yv();

  int start = 0;
  for (int pIter = 0; pIter < m_numPolys; pIter++) {

    if (pIter > 0) start += m_numVerts[pIter - 1];

    for (int v = 0; v < m_numVerts[pIter]; v++) {

      anno A;
      A.x     = xv[start + v];
      A.y     = yv[start + v];
      A.label = num2str(start +v);
      m_vertIndexAnno.push_back(A);
    }

  }

  return;
}


void dPoly::compVertIndexAnno() {

  m_vertIndexAnno.clear();

  const double * xv = get_xv();
  const double * yv = get_yv();
  
  int start = 0;
  for (int pIter = 0; pIter < m_numPolys; pIter++) {

    if (pIter > 0) start += m_numVerts[pIter - 1];

    for (int v = 0; v < m_numVerts[pIter]; v++) {

      anno A;
      A.x     = xv[start + v];
      A.y     = yv[start + v];
      A.label = num2str(v);
      m_vertIndexAnno.push_back(A);
    }

  }

  return;
}

// Place the poly index annotation at the poly center of graviity
void dPoly::compPolyIndexAnno() {

  m_polyIndexAnno.clear();

  const double * xv = get_xv();
  const double * yv = get_yv();
  
  int start = 0;
  for (int pIter = 0; pIter < m_numPolys; pIter++) {

    if (pIter > 0)
      start += m_numVerts[pIter - 1];
    
    double x = 0.0, y = 0.0, num = 0.0;
    for (int v = 0; v < m_numVerts[pIter]; v++) {
      x   += xv[start + v];
      y   += yv[start + v];
      num += 1.0;
    }

    if (num > 0) {
      anno A;
      A.x     = x/num;
      A.y     = y/num;
      A.label = num2str(pIter);
      m_polyIndexAnno.push_back(A);
    }
  }
  
  return;
}

void dPoly::compLayerAnno() {

  m_layerAnno.clear();

  const double * xv = get_xv();
  const double * yv = get_yv();

  int start = 0;
  for (int pIter = 0; pIter < m_numPolys; pIter++) {

    if (pIter > 0) start += m_numVerts[pIter - 1];

    for (int v = 0; v < m_numVerts[pIter]; v++) {

      anno A;
      int vn = (v+1)%m_numVerts[pIter];

      A.x     = (xv[start + v] + xv[start + vn])/2.0; // put anno at midpt
      A.y     = (yv[start + v] + yv[start + vn])/2.0; // put anno at midpt
      A.label = m_layers[pIter];
      m_layerAnno.push_back(A);
    }

  }

  return;
}

void dPoly::eraseAnno(int annoIndex) {

  assert(0 <= annoIndex && annoIndex < (int)m_annotations.size());

  m_annotations.erase(m_annotations.begin() + annoIndex,
                      m_annotations.begin() + annoIndex + 1);
  return;
}

void dPoly::findClosestAnnotation(// inputs
                                  double x0, double y0,
                                  // outputs
                                  int & annoIndex,
                                  double & min_dist
                                  ) const {

  // Given a point and a set of polygons, find the annotation
  // closest to the given point. Return the closest index and the
  // distance to it from the given point. Return DBL_MAX if
  // there are no annotations

  min_dist = DBL_MAX;
  annoIndex = -1;

  for (size_t aIter = 0; aIter < m_annotations.size(); aIter++) {
    double dist = distance(x0, y0, m_annotations[aIter].x, m_annotations[aIter].y);
    if (dist <= min_dist) {
      annoIndex = aIter;
      min_dist  = dist;
    }
  }

  return;
}


// Given a point and a set of polygons, find the polygon vertex
// closest to the given point. Return the closest vertex and the
// distance to it from the given point. Return DBL_MAX if the
// polygon is empty.
void dPoly::findClosestPolyVertex(// inputs
                                  double x0, double y0,
                                  // outputs
                                  int & polyIndex,
                                  int & vertIndex,
                                  double & min_x, double & min_y,
                                  double & min_dist
                                  ) const{

  min_x = x0; min_y = y0; min_dist = DBL_MAX;
  polyIndex = -1; vertIndex = -1;

  int start = 0;
  for (int pIter = 0; pIter < m_numPolys; pIter++) {

    if (pIter > 0) start += m_numVerts[pIter - 1];

    for (int vIter = 0; vIter < m_numVerts[pIter]; vIter++) {
      double dist = distance(x0, y0, m_xv[start + vIter], m_yv[start + vIter]);
      if (dist <= min_dist) {
        polyIndex = pIter;
        vertIndex = vIter;
        min_dist  = dist;
        min_x     = m_xv[start + vIter];
        min_y     = m_yv[start + vIter];
      }
    }

  }

  return;
}

// Given a point and a set of polygons, find the polygon edge
// closest to the given point and the location on the edge where the
// smallest distance is achieved. Return the index of the polygon
// where the closest distance is achieved, as well as the point at
// which that distance is achieved and the smallest distance itself.
void dPoly::findClosestPolyEdge(//inputs
                                 double x0, double y0,
                                 // outputs
                                 int & polyIndex, int & vertIndex,
                                 double & minX, double & minY, double & minDist
                                 ) const{

  polyIndex = -1;
  vertIndex = -1;
  minX     = DBL_MAX, minY = DBL_MAX;
  minDist  = DBL_MAX;

  int start = 0;
  double xval, yval;
  for (int pIter = 0; pIter < m_numPolys; pIter++) {

    if (pIter > 0) start += m_numVerts[pIter - 1];

    for (int vIter = 0; vIter < m_numVerts[pIter]; vIter++) {

      int beg = start + vIter;
      int end = start + (vIter + 1)%m_numVerts[pIter];

      double dist = DBL_MAX;
      minDistFromPtToSeg(// inputs
                         x0, y0, m_xv[beg], m_yv[beg], m_xv[end], m_yv[end],
                         // outputs
                         xval, yval, dist
                         );

      if (dist <= minDist) {
        polyIndex = pIter;
        vertIndex = vIter;
        minX      = xval;
        minY      = yval;
        minDist   = dist;
      }

    }

  }

  return;
}

void dPoly::eraseOnePoly(int polyIndex) {

  assert(0 <= polyIndex && polyIndex < m_numPolys);

  map<int, int> mark;
  mark[polyIndex] = 1;
  eraseMarkedPolys(mark);

  return;
}

void dPoly::insertVertex(int polyIndex, int vertIndex,
                         double x, double y) {

  assert(0 <= polyIndex && polyIndex < m_numPolys);
  assert(0 <= vertIndex && vertIndex < m_numVerts[polyIndex] + 1);

  int start = 0;
  for (int pIter = 0; pIter < polyIndex; pIter++) {
    start += m_numVerts[pIter];
  }

  int iv = start + vertIndex;
  m_xv.insert(m_xv.begin() + iv, x);
  m_yv.insert(m_yv.begin() + iv, y);

  m_totalNumVerts++;
  m_numVerts[polyIndex]++;

  m_vertIndexAnno.clear();

  return;
}

void dPoly::eraseVertex(int polyIndex, int vertIndex) {

  assert(0 <= polyIndex && polyIndex < m_numPolys);
  assert(0 <= vertIndex && vertIndex < m_numVerts[polyIndex]);

  int start = 0;
  for (int pIter = 0; pIter < polyIndex; pIter++) {
    start += m_numVerts[pIter];
  }

  int iv = start + vertIndex;
  m_xv.erase(m_xv.begin() + iv, m_xv.begin() + iv + 1);
  m_yv.erase(m_yv.begin() + iv, m_yv.begin() + iv + 1);

  m_totalNumVerts--;
  m_numVerts[polyIndex]--;

  m_vertIndexAnno.clear();

  return;
}

void dPoly::changeVertexValue(int polyIndex, int vertIndex, double x, double y) {

  assert(0 <= polyIndex && polyIndex < m_numPolys);
  assert(0 <= vertIndex && vertIndex < m_numVerts[polyIndex]);

  int start = 0;
  for (int pIter = 0; pIter < polyIndex; pIter++) {
    start += m_numVerts[pIter];
  }

  m_xv[start + vertIndex] = x;
  m_yv[start + vertIndex] = y;

  m_vertIndexAnno.clear();

  return;
}

void dPoly::shiftEdge(int polyIndex, int vertIndex, double shift_x, double shift_y) {

  assert(0 <= polyIndex && polyIndex < m_numPolys);
  assert(0 <= vertIndex && vertIndex < m_numVerts[polyIndex]);

  int start = 0;
  for (int pIter = 0; pIter < polyIndex; pIter++) {
    start += m_numVerts[pIter];
  }

  // Beginning point of the edge
  m_xv[start + vertIndex] += shift_x;
  m_yv[start + vertIndex] += shift_y;

  m_vertIndexAnno.clear();

  // End point of the edge
  if (m_numVerts[polyIndex] <= 1) return;
  int vertIndexEnd = (vertIndex + 1)%m_numVerts[polyIndex];
  m_xv[start + vertIndexEnd] += shift_x;
  m_yv[start + vertIndexEnd] += shift_y;

  return;
}

void dPoly::shiftOnePoly(int polyIndex, double shift_x, double shift_y) {

  assert(0 <= polyIndex && polyIndex < m_numPolys);

  int start = 0;
  for (int pIter = 0; pIter < polyIndex; pIter++) {
    start += m_numVerts[pIter];
  }

  for (int vIter = 0; vIter < m_numVerts[polyIndex]; vIter++) {
    m_xv[start + vIter] += shift_x;
    m_yv[start + vIter] += shift_y;
  }

  m_vertIndexAnno.clear();

  return;
}

void dPoly::shiftOneAnno(int index, double shift_x, double shift_y) {

  assert(0 <= index && index < (int)m_annotations.size());

  m_annotations[index].x += shift_x;
  m_annotations[index].y += shift_y;

  m_vertIndexAnno.clear();

  return;
}

void dPoly::shiftMarkedPolys(std::map<int, int> const & mark, double shift_x, double shift_y) {

  for (auto it = mark.begin(); it != mark.end(); it++)
    shiftOnePoly(it->first, shift_x, shift_y);

  return;
}

void dPoly::shiftMarkedAnnos(std::map<int, int> const & amark, double shift_x, double shift_y) {

  for (auto it = amark.begin(); it != amark.end(); it++)
    shiftOneAnno(it->first, shift_x, shift_y);

  return;
}

void dPoly::reverseMarkedPolys(std::map<int, int> const & mark) {

  for (auto it = mark.begin(); it != mark.end(); it++)
    reverseOnePoly(it->first);

  return;
}

std::vector<int> dPoly::get_startingIndices() const{
	std::vector<int> start_indices;

	int start = 0;
	start_indices.push_back(start);
	for (int pIter = 0; pIter < m_numPolys; pIter++) {
		start += m_numVerts[pIter];
		start_indices.push_back(start);
	}
	return start_indices;
}

void dPoly::extractOnePoly(int polyIndex, // input
                           dPoly & poly,   // output
						   int start_index) const {

  assert(0 <= polyIndex && polyIndex < m_numPolys);

  int start = 0;
  if (start_index >= 0)
	  start = start_index;
  else {

	  for (int pIter = 0; pIter < polyIndex; pIter++) {
		  start += m_numVerts[pIter];
	  }
  }

  poly.setPolygon(m_numVerts[polyIndex],
                  vecPtr(m_xv) + start, vecPtr(m_yv) + start,
                  m_isPolyClosed[polyIndex],
                  m_colors[polyIndex],
                  m_layers[polyIndex]);

  return;
}

void dPoly::extractMarkedPolys(std::map<int, int> const& mark, // input
                               dPoly & polys) const {          // output 

  const auto &start_ids = get_startingIndices();

  polys.reset();
  for (int pIter = 0; pIter < m_numPolys; pIter++) {
    if (mark.find(pIter) == mark.end()) continue;
    dPoly lPoly;
    extractOnePoly(pIter, // input
                   lPoly,  // output
				   start_ids[pIter]
                   );
    polys.appendPolygons(lPoly);
  }

  return;
}

// Reverse orientation of all polygons
void dPoly::reverse(){

  int start = 0;
  for (int pIter = 0; pIter < m_numPolys; pIter++){
    if (pIter > 0) start += m_numVerts[pIter - 1];
    std::reverse(vecPtr(m_xv) + start, vecPtr(m_xv) + start + m_numVerts[pIter]);
    std::reverse(vecPtr(m_yv) + start, vecPtr(m_yv) + start + m_numVerts[pIter]);
  }
  
  return;
}
  
void dPoly::reverseOnePoly(int polyIndex){

  assert(0 <= polyIndex && polyIndex < m_numPolys);

  int start = 0;
  for (int pIter = 0; pIter < polyIndex; pIter++){
    start += m_numVerts[pIter];
  }

  std::reverse(vecPtr(m_xv) + start, vecPtr(m_xv) + start + m_numVerts[polyIndex]);
  std::reverse(vecPtr(m_yv) + start, vecPtr(m_yv) + start + m_numVerts[polyIndex]);

  return;
}

namespace dPoly_local_functions{
  struct ptAndIndex{
    dPoint point;
    double area;
    int    index;
  };
  bool greaterThanPtIndex (ptAndIndex P, ptAndIndex Q) {

    if (greaterThan(P.point, Q.point) ) return true;
    if (greaterThan(Q.point, P.point) ) return false;

    return (P.area > Q.area);

  }
}

void dPoly::sortFromLargestToSmallest(bool counter_cc) {

  // Sort the polygons so that if polygon A is inside of polygon B, then
  // polygon B shows up before polygon A after sorting.

  // Use the fact that if A is inside of B, then the dimensions/area of A
  // are no more than those of B.

  using namespace dPoly_local_functions;

  // Find the bounding boxes of polygons
  vector<double> xll, yll, xur, yur;
  bdBoxes(xll, yll, xur, yur);

  int numPolys = xll.size();

  vector<ptAndIndex> boxDims;
  int start = 0;
  boxDims.resize(numPolys);
  for (int s = 0; s < numPolys; s++) {

    if (s > 0) start += m_numVerts[s - 1];
    int numV = m_numVerts[s];

    boxDims[s].point = dPoint( xur[s] - xll[s], yur[s] - yll[s] );
    boxDims[s].area  = abs(signedPolyArea(numV, vecPtr(m_xv) + start, vecPtr(m_yv) + start,
                                          counter_cc));
    boxDims[s].index = s;
  }

  // Sort the bounding boxes, this will tell us how to sort the polygons
  sort(boxDims.begin(), boxDims.end(), greaterThanPtIndex);

  // Sort the polygons using auxiliary storage

  vector<double> l_xv           = m_xv;
  vector<double> l_yv           = m_yv;
  vector<int>    l_numVerts     = m_numVerts;
  vector<char>   l_isPolyClosed = m_isPolyClosed;
  vector<string> l_colors       = m_colors;
  vector<string> l_layers       = m_layers;

  for (int s = 0; s < numPolys; s++) {
    int index          = boxDims[s].index;
    m_numVerts     [s] = l_numVerts     [index];
    m_isPolyClosed [s] = l_isPolyClosed [index];
    m_colors       [s] = l_colors       [index];
    m_layers       [s] = l_layers       [index];
  }

  start = 0;
  for (int s = 0; s < numPolys; s++) {

    if (s > 0) start += m_numVerts[s - 1];

    int index = boxDims[s].index;
    int start2 = 0;
    for (int t = 0; t < index; t++) start2 += l_numVerts[t];

    for (int t = 0; t < m_numVerts[s]; t++) {
      m_xv[start + t] = l_xv[start2 + t];
      m_yv[start + t] = l_yv[start2 + t];
    }

  }

}

void dPoly::sortBySizeAndMaybeAddBigContainingRect(// inputs
                                                   double bigXll, double bigYll,
                                                   double bigXur, double bigYur,
                                                   bool counter_cc) {

  // Sort the polygons from largest to smallest by size. If the
  // largest one is going clockwise, it is a hole. In this case, add a
  // large box going counter-clockwise so that we see the hole as a
  // hole in this box. This is important only when polygons are filled
  // (and holes are of background color).

  sortFromLargestToSmallest(counter_cc);

  if (get_numPolys() <= 0) return;

  const double         * xv       = get_xv();
  const double         * yv       = get_yv();
  const int            * numVerts = get_numVerts();
  const vector<string> & colors   = get_colors();
  const vector<string> & layers   = get_layers();

  double signedArea = signedPolyArea(numVerts[0], xv, yv, counter_cc);

  if (signedArea >= 0) return; // Outer poly is correctly oriented

  double xll, yll, xur, yur;
  bdBox(xll, yll, xur, yur);
  bigXll = min(xll, bigXll); bigXur = max(xur, bigXur);
  bigYll = min(yll, bigYll); bigYur = max(yur, bigYur);

  // Add the extra term to ensure we get a box strictly bigger than
  // the polygons and the big window (this will help with sorting below).
  double extra = max(abs(bigXur - bigXll), abs(bigYur - bigYll)) + 10.0;
  bigXll -= extra; bigXur += extra;
  bigYll -= extra; bigYur += extra;

  bool isPolyClosed = true;
  appendRectangle(bigXll, bigYll, bigXur, bigYur, isPolyClosed, colors[0], layers[0]);

  // Reorder the updated set of polygons
  sortFromLargestToSmallest(counter_cc);

  return;
}


void dPoly::enforce45() {

  // Enforce that polygon vertices are integers and the angles are 45x.

  int start = 0;
  for (int pIter = 0; pIter < m_numPolys; pIter++) {

    if (pIter > 0) start += m_numVerts[pIter - 1];

    bool isClosedPolyLine = true;
    int numV    = m_numVerts[pIter];
    double * px = vecPtr(m_xv) + start;
    double * py = vecPtr(m_yv) + start;
    snapPolyLineTo45DegAngles(isClosedPolyLine, numV, px, py);

  }

  return;
};

// This is not used in stereo_gui as it cannot handle mutiple columns.
// See instead read_poly_file() in ASP.
bool dPoly::readPoly(std::string filename,
                     // If isPointCloud is true, treat each point as a
                     // singleton polygon
                     bool isPointCloud) {

  reset();

  m_isPointCloud = isPointCloud;

  // To do: The test below will succeed if filename is a directory.
  // This is probably not right.
  ifstream fh(filename.c_str());
  if (!fh) {
    cerr << "Error: Could not open " << filename << endl;
    return false;
  }

  // The current polygon has vertices in the range [beg, end)
  int beg = 0, end = 0;

  anno annotation;
  string layer, line;
  string color = "yellow"; // default color for polygons

  while (getline(fh, line)) {

    bool isLastLine = (fh.peek() == EOF);
    
    // If line starts with "anno "
    if (line.find("anno ") == 0 && searchForAnnotation(line, annotation)) {
      m_annotations.push_back(annotation);
      continue;
    }

    // Convert to lowercase, after extracting the annotation
    transform(line.begin(), line.end(), line.begin(), ::tolower);

    // If the current line has a color, store it in 'color'. Such a line can
    // start with 'color='. Here we do not just skip the remaining logic below
    // as would like to record the color later.
    if (line.find("color") == 0)
      searchForColor(line, color);

    char * linePtr = (char*)line.c_str(); // To do: Avoid this casting hack.

    // Replace comma with space, to be able to use comma as separator
    for (int s = 0; s < (int)strlen(linePtr); s++) {
      if (linePtr[s] == ',') linePtr[s] = ' ';
    }

    // Ignore any text after the comment character, which is '#' or '!'
    for (int s = 0; s < (int)strlen(linePtr); s++) {
      if (linePtr[s] == '#' || linePtr[s] == '!') {
        for (int t = s; t < (int)strlen(linePtr); t++) {
          linePtr[t] = '\0';
        }
        break;
      }
    }

    // Extract the coordinates of the current vertex and the layer
    // The format we expect is: x y ; layerNo (e.g., 2.3 -1.2 ; 5:16)
    istringstream iss_xy (line);
    double x, y;
    if ( iss_xy >> x >> y ) {

      // This line has valid coordinates, which we read in x and y.
      m_xv.push_back(x);
      m_yv.push_back(y);
      end++;

      if (end == beg + 1) {
        // Find the layer for the current point only if this point
        // is the first point in the polygon
        searchForLayer(line, layer);
      }

    }

    // If this is the last line in the file, or if we encountered a
    // "next" statement, or if we treat a polygon as just a set of
    // points (point cloud) then close the current polygon and start a
    // new one.
    istringstream iss_next(line);
    string val;
    bool isLastVertexOfCurrPoly = ( isLastLine                               ||
                                    ( (iss_next >> val) && (val == "next") ) ||
                                    isPointCloud
                                    );
    bool isCurrPolyNonEmpty = (beg < end);
    if (isLastVertexOfCurrPoly && isCurrPolyNonEmpty) {

      assert( end == (int)m_xv.size() && end == (int)m_yv.size() );

      if (beg < end - 1              &&
          m_xv[beg] == m_xv[end - 1] &&
          m_yv[beg] == m_yv[end - 1]) {
        // The first vertex equals to the last vertex in the current
        // polygon. That means that this is a true polygon rather
        // than a polygonal line. Don't store the last
        // vertex.
        end--;
        m_xv.resize(end);
        m_yv.resize(end);
        m_isPolyClosed.push_back(true);
      }else{
        m_isPolyClosed.push_back(false);
      }

      m_layers.push_back(layer);
      m_colors.push_back(color);

      m_numPolys++;
      m_numVerts.push_back(end - beg);
      m_totalNumVerts = end;

      // Start a new polygon
      beg = end;

    } // End processing the current polygon in the list of polygons

  } // End reading the file and processing all polygons

  return true; // success

}

void dPoly::writePoly(std::string const& filename, std::string const& defaultColor, 
                      bool emptyLineAsSeparator) {

  std::ofstream out(filename.c_str());
  if (!out.is_open()) {
    cerr << "Error: Could not write: " << filename << endl;
    return;
  }

  this->writePoly(out, defaultColor, emptyLineAsSeparator);
}

void dPoly::writePoly(std::ofstream & out, std::string const& defaultColor,
                      bool emptyLineAsSeparator) {

  out.precision(17);
  
  string color = defaultColor, prevColor = defaultColor;

  int start = 0, annoCount = 0;
  for (int pIter = 0; pIter < m_numPolys; pIter++) { // iterate over polygons

    if (m_numVerts[pIter] <= 0) continue; // skip empty polygons
    if (pIter > 0) start += m_numVerts[pIter - 1];

    if (pIter < (int)m_colors.size())
      color = m_colors[pIter];
    if (color != prevColor || pIter == 0)
      out << "color = " << color << endl;
    prevColor = color;

    string layer = "";
    if (pIter < (int)m_layers.size()) layer = m_layers[pIter];
    if (!layer.empty()) layer = " ; " + layer;

    bool isPolyClosed = true;
    if ( pIter < (int)m_colors.size() ) isPolyClosed = m_isPolyClosed[pIter];

    for (int vIter = 0; vIter < m_numVerts[pIter]; vIter++) { // iterate over vertices

      int pos = start + vIter;
      out << m_xv[pos] << ' ' << m_yv[pos] << layer << endl;

      // Put one annotation for each vertex, if possible
      if (annoCount < (int)m_annotations.size() ) {
        m_annotations[annoCount].appendTo(out);
        annoCount++;
      }

    }

    // Turn off the logic below. Do not repeat first vertex at the end. 
    // It is a display preference if the polygon is closed or not.
#if 0 
    if (!m_isPointCloud && isPolyClosed) {
      // Repeat last vertex for closed poly
    
      assert(m_numVerts[pIter] > 0);
     out << m_xv[start] << ' ' << m_yv[start] << layer << endl;
    }
#endif
    
    // Give a choice to use just an empty line as a separator. Note
    // that the reader can only handle the string "NEXT" being a
    // separator.
    if (!m_isPointCloud) {
      if (emptyLineAsSeparator)
        out <<"\n";
      else
        out << "NEXT" << endl;
    }
  }

  // Write the remaining annotations
  for (int a = annoCount; a < (int)m_annotations.size(); a++) {
    m_annotations[a].appendTo(out);
  }

  out.close();

  return;
}

bool dPoly::getColorInCntFile(const std::string & line, std::string & color) {

  // Minor function. Out of the line: "#Color = #ff00" return the string "#ff00"
  // and append "00" to it to make it a valid color.

  istringstream iss(line);
  string equal, colorTag;
  if ( ! (iss >> colorTag >> equal >> color)        ) return false;
  if ( colorTag != "#Color" && colorTag != "#color" ) return false;
  while (color[0] == '#' && color.size() < 7 ) color += "0";

  return true;
}

bool dPoly::read_pol_or_cnt_format(std::string filename,
                                   std::string type,
                                   bool isPointCloud){

  // Read in two very simple and related polygon formats, named pol and cnt.

  assert(type == "pol" || type == "cnt");

  string color = "yellow";
  string layer = "";

  reset();
  m_isPointCloud = isPointCloud;

  ifstream fh(filename.c_str());
  if( !fh ) {
    cerr << "Could not open " << filename << endl;
    return false;
  }

  // Bypass the lines starting with comments. Extract the color.
  while (1) {
    char c = fh.peek();
    if (c != '#' && c != '!') break;
    string line;
    getline(fh, line);
    string lColor;
    if (getColorInCntFile(line, lColor)) color = lColor;
  }

  // Parse the header for pol files.
  double tmp;
  if (type == "pol" && !(fh >> tmp >> tmp >> tmp >> tmp) ) return false;

  while (1) {

    int numVerts = 0;

    if (type == "pol") {
      // Extract the number of vertices for pol files
     if (! (fh >> tmp >> numVerts) ) return true;  // no more vertices
     if (! (fh >> tmp >> tmp) )      return false; // invalid format
    }else{
      // Extract the number of vertices and/or color for cnt file.
      // Skip lines with comments.
      string line;
      if (!getline(fh, line)) return true;
      string lColor; if (getColorInCntFile(line, lColor)) color = lColor;
      if (!line.empty() && line[0] == '#') continue;
      numVerts = int(atof(line.c_str()));
    }

    // Now that we know how many vertices to expect, try reading them
    // from the file. Stop if we are unable to find the expected vertices.

    m_numPolys++;
    m_numVerts.push_back(numVerts);
    m_totalNumVerts += numVerts;
    m_colors.push_back(color);
    m_layers.push_back(layer);

    for (int s = 0; s < numVerts; s++) {
      double x, y;
      if (! (fh >> x >> y) ) return false;
      m_xv.push_back(x);
      m_yv.push_back(y);
    }

    int l = m_totalNumVerts;
    if (l > 0 && numVerts >= 2            &&
        m_xv[l - 1] == m_xv[l - numVerts] &&
        m_yv[l - 1] == m_yv[l - numVerts]
	) {
	// Remove last repeated vertex
	m_xv.pop_back();
	m_yv.pop_back();
	m_totalNumVerts--;
        m_numVerts[m_numVerts.size() - 1]--;
        m_isPolyClosed.push_back(true);
    }else if (type == "pol") {
      // pol files are always closed
      m_isPolyClosed.push_back(true);
    }else{
      m_isPolyClosed.push_back(false);
    }

  }

   return true;
}

void dPoly::set_pointCloud(const std::vector<dPoint> & P, std::string color,
                           std::string layer) {

  // Form a dPoly structure from a set of points
  reset();
  m_isPointCloud = true;
  for (int s = 0; s < (int)P.size(); s++) {
    m_totalNumVerts++;
    m_numPolys++;
    m_numVerts.push_back(1);
    m_layers.push_back(layer);
    m_isPolyClosed.push_back(true);
    m_colors.push_back(color);
    m_xv.push_back(P[s].x);
    m_yv.push_back(P[s].y);
  }

  return;
}

void dPoly::buildGrid(double xl, double yl, double xh, double yh,
                      double gridSize, std::string gridColor) {

  reset();

  if (gridSize <= 0) {
    cerr << "Must have positive grid size." << endl;
    return;
  }

  xl = gridSize*floor(xl/gridSize); xh = gridSize*ceil(xh/gridSize);
  yl = gridSize*floor(yl/gridSize); yh = gridSize*ceil(yh/gridSize);

  int nx = (int)ceil((xh - xl)/gridSize);
  int ny = (int)ceil((yh - yl)/gridSize);

  string layer = "";
  bool isPolyClosed = true;

  for (int i = 0; i <= nx; i++) {
    // Build vertical lines from left to right
    appendRectangle(xl + i*gridSize, yl, xl + i*gridSize, yh,
                    isPolyClosed,
                    gridColor, layer
                    );
  }
  for (int i = 0; i <= ny; i++) {
    // Build horizontal lines from bottom to top
    appendRectangle(xl, yl + i*gridSize, xh, yl + i*gridSize,
                    isPolyClosed,
                    gridColor, layer
                    );
  }

  return;
}

void dPoly::markPolysIntersectingBox(// Inputs
                                     double xll, double yll,
                                     double xur, double yur,
                                     // Outputs
                                     std::map<int, int> & mark) const {

  mark.clear();

  const auto &start_ids = get_startingIndices();

  dRect box(xll, yll, xur, yur);
  int beg_inex = 0;

  dPoly onePoly, clippedPoly;
  for (int polyIndex = 0; polyIndex < m_numPolys; polyIndex++) {

	  if (m_numVerts[polyIndex] == 1){

		  if (box.isInSide(m_xv[beg_inex], m_yv[beg_inex])){
			  mark[polyIndex] = 1;
		  }
	  } else {
		  extractOnePoly(polyIndex, // input
				  onePoly,
				  start_ids[polyIndex]);  // output

		  // A polygon intersects a rectangle if cut to the rectangle
		  // it returns a non-empty polygon
		  onePoly.clipPoly(xll, yll, xur, yur, // inputs
				  clippedPoly);       // outputs
		  if (clippedPoly.get_totalNumVerts() != 0){
			  mark[polyIndex] = 1;
		  }
	  }

	  beg_inex += m_numVerts[polyIndex];

  }

  return;
}

void dPoly::markAnnosIntersectingBox(// Inputs
                                     double xll, double yll,
                                     double xur, double yur,
                                     // Outputs
                                     std::map<int, int> & mark) const {
  mark.clear();

  for (size_t aIter = 0; aIter < m_annotations.size(); aIter++) {
    double x = m_annotations[aIter].x;
    double y = m_annotations[aIter].y;

    if (x >= xll && x <= xur && y >= yll && y <= yur)
      mark[aIter] = 1;
  }  
}

void dPoly::replaceOnePoly(int polyIndex, int numV, const double* x, const double* y) {

  assert(0 <= polyIndex && polyIndex < m_numPolys);

  int start = 0;
  for (int pIter = 0; pIter < polyIndex; pIter++) start += m_numVerts[pIter];

  m_xv.erase(m_xv.begin() + start, m_xv.begin() + start + m_numVerts[polyIndex]);
  m_yv.erase(m_yv.begin() + start, m_yv.begin() + start + m_numVerts[polyIndex]);

  m_xv.insert(m_xv.begin() + start, x, x + numV);
  m_yv.insert(m_yv.begin() + start, y, y + numV);

  m_numVerts[polyIndex] = numV;
  m_totalNumVerts = m_xv.size();

  m_vertIndexAnno.clear();
  m_layerAnno.clear();

  return;
}

void dPoly::eraseMarkedPolys(std::map<int, int> const& mark) {

  // Erase the polygons matching the given mark.
  // See also the function named eraseOnePoly().

  vector<char> dmark, imark;
  dmark.assign(m_totalNumVerts, 0);
  imark.assign(m_numPolys, 0);

  int start = 0;
  for (int pIter = 0; pIter < m_numPolys; pIter++) {
    if (pIter > 0) start += m_numVerts[pIter - 1];
    if (mark.find(pIter) == mark.end()) continue;

    imark[pIter] = 1;
    for (int vIter = 0; vIter < m_numVerts[pIter]; vIter++) dmark[start + vIter] = 1;
  }

  eraseMarkedElements(m_xv, dmark);
  eraseMarkedElements(m_yv, dmark);

  eraseMarkedElements(m_isPolyClosed,  imark);
  eraseMarkedElements(m_colors,        imark);
  eraseMarkedElements(m_layers,        imark);
  eraseMarkedElements(m_numVerts,      imark);

  m_totalNumVerts = m_xv.size();
  m_numPolys      = m_numVerts.size();

  m_vertIndexAnno.clear();
  m_layerAnno.clear();

  return;
}

void dPoly::eraseMarkedAnnos(std::map<int, int> const& mark) {
  vector<char> amark;
  amark.assign(m_annotations.size(), 0);
  for (size_t it = 0; it < m_annotations.size(); it++) {
    if (mark.find(it) == mark.end()) continue;
    amark[it] = 1;
  }

  eraseMarkedElements(m_annotations, amark);
}

void dPoly::erasePolysIntersectingBox(double xll, double yll, double xur, double yur) {

  map<int, int> mark;
  markPolysIntersectingBox(xll, yll, xur, yur, // Inputs
                           mark);              // Outputs
  eraseMarkedPolys(mark);
}

void dPoly::eraseAnnosIntersectingBox(double xll, double yll, double xur, double yur) {
  map<int, int> amark;
  markAnnosIntersectingBox(xll, yll, xur, yur, // Inputs
                           amark);             // Outputs
  eraseMarkedAnnos(amark);

  return;
}

void dPoly::appendAndShiftMarkedPolys(// Inputs
                                      std::map<int, int> & mark,
                                      double shift_x, double shift_y) {

  dPoly polys;

  extractMarkedPolys(mark,    // input
                     polys); // output
  polys.shift(shift_x, shift_y);
  int start = m_numPolys;
  appendPolygons(polys);
  int end = m_numPolys;

  // Remove the mark from the original polygons and mark the newly appended polygons.
  mark.clear();
  for (int s = start; s < end; s++) mark[s] = 1;

  return;
}

}}
