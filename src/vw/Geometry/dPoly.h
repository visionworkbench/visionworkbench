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

#ifndef VW_GEOMETRY_DPOLY_H
#define VW_GEOMETRY_DPOLY_H

#include <vector>
#include <map>
#include <vw/Geometry/baseUtils.h>
#include <vw/Geometry/geomUtils.h>

namespace vw { namespace geometry {

// A class holding a set of polygons in double precision
class dPoly{

public:

  dPoly(){
    reset();
  }

  void reset();

  bool read_pol_or_cnt_format(std::string filename,
                              std::string type,
                              bool isPointCloud = false
                              );

  bool readPoly(std::string filename,
                bool isPointCloud = false
                );

  void writePoly(std::string filename, std::string defaultColor = "yellow");
  void bdBoxCenter(double & mx, double & my) const;

  void appendPolygon(int numVerts,
                     const double * xv,
                     const double * yv,
                     bool isPolyClosed,
                     const std::string & color,
                     const std::string & layer
                     );

  void appendPolygons(const dPoly & poly);

  void appendRectangle(double xl, double yl, double xh, double yh,
                       bool isPolyClosed,
                       const std::string & color, const std::string & layer
                       );

  void setRectangle(double xl, double yl, double xh, double yh,
                    bool isPolyClosed,
                    const std::string & color, const std::string & layer
                    );

  bool isXYRect();

  void clipPoly(// inputs
                double clip_xll, double clip_yll,
                double clip_xur, double clip_yur,
                dPoly & clippedPoly // output
                );

  void shift(double shift_x, double shift_y);
  void rotate(double angle);
  void scale(double scale);

  const int    * get_numVerts         () const { return vecPtr(m_numVerts); }
  const double * get_xv               () const { return vecPtr(m_xv);       }
  const double * get_yv               () const { return vecPtr(m_yv);       }
  double * get_xv                     () { return vecPtr(m_xv);             }  // non-const
  double * get_yv                     () { return vecPtr(m_yv);             }
  int get_numPolys                    () const { return m_numPolys;         }
  int get_totalNumVerts               () const { return m_totalNumVerts;    }
  std::vector<char> get_isPolyClosed  () const { return m_isPolyClosed;     }
  std::vector<std::string> get_colors () const { return m_colors;           }
  std::vector<std::string> get_layers () const { return m_layers;           }

  void set_color(std::string color);

  void set_isPolyClosed(bool isPolyClosed);

  void eraseMarkedPolys(// Inputs
                        std::map<int, int> & mark
                        );
  void erasePolysIntersectingBox(double xll, double yll, double xur, double yur);
  void appendAndShiftMarkedPolys(// Inputs
                                 std::map<int, int> & mark,
                                 double shift_x, double shift_y
                                 );
  void set_isPointCloud(bool isPointCloud){ m_isPointCloud = isPointCloud; }
  bool isPointCloud() { return m_isPointCloud;}

  void set_pointCloud(const std::vector<dPoint> & P, std::string color,
                      std::string layer);
  void buildGrid(double xl, double yl, double xh, double yh,
                 double gridSize, std::string gridColor);
  void markPolysIntersectingBox(// Inputs
                                double xll, double yll,
                                double xur, double yur,
                                // Outputs
                                std::map<int, int> & mark
                                ) const;
  void replaceOnePoly(int polyIndex, int numV, const double* x, const double* y);
  // Annotations
  void get_annotations (std::vector<anno> & annotations) const;
  void get_layerAnno(std::vector<anno> & annotations) const;
  void get_vertIndexAnno(std::vector<anno> & annotations) const;
  void set_annotations(const std::vector<anno> & A);
  void set_layerAnno(const std::vector<anno> & annotations);
  void set_vertIndexAnno(const std::vector<anno> & annotations);

  void addAnno(const anno & A){m_annotations.push_back(A); }
  void compVertIndexAnno();
  void compLayerAnno();

  void bdBox(double & xll, double & yll, double & xur, double & yur) const;

  void bdBoxes(std::vector<double> & xll, std::vector<double> & yll,
               std::vector<double> & xur, std::vector<double> & yur) const;

  void setPolygon(int numVerts,
                  const double * xv,
                  const double * yv,
                  bool isPolyClosed,
                  const std::string & color,
                  const std::string & layer
                  );

  void eraseAnno(int annoIndex);

  void findClosestAnnotation(// inputs
                             double x0, double y0,
                             // outputs
                             int & annoIndex,
                             double & min_dist
                             ) const;

  void findClosestPolyVertex(// inputs
                             double x0, double y0,
                             // outputs
                             int & polyIndex,
                             int & vertIndex,
                             double & min_x, double & min_y,
                             double & min_dist
                             ) const;

  void findClosestPolyEdge(//inputs
                           double x0, double y0,
                           // outputs
                           int & polyIndex, int & vertIndex,
                           double & minX, double & minY, double & minDist
                           ) const;

  void eraseOnePoly(int polyIndex);
  void insertVertex(int polyIndex, int vertIndex,
                    double x, double y);
  void eraseVertex(int polyIndex, int vertIndex);
  void changeVertexValue(int polyIndex, int vertIndex, double x, double y);
  void shiftEdge(int polyIndex, int vertIndex, double shift_x, double shift_y);
  void shiftOnePoly(int polyIndex, double shift_x, double shift_y);
  void shiftMarkedPolys(const std::map<int, int> & mark, double shift_x, double shift_y);
  void extractOnePoly(int polyIndex, // input
                      dPoly & poly   // output
                      ) const;
  void extractMarkedPolys(std::map<int, int> & mark, // input
                          dPoly & polys              // output
                          ) const;
  void reverse(); // reverse orientation
  void reverseOnePoly(int polyIndex);
  void sortFromLargestToSmallest();

  void sortBySizeAndMaybeAddBigContainingRect(// inputs
                                              double bigXll, double bigYll,
                                              double bigXur, double bigYur);

  void enforce45();

private:

  bool getColorInCntFile(const std::string & line, std::string & color);
  void get_annoByType(std::vector<anno> & annotations, int annoType);
  void set_annoByType(const std::vector<anno> & annotations, int annoType);

  // If isPointCloud is true, treat each point as a set of unconnected points
  bool                     m_isPointCloud;

  std::vector<double>      m_xv;
  std::vector<double>      m_yv;
  std::vector<int>         m_numVerts;
  int                      m_numPolys;
  int                      m_totalNumVerts;
  std::vector<char>        m_isPolyClosed;
  std::vector<std::string> m_colors;
  std::vector<std::string> m_layers;
  std::vector<anno>        m_annotations;
  std::vector<anno>        m_vertIndexAnno; // Anno showing vertex index
  std::vector<anno>        m_layerAnno;     // Anno showing layer number

};

}}

#endif // VW_GEOMETRY_DPOLY_H
