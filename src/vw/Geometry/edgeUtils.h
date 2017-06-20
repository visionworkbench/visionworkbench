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

#ifndef VW_GEOMETRY_EDGEUTILS_H
#define VW_GEOMETRY_EDGEUTILS_H

namespace vw { namespace geometry {
  
  bool edgeIntersectsBox(// Input: arbitrary edge
                         double bx, double by,
                         double ex, double ey,
                         // Input: Box
                         double xl, double yl,
                         double xh, double yh
                         );
  bool edgeIntersectsHorizontalEdge(// Input: arbitrary edge
                                    double x0, double y0,
                                    double x1, double y1,
                                    // Input: horizontal edge
                                    double begx, double endx,
                                    double yval
                                    );
  bool isPointOnEdge(double x0, double y0, double x1, double y1,
                     double x, double y);
  bool collinearEdgesIntersect(// Input: first edge
                                    double ax0, double ay0,
                                    double ax1, double ay1,
                                    // Input: second edge
                                    double bx0, double by0,
                                    double bx1, double by1,
                                    // Output: intersection
                                    // if it exists
                                    double & x, double & y
                                    );
  bool edgesIntersect(// Input: first edge
                      double ax0, double ay0,
                      double ax1, double ay1,
                      // Input: second edge
                      double bx0, double by0,
                      double bx1, double by1,
                      // Output: intersection if it exists
                      double & x, double & y
                      );
  void cutEdge(double x0, double y0, double x1, double y1,
               double nx, double ny, double H,
               double & cutx, double & cuty);
  
}}

#endif // VW_GEOMETRY_EDGEUTILS_H

