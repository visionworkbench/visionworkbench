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


#pragma once

#include <list>
#include <map>
#include <deque>
#include <algorithm>
#include <utility>
#include <vw/Image.h>

/*
 * Point/Vector types and operators
 */

typedef vw::Vector2 ContourPoint;
//typedef vw::Vector3 ContourPoint; // (x,y,t)

/*
inline bool         operator==(ContourPoint a, ContourPoint b);
inline ContourPoint operator+ (ContourPoint a, ContourPoint b);
inline ContourPoint operator- (ContourPoint a, ContourPoint b);
inline float        operator* (ContourPoint a, ContourPoint b);
inline ContourPoint operator* (ContourPoint a, float s);
inline ContourPoint operator* (float s, ContourPoint a);
inline float norm_2(ContourPoint);
inline void normalize(ContourPoint &p);
*/

/*
 * Contour types and functions
 */
struct ContourSegment {
    ContourPoint a;
    ContourPoint b;
    int level;
};

typedef std::list<ContourPoint> PointContour;
typedef std::multimap<int, PointContour > PointContourSet;
typedef std::list<ContourSegment > SegmentList;
typedef vw::math::Vector<vw::Vector2, 4> VWBezierCurve;
typedef std::list<VWBezierCurve> BezierContour;
typedef std::multimap<int, BezierContour > BezierContourSet;

void add_segment(PointContourSet& cset, ContourSegment& s);
bool join(PointContour& c1, PointContour& c2);
vw::Vector2 tangent(const PointContour c, PointContour::iterator p);
void chord_length_parameterize(PointContour c,
        PointContour::iterator first,
        PointContour::iterator last);

/*
 * Contouring functions
 */
void conrec(vw::ImageView<float>& dem, PointContourSet& cset,
            int cint, float nodataval, std::list<ContourSegment>& seglist);


