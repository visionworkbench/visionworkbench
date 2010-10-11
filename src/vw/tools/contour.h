// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
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


