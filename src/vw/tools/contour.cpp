// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/*
 * The conrec function in this file implements the CONREC algorithm,
 * published by Paul Bourke: http://local.wasp.uwa.edu.au/~pbourke/papers/conrec/
 *
 * This implementation of CONREC is based on the C++ implementation by Nicholas Yue:
 * http://local.wasp.uwa.edu.au/~pbourke/papers/conrec/conrec_cxx.txt
 *
 */

#include <list>
#include <map>
#include <algorithm>
#include <utility>
#include <vw/Image.h>
#include <cmath>

#include "contour.h"

#define xsect(p1,p2) (h[p2]*xh[p1]-h[p1]*xh[p2])/(h[p2]-h[p1])
#define ysect(p1,p2) (h[p2]*yh[p1]-h[p1]*yh[p2])/(h[p2]-h[p1])

/*
// Point (and Vector) operators
inline bool         operator==(ContourPoint a, ContourPoint b) {
    return (a[0] == b[0] && a[1] == b[1]) ? true : false;
}
inline ContourPoint operator+ (ContourPoint a, ContourPoint b) {
    return ContourPoint(a[0] + b[0], a[1] + b[1]);
}
inline ContourPoint operator- (ContourPoint a, ContourPoint b) {
    return ContourPoint(a[0] - b[0], a[1] - b[1]);
}
inline float operator* (ContourPoint a, ContourPoint b) {
    return a[0] * b[0] + a[1] * b[1];
}
inline ContourPoint operator* (ContourPoint a, float s) {
    return ContourPoint(a[0] * s, a[1] * s);
}
inline ContourPoint operator* (float s, ContourPoint a) {
    return ContourPoint(a[0] * s, a[1] * s);
}
inline float norm(ContourPoint p) {
    return sqrt(pow(p[0],2) + pow(p[1],2));
}

inline void normalize(ContourPoint &p) {
    float n = norm(p);
    p[0] /= n; p[1] /= n;
}
*/

inline float min_nodata(float a, float b, float nodata) {
    if (a == nodata && b == nodata)
        return nodata;
    else if (a == nodata && b != nodata)
        return b;
    else if (a != nodata && b == nodata)
        return a;
    else return std::min(a,b);
}

inline float max_nodata(float a, float b, float nodata) {
    if (a == nodata && b == nodata)
        return nodata;
    else if (a == nodata && b != nodata)
        return b;
    else if (a != nodata && b == nodata)
        return a;
    else return std::max(a,b);
}

inline bool closed(PointContour& c) {
    return c.front() == c.back();
}

/*
void PointContour::add_point(Point2 p) {
    d.push_back(p);
    u.push_back(0.0);
}
*/

/*
void PointContour::splice(std::deque<Point2>::iterator position, PointContour& c2, bool reverse) {
    assert(position == d.begin() || position == d.end());

    if (reverse) {
        std::deque<Point2>::reverse_iterator c2d_riter;
        std::deque<double>::reverse_iterator c2u_riter;
        if (position == d.begin()) {
            for (c2d_riter = c2.d.rbegin(), c2u_riter = c2.u.rbegin();
                 c2d_riter != c2.d.rend();
                 ++c2d_riter, ++c2u_riter)
            {
                d.push_front(*c2d_riter);
                u.push_front(*c2u_riter);
            }
        } else if (position == d.end()) {
            for (c2d_riter = c2.d.rbegin(), c2u_riter = c2.u.rbegin();
                 c2d_riter != c2.d.rend();
                 ++c2d_riter, ++c2u_riter)
            {
                d.push_back(*c2d_riter);
                u.push_back(*c2u_riter);
            }
        }
    } else {
        std::deque<Point2>::iterator c2d_iter;
        std::deque<double>::iterator c2u_iter;
        if (position == d.begin()) {
            for (c2d_iter = c2.d.begin(), c2u_iter = c2.u.begin();
                 c2d_iter != c2.d.end();
                 ++c2d_iter, ++c2u_iter)
            {
                d.push_front(*c2d_iter);
                u.push_front(*c2u_iter);
            }
        } else if (position == d.end()) {
            for (c2d_iter = c2.d.begin(), c2u_iter = c2.u.begin();
                 c2d_iter != c2.d.end();
                 ++c2d_iter, ++c2u_iter)
            {
                d.push_back(*c2d_iter);
                u.push_back(*c2u_iter);
            }
        }
    }
}
*/

bool join(PointContour& c1, PointContour& c2) {
    ContourPoint c1_front = c1.front();
    ContourPoint c1_back = c1.back();
    ContourPoint c2_front = c2.front();
    ContourPoint c2_back = c2.back();

    if (c1_back == c2_front) {
        c1.splice(c1.end(), c2);
    }
    else if (c1_front == c2_back) {
        c1.splice(c1.begin(), c2);
    }
    else if (c1_back == c2_back) {
        if (c2.size() < c1.size()) {
            c2.reverse();
            c1.splice(c1.end(), c2);
        }
        else {
            c1.reverse();
            c1.splice(c1.begin(), c2);
        }
    }
    else if (c1_front == c2_front) {
        if (c2.size() < c1.size()) {
            c2.reverse();
            c1.splice(c1.begin(), c2);
        }
        else {
            c1.reverse();
            c1.splice(c1.end(), c2);
        }
    }
    else
        return false;

    return true;
}

/*
bool join(PointContour& c1, PointContour& c2) {
    Point2 c1_front = c1.d.front();
    Point2 c1_back = c1.d.back();
    Point2 c2_front = c2.d.front();
    Point2 c2_back = c2.d.back();

    if (c1_back == c2_front)
        c1.splice(c1.d.end(), c2, false);
    else if (c1_front == c2_back)
        c1.splice(c1.d.begin(), c2, false);
    else if (c1_back == c2_back)
        c1.splice(c1.d.end(), c2, true);
    else if (c1_front == c2_front)
        c1.splice(c1.d.begin(), c2, true);
    else
        return false;

    return true;
}
*/

void add_segment(PointContourSet& cset, ContourSegment& s) {
    PointContourSet::iterator iter, iter2, ub_iter;
    int level = s.level;

    PointContour c_new;
    c_new.push_back(s.a);
    c_new.push_back(s.b);

    bool matched = false;

    // for each contour with the same level
    ub_iter = cset.upper_bound(level);
    for (iter = cset.lower_bound(level); iter != ub_iter; iter++) {

        if (closed((*iter).second)) continue;

        if (join((*iter).second, c_new)) {
            matched = true;

            // check for another contour that matches the joined contour
            for (iter2 = cset.lower_bound(level); iter2 != ub_iter; iter2++) {
                // ignore first contour
                if (iter2 == iter) continue;

                //PointContour c2 = (*iter2).second;
                // skip closed contours
                if (closed((*iter2).second)) continue;

                if (join((*iter).second, (*iter2).second)) {
                    cset.erase(iter2);
                    break;
                }
            }
            break;
        }
    }

    if (!matched) {
        // either no contours for this level, or this segment doesn't
        // match any of the existing contours. In either case add a
        // new contour for this level
        cset.insert(make_pair(level, c_new));
    }
}


vw::Vector2 tangent(PointContour c, PointContour::iterator p) {
    vw::Vector2 t;
    PointContour::iterator b = c.begin();
    PointContour::iterator e = --(c.end());
    if (p == b) {
        t = *(++p) - *(--p);
    } else if (p == e) {
        t = *(--p) - *(++p);
    } else {
        vw::Vector2 v1, v2;
        v1 = *(--p) - *(++p);
        v2 = *p - *(++p);
        t = vw::Vector2((v1[0] + v2[0])/2.0, (v1[1] + v2[1])/2.0);
    }
    normalize(t);
    return t;
}

void chord_length_parameterize(PointContour c, PointContour::iterator first, PointContour::iterator last) {
    PointContour::iterator iter;

    (*first)[2] = 0.0;
    iter = first;
    for (++iter; iter != last; ++iter) {
        ContourPoint a = *iter;
        ContourPoint b = *(--iter);

        (*(++iter))[2] = b[2] + norm_2(a - b);
    }

    double len = (*(--last))[2];
    last++;
    iter = first;
    for (++iter; iter != last; ++iter) {
        (*iter)[2] = (*iter)[2]/len;
    }
}

void conrec(vw::ImageView<float>& dem, PointContourSet& cset,
                    int cint, float nodataval,
                    std::list<ContourSegment>& seglist) {
    int m1,m2,m3,case_value;
    double zmin,zmax;
    register int c,i,j,m;
    double h[5];
    int sh[5];
    double xh[5], yh[5];
    int im[4] = {0,1,1,0}, jm[4] = {0,0,1,1};
    int castab[3][3][3] =
        { { {0,0,8},{0,2,5},{7,6,9} },
          { {0,3,4},{1,3,1},{4,3.0} },
          { {9,6,7},{5,2,0},{8,0,0} } };
    ContourSegment seg;

    vw::vw_out(vw::InfoMessage, "console") << "Running CONREC\n";
    vw::vw_out(vw::DebugMessage, "console") << "\tFinding contours\n";
    for (i=0; i < dem.cols()-1; i++) {
        for (j=0; j < dem.rows()-1; j++) {
            zmin = min_nodata( min_nodata(dem(i,j),   dem(i,j+1),   nodataval),
                               min_nodata(dem(i+1,j), dem(i+1,j+1), nodataval),
                               nodataval);
            if (zmin == nodataval) continue;

            zmax = max_nodata( max_nodata(dem(i,j),   dem(i,j+1),   nodataval),
                               max_nodata(dem(i+1,j), dem(i+1,j+1), nodataval),
                               nodataval);

            int cmin = ceil(zmin / cint) * cint;
            int cmax = floor(zmax / cint) * cint;
            for (c = cmin; c <= cmax; c += cint) {
                //printf("(%d,%d) c: %d\n",i,j,c);
                int goodvals = 0;
                h[0] = 0;
                for (m = 4; m >= 0; m--) {
                    if (m > 0) {
                        if (dem(i+im[m-1], j+jm[m-1]) == nodataval)
                            h[m] = nodataval;
                        else {
                            h[m] = dem(i+im[m-1], j+jm[m-1]) - c;
                            h[0] += h[m];
                            goodvals++;
                        }
                        xh[m] = i + im[m-1];
                        yh[m] = j + jm[m-1];
                        //printf("h[%d]: %0.2f (%0.2f)\n",m,dem(i+im[m-1],j+jm[m-1]),h[m]);
                    } else {
                        h[0] /= goodvals;
                        xh[0] = i + 0.5;
                        yh[0] = j + 0.5;
                        //printf("h[%d]: ... (%0.2f)\n",m,h[m]);
                    }


                    if (h[m] > 0.0)
                        sh[m] = 1;
                    else if (h[m] < 0.0)
                        sh[m] = -1;
                    else
                        sh[m] = 0;
                }
                //=================================================================
                //
                // Note: at this stage the relative heights of the corners and the
                // centre are in the h array, and the corresponding coordinates are
                // in the xh and yh arrays. The centre of the box is indexed by 0
                // and the 4 corners by 1 to 4 as shown below.
                // Each triangle is then indexed by the parameter m, and the 3
                // vertices of each triangle are indexed by parameters m1,m2,and
                // m3.
                // It is assumed that the centre of the box is always vertex 2
                // though this isimportant only when all 3 vertices lie exactly on
                // the same contour level, in which case only the side of the box
                // is drawn.
                //
                //
                //      vertex 4 +-------------------+ vertex 3
                //               | \               / |
                //               |   \    m-3    /   |
                //               |     \       /     |
                //               |       \   /       |
                //               |  m=4    X   m=2   |       the centre is vertex 0
                //               |       /   \       |
                //               |     /       \     |
                //               |   /    m=1    \   |
                //               | /               \ |
                //      vertex 1 +-------------------+ vertex 2
                //
                //
                //
                //               Scan each triangle in the box
                //
                //=================================================================
                for (m=1;m<=4;m++) {
                    m1 = m;
                    m2 = 0;
                    m3 = (m==4) ? 1 : m+1;

                    if (h[m1] == nodataval || h[m3] == nodataval)
                        continue;

                    case_value = castab[sh[m1]+1][sh[m2]+1][sh[m3]+1];
                    if (case_value!=0) {
                        ContourSegment seg;
                        seg.level = c;
                        switch (case_value) {
                            //===========================================================
                            //     Case 1 - Line between vertices 1 and 2
                            //===========================================================
                            case 1:
                                seg.a = ContourPoint(xh[m1], yh[m2]);
                                seg.b = ContourPoint(xh[m2], yh[m2]);
                                break;
                            //===========================================================
                            //     Case 2 - Line between vertices 2 and 3
                            //===========================================================
                            case 2:
                                seg.a = ContourPoint(xh[m2], yh[m2]);
                                seg.b = ContourPoint(xh[m3], yh[m3]);
                                break;
                            //===========================================================
                            //     Case 3 - Line between vertices 3 and 1
                            //===========================================================
                            case 3:
                                seg.a = ContourPoint(xh[m3], yh[m3]);
                                seg.b = ContourPoint(xh[m1], yh[m1]);
                                break;
                            //===========================================================
                            //     Case 4 - Line between vertex 1 and side 2-3
                            //===========================================================
                            case 4:
                                seg.a = ContourPoint(xh[m1], yh[m1]);
                                seg.b = ContourPoint(xsect(m2,m3), ysect(m2,m3));
                                break;
                            //===========================================================
                            //     Case 5 - Line between vertex 2 and side 3-1
                            //===========================================================
                            case 5:
                                seg.a = ContourPoint(xh[m2], yh[m2]);
                                seg.b = ContourPoint(xsect(m3,m1), ysect(m3,m1));
                                break;
                            //===========================================================
                            //     Case 6 - Line between vertex 3 and side 1-2
                            //===========================================================
                            case 6:
                                seg.a = ContourPoint(xh[m3], yh[m3]);
                                seg.b = ContourPoint(xsect(m1,m2), ysect(m1,m2));
                                break;
                            //===========================================================
                            //     Case 7 - Line between sides 1-2 and 2-3
                            //===========================================================
                            case 7:
                                seg.a = ContourPoint(xsect(m1,m2), ysect(m1,m2));
                                seg.b = ContourPoint(xsect(m2,m3), ysect(m2,m3));
                                break;
                            //===========================================================
                            //     Case 8 - Line between sides 2-3 and 3-1
                            //===========================================================
                            case 8:
                                seg.a = ContourPoint(xsect(m2,m3), ysect(m2,m3));
                                seg.b = ContourPoint(xsect(m3,m1), ysect(m3,m1));
                                break;
                            //===========================================================
                            //     Case 9 - Line between sides 3-1 and 1-2
                            //===========================================================
                            case 9:
                                seg.a = ContourPoint(xsect(m3,m1), ysect(m3,m1));
                                seg.b = ContourPoint(xsect(m1,m2), ysect(m1,m2));
                                break;
                            default:
                                break;
                        }
                        seglist.push_back(seg);
                    }
                }
            }
        }
    }
    vw::vw_out(vw::DebugMessage, "console")
        << "\tCONREC found " << seglist.size() << " segments" << std::endl;
}


