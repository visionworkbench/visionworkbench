// __BEGIN_LICENSE__
//  Copyright (c) 2006-2025, United States Government as represented by the
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

#ifndef VW_CARTOGRAPHY_SHAPEFILE_H
#define VW_CARTOGRAPHY_SHAPEFILE_H

#include <ogrsf_frmts.h>

#include <vw/Math/BBox.h>
#include <vw/Math/Vector.h>
#include <vw/Cartography/GeoReference.h>
#include <vw/Cartography/GeoReferenceUtils.h>
#include <vw/Geometry/dPoly.h>

#include <vector>
#include <map>

// Shape file functions

namespace vw { namespace geometry {

// Convert a single point to OGRPoint
void toOGR(double x, double y, OGRPoint & P);

// Convert a polygonal line to OGRLineString
void toOGR(const double * xv, const double * yv, int startPos, int numVerts,
      OGRLineString & L);

// Convert a single polygon to OGRLinearRing
void toOGR(const double * xv, const double * yv, int startPos, int numVerts,
      OGRLinearRing & R);

// Convert a polygon set to OGRPolygon
void toOGR(vw::geometry::dPoly const& poly, OGRPolygon & P);

// Convert a set of 3D points to OGRPolygon. Ignore the third coordinate.
void toOGR(std::vector<vw::Vector3> const& points, OGRPolygon & P);

// Read a polygon set from OGR
void fromOGR(OGRPolygon *poPolygon, std::string const& poly_color,
              std::string const& layer_str, vw::geometry::dPoly & poly);

// Read an OGR polygonal line (line string)
void fromOGR(OGRLineString *poLineString, std::string const& poly_color,
              std::string const& layer_str, vw::geometry::dPoly & poly);

// Read polygons from OGR multipolygon
void fromOGR(OGRMultiPolygon *poMultiPolygon, std::string const& poly_color,
              std::string const& layer_str, std::vector<vw::geometry::dPoly> & polyVec,
              bool append);

// Read polygons from OGR geometry
void fromOGR(OGRGeometry *poGeometry, std::string const& poly_color,
              std::string const& layer_str, std::vector<vw::geometry::dPoly> & polyVec,
              bool append);

// Read a shapefile
void read_shapefile(std::string const& file,
        std::string const& poly_color,
        bool & has_geo,
        vw::cartography::GeoReference & geo,
        std::vector<vw::geometry::dPoly> & polyVec);

// Write a set of dPolys to a shapefile
void write_shapefile(std::string const& file,
          bool has_geo,
          vw::cartography::GeoReference const& geo,
          std::vector<vw::geometry::dPoly> const& polyVec);

// Write a single dPoly to a shapefile (this can still have multiple polygons)
void write_shapefile(std::string const& file,
          bool has_geo,
          vw::cartography::GeoReference const& geo,
          vw::geometry::dPoly const& poly);

// Bounding box of a shapefile
void shapefile_bdbox(const std::vector<vw::geometry::dPoly> & polyVec,
          // outputs
          double & xll, double & yll,
          double & xur, double & yur);

// Merge polygons into polyVec. The inputs are in ogr_polys, but they be
// changed by this function. The output is in polyVec.
void mergeOGRPolygons(std::string const& poly_color,
                      std::string const& layer_str,
                      std::vector<OGRGeometry*>& ogr_polys,
                      std::vector<vw::geometry::dPoly>& polyVec); // output

// Find the convex hull of a set of 3D points to OGRPolygon. Ignore the third coordinate.
void convexHull(std::vector<vw::Vector3> const& points,
                std::vector<vw::geometry::dPoly> & polyVec); // output

}} // end namespace vw::geometry

#endif // VW_CARTOGRAPHY_SHAPEFILE_H
