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

#ifndef VW_GEOMETRY_SHAPEFILE_H
#define VW_GEOMETRY_SHAPEFILE_H

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


  void toOGR(const double * xv, const double * yv, int startPos, int numVerts,
	     OGRLinearRing & R);
  
  void toOGR(vw::geometry::dPoly const& poly, OGRPolygon & P);
  
  void fromOGR(OGRPolygon *poPolygon, std::string const& poly_color,
               std::string const& layer_str, vw::geometry::dPoly & poly);
  
  void fromOGR(OGRLineString *poLineString, std::string const& poly_color,
               std::string const& layer_str, vw::geometry::dPoly & poly);

  void fromOGR(OGRMultiPolygon *poMultiPolygon, std::string const& poly_color,
               std::string const& layer_str, std::vector<vw::geometry::dPoly> & polyVec,
               bool append);

  void fromOGR(OGRGeometry *poGeometry, std::string const& poly_color,
               std::string const& layer_str, std::vector<vw::geometry::dPoly> & polyVec,
               bool append);
  
  void read_shapefile(std::string const& file,
		      std::string const& poly_color,
		      bool & has_geo, 
		      vw::cartography::GeoReference & geo,
		      std::vector<vw::geometry::dPoly> & polyVec);

  void write_shapefile(std::string const& file,
		       bool has_geo,
		       vw::cartography::GeoReference const& geo, 
		       std::vector<vw::geometry::dPoly> const& polyVec);
  
  void shapefile_bdbox(const std::vector<vw::geometry::dPoly> & polyVec,
		       // outputs
		       double & xll, double & yll,
		       double & xur, double & yur);
  
}}

#endif // VW_GEOMETRY_SHAPEFILE_H
