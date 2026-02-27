// __BEGIN_LICENSE__
//  Copyright (c) 2006-2026, United States Government as represented by the
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

#include <vw/BundleAdjustment/ControlNetwork.h>
#include <vw/Core/Log.h>

#include <boost/foreach.hpp>

namespace vw {
namespace ba {

  // ControlMeasure

  ControlMeasure::ControlMeasure(float col, float row,
                                  float col_sigma, float row_sigma,
                                  uint64 image_id,
                                  ControlMeasureType type):
    m_col(col), m_row(row), m_col_sigma(col_sigma),
    m_row_sigma(row_sigma), m_image_id(image_id),
    m_ignore(false), m_type(type) {
  }

  ControlMeasure::ControlMeasure(ControlMeasureType type):
    m_col(0), m_row(0), m_col_sigma(0), m_row_sigma(0),
    m_image_id(0), m_ignore(false), m_type(type) {
  }

  std::string ControlMeasure::get_image_name(ControlNetwork const& net) const {
    if (net.get_image_list().size() <= m_image_id)
      return "";
    return net.get_image_list()[m_image_id];
  }

  // ControlPoint

  ControlPoint::ControlPoint(ControlPointType type): m_type(type) {
    m_ignore = false;
    m_id = "Null";
    m_sigma = Vector3(1.0, 1.0, 1.0);
  }

  void ControlPoint::add_measure(ControlMeasure const& measure) {
    m_measures.push_back(measure);
  }

  void ControlPoint::add_measures(std::vector<ControlMeasure> const& measures) {
    m_measures.insert(m_measures.end(), measures.begin(), measures.end());
  }

  void ControlPoint::delete_measure(size_t index) {
    if (index >= this->size())
      vw_throw(LogicErr() << "ControlPoint::delete_measure -- index "
               << index << " exceeds control point dimensions.");

    iterator iter = this->begin();
    iter += index;
    m_measures.erase(iter);
  }

  size_t ControlPoint::find(ControlMeasure const& query) {
    for (size_t i = 0; i < m_measures.size(); ++i)
      if (m_measures[i] == query)
        return i;
    return m_measures.size();
  }

  // ControlNetwork

  ControlNetwork::ControlNetwork(std::string id,
                                 ControlNetworkType type,
                                 std::string target_name,
                                 std::string descrip,
                                 std::string user_name):
    m_targetName(target_name), m_networkId(id),
    m_description(descrip), m_userName(user_name), m_type(type) {
  }

  void ControlNetwork::add_control_point(ControlPoint const& point) {
    if (m_type != ControlNetwork::ImageToGround &&
         point.type() == ControlPoint::GroundControlPoint)
      m_type = ControlNetwork::ImageToGround;

    m_control_points.push_back(point);
  }

  void ControlNetwork::add_control_points(std::vector<ControlPoint> const& points) {
    if (m_type != ControlNetwork::ImageToGround) {
      BOOST_FOREACH(ControlPoint const& cp, points)
        if (cp.type() == ControlPoint::GroundControlPoint) {
          m_type = ControlNetwork::ImageToGround;
          break;
        }
    }

    m_control_points.insert(m_control_points.end(),
                            points.begin(), points.end());
  }

  void ControlNetwork::delete_control_point(size_t index) {
    if (index >= this->size())
      vw_throw(LogicErr()
               << "ControlNetwork::delete_control_point -- index "
               << index << " exceeds control network dimensions.");

    iterator iter = this->begin();
    iter += index;
    m_control_points.erase(iter);
  }

  size_t ControlNetwork::find_measure(ControlMeasure const& query) {
    for (size_t i = 0; i < m_control_points.size(); ++i)
      if (m_control_points[i].find(query) != m_control_points[i].size())
        return i;
    return m_control_points.size();
  }

  void ControlNetwork::write_in_gcp_format(std::string const& filename,
                                           cartography::Datum const& d) const {

    const std::string UNSPECIFIED_DATUM = "unspecified_datum";
    if (d.name() == UNSPECIFIED_DATUM)
      vw_throw(ArgumentErr()
               << "FATAL: No datum was specified. "
               << "Cannot save control network as csv.\n");

    std::ofstream ofs(filename.c_str());
    ofs.precision(17);

    // It is important to keep track of the datum
    // because the elevations are relative to it
    ofs << "# WKT: " << d.get_wkt() << "\n";

    int count = 0;
    for (const_iterator iter = begin(); iter != end(); iter++) {
      count++;

      // lon,lat,height -> lat,lon,height
      Vector3 llr = d.cartesian_to_geodetic((*iter).position());
      std::swap(llr[0], llr[1]);

      Vector3 sigma = (*iter).sigma();
      for (size_t ipt = 0; ipt < sigma.size(); ipt++)
        if (sigma[ipt] <= 0)
          sigma[ipt] = 1;

      ofs << count    << ' '
          << llr  [0] << ' ' << llr  [1] << ' ' << llr[2] << ' ';
      ofs << sigma[0] << ' ' << sigma[1] << ' ' << sigma[2] << ' ';

      for (ControlPoint::const_iterator measure = (*iter).begin();
            measure != (*iter).end(); ++measure) {

        if (measure->image_id() >= m_image_names.size())
          vw_throw(ArgumentErr()
                   << "ControlNetwork::write_in_gcp_format: "
                   << "Measure image ID " << measure->image_id()
                   << " exceeds image count: "
                   << m_image_names.size());

        std::string image_name = m_image_names[measure->image_id()];

        ofs << image_name << ' '
            << measure->position()[0] << ' '
            << measure->position()[1] << ' '
            << measure->sigma()[0]    << ' '
            << measure->sigma()[1];

        if (measure+1 != (*iter).end())
          ofs << ' ';
        else
          ofs << "\n";
      }
    }
    ofs.close();
  }

  // ostream operators

  std::ostream& operator<<(std::ostream& os,
                            ControlMeasure const& measure) {
    os << measure.image_id() << ":" << measure.position();
    return os;
  }

  std::ostream& operator<<(std::ostream& os,
                            ControlPoint const& point) {
    os << "[Control Point: " << point.position() << "] ";
    BOOST_FOREACH(ControlMeasure const& cm, point)
      os << cm << " ";
    os << "\n";
    return os;
  }

  std::ostream& operator<<(std::ostream& os,
                            ControlNetwork const& cnet) {
    os << "Control Network: " << cnet.size() << " points.\n";
    BOOST_FOREACH(ControlPoint const& cp, cnet)
      os << "\t" << cp;
    os << "\n";
    return os;
  }

}} // namespace vw::ba
