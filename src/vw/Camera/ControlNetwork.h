// __BEGIN_LICENSE__
// 
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
// 
// Copyright 2006 Carnegie Mellon University. All rights reserved.
// 
// This software is distributed under the NASA Open Source Agreement
// (NOSA), version 1.3.  The NOSA has been approved by the Open Source
// Initiative.  See the file COPYING at the top of the distribution
// directory tree for the complete NOSA document.
// 
// THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
// KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
// LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
// SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
// A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
// THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
// DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
// 
// __END_LICENSE__

/// \file BundleAdjust.h
/// 
/// Optimization classes for carrying out bundle adjustment of many
/// camera images.

#ifndef __VW_CAMERA_CONTROL_NETWORK_H__
#define __VW_CAMERA_CONTROL_NETWORK_H__

// STL
#include <string>
#include <vector>
#include <fstream>
#include <iomanip>

#include <vw/Math/Vector.h>

namespace vw {
namespace camera {

  /// A ControlMeasure identifies a pixel in an image that corresponds
  /// to a control point.  In addition to the location of the pixel, the
  /// control measure also stores the uncertainty of the measurement,
  /// and a identifier for the image from which it was derived.
  class ControlMeasure {
    float m_col, m_row, m_col_sigma, m_row_sigma;
    std::string m_date_time, m_description;
    double m_ephemeris_time;
    int m_image_id;
    bool m_ignore;

  public:

    // Constructor
    ControlMeasure( float col, float row, float col_sigma, float row_sigma, int image_id) :
      m_col(col), m_row(row), m_col_sigma(col_sigma), m_row_sigma(row_sigma), m_image_id(image_id) {

      m_date_time = "not recorded";
      m_description = "none";
      m_ignore = false;

      m_ephemeris_time = 0;
    }

    ControlMeasure() : 
      m_col(0), m_row(0), m_col_sigma(0), m_row_sigma(0), m_image_id(0) {

      m_date_time = "not recorded";
      m_description = "none";
      m_ignore = false;

      m_ephemeris_time = 0;
    }

    /// Return the pixel location for this measurement.
    Vector2 position() const { return Vector2(m_col, m_row); }

    /// Set the pixel location for this measurement.
    void set_position(float col, float row) { 
      m_col = col;
      m_row = row;
    }

    /// Set the pixel location for this measurement.
    void set_position(Vector2 position) { 
      m_col = position[0];
      m_row = position[1];
    }

    /// Return the measurement error for this point.
    Vector2 sigma() const { return Vector2(m_col_sigma, m_row_sigma); }

    /// Set the measurement error for this point.
    void set_sigma(float col_sigma, float row_sigma) {
      m_col_sigma = col_sigma;
      m_row_sigma = row_sigma;
    }

    /// Set the measurement error for this point.
    void set_sigma(Vector2 sigma) { 
      m_col_sigma = sigma[0]; 
      m_row_sigma = sigma[1];
    }
 
    /// Return the identifier for the image from which this control
    /// point was derived.  
    int image_id() const { return m_image_id; }

    /// Set the identifier for the image from which this control
    /// point was derived.  Any string value can be used here, but it
    /// should be something that can be used to locate the original
    /// image.
    void set_image_id(int image_id) { m_image_id = image_id; }


    /// Return the description of this control measure.  
    std::string description() const { return m_description; }

    /// Set the description of this control measure.  This typically
    /// includes a description of the source imagery used and the person
    /// or software that was used to make the measurement.
    void set_description(std::string const& description) { m_description = description; }


    /// Query when this control measurement was made.
    std::string date_time() const { return m_date_time; }

    /// Set when this control measurement was made.
    void set_date_time(std::string const& date_time) { m_date_time = date_time; }

    /// Query whether this control measurement should be ignored in a bundle adjustment.
    bool ignore() const { return m_ignore; }

    /// Set whether this control measurement should be ignored in a bundle adjustment.
    void set_ignore(bool state) { m_ignore = state; }
     
    /// Query what was the ephemeris time of when this pixel was taken
    double ephemeris_time() const { return m_ephemeris_time; }
    
    /// Set the ephemeris time of when this pixel was taken
    void set_ephemeris_time( double const& time ) { m_ephemeris_time = time; }

  };

  /// Two control meaures are considered equal if their position,
  /// sigma, and image_id are equal.
  inline bool operator==( ControlMeasure const& m1, ControlMeasure const& m2 ) {
    if (m1.position() == m2.position() &&
        m1.sigma() == m2.sigma() &&
        m1.image_id() == m2.image_id() && 
	m1.ephemeris_time() == m2.ephemeris_time() ) 
      return true;
    //else
    return false;
  }

  std::ostream& operator<<( std::ostream& os, ControlMeasure const& measure);

  /// ControlPoints are 3D locations in geographic coordinates (lon,
  /// lat, radius) that are associated with a certain number of
  /// ControlMeasures.  Each ControlMeasure is an observation where
  /// this ControlPoint was located in an image.
  class ControlPoint {

    std::string m_id;
    std::vector<ControlMeasure> m_measures;
    bool m_ignore;
    Vector3 m_position;
    Vector3 m_sigma;

  public:
    typedef std::vector<ControlMeasure>::iterator iterator;
    typedef std::vector<ControlMeasure>::const_iterator const_iterator;
    
    enum ControlPointType { GroundControlPoint, TiePoint };
    ControlPointType m_type;

    /// Constructor 
    ControlPoint(ControlPointType type = ControlPoint::TiePoint) : m_type(type) {
      m_ignore = false;
      m_id = "not set";
    }

    ControlPointType type() const { return m_type; }
    void set_type(ControlPointType type) { m_type = type; }
    
    // Iterators
    iterator begin() { return m_measures.begin(); }
    iterator end() { return m_measures.end(); }
    const_iterator begin() const { return m_measures.begin(); }
    const_iterator end() const { return m_measures.end(); }

    std::string id() const { return m_id; }
    void set_id(std::string id) { m_id = id; }
    
    bool ignore() const { return m_ignore; }
    void set_ignore(bool ignore) { m_ignore = ignore; }

    /// Returns the number of control measures associated with this
    /// control point.
    unsigned size() const { return m_measures.size(); }

    /// Associate a single control measure with this ControlPoint
    void add_measure(ControlMeasure const& measure) { m_measures.push_back(measure); }

    /// Associate multiple control measures with this ControlPoint
    void add_measures(std::vector<ControlMeasure> const& measures) { 
      m_measures.insert(m_measures.end(), measures.begin(), measures.end());
    }

    /// Remove the control point at the specified index.
    void delete_measure(unsigned index) {
      if (index >= this->size())
        vw_throw(LogicErr() << "ControlPoint::delete_control_point -- index " << index << " exceeds control point dimensions.");

      iterator iter = this->begin();
      for (unsigned i=0; i<index; ++i)
        ++iter;
      
      m_measures.erase(iter);
    }
    
    /// Access a specific control measure that is associated with this control point.
    ControlMeasure& operator[] (int index) { return m_measures[index]; }
    const ControlMeasure& operator[] (int index) const { return m_measures[index]; }
    
    /// Clear all the control measures associated with this control point.
    void clear() { m_measures.clear(); }

    /// Locate a control measure that is equal to the query.  Returns
    /// this->size() if no match is found.
    unsigned find(ControlMeasure const& query) {
      for (unsigned i = 0; i < m_measures.size(); ++i) 
        if (m_measures[i] == query)  
          return i;
      // If no matches are found, return m_measures.size() (the last index + 1)
      return m_measures.size();
    }

    /// Set the position of the control point
    void set_position(double lon, double lat, double radius) {
      m_position = Vector3(lon,lat,radius);
    }

    /// Set the position of the control point as [lon, lat, radius]
    void set_position(Vector3 position) { m_position = position; }
    
    /// Returns the position of the control point as [longitude,
    /// latitude, radius]
    Vector3 position() const { return m_position; }

    /// Set the uncertainty of the control point
    void set_sigma(double lon_sigma, double lat_sigma, double radius_sigma) {
      m_sigma = Vector3(lon_sigma, lat_sigma, radius_sigma);
    }

    /// Set the uncertainty of the control point
    void set_sigma(Vector3 sigma) { m_sigma = sigma; }

    /// Returns the uncertainty of the control point as [longitude_sigma,
    /// latitude_sigma, radius_sigma]
    Vector3 sigma() const { return m_sigma; }

  };

  std::ostream& operator<<( std::ostream& os, ControlPoint const& point);

  /// The control network contains a list of control points (either ground control
  /// points or tie points).
  ///
  /// Things to add:
  /// - read/write controlnet file
  /// - various methods for computing error
  /// - assoc. with image list/serial number
  ///
  class ControlNetwork {
    std::string m_description;
    std::string m_target_body;
    std::vector<ControlPoint> m_control_points;

  public:

    typedef std::vector<ControlPoint>::iterator iterator;
    typedef std::vector<ControlPoint>::const_iterator const_iterator;

    /// Constructor
    ControlNetwork(std::string description, std::string target_body = "Unknown") :
      m_description(description), m_target_body(target_body) {}

    /// Returns the number of control measures associated with this
    /// control point.
    unsigned size() const { return m_control_points.size(); }

    /// Return the number of Control Points that are Ground Control Points (GCPs)
    unsigned num_ground_control_points() const { 
      unsigned count=0;
      for (unsigned i=0; i<this->size(); ++i) {
        if ((*this)[i].type() == ControlPoint::GroundControlPoint)
          ++count;
      }
      return count;
    }

    unsigned num_3d_tie_points() const { 
      unsigned count=0;
      for (unsigned i=0; i<this->size(); ++i) {
        if ((*this)[i].type() == ControlPoint::TiePoint)
          ++count;
      }
      return count;
    }

    // Iterators
    iterator begin() { return m_control_points.begin(); }
    iterator end() { return m_control_points.end(); }
    const_iterator begin() const { return m_control_points.begin(); }
    const_iterator end() const { return m_control_points.end(); }

    /// Associate a single control measure with this ControlPoint
    void add_control_point(ControlPoint const& point) { m_control_points.push_back(point); }

    /// Associate multiple control measures with this ControlPoint
    void add_control_points(std::vector<ControlPoint> const& points) { 
      m_control_points.insert(m_control_points.end(), points.begin(), points.end());
    }

    /// Remove the control point at the specified index.
    void delete_control_point(unsigned index) {
      if (index >= this->size())
        vw_throw(LogicErr() << "ControlNetwork::delete_control_point -- index " << index << " exceeds control network dimensions.");

      iterator iter = this->begin();
      for (unsigned i=0; i<index; ++i)
        ++iter;
      
      m_control_points.erase(iter);
    }
    
    /// Access a specific control measure that is associated with this control point.
    ControlPoint& operator[] (int index) { return m_control_points[index]; }
    const ControlPoint& operator[] (int index) const { return m_control_points[index]; }

    /// Clear all the control points associated with this control network.
    void clear() { m_control_points.clear(); }

    /// Locate a control point that contains the control measure that
    /// is equal to the query.  Returns this->size() if no match is
    /// found.
    unsigned find_measure(ControlMeasure const& query) {
      for (unsigned i = 0; i < m_control_points.size(); ++i) 
        if (m_control_points[i].find(query) != m_control_points[i].size()) 
          return i;
      // Otherwise...
      return m_control_points.size();
    }

    /// Write the current control network to a file on disk.
    void write_control_network(std::string filename) {
      std::ofstream ofile(filename.c_str());
      ofile << this->size() << "\n";
      
      for (unsigned c=0; c < this->size(); ++c) {
        ControlPoint& cpoint = (*this)[c];
        ofile << std::setprecision(18) << cpoint.size() << " " << cpoint.position()[0] << " " << cpoint.position()[1] << " " << cpoint.position()[2] << " " << cpoint.sigma()[0] << " " << cpoint.sigma()[1]<< " " << cpoint.sigma()[2] << " " << cpoint.ignore() << "\n";
        for (unsigned m = 0; m < cpoint.size(); ++m) {
          ControlMeasure& cmeasure = cpoint[m];
          ofile << std::setprecision(18) << cmeasure.image_id() << " " << cmeasure.position()[0] << " " << cmeasure.position()[1] << " " << cmeasure.sigma()[0] << " " << cmeasure.sigma()[1] << " " << cmeasure.ignore() << " " <<  cmeasure.ephemeris_time() << "\n";
        }
      }
      ofile.close();
    }

    /// Add the contents of the file on disk to the existing data in
    /// this control network.
    void read_control_network(std::string filename) {
      std::ifstream ifile(filename.c_str());

      int total_control_points;
      ifile >> total_control_points;

      std::vector<ControlPoint> cpoints(total_control_points);
      
      unsigned num_measures;
      Vector3 pos;
      Vector3 pos_sigma;
      bool ignore_point;
      for (int c=0; c < total_control_points; ++c) {
        ifile >> num_measures >> pos[0] >> pos[1] >> pos[2] >> pos_sigma[0] >> pos_sigma[1] >> pos_sigma[2] >> ignore_point;

        cpoints[c].set_position(pos);
        cpoints[c].set_sigma(pos_sigma);
        cpoints[c].set_ignore(ignore_point);

        std::vector<ControlMeasure> cmeasures(num_measures);
        int image_id;
        Vector2 measure_pos, measure_sigma;
	double ephemeris;
        bool ignore_measure;
        for (unsigned m=0; m < num_measures; ++m) {
          ifile >> image_id >> measure_pos[0] >> measure_pos[1] >> measure_sigma[0] >> measure_sigma[1] >> ignore_measure >> ephemeris;
          cmeasures[m].set_image_id(image_id);
          cmeasures[m].set_position(measure_pos);
          cmeasures[m].set_sigma(measure_sigma);
          cmeasures[m].set_ignore(ignore_measure);
	  cmeasures[m].set_ephemeris_time(ephemeris);
        }
        cpoints[c].add_measures(cmeasures);
      }
      this->add_control_points(cpoints);
    }

  };

  std::ostream& operator<<( std::ostream& os, ControlNetwork const& cnet);

}} // namespace vw::camera

#endif // __VW_CAMERA_CONTROL_NETWORK_H__
