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


/// \file KML.h
///
/// Provides support for writing Google Earth KMLs. At the moment the
/// focus is support for placemarks.
///
#ifndef __VW_FILEIO_KML_H__
#define __VW_FILEIO_KML_H__

// STL
#include <string>
#include <fstream>
#include <stack>

// Vision Workbench
#include <vw/Math/Vector.h>
#include <vw/Math/Quaternion.h>

// Boost
#include <boost/algorithm/string.hpp>
#include <boost/thread/xtime.hpp>
#include <boost/filesystem/fstream.hpp>

namespace vw {

  // TabCount:
  // Used to keep track of tabbing in KML File
  struct TabCount {
    int count;
  };

  std::ostream& operator<<( std::ostream& os, TabCount const& tab);

  // KMLFile:
  // Class wrapper for KML file.
  class KMLFile {
    boost::filesystem::ofstream m_output_file;
    TabCount m_tab;
    std::string m_filename;   // Output filename
    std::string m_name;       // What the file refers to itself as
    std::string m_directory;  // Directory location (used for trees)
    std::stack<std::string> m_bracket_names;
  public:
    KMLFile( std::string filename,
             std::string name,
             std::string directory="" );
    ~KMLFile(); // Closes file out

    // Access internal
    std::string filename () const { return m_filename;  }
    std::string name     () const { return m_name;      }
    std::string directory() const { return m_directory; }

    // Lower Level Writing Functions
    void open_bracket  ( std::string name );
    void close_bracket ();
    void close_brackets( int );
    void close_all_brackets();

    // Low Level Functions
    void enter_folder( std::string name="",
                       std::string description="" );
    void exit_folder();
    
    /// Add a single placemark to the KML file
    void append_placemark( double lon, double lat,
                           std::string name="",
                           std::string description="",
                           std::string style="",
                           double altitude=0,
                           bool extrude=false );
                           
    /// Add a line segment to the KML file
    void append_line( std::vector<Vector3> coordinates,
                      std::string name="",
                      std::string style="");
                           
    /// Add a 3d model to the KML file
    void append_model( std::string path_to_model,
                       double lon, double lat,
                       vw::Quaternion<double> pose,
                       std::string name,
                       std::string description,
                       double altitude,
                       double scale );
                       
    /// LatLonAltBox: This is a bounding box, that only displays contents
    /// when viewer is inside box.
    void append_latlonaltbox( float north,
                              float south,
                              float east,
                              float west );
                              
    /// Lod: This sets the min max pixel viewing range for an object
    void append_lod( float min, float max );
    
    /// Style: Defines an Icon to use later
    /// - If hide_label is set, the placemark labels will only be shown if highlighted.
    void append_style( std::string id, std::string color_hex,
                       float scale, std::string image_url,
                       bool hide_label=false );

    /// Style: Defines a line style to use later
    void append_line_style( std::string id, std::string color_hex,
                            float width=1.0);
   
    /// StyleMap: Maps two styles together to create a bipolar icon
    void append_stylemap( std::string id,
                          std::string style_normal,
                          std::string style_highlight );
    void append_network( std::string link,
                         double north, double south,
                         double east, double west );

    void close_kml(); // If it seems the file wasn't finished, try this.
  protected:
    void open_kml();

  };

  // High Level Tools!
  // Not Existing Yet!
}

#endif//__VW_FILEIO_KML_H__
