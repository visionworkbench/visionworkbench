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

/// \file BundleAdjustReport.h
/// 
/// Report option for Bundle Adjustment

#ifndef __VW_CAMERA_BUNDLE_ADJUST_REPORT_H_
#define __VW_CAMERA_BUNDLE_ADJUST_REPORT_H_

// STL
#include <string>
#include <iostream>
#include <iomanip>
#include <time.h>

// Vision Workbench
#include <vw/Core/Log.h>
#include <vw/Camera/BundleAdjust.h>

// Boost
#include <boost/algorithm/string.hpp>
#include <boost/thread/xtime.hpp>
// Posix time is not fully supported in the version of Boost for RHEL
// Workstation 4
#ifdef __APPLE__
#include <boost/date_time/posix_time/posix_time.hpp>
#else
#include <ctime>
#endif

// Bundle Adjust Report has many different reporting levels much like
// VW's log. Higher number options include all of the options of the
// numbers below them. Below is the report level numbers:
//    0   - CommandLine Error & Final report
//    10  - CommandLine Iteration Error (classic)
//    20  - Write Report file
//    25  - Write Bundlevis file (binary state information)
//    30  - Write Stereo Triangulation Error Max/Mean/Min
//    35  - Write Bundlevis Stereo Triangulation Error (binary state 
//          information)
//    100 - Write Debug Error Vectors   (big human readable)
//    110 - Write Debug Jacobian Matrix (massive human readable)

namespace vw {
namespace camera { 

  // Posix time is not fully supported in the version of Boost for RHEL
  // Workstation 4

  enum ReportLevel {
    MinimalReport = 0,
    ClassicReport = 10,
    ReportFile = 20,
    BundlevisFile = 25,
    TriangulationReport = 30,
    BundlevisTriangulation = 35,
    DebugErrorReport = 100,
    DebugJacobianReport = 110
  };

  template <class BundleAdjustModelT, class BundleAdjusterT>
  class BundleAdjustReport {
  private:
    multi_ostream m_human_both;
    std::ofstream m_human_report;
    std::ofstream m_bundlevis_binary;
    std::string   m_file_prefix;
    
    BundleAdjustModelT& m_model;
    BundleAdjusterT& m_adjuster;

    #ifndef __APPLE__
    inline std::string current_posix_time_string()  {
      char time_string[2048];
      time_t t = time(0);
      struct tm* time_struct = localtime(&t);
      strftime(time_string, 2048, "%F %T", time_struct);
      return std::string(time_string);
    }
#else
    inline std::string current_posix_time_string() {
      std::ostringstream time_string_stream;
      time_string_stream << boost::posix_time::second_clock::local_time();
      return time_string_stream.str();
    }
#endif
    
  public:
    BundleAdjustReport(std::string const& name, 
		       BundleAdjustModelT& model,
		       BundleAdjusterT& adjuster,
		       int report_lvl=10)  : bundleadjust_name(name), m_model(model), m_adjuster(adjuster), report_level(10) {
      
      m_human_both.add(std::cout);
      
      m_file_prefix = bundleadjust_name;
      boost::to_lower( m_file_prefix );
      boost::replace_all( m_file_prefix, " ", "_" );
      
      // Loading up files
      if ( report_level >= ReportFile ) { 
	// Report file
	std::ostringstream human_file;
	human_file << m_file_prefix + ".report";
	m_human_report.open( human_file.str().c_str(), 
			     std::ios::out );
	m_human_both.add( m_human_report );
	if ( report_level >= BundlevisFile ) {
	  // Bundlevis file
	  std::ostringstream bundlevis_file;
	  bundlevis_file << m_file_prefix + ".bvis";
	  m_bundlevis_binary.open( bundlevis_file.str().c_str(), 
				   std::ios::out );
	}
      }
      
      // Writing Report Intro
      m_human_both << "[" << current_posix_time_string() << "]\tStarted "
		   << bundleadjust_name << " Bundle Adjustment.\n";
      if ( report_level >= ClassicReport ) {
	m_human_both << "\tNumber Camera params: " << m_model.camera_params_n
		 << std::endl;
	m_human_both << "\tNumber Point params:  " << m_model.point_params_n
		     << std::endl;
	m_human_both << "\tCost Function:        " << m_adjuster.costfunction() 
		     << std::endl;
	m_human_both << "\tInitial Lambda:       " << m_adjuster.lambda()
		     << std::endl;
	generic_readings();
      }
      
      // Done
      m_human_both << "\n";
    }
    
    // This is a callback from inside the loop of iterations
    void loop_tie_in( void ) {
      m_human_both << "[" << current_posix_time_string() << "]\tFinished Iteration "
		   << m_adjuster.iterations() << std::endl;
      
      if ( report_level >= ClassicReport ) {
	generic_readings();
	
	if (report_level >= TriangulationReport)
	  triangulation_readings();
      }
  
      m_human_both << "\n";
    }
    // This is a callback for just exit the loop of iterations
    void end_tie_in( void ) { 
      m_human_both << "[" << current_posix_time_string() << "]\tFinished Bundle Adjustment\n";
      m_human_both << "\tNumber of Iterations: " << m_adjuster.iterations() - 1 << std::endl;
      
      // Closing all files out
      m_human_both.remove( m_human_report );
      m_human_both.remove( std::cout );
      
      if (m_human_report.is_open())
	m_human_report.close();
      if (m_bundlevis_binary.is_open())
	m_bundlevis_binary.close();
    }
    // This will display the current error statistics
    std::ostream& operator() ( void ) {
      // These are for special comments to the report
      return m_human_both << "[" << current_posix_time_string() << "] : ";
    }
    
    // Repeated sections of reporting
    void generic_readings( void ) {
      std::vector<double> image_errors, camera_position_errors, camera_pose_errors, gcp_errors;
  
      m_model.image_errors(image_errors);
      m_model.camera_position_errors(camera_position_errors);
      m_model.camera_pose_errors(camera_pose_errors);
      m_model.gcp_errors(gcp_errors);

      // Statistics gathering for image errors
      double min_image = *(std::min_element(image_errors.begin(),
					    image_errors.end()));
      double max_image = *(std::max_element(image_errors.begin(),
					    image_errors.end()));
      double mean_image=0, stddev_image=0;
      for (unsigned i = 0; i < image_errors.size(); i++ ) 
	mean_image += image_errors[i];
      mean_image /= image_errors.size();
      for (unsigned i = 0; i < image_errors.size(); i++ )
	stddev_image += image_errors[i]*image_errors[i];
      stddev_image /= image_errors.size();
      stddev_image = sqrt( stddev_image - mean_image*mean_image );
      
      // Statistics gathering for camera_position_errors
      double min_cam_position = *(std::min_element(camera_position_errors.begin(),
						   camera_position_errors.end()));
      double max_cam_position = *(std::max_element(camera_position_errors.begin(),
						   camera_position_errors.end()));
      double mean_cam_position=0, stddev_cam_position=0;
      for (unsigned i = 0; i < camera_position_errors.size(); i++ ) 
	mean_cam_position += camera_position_errors[i];
      mean_cam_position /= camera_position_errors.size();
      for (unsigned i = 0; i < camera_position_errors.size(); i++ )
	stddev_cam_position += camera_position_errors[i]*camera_position_errors[i];
      stddev_cam_position /= camera_position_errors.size();
      stddev_cam_position = sqrt( stddev_cam_position - mean_cam_position*mean_cam_position );
      
      // Statistics gathering for camera_pose_errors
      double min_cam_pose = *(std::min_element(camera_pose_errors.begin(),
					       camera_pose_errors.end()));
      double max_cam_pose = *(std::max_element(camera_pose_errors.begin(),
					       camera_pose_errors.end()));
      double mean_cam_pose=0, stddev_cam_pose=0;
      for (unsigned i = 0; i < camera_pose_errors.size(); i++ ) 
	mean_cam_pose += camera_pose_errors[i];
      mean_cam_pose /= camera_pose_errors.size();
      for (unsigned i = 0; i < camera_pose_errors.size(); i++ )
	stddev_cam_pose += camera_pose_errors[i]*camera_pose_errors[i];
      stddev_cam_pose /= camera_pose_errors.size();
      stddev_cam_pose = sqrt( stddev_cam_pose - mean_cam_pose*mean_cam_pose );

      // Statistices gathering for ground control points
      double min_gcp=0,max_gcp,mean_gcp=0,stddev_gcp=0;
      if ( gcp_errors.size() != 0 ) {
	min_gcp = *(std::min_element(gcp_errors.begin(),
				     gcp_errors.end()));
	max_gcp = *(std::max_element(gcp_errors.begin(),
				     gcp_errors.end()));
	for (unsigned i = 0; i < gcp_errors.size(); i++ ) 
	  mean_gcp += gcp_errors[i];
	mean_gcp /= gcp_errors.size();
	for (unsigned i = 0; i < gcp_errors.size(); i++ )
	  stddev_gcp += gcp_errors[i]*gcp_errors[i];
	stddev_gcp /= gcp_errors.size();
	stddev_gcp = sqrt( stddev_gcp - mean_gcp*mean_gcp );
      }
      
      // Sharing now
      m_human_both << "\tImage [min: " << min_image << " mean: " << mean_image
		   << " max: " << max_image << " dev: " << stddev_image << "] "
		   << m_model.image_unit() << std::endl;
      m_human_both << "\tCam Position [min: " << min_cam_position << " mean: "
		   << mean_cam_position << " max: " << max_cam_position
		   << " dev: " << stddev_cam_position << "] "
		   << m_model.camera_position_unit() << std::endl;
      m_human_both << "\tCam Pose [min: " << min_cam_pose << " mean: "
		   << mean_cam_pose << " max: " << max_cam_pose << " dev: "
		   << stddev_cam_pose << "] " << m_model.camera_pose_unit()
		   << std::endl;
      if (gcp_errors.size() != 0 )
	m_human_both << "\tGCP [min: " << min_gcp << " mean: " << mean_gcp
		     << " max: " << max_gcp << " dev: " << stddev_gcp << "] "
		     << m_model.gcp_unit() << std::endl;
      else
	m_human_both << "\tGCP [N/A]\n";
    }
    void triangulation_readings( void ) {
      std::cout << "Not implemented";
    }

    // Public variables
    std::string bundleadjust_name;
    int report_level;
  };

}} // End namespace

#endif//__VW_CAMERA_BUNDLE_ADJUST_REPORT_H_
