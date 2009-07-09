
// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
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
#include <fstream>

// Vision Workbench
#include <vw/Core/Log.h>
#include <vw/Camera/BundleAdjust.h>
#include <vw/Camera/CameraModel.h>
#include <vw/Stereo/StereoModel.h>
#include <vw/Cartography/PointImageManipulation.h>
#include <vw/FileIO/KML.h>

// Boost
#include <boost/algorithm/string.hpp>
#include <boost/thread/xtime.hpp>
#include <boost/bind.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/fstream.hpp>

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

  // Small Tools
  enum ReportLevel {
    MinimalReport = 0,
    ClassicReport = 10,
    ReportFile = 20,
    BundlevisFile = 25,
    TriangulationReportAtEnds = 27,
    TriangulationReport = 30,
    BundlevisTriangulation = 35,
    DebugErrorReport = 100,
    DebugJacobianReport = 110
  };

  bool vector_sorting( Vector3 i, Vector3 j);

  void write_kml_styles( KMLFile& kml );
  void write_gcps_kml( KMLFile& kml,
		       ControlNetwork const& cnet );
  void write_3d_est_kml( KMLFile& kml,
			 ControlNetwork const& cnet,
			 std::vector<double>& image_errors );
  void recursive_kml_placemark( KMLFile& kml,
				std::list<Vector3>& list,
				std::string const& name,
				double& min, double& max,
				float& north, float& south,
				float& east, float& west,
				int recursive_lvl );

  // Bundle Adjust Report Code
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
		       int report_lvl=10)  :  m_model(model), m_adjuster(adjuster), bundleadjust_name(name), report_level(report_lvl) {
      
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
	m_human_both << "Bundle Adjustment Parameters:\n";
	m_human_both << "\tNumber Camera params: " << m_model.camera_params_n
		 << std::endl;
	m_human_both << "\tNumber Point params:  " << m_model.point_params_n
		     << std::endl;
	m_human_both << "\tCost Function:        " 
		     << m_adjuster.costfunction().name_tag()
		     << "\t" << m_adjuster.costfunction().threshold() 
		     << std::endl;
	m_human_both << "\tInitial Lambda:       " << m_adjuster.lambda()
		     << std::endl;
	m_human_report << "\tA Inverse Covariance: " 
		       << m_model.A_inverse_covariance(0) << std::endl;
	m_human_report << "\tB Inverse Covariance: " 
		       << m_model.B_inverse_covariance(0) << std::endl;
	if (!m_adjuster.camera_constraint())
	  m_human_report << "\tCamera Constraint shut off!\n"; 
	if (!m_adjuster.gcp_constraint())
	  m_human_report << "\tGCP Constraint shut off!\n";
	m_human_both << "\nStarting Error:\n";
	generic_readings();
      }

      if ( report_level >= TriangulationReportAtEnds )
	triangulation_readings();
      
      // Done
      m_human_both << "\n\n";
    }
    
    // This is a callback from inside the loop of iterations
    void loop_tie_in( void ) {
      m_human_both << "[" << current_posix_time_string() << "]\tFinished Iteration "
		   << m_adjuster.iterations() << std::endl;
      
      if ( report_level >= ClassicReport )
	generic_readings();
	
      if ( report_level >= TriangulationReport )
	triangulation_readings();
  
      m_human_both << "\n";
    }
    // This is a callback for just exit the loop of iterations
    void end_tie_in( void ) { 
      m_human_both << "[" << current_posix_time_string() 
		   << "]\tFinished Bundle Adjustment\n";
      m_human_both << "\tNumber of Iterations: " 
		   << m_adjuster.iterations() << "\n\n";
      if ( report_level >= ClassicReport )
	generic_readings();
      if ( report_level >= TriangulationReportAtEnds )
	triangulation_readings();
      
      // Closing all files out
      m_human_both.remove( m_human_report );
      m_human_both.remove( std::cout );
      
      if (m_human_report.is_open())
	m_human_report.close();
      if (m_bundlevis_binary.is_open())
	m_bundlevis_binary.close();

      if ( report_level >= BundlevisTriangulation ) {
	// Writing an image mean file (TEMPORARY)
	std::cout << "Writing image mean\n";
	std::ofstream image_mean_file;
	image_mean_file.open( "image_mean.err", std::ios::binary | std::ios::out );
	boost::shared_ptr<ControlNetwork> network = m_model.control_network();
	int size1 = network->size();
	image_mean_file.write((char*)&size1, sizeof(int));

	// Right Now getting the data I need
	std::vector<double> image_errors;
	m_model.image_errors(image_errors);
	int size2 = image_errors.size();
	image_mean_file.write((char*)&size2, sizeof(int));

	for ( unsigned i = 0; i < image_errors.size(); i++ )
	  image_mean_file.write((char*)&(image_errors[i]),sizeof(double));
	
	image_mean_file.close();
      }
      
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
      for (unsigned i = 0; i < image_errors.size(); i++ ) {
	mean_image += image_errors[i]/double(image_errors.size());
	stddev_image += image_errors[i]*image_errors[i]/double(image_errors.size());
      }
      //mean_image /= image_errors.size();
      //stddev_image /= image_errors.size();
      stddev_image = sqrt( stddev_image - mean_image*mean_image );
      
      // Statistics gathering for camera_position_errors
      double min_cam_position = *(std::min_element(camera_position_errors.begin(),
						   camera_position_errors.end()));
      double max_cam_position = *(std::max_element(camera_position_errors.begin(),
						   camera_position_errors.end()));
      double mean_cam_position=0, stddev_cam_position=0;
      for (unsigned i = 0; i < camera_position_errors.size(); i++ ) {
	mean_cam_position += camera_position_errors[i];
	stddev_cam_position += camera_position_errors[i]*camera_position_errors[i];
      }
      mean_cam_position /= camera_position_errors.size();
      stddev_cam_position /= camera_position_errors.size();
      stddev_cam_position = sqrt( stddev_cam_position - mean_cam_position*mean_cam_position );
      
      // Statistics gathering for camera_pose_errors
      double min_cam_pose = *(std::min_element(camera_pose_errors.begin(),
					       camera_pose_errors.end()));
      double max_cam_pose = *(std::max_element(camera_pose_errors.begin(),
					       camera_pose_errors.end()));
      double mean_cam_pose=0, stddev_cam_pose=0;
      for (unsigned i = 0; i < camera_pose_errors.size(); i++ ) {
	mean_cam_pose += camera_pose_errors[i];
	stddev_cam_pose += camera_pose_errors[i]*camera_pose_errors[i];
      }
      mean_cam_pose /= camera_pose_errors.size();
      stddev_cam_pose /= camera_pose_errors.size();
      stddev_cam_pose = sqrt( stddev_cam_pose - mean_cam_pose*mean_cam_pose );

      // Statistics gathering for ground control points
      double min_gcp=0,max_gcp=0,mean_gcp=0,stddev_gcp=0;
      if ( gcp_errors.size() != 0 ) {
	min_gcp = *(std::min_element(gcp_errors.begin(),
				     gcp_errors.end()));
	max_gcp = *(std::max_element(gcp_errors.begin(),
				     gcp_errors.end()));
	for (unsigned i = 0; i < gcp_errors.size(); i++ ) {
	  mean_gcp += gcp_errors[i];
	  stddev_gcp += gcp_errors[i]*gcp_errors[i];
	}
	mean_gcp /= gcp_errors.size();
	stddev_gcp /= gcp_errors.size();
	stddev_gcp = sqrt( stddev_gcp - mean_gcp*mean_gcp );
      }
      
      // Sharing now
      m_human_both << "\tLambda: " << m_adjuster.lambda() << std::endl;
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
      if (!m_adjuster.camera_constraint())
	m_human_both << "\t[NOTE: Camera Constraint Shutoff!]\n";
      if (gcp_errors.size() != 0 )
	m_human_both << "\tGCP [min: " << min_gcp << " mean: " << mean_gcp
		     << " max: " << max_gcp << " dev: " << stddev_gcp << "] "
		     << m_model.gcp_unit() << std::endl;
      else
	m_human_both << "\tGCP [N/A]\n";
      if (!m_adjuster.gcp_constraint())
	m_human_both << "\t[NOTE: GCP Constraint Shutoff!]\n";
    }
    
    void stereo_errors( std::vector<double>& stereo_errors ) {
      // Where all the measurement errors will go
      stereo_errors.clear();
      // Pulling out the required models and control network
      std::vector<boost::shared_ptr<camera::CameraModel> > camera_models = m_model.adjusted_cameras();
      boost::shared_ptr<ControlNetwork> network = m_model.control_network();
      
      for (unsigned i = 0; i < network->size(); ++i) {
	for (unsigned j = 0; j+1 < (*network)[i].size(); ++j) {
	  int cam1_index = (*network)[i][j].image_id();
	  int cam2_index = (*network)[i][j+1].image_id();

	  stereo::StereoModel sm( *camera_models[cam1_index],
				  *camera_models[cam2_index] );

	  double error;
	  Vector3 pos = sm( (*network)[i][j].position(),
			    (*network)[i][j+1].position(),
			    error );
	  stereo_errors.push_back(error);
	}
      }
    }
    void stereo_gcp_errors( std::vector<double>& stereo_gcp_errors ) {
      // Where all the measurement errors will go
      stereo_gcp_errors.clear();

      // Pulling out the required models and control network
      std::vector<boost::shared_ptr<camera::CameraModel> > camera_models = m_model.adjusted_cameras();
      boost::shared_ptr<ControlNetwork> network = m_model.control_network();

      for (unsigned i = 0; i < network->size(); ++i) {
	if ((*network)[i].type() == ControlPoint::GroundControlPoint) {
	  for (unsigned j = 0; j+1 < (*network)[i].size(); ++j) {
	    int cam1_index = (*network)[i][j].image_id();
	    int cam2_index = (*network)[i][j+1].image_id();
	    
	    stereo::StereoModel sm( *camera_models[cam1_index],
				    *camera_models[cam2_index] );
	    double error;
	    Vector3 pos = sm( (*network)[i][j].position(),
			      (*network)[i][j+1].position(),
			      error );
	    stereo_gcp_errors.push_back(error);
	  }
	}
      }
    }
    void triangulation_readings( void ) {
      std::vector<double> errors, gcp_errors;
      stereo_errors( errors );
      stereo_gcp_errors( gcp_errors );
      
      // Error for Statistics
      // All points:
      double min_tri=0,max_tri=0,mean_tri=0,stddev_tri=0;
      min_tri = *(std::min_element(errors.begin(),
				   errors.end()));
      max_tri = *(std::max_element(errors.begin(),
				   errors.end()));
      for (unsigned i = 0; i < errors.size(); i++ )
	mean_tri += errors[i];
      mean_tri /= errors.size();
      for (unsigned i = 0; i < errors.size(); i++ )
	stddev_tri += errors[i]*errors[i];
      stddev_tri /= errors.size();
      stddev_tri = sqrt( stddev_tri - mean_tri*mean_tri );
      // GCPs only"
      double min_gcp=0,max_gcp=0,mean_gcp=0,stddev_gcp=0;
      if ( gcp_errors.size() > 0 ) {
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

      // Sharing the information now
      m_human_both << "\tStereo Tri Error [min: " << min_tri 
		   << " mean: " << mean_tri
		   << "\n\t                 max: " << max_tri 
		   << " dev: " << stddev_tri << "]\n";
      if ( gcp_errors.size() > 0 ) 
	m_human_both << "\tStereo GCP Error [min: " << min_gcp 
		     << " mean: " << mean_gcp
		     << "\n\t                 max: " << max_gcp 
		     << " dev: " << stddev_gcp << "]\n";
      else
	m_human_both << "\tStereo GCP Error [N/A]\n";

      if ( report_level >= DebugErrorReport ) {
	m_human_report << "Stereo GCP Debug Error:\n[";
	for ( unsigned i = 0; i < gcp_errors.size(); i++ ) {
	  m_human_report << gcp_errors[i] << ",";
	}
	m_human_report << "]\n";
      }
    }
    // Write KML?
    void write_control_network_kml( bool gcps_only=true ) {
      m_file_prefix = bundleadjust_name;
      boost::to_lower( m_file_prefix );
      boost::replace_all( m_file_prefix, " ", "_" );

      std::ostringstream kml_filename;
      kml_filename << m_file_prefix + ".kml";
      KMLFile kml( kml_filename.str(), 
		   bundleadjust_name,
		   m_file_prefix );
      write_kml_styles( kml );

      boost::shared_ptr<ControlNetwork> network = m_model.control_network();
      
      write_gcps_kml( kml, *network );
      if ( !gcps_only ) {
	std::vector<double> image_errors;
	m_model.image_errors( image_errors );
	write_3d_est_kml( kml, *network, image_errors );
      }
    }

    // Public variables
    std::string bundleadjust_name;
    int report_level;
  };

}} // End namespace

#endif//__VW_CAMERA_BUNDLE_ADJUST_REPORT_H_
