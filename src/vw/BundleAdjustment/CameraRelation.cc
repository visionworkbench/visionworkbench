// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file CameraRelation.cc
///

#include <vw/BundleAdjustment/CameraRelation.h>

namespace vw {
namespace ba {

  // Implementations for specific features
  ControlMeasure
  IPFeature::control_measure() const {
    return ControlMeasure( m_ip.x, m_ip.y,
                           m_ip.scale, m_ip.scale,
                           m_camera_id );
  }

  ControlMeasure
  JFeature::control_measure() const {
    return ControlMeasure( m_location[0], m_location[1],
                           m_scale[0],    m_scale[1],
                           m_point_id );
  }

  // Camera Relation Network
  template <class FeatureT>
  void CameraRelationNetwork<FeatureT>::add_node( cnode const& node ) {
    m_nodes.push_back( node );
  }

  template <class FeatureT>
  void CameraRelationNetwork<FeatureT>::read_controlnetwork( ControlNetwork const& /*cnet*/ ) {
    vw_throw( NoImplErr() << "Not Finished Yet!" );
  }

  template <class FeatureT>
  void CameraRelationNetwork<FeatureT>::write_controlnetwork( ControlNetwork & cnet ) const {

    if ( this->size() == 0 )
      vw_throw( ArgumentErr() << "CameraRelation network is empty." );
    cnet.clear();

    // This process is destructive, making a copy of self.
    CameraRelationNetwork<FeatureT> crn = (*this);

    // On top of building the control network, we're also going filter
    // out 'spiral' type errors. Features that managed to link to the
    // image twice.
    int spiral_error_count = 0;
    {
      TerminalProgressCallback progress("ba","Assembly:  ");
      progress.report_progress(0);
      for ( uint32 i = 0; i < crn.size() - 1; i++ ) {
        progress.report_progress(float(i)/float(crn.size()-1));
        typedef boost::shared_ptr<FeatureT> f_ptr;
        typedef typename std::list<f_ptr>::iterator f_list_iter;

        // Iterating over matched relations inside a camera node and
        // building control points for them.
        for ( f_list_iter iter = crn[i].begin();
              iter != crn[i].end(); iter++ ) {
          // 1.) Building a listing of interest point for a control
          // point
          std::list<f_ptr> interestpts;
          interestpts.push_back(*iter);
          (*iter)->list_connections( interestpts );

          ControlPoint cpoint( ControlPoint::TiePoint );

          // 2.) Adding this location
          cpoint.add_measure( (*iter)->control_measure() );

          // 3.) Adding and Removing measures in all other locations
          f_list_iter measure = interestpts.begin();
          measure++;
          for ( ; measure != interestpts.end(); measure++ ) {
            crn[(*measure)->m_camera_id].relations.remove(*measure);
            cpoint.add_measure( (*measure)->control_measure() );
          }

          // 4.) Removing this location finally
          iter = crn[i].relations.erase(iter);
          iter--;

          // 5.) Checking for spiral error
          {
            std::list<unsigned> previous_camera;
            bool match = false;
            for ( f_list_iter interest = interestpts.begin();
                  interest != interestpts.end(); interest++ ) {
              for ( std::list<unsigned>::iterator previous = previous_camera.begin();
                    previous != previous_camera.end(); previous++ ) {
                if ((*previous) == (*interest)->m_camera_id) {
                  match = true;
                  break;
                }
              }
              previous_camera.push_back( (*interest)->m_camera_id );
              if ( match ) {
                continue;
              }
            }
            if ( match ) {
              spiral_error_count++;
              continue;
            }
          }

          // 6.) Did you pass? Sweet you're in the gang!
          if ( cpoint.size() > 1 )
            cnet.add_control_point( cpoint );
        } // end of iteration over relations
      } // end of iteration over camera nodes
      progress.report_finished();
      if ( spiral_error_count != 0 )
        vw_out(WarningMessage,"ba") << "\t"
                                        << spiral_error_count
                                        << " control points removed due to spiral errors.\n";
    }
    VW_ASSERT( cnet.size() != 0,
               Aborted() << "Failed to load any points, Control Network empty\n" );
  }

  // Explicit template Instantiation
#define VW_INSTANTIATE_CAMERA_RELATION_TYPES(FEATURET) \
  template class CameraRelationNetwork<FEATURET >; \

  VW_INSTANTIATE_CAMERA_RELATION_TYPES(IPFeature)
  VW_INSTANTIATE_CAMERA_RELATION_TYPES(JFeature)

}}
