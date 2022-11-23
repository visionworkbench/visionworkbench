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


/// \file CameraRelation.cc
///

#include <vw/BundleAdjustment/CameraRelation.h>
#include <boost/foreach.hpp>

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
                           m_camera_id );
  }

  std::ostream& operator<<( std::ostream& os, IPFeature const& feat ) {
    os << "IPFeature( (" << feat.m_ip.x << "," << feat.m_ip.y << ")@"
       << feat.m_camera_id << " links "
       << feat.m_connections.size() << " )";
    return os;
  }

  std::ostream& operator<<( std::ostream& os, JFeature const& feat ) {
    os << "JFeature( " << feat.m_point_id << " " << feat.m_location
       << "@" << feat.m_camera_id << " links "
       << feat.m_connections.size() << " )";
    return os;
  }

  // Camera Relation Network
  template <class FeatureT>
  void CameraRelationNetwork<FeatureT>::add_node( cnode const& node ) {
    m_nodes.push_back( node );
  }

  template <class FeatureT>
  void CameraRelationNetwork<FeatureT>::build_map() {
    typedef          boost::shared_ptr<FeatureT>            f_ptr;
    typedef          boost::weak_ptr<FeatureT>              w_ptr;
    typedef typename std::list<f_ptr>::iterator             list_it;
    typedef typename std::list<f_ptr>::const_iterator       list_cit;
    typedef typename std::map<size_t,w_ptr>::const_iterator map_cit;

    for ( size_t j = 0; j < m_nodes.size(); j++ ) {
      // Building maps for all features
      for ( list_it fiter = m_nodes[j].relations.begin();
            fiter != m_nodes[j].relations.end(); fiter++ ) {
        (**fiter).build_map();
      }

      // Finally building multimap
      m_nodes[j].map.clear();
      // Iterating through all features within this camera
      for ( list_cit fiter = m_nodes[j].relations.begin();
            fiter != m_nodes[j].relations.end(); fiter++ ) {
        // Iterating through all features that our feature connects to
        for ( map_cit miter = (**fiter).m_map.begin();
              miter != (**fiter).m_map.end(); miter++ ) {
          m_nodes[j].map.insert( std::make_pair( (*miter).first, *fiter  ) );
        }
      } // end iterating through this camera's features
    }   // end iterating through cameras
  }

  template <class FeatureT>
  void CameraRelationNetwork<FeatureT>::read_controlnetwork( ControlNetwork const& cnet ) {
    typedef boost::shared_ptr<FeatureT> f_ptr;
    typedef boost::weak_ptr<FeatureT>   w_ptr;
    m_nodes.clear();

    // Add a node for each image in the list
    for ( size_t i = 0; i < cnet.get_image_list().size(); i++ )
      this->add_node( CameraNode<FeatureT>( i, "" ));
    
    for ( size_t point_id = 0; point_id < cnet.size(); point_id++ ) {
      std::vector<f_ptr> features_added;
      // Building up features to be added and linking to camera nodes
      BOOST_FOREACH( ControlMeasure const& cm, cnet[point_id] ) {
        
        // Extra robustness, ensure nodes exist for all encountered
        // indices. This is necessary since the cnet is sometimes
        // created without populating the image list.
        if ( cm.image_id() >= this->size() ) {
          for ( size_t i = this->size(); i <= cm.image_id(); i++ ) {
            this->add_node( CameraNode<FeatureT>( i, "" ) );
          }
        }
        
        // Appending to list
        features_added.push_back( f_ptr( new FeatureT(cm, point_id) ) );
        
        // Attaching to camera node
        m_nodes[cm.image_id()].relations.push_back( features_added.back() );
      } // End loop through ControlMeasures
      
      // Doubly-linking features together
      typedef typename std::vector<f_ptr>::iterator fvi_ptr;
      for ( fvi_ptr first = features_added.begin();
            first < features_added.end() - 1; first++ ) {
        for ( fvi_ptr second = first + 1;
              second < features_added.end(); second++ ) {
          (*first)->connection( w_ptr( *second ), false );
          (*second)->connection( w_ptr( *first ), false );
        }
      }

    } // end for through control points

    // setting up maps
    this->build_map();
  }

  template <class FeatureT>
  bool CameraRelationNetwork<FeatureT>::assemble_cnet(ControlNetwork & cnet) const {

    if ( this->size() == 0 )
      vw_throw( ArgumentErr() << "CameraRelation network is empty." );
    cnet.clear();

    // This process is destructive, so start by making a copy of self
    // that we will modify.
    CameraRelationNetwork<FeatureT> crn = (*this);

    // On top of building the control network, we're also going filter
    // out 'spiral' type errors. Features that managed to link to the
    // image twice.
    int spiral_error_count = 0;
    {
      TerminalProgressCallback progress("ba","Assembly:  ");
      progress.report_progress(0);
      for ( size_t i = 0; i < crn.size() - 1; i++ ) {
        progress.report_progress(float(i)/float(crn.size()-1));
        typedef          boost::weak_ptr<FeatureT>   w_ptr;
        typedef          boost::shared_ptr<FeatureT> f_ptr;
        typedef typename std::list<f_ptr>::iterator  f_list_iter;
        typedef typename std::list<w_ptr>::iterator  w_list_iter;

        // Iterating over matched relations inside a camera node and
        // building control points for them.
        for ( f_list_iter iter = crn[i].begin();
              iter != crn[i].end(); iter++ ) {
          // 1.) Building a listing of interest point for a control
          // point
          std::list<w_ptr> interestpts;
          interestpts.push_back(*iter);
          (*iter)->list_connections( interestpts );

          ControlPoint cpoint( ControlPoint::TiePoint );

          // 2.) Adding this location
          cpoint.add_measure( (*iter)->control_measure() );

          // 3.) Adding and Removing measures in all other locations
          w_list_iter measure = interestpts.begin();
          measure++;
          for ( ; measure != interestpts.end(); measure++ ) {
            crn[(*measure).lock()->m_camera_id].relations.remove((*measure).lock());
            cpoint.add_measure( (*measure).lock()->control_measure() );
          }

          // 4.) Removing this location finally
          iter = crn[i].relations.erase(iter);
          iter--;

          // 5.) Checking for spiral error
          {
            std::list<size_t> previous_camera;
            bool match = false;
            for ( w_list_iter interest = interestpts.begin();
                  interest != interestpts.end(); interest++ ) {
              BOOST_FOREACH( size_t previous, previous_camera ) {
                if ( previous == interest->lock()->m_camera_id ) {
                  match = true;
                  break;
                }
              }
              previous_camera.push_back( (*interest).lock()->m_camera_id );
              if ( match )
                continue;
            }
            if ( match ) {
              spiral_error_count++;
              continue;
            }
          }

          // 6.) Did you pass? 
          //     - CPoints with only a single measure are GCPs
          if ( cpoint.size() > 0 )
            cnet.add_control_point( cpoint );
        } // end of iteration over relations
      } // end of iteration over camera nodes
      progress.report_finished();
      if ( spiral_error_count != 0 )
        vw_out(WarningMessage,"ba") << "\t" << spiral_error_count
                                    << " control points removed due to spiral errors.\n";
    }
    if ( cnet.size() == 0) {
      vw_out(WarningMessage,"ba")
        << "Failed to load any points, control network is empty.";
      return false;
    }

    return true;
  }

  // Explicit template instantiation.
  // TODO(oalexan1): Use this kind of logic more often to avoid long
  // compile times, if it is known which classes will be used as
  // template arguments.
  template class CameraRelationNetwork<IPFeature>;
  template class CameraRelationNetwork<JFeature>;

}}
