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
  void CameraRelationNetwork<FeatureT>::from_cnet( ControlNetwork const& cnet ) {
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

  // Explicit template instantiation.
  // TODO(oalexan1): Use this kind of logic more often to avoid long
  // compile times, if it is known which classes will be used as
  // template arguments.
  template class CameraRelationNetwork<IPFeature>;
  template class CameraRelationNetwork<JFeature>;

}}
