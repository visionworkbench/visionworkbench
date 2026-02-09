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


/// \file CameraRelation.h Camera Relations Network
///
/// This is not a replacement for the Control Network, but instead a
/// different approach of viewing the same data. It's optimized in a way
/// to allow for fast insertions of match data, it's great for organizing
/// jacobians right before insertion into hessians, and it also provides a
/// constant structure for handling internal variables inside of bundle
/// adjustment.
///
/// Control Networks are still useful for providing an organized
/// manner for saving data and for interfacing with ISIS. Control Networks
/// are also easier to understand at first glance. A control network can
/// be created from a camera relations network and the opposite is also
/// true.

#ifndef __VW_BUNDLEADJUSTMENT_CAMERA_RELATIONS_H__
#define __VW_BUNDLEADJUSTMENT_CAMERA_RELATIONS_H__

#include <vw/Math/Matrix.h>
#include <vw/BundleAdjustment/ControlNetwork.h>
#include <vw/InterestPoint/InterestPoint.h>
#include <boost/foreach.hpp>
#include <map>

namespace vw {
namespace ba {

/** \addtogroup BundleAdjustment
 *  @{
 */

  // Feature Base
  template <class ImplT>
  struct FeatureBase {

    FeatureBase( size_t id ) : m_camera_id(id) {}

    // Access to derived type
    inline ImplT      & impl()       { return static_cast<ImplT      &>(*this); }
    inline ImplT const& impl() const { return static_cast<ImplT const&>(*this); }

    // Features use weak points as they can't be allowed to have
    // ownership. The reason is that 2 features can keep themselves from
    // being deleted if they have shared_ptrs of each
    // other. Shared_ptrs and the control of deletion is handled by Camera
    // Nodes below.
    typedef boost::weak_ptr<ImplT> w_ptr;
    std::list<w_ptr> m_connections;
    std::map<size_t,w_ptr> m_map; // Alternative access compared to
    // connections. This should only be used after spiral error detection has
    // been performed. Map can used to find quickly the feature that belongs
    // to a camera specific camera.

    size_t m_camera_id;

    /// Returns a string identifying feature type
    std::string type() { return impl().type(); }

    /// Connects this feature to another
    void connection( w_ptr con, bool check=true ) {
      if (check) {
        BOOST_FOREACH( w_ptr connection, m_connections )
          if ( connection.lock() == con.lock() )
            return;
      }
      m_connections.push_back( con );
    }
    
    // List all features connected to this one. This will traverse all
    // connections and append them to the provided list.
    void list_connections( std::list<w_ptr>& listing ) {
      BOOST_FOREACH( w_ptr connection, m_connections ) {
        bool contains = false;
        BOOST_FOREACH( w_ptr pconnection, listing )
          if ( pconnection.lock() == connection.lock() ) {
            contains = true;
            break;
          }

        if ( !contains ) {
          listing.push_back(connection);
          connection.lock()->list_connections( listing );
        }
      }
    }

    void build_map() {
      m_map.clear();
      BOOST_FOREACH( w_ptr connection, m_connections )
        m_map[ connection.lock()->m_camera_id ] = connection;
    }

    // Interface to aid conversion with Control Network
    ControlMeasure control_measure() const { return impl().control_measure(); }
  }; // End class FeatureBase

  // Interest Point Feature
  // - Intended for fast insertion of IP matches
  struct IPFeature : public FeatureBase<IPFeature> {
    ip::InterestPoint m_ip;

    // Standard Constructor
    IPFeature ( ip::InterestPoint const& ip, size_t const& id ) :
          FeatureBase<IPFeature>(id), m_ip(ip) {}

    // For building from control networks
    IPFeature( ControlMeasure const& cmeas, size_t const& /*point_id*/ ) :
          FeatureBase<IPFeature>( cmeas.image_id() ) {
      m_ip = ip::InterestPoint( cmeas.position()[0], cmeas.position()[1],
                                cmeas.sigma()[0] );
    }

    IPFeature( ControlMeasure const& cmeas,
               size_t const& /*point_id*/,
               size_t image_id ) :
          FeatureBase<IPFeature>( image_id ) {
      m_ip = ip::InterestPoint( cmeas.position()[0], cmeas.position()[1],
                                cmeas.sigma()[0] );
    }

    std::string type() { return "IP"; }
    ControlMeasure control_measure() const;
  };

  std::ostream& operator<<( std::ostream& os, IPFeature const& feat );

  // Jacobian Feature
  // - Intended for internal use in BA
  struct JFeature : public FeatureBase<JFeature> {
    Matrix<double> m_w, m_y; // W = product of Jacobians, Y = product of W

    size_t   m_point_id;
    Vector2f m_location;
    Vector2f m_scale;

    // Standard Constructor
    JFeature ( size_t const& point_id, size_t const& camera_id ) :
        FeatureBase<JFeature>(camera_id), m_point_id(point_id) {}

    // For building from control networks
    JFeature ( ControlMeasure const& cmeas, size_t const& point_id ) :
        FeatureBase<JFeature>( cmeas.image_id() ), m_point_id(point_id) {
      m_location = cmeas.position();
      m_scale    = cmeas.sigma();
    }

    std::string type() { return "J"; }
    ControlMeasure control_measure() const;
  };

  std::ostream& operator<<( std::ostream& os, JFeature const& feat );

  // Camera Relation Class
  template <class FeatureT>
  struct CameraNode {
    typedef boost::shared_ptr<FeatureT> f_ptr;
    size_t      id;
    std::string name;
    std::list<f_ptr> relations;
    std::multimap< size_t, f_ptr> map; // Provides alternative access
    // to links inside the elements in relations. It's keyed by all
    // cameras. Example for using this map is for looking via equal_range(x)
    // to find all features in this camera that connect to camera x.

    CameraNode ( size_t tid, std::string tname ) :
      id(tid), name(tname) {}

    // Iterator access (saves a section for me the monkey)
    typedef typename std::list<f_ptr>::iterator iterator;
    typedef typename std::list<f_ptr>::const_iterator const_iterator;
    iterator       begin()       { return relations.begin(); }
    iterator       end  ()       { return relations.end();   }
    const_iterator begin() const { return relations.begin(); }
    const_iterator end  () const { return relations.end();   }
  };

  /// Alternate storage method of the information in a ControlNetwork object.
  /// - See the description at the top of the file.
  template <class FeatureT>
  class CameraRelationNetwork {
    typedef CameraNode<FeatureT> cnode;
    std::vector<cnode> m_nodes; // One for each image/camera instance.

  public:

    // Iterator access
    typedef typename std::vector<cnode>::iterator       iterator;
    typedef typename std::vector<cnode>::const_iterator const_iterator;

    size_t size() const { return m_nodes.size(); }

    cnode& operator[]( int32 const& i ) { return m_nodes[i]; }
    cnode const& operator[]( int32 const& i ) const { return m_nodes[i]; }

    iterator       begin()       { return m_nodes.begin(); }
    iterator       end  ()       { return m_nodes.end();   }
    const_iterator begin() const { return m_nodes.begin(); }
    const_iterator end  () const { return m_nodes.end();   }

    // Complex functions
    void add_node( cnode const& node );
    void build_map();
    void from_cnet ( ControlNetwork const& cnet );
  };


/** @}  End group BundleAdjustment*/

}}

#endif//__VW_BUNDLEADJUSTMENT_CAMERA_RELATIONS_H__
