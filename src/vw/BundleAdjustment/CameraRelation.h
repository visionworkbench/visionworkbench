// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// Camera Relations Network
//
// This is not a replacement for the Control Network, but instead a
// different approach of viewing the same data. It's optimized in a way
// to allow for fast insertions of match data, it's great for organizing
// jacobians right before insertion into hessians, and it also provides a
// constant structure for handling internal variables inside of bundle
// adjustment.
//
// Control Networks are still useful for providing an organized
// manner for saving data and for interfacing with ISIS. Control Networks
// are also easier to understand at first glance. A control network can
// be created from a camera relations network and the opposite is also
// true.

#ifndef __VW_BUNDLEADJUSTMENT_CAMERA_RELATIONS_H__
#define __VW_BUNDLEADJUSTMENT_CAMERA_RELATIONS_H__

#include <vw/Math/Matrix.h>
#include <vw/BundleAdjustment/ControlNetwork.h>
#include <vw/InterestPoint/InterestData.h>
#include <map>

namespace vw {
namespace ba {

  // Feature Base
  template <class ImplT>
  struct FeatureBase {

    FeatureBase( uint32 id ) : m_camera_id(id) {}

    // Access to derived type
    inline ImplT& impl() { return static_cast<ImplT&>(*this); }
    inline ImplT const& impl() const { return static_cast<ImplT const&>(*this); }

    // Features use weak points as they can't be allowed to have
    // ownership. The reason is that 2 features can keep themselves from
    // being deleted if they have shared_ptrs of each
    // other. Shared_ptrs and the control of deletion is handled by Camera
    // Nodes below.
    typedef boost::weak_ptr<ImplT> w_ptr;
    std::list<w_ptr> m_connections;
    std::map<uint32,w_ptr> m_map; // Alternative access compared to
    // connections. This should only be used after spiral error detection has
    // been performed. Map can used to find quickly the feature that belongs
    // to a camera specific camera.

    uint32 m_camera_id;

    // Returns a string identifying feature type
    std::string type() { return impl().type(); }
    // Connects this feature to another
    void connection( w_ptr con, bool check=true ) {
      if (check) {
        typedef typename std::list<w_ptr>::iterator list_it;
        for ( list_it iter = m_connections.begin();
              iter != m_connections.end(); iter++ ) {
          if ( (*iter).lock() == con.lock() ) // Compare pointers only
            return;
        }
      }
      m_connections.push_back( con );
    }
    // List all features connected to this one. This will traverse all
    // connections and append them to the provided list.
    void list_connections( std::list<w_ptr>& listing ) {
      typedef typename std::list<w_ptr>::iterator list_it;
      for ( list_it iter = m_connections.begin();
            iter != m_connections.end(); iter++ ) {
        bool contains = false;
        for ( list_it diter = listing.begin();
              diter != listing.end(); diter++ )
          if ( (*diter).lock() == (*iter).lock() ) {
            contains = true;
            break;
          }

        if ( !contains ) {
          listing.push_back(*iter);
          (*iter).lock()->list_connections( listing );
        }
      }
    }
    void build_map() {
      m_map.clear();
      typedef typename std::list<w_ptr>::const_iterator list_it;
      for ( list_it fiter = m_connections.begin();
            fiter != m_connections.end(); fiter++ )
        m_map[ (*fiter).lock()->m_camera_id ] = *fiter;
    }

    // Interface to aid conversion with Control Network
    ControlMeasure control_measure() const { return impl().control_measure(); }
  };

  // Interest Point Feature
  // - Intended for fast insertion of IP matches
  struct IPFeature : public FeatureBase<IPFeature> {
    ip::InterestPoint m_ip;

    // Standard Constructor
    IPFeature ( ip::InterestPoint const& ip,
                uint32 const& id ) :
    FeatureBase<IPFeature>(id), m_ip(ip) {}
    // For building from control networks
    IPFeature ( ControlMeasure const& cmeas,
                uint32 const& /*point_id*/,
                int32 image_id = -1 ) :
    FeatureBase<IPFeature>( image_id != -1 ? image_id : cmeas.image_id() ) {
      m_ip = ip::InterestPoint( cmeas.position()[0],
                                cmeas.position()[1],
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

    uint32 m_point_id;
    Vector2f m_location;
    Vector2f m_scale;

    // Standard Constructor
    JFeature ( uint32 const& point_id, uint32 const& camera_id ) :
    FeatureBase<JFeature>(camera_id), m_point_id(point_id) {}
    // For building from control networks
    JFeature ( ControlMeasure const& cmeas,
               uint32 const& point_id ) :
    FeatureBase<JFeature>( cmeas.image_id() ), m_point_id(point_id) {
      m_location = cmeas.position();
      m_scale = cmeas.sigma();
    }

    std::string type() { return "J"; }
    ControlMeasure control_measure() const;
  };

  std::ostream& operator<<( std::ostream& os, JFeature const& feat );

  // Camera Relation Class
  template <class FeatureT>
  struct CameraNode {
    typedef boost::shared_ptr<FeatureT> f_ptr;
    uint32 id;
    std::string name;
    std::list<f_ptr> relations;
    std::multimap< uint32, f_ptr> map; // Provides alternative access
    // to links inside the elements in relations. It's keyed by all
    // cameras. Example for using this map is for looking via equal_range(x)
    // to find all features in this camera that connect to camera x.

    CameraNode ( uint32 tid, std::string tname ) :
    id(tid), name(tname) {}

    // Iterator access (saves a section for me the monkey)
    typedef typename std::list<f_ptr>::iterator iterator;
    typedef typename std::list<f_ptr>::const_iterator const_iterator;
    iterator begin() { return relations.begin(); }
    iterator end()   { return relations.end();   }
    const_iterator begin() const { return relations.begin(); }
    const_iterator end()   const { return relations.end();   }
  };

  // Base of the Camera Relation Network
  template <class FeatureT>
  class CameraRelationNetwork {
    typedef CameraNode<FeatureT> cnode;
    std::vector<cnode> m_nodes;

  public:

    // Iterator access
    typedef typename std::vector<cnode>::iterator iterator;
    typedef typename std::vector<cnode>::const_iterator const_iterator;
    uint32 size() const { return m_nodes.size(); }
    cnode& operator[]( int32 const& i ) { return m_nodes[i]; }
    iterator begin() { return m_nodes.begin(); }
    iterator end()   { return m_nodes.end();   }
    const_iterator begin() const { return m_nodes.begin(); }
    const_iterator end()   const { return m_nodes.end();   }

    // Complex functions
    void add_node( cnode const& node );
    void build_map();
    void read_controlnetwork( ControlNetwork const& cnet );
    void write_controlnetwork( ControlNetwork & cnet ) const;

  };

}}

#endif//__VW_BUNDLEADJUSTMENT_CAMERA_RELATIONS_H__
