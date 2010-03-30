// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
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

namespace vw {
namespace ba {

  // Feature Base
  template <class ImplT>
  struct FeatureBase {

    FeatureBase( uint32 id ) : m_camera_id(id) {}

    // Access to derived type
    inline ImplT& impl() { return static_cast<ImplT&>(*this); }
    inline ImplT const& impl() const { return static_cast<ImplT const&>(*this); }

    typedef boost::shared_ptr<ImplT> f_ptr;
    std::list<f_ptr> m_connections;
    uint32 m_camera_id;

    // Returns a string identifying feature type
    std::string type() { return impl().type(); }
    // Connects this feature to another
    void connection( f_ptr& con, bool check=true ) {
      if (check) {
        typedef typename std::list<f_ptr>::iterator list_it;
        for ( list_it iter = m_connections.begin();
              iter != m_connections.end(); iter++ ) {
          if ( (*iter).get() == con.get() ) // Compare pointers only
            return;
        }
      }
      m_connections.push_back( con );
    }
    // List all features connected to this one. This will traverse all
    // connections and append them to the provided list.
    void list_connections( std::list<f_ptr>& listing ) {
      typedef typename std::list<f_ptr>::iterator list_it;
      for ( list_it iter = m_connections.begin();
            iter != m_connections.end(); iter++ ) {
        bool contains = false;
        for ( list_it diter = listing.begin();
              diter != listing.end(); diter++ )
          if ( *diter == *iter ) {
            contains = true;
            break;
          }

        if ( !contains ) {
          listing.push_back(*iter);
          (*iter)->list_connections( listing );
        }
      }
    }
    // Interface to create Control Measure
    ControlMeasure control_measure() const { return impl().control_measure(); }
  };

  // Interest Point Feature
  // - Intended for fast insertion of IP matches
  struct IPFeature : public FeatureBase<IPFeature> {
    typedef boost::shared_ptr<IPFeature> f_ptr;
    ip::InterestPoint m_ip;

    IPFeature ( ip::InterestPoint const& ip,
                uint32 const& id ) :
    FeatureBase<IPFeature>(id), m_ip(ip) {}
    std::string type() { return "IP"; }
    ControlMeasure control_measure() const;
  };

  // Jacobian Feature
  // - Intended for internal use in BA
  struct JFeature : public FeatureBase<JFeature> {
    typedef boost::shared_ptr<JFeature> f_ptr;
    Matrix<double> m_a_jacobian;
    Matrix<double> m_b_jacobian;
    uint32 m_point_id;
    Vector2f m_location;
    Vector2f m_scale;

    JFeature ( uint32 const& i, uint32 const& j ) :
    FeatureBase<JFeature>(j), m_point_id(i) {}
    std::string type() { return "J"; }
    ControlMeasure control_measure() const;
  };

  // Camera Relation Class
  template <class FeatureT>
  struct CameraNode {
    typedef boost::shared_ptr<FeatureT> f_ptr;
    int id;
    std::string name;
    std::list<f_ptr> relations;

    CameraNode ( int tid, std::string tname ) :
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
    uint32 size( void ) const { return m_nodes.size(); }
    cnode& operator[]( int32 const& i ) { return m_nodes[i]; }
    iterator begin() { return m_nodes.begin(); }
    iterator end()   { return m_nodes.end();   }
    const_iterator begin() const { return m_nodes.begin(); }
    const_iterator end()   const { return m_nodes.end();   }

    // Complex functions
    void add_node( cnode const& node );
    void read_controlnetwork( ControlNetwork const& cnet );
    void write_controlnetwork( ControlNetwork & cnet ) const;

  };

}}

#endif//__VW_BUNDLEADJUSTMENT_CAMERA_RELATIONS_H__
