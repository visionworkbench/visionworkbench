// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_PLATE_TREE_H__
#define __VW_PLATE_TREE_H__

#include <vw/Core/Log.h>
#include <vw/Core/Exception.h>
#include <vw/Core/FundamentalTypes.h>

// TODO: Remove this once we fix invalidate()
#include <vw/Plate/Exception.h>
#include <vw/Plate/ProtoBuffers.pb.h>

#include <boost/shared_array.hpp>
#include <map>

namespace vw {
namespace platefile {

  // A subclass of TreeMapFunc can be used to iterate over a tree
  // using the TreeNode::map() function below.
  struct TreeMapFunc {
    virtual ~TreeMapFunc() {}
    virtual void operator()(int32 col, int32 row, int32 depth) = 0;
  };

  template <class ElementT>
  class TreeNode {
    typedef std::map<vw::int32,ElementT,std::greater<vw::int32> > record_type;

    TreeNode<ElementT> *m_parent;
    std::vector<boost::shared_ptr<TreeNode<ElementT> > > m_children;
    record_type m_records;
    int m_max_depth;

    // ------------------------ Private Methods -----------------------------

    // Do some quick arithmetic here to determine which child to follow
    // in the quad tree for a given col, row, level request.
    int compute_child_id(int col, int row, int level, int current_level) {

      int tile_x = col / int(pow(2,level-current_level));
      int tile_y = row / int(pow(2,level-current_level));

      if (tile_x >= pow(2, current_level) ||
          tile_y >= pow(2, current_level))
        vw_throw(TileNotFoundErr() << "TreeNode: invalid index (" 
                 << col << " " << row << " at level " << level << ").");
      tile_x %= 2;
      tile_y %= 2;

      // For debugging
      // std::cout << "Adding at " << tile_x << " " << tile_y << "   " 
      //           << col << " " << row << "  " << level << " " << current_level << "\n";
      if (tile_x == 0 && tile_y == 0) 
        return 0;
      else if (tile_x == 1 && tile_y == 0) 
        return 1;
      else if (tile_x == 0 && tile_y == 1) 
        return 2;
      else
        return 3;
    }

    // Sets the child of this node with the 'id' according to the above
    // index scheme.
    //
    void set_child(int id, boost::shared_ptr<TreeNode> node) {
      m_children[id] = node;
    }

    // Insert the child of this node, but preserves the previous child's descendents.
    //
    void insert_child(int id, boost::shared_ptr<TreeNode> node) {
      
      // First we save the old child and replace it with the new one.
      boost::shared_ptr<TreeNode> old_child = m_children[id];
      m_children[id] = node;

      // Then, if the old child existed, we transfer the old child's
      // children (our grandchildren) to the new child node.
      if (old_child)
        for (int i = 0; i < 4 ; ++i) 
          m_children[id]->set_child(i, old_child->child(i));
    }

    // Invalidate all of the nodes leading up to the specified node so
    // that they can be regenerated on the next mipmapping pass.
    // 
    // TODO: I'm violating the abstraction barrier here for the sake
    // of expediency.  However, the tree structure really ought to be
    // kept seperate from the set_status() method of the record.
    // Maybe we need to feed in some sort of invalidation functor to
    // the tree structure to keep the barrier intact?
    bool invalidate_records_helper(int col, int row, int level, int current_level) {
      // If we have reached the requested depth, then we must be at
      // the node we want!  Return the ElementT!
      if (current_level == level) {

        // Handle the edge case where the user has requested a tile
        // outside of the 1x1 bounds of the root level.
        if ( level == 0 && (col !=0 || row != 0) )
          vw_throw(TileNotFoundErr() << "TreeNode: invalid index (" << col << " " << row << ").");
        return true;
        
      // Otherwise, we go recurse deeper into the tree....
      } else {        
        int child_id = compute_child_id(col, row, level, current_level + 1);
        if (m_children[child_id]) {
          
          // If a branch of the tree is found, we dive deeper.
          bool success = m_children[child_id]->invalidate_records_helper(col, row, level, 
                                                                         current_level + 1);
          if (success && m_records.size() > 0 &&
              (*(m_records.begin())).second.status() == INDEX_RECORD_VALID || 
              (*(m_records.begin())).second.status() == INDEX_RECORD_PRIMED )
            (*(m_records.begin())).second.set_status(INDEX_RECORD_STALE);
          
          return success;
          
        } else {
          return false;
        }
      }
    }

    // Search for a node at a given col, row, and level.
    ElementT search_helper(int col, int row, int level, int transaction_id, int current_level) {
      // For debugging:
      // std::cout << "Call to search_helper(" << col << " " << row << " " << level << " " 
      //           << transaction_id << " " << current_level << ")\n";

      // If we have reached the requested depth, then we must be at
      // the node we want!  Return the ElementT!
      if (current_level == level) {

        // Handle the edge case where the user has requested a tile
        // outside of the 1x1 bounds of the root level.
        if ( level == 0 && (col !=0 || row != 0) )
          vw_throw(TileNotFoundErr() << "TreeNode: invalid index (" << col << " " << row << ").");

        return this->value(transaction_id);

      // Otherwise, we go recurse deeper into the tree....
      } else {
        
        int child_id = compute_child_id(col, row, level, current_level + 1);
        
        if (m_children[child_id]) {

          // If a branch of the tree is found, we dive deeper.
          return m_children[child_id]->search_helper(col, row, level, 
                                                     transaction_id, current_level + 1);

        } else {
          // If not, we throw an exception.
          vw_throw(TileNotFoundErr() << "Tile search [" << col << " " << row << " " 
                   << current_level << "] failed at depth " << current_level << "\n");
        }
      }
      // If not, we throw an exception.
      vw_throw(TileNotFoundErr() << "Tile search [" << col << " " << row << " " 
               << current_level << "] failed at depth " << current_level << "\n");
    }

    // Recursively call a function with valid [col, row, level] entries.
    //
    //    |---|---|
    //    | 0 | 1 |
    //    |---+---|
    //    | 2 | 3 |
    //    |---|---|
    //
    void map_helper(boost::shared_ptr<TreeMapFunc> func, int col, int row, int level) {

      // Call the function for the current level.
      (*func)(col, row, level);

      // Call the function for future levels.
      if ( this->child(0) ) 
        this->child(0)->map_helper(func, col*2, row*2, level + 1);

      if ( this->child(1) ) 
        this->child(1)->map_helper(func, col*2+1, row*2, level + 1);

      if ( this->child(2) ) 
        this->child(2)->map_helper(func, col*2, row*2+1, level + 1);

      if ( this->child(3) ) 
        this->child(3)->map_helper(func, col*2+1, row*2+1, level + 1);
    }

    void insert_helper(ElementT const& record, 
                       int col, int row, int level, 
                       int transaction_id, int current_level,
                       bool insert_at_all_levels) {

      // If we have reached the requested depth, then we must be at
      // the node we want!  Return the ElementT!
      if (current_level == level) {

        // Handle the edge case where the user has requested a tile
        // outside of the 1x1 bounds of the root level.
        if ( level == 0 && (col !=0 || row != 0) )
          vw_throw(TileNotFoundErr() << "TreeNode: invalid index (" << col << " " << row << ").");
        // Add the record for the given transaction ID.
        m_records[transaction_id] = record;

      // Otherwise, we need to recurse further into the tree....
      } else {

        int child_id = this->compute_child_id(col, row, level, current_level+1);
        
        // If the child we need is not yet created, we create it, add
        // it as our child, and recurse down it.
        if (!m_children[child_id]) {
          boost::shared_ptr<TreeNode> node( new TreeNode(this, ElementT(), transaction_id ) );
          this->set_child(child_id, node);
        }
 
        // Insert the record at this level if requested.  Right now this is mostly a
        // special feature used by the transaction_request() code in
        // LocalIndex to "prime" the tree with empty index entries.
        if (insert_at_all_levels)
          m_records[transaction_id] = record;
        
        m_children[child_id]->insert_helper(record, col, row, level, 
                                            transaction_id, current_level+1,
                                            insert_at_all_levels);
      }

    }


    /// XXX : This is an abstraction violation!!
    void clean_branch_helper(int col, int row, int level, int transaction_id, int current_level) {

      // If we have reached the requested depth, then we must be at
      // the node we want!  Return the ElementT!
      if (current_level == level) {

        // Handle the edge case where the user has requested a tile
        // outside of the 1x1 bounds of the root level.
        if ( level == 0 && (col !=0 || row != 0) )
          vw_throw(TileNotFoundErr() << "TreeNode: invalid index (" << col << " " << row << ").");

        try {
          ElementT val = this->value(transaction_id);
          if (val.status() == INDEX_RECORD_PRIMED) 
            m_records.erase(transaction_id);
        } catch (TileNotFoundErr &e) {
          // do nothing if the value is not found... there is no need to erase it.
        }

      // Otherwise, we go recurse deeper into the tree....
      } else {
        
        int child_id = compute_child_id(col, row, level, current_level + 1);
        
        if (m_children[child_id]) {

          // If a branch of the tree is found, we dive deeper.
          m_children[child_id]->clean_branch_helper(col, row, level, 
                                                    transaction_id, current_level + 1);

          try {
            ElementT val = this->value(transaction_id);
            if (val.status() == INDEX_RECORD_PRIMED) 
              m_records.erase(transaction_id);
          } catch (TileNotFoundErr &e) {
            // do nothing if the value is not found
          }

        } 
      }
    }

    void print_helper(int current_level) const {
      vw_out(0) << (*(m_records.begin())).second.status() << "\n";
      for (int i = 0; i < 4; ++i) 
        if ( this->child(i) ) {
          for (int l = 0; l < current_level+1; ++l) 
            vw_out(0) << "  ";
          vw_out(0) << "[ " << (current_level+1) 
                     << " -- Child " << i << " ] - " ;
          this->child(i)->print_helper(current_level + 1);
        }
    }

  public: 

    // ------------------------ Public Methods -----------------------------

    /// Use this constructor for the root of the tree....
    TreeNode() : m_parent(0), m_max_depth(0) {
      m_records[0] = ElementT();
      m_children.resize(4);
    }

    /// ... or use this contructor for the root of the tree.
    TreeNode(ElementT const& record, int transaction_id) : m_parent(0), m_max_depth(0) {
      m_records[transaction_id] = record;
      m_children.resize(4);
    }

    /// Use this contstructor to add record data immediately.
    TreeNode(TreeNode *parent, ElementT const& record, int transaction_id) :
      m_parent(parent), m_max_depth(0) {
      m_records[transaction_id] = record;
      m_children.resize(4);
    }
    
    // Return the child of this node with the 'id' according to the
    // following index scheme:
    //
    //    |---|---|
    //    | 0 | 1 |
    //    |---+---|
    //    | 2 | 3 |
    //    |---|---|
    //
    boost::shared_ptr<TreeNode> child(int id) const { return m_children[id]; }

    int num_children() const {
      int n = 0;
      for (int i = 0; i < 4; ++ i) 
        if (m_children[i])
          ++n;
      return n;
    }

    ElementT value_helper(int transaction_id) { 
      typename record_type::iterator it = m_records.begin();

      // A transaction ID of -1 indicates that we should return the
      // most recent tile, regardless of its transaction id.
      if (transaction_id == -1 && it != m_records.end())
        return (*it).second;

      // Otherwise, we search for the most recent record that happened
      // before on on the queried ephoch.
      while (it != m_records.end()) {
        if ((*it).first <= transaction_id)
          return (*it).second;
        ++it;
      }

      // If we reach this point, then there are no entries before
      // the given transaction_id, so we return an empty (and invalid) record.
      vw_throw(TileNotFoundErr() << "Tiles exist at this location, but none before transaction_id=" 
               << transaction_id << "\n");
      return ElementT(); // never reached
    }

    /// Access the data member of this node.
    ElementT value(int transaction_id) { 
      return value_helper(transaction_id);
    }

    const ElementT value(int transaction_id) const { 
      return value_helper(transaction_id);
    }

    /// Query max depth of tree.
    int max_depth() const { return m_max_depth; }


    // Search for a node at a given col, row, and level.  Note: a
    // transaction ID of -1 indicates that we should return the
    // most recent tile, regardless of its transaction id.
    ElementT search(int col, int row, int level, int transaction_id) {
      return search_helper(col, row, level, transaction_id, 0);
    }

    // Search for a node at a given col, row, and level.  Note: a
    // transaction ID of -1 indicates that we should return the
    // most recent tile, regardless of its transaction id.
    void clean_branch(int col, int row, int level, int transaction_id) {
      clean_branch_helper(col, row, level, transaction_id, 0);
    }

    // Search for a node at a given col, row, and level.
    bool invalidate_records(int col, int row, int level) {
      return invalidate_records_helper(col, row, level, 0);
    }

    // Insert an ElementT at a given position.  Intermediate nodes
    // in the tree are created (with empty ElementTs) in the tree
    // along the way, as needed.
    void insert(ElementT const& record, int col, int row, 
                int level, int transaction_id,
                bool insert_at_all_levels = false) {
      this->insert_helper(record, col, row, level, transaction_id, 0, insert_at_all_levels);
      if (level > m_max_depth)
        m_max_depth = level;
    }

    /// Print the tree.  (Use only for debugging small trees....)
    void print() const { 
      vw_out(0) << "[ 0 -- Child 0 ] - ";
      this->print_helper(0); 
    }

    void map(boost::shared_ptr<TreeMapFunc> func) {
      this->map_helper(func, 0, 0, 0);
    }

  };



}} // namespace vw

#endif // __VW_PLATE_TREE_H__
