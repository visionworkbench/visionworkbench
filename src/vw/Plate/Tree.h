// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#ifndef __VW_PLATE_TREE_H__
#define __VW_PLATE_TREE_H__

#include <vw/Core/Log.h>
#include <vw/Core/Exception.h>
#include <vw/Core/FundamentalTypes.h>

#include <boost/shared_array.hpp>

namespace vw {
namespace platefile {

  /// TileNotFound exception
  VW_DEFINE_EXCEPTION(TileNotFoundErr, Exception);

  /// IndexError exception
  VW_DEFINE_EXCEPTION(IndexErr, Exception);


  template <class ElementT>
  class TreeNode {

    TreeNode<ElementT> *m_parent;
    std::vector<boost::shared_ptr<TreeNode<ElementT> > > m_children;
    ElementT m_record;

    // ------------------------ Private Methods -----------------------------

    // Do some quick arithmetic here to determine which child to follow
    // in the quad tree for a given col, row, level request.
    int compute_child_id(int col, int row, int level, int current_level) {

      int tile_x = col / int(pow(2,level-current_level));
      int tile_y = row / int(pow(2,level-current_level));

      if (tile_x >= pow(2, current_level) ||
          tile_y >= pow(2, current_level))
        vw_throw(IndexErr() << "TreeNode: invalid index (" << col << " " << row << " at level " << level << ").");
      tile_x %= 2;
      tile_y %= 2;

      int child_id;

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

    // Search for a node at a given col, row, and level.
    ElementT search_helper(int col, int row, int level, int current_level) {

      // If we have reached the requested depth, then we must be at
      // the node we want!  Return the ElementT!
      if (current_level == level) {

        // Handle the edge case where the user has requested a tile
        // outside of the 1x1 bounds of the root level.
        if ( level == 0 && (col !=0 || row != 0) )
          vw_throw(IndexErr() << "TreeNode: invalid index (" << col << " " << row << ").");

        return m_record;
        
      // Otherwise, we go recurse deeper into the tree....
      } else {
        
        int child_id = compute_child_id(col, row, level, current_level + 1);
        
        if (m_children[child_id]) {

          // If a branch of the tree is found, we dive deeper.
          return m_children[child_id]->search_helper(col, row, level, current_level + 1);

        } else {
          // If not, we throw an exception.
          vw_throw(TileNotFoundErr() << "Tile search [" << col << " " << row << " " 
                   << current_level << "] failed at depth " << current_level << "\n");

        }
      }
    }

    void insert_helper(ElementT const& record, int col, int row, int level, int current_level) {

      // If we have reached the requested depth, then we must be at
      // the node we want!  Return the ElementT!
      if (current_level == level) {

        // Handle the edge case where the user has requested a tile
        // outside of the 1x1 bounds of the root level.
        if ( level == 0 && (col !=0 || row != 0) )
          vw_throw(IndexErr() << "TreeNode: invalid index (" << col << " " << row << ").");

        // TODO: This is where we could keep the history of ElementTs
        // that were used for this node.  This would allow us to go
        // "back in time" in the journal of the blob.  Let's add that
        // feature someday!!
        m_record = record;

      // Otherwise, we need to recurse further into the tree....
      } else {

        int child_id = this->compute_child_id(col, row, level, current_level+1);
        
        // If the child we need is not yet created, we create it, add
        // it as our child, and recurse down it.
        if (!m_children[child_id]) {
          boost::shared_ptr<TreeNode> node( new TreeNode(this, ElementT() ) );
          this->set_child(child_id, node);
        }
         
        m_children[child_id]->insert_helper(record, col, row, level, current_level+1);
      }

    }

    void print_helper(int current_level) const {
      vw_out(0) << m_record.valid() << "\n";
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

    /// Use this contructor for the root of the tree....
    TreeNode() {
      m_children.resize(4);
    }

    /// ... or use this contructor for the root of the tree.
    TreeNode(ElementT const& record) : m_record(record) {
      m_children.resize(4);
    }

    /// Use this contstructor to add record data immediately.
    TreeNode(TreeNode *parent, ElementT const& record) :
      m_parent(parent), m_record(record) {
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

    /// Access the data member of this node.
    ElementT& value() { return m_record; }
    const ElementT& value() const { return m_record; }



    // Search for a node at a given col, row, and level.
    ElementT search(int col, int row, int level) {
      return search_helper(col, row, level, 0);
    }

    // Insert an ElementT at a given position.  Intermediate nodes
    // in the tree are created (with empty ElementTs) in the tree
    // along the way, as needed.
    void insert(ElementT const& record, int col, int row, int level) {
      this->insert_helper(record, col, row, level, 0);
    }

    /// Print the tree.  (Use only for debugging small trees....)
    void print() const { 
      vw_out(0) << "[ 0 -- Child 0 ] - ";
      this->print_helper(0); 
    }

    
    /// Serialize a tree as binary data for storage on disk.
    void serialize(std::ostream &ostr) {

      // First we serialize our record data.
      int ob_size = m_record.ByteSize();
      ostr.write((char*)(&ob_size), sizeof(ob_size));      
      m_record.SerializeToOstream(&ostr);

      std::string debug;
      m_record.SerializeToString(&debug);

      // Then we serialize our children.  First we write the number of
      // children at this node.  If a child exists, we then write
      // its ID, then we serialize it.  This will cause the entire
      // child branch of the tree to be serialized in a depth-first
      // manner.
      uint8 num_children = uint8(this->num_children());
      ostr.write( (char*)&num_children, sizeof(num_children) );
          
      for (uint8 i=0; i < 4; ++i) {
        if (m_children[i]) {
          ostr.write( (char*)&i, sizeof(i) );
          m_children[i]->serialize(ostr);
        }
      }
    }

    /// Deserialize a tree as binary data for storage on disk.
    void deserialize(std::istream &istr) {

      if (istr.eof())
        vw_throw(IOErr() << "Tree::deserialize() reached the end of the index file prematurely!");

      // First we de-serialize our record data.
      int ib_size;
      istr.read((char*)(&ib_size), sizeof(ib_size));
      boost::shared_array<char> buffer = boost::shared_array<char>(new char[ib_size]);
      istr.read(buffer.get(), ib_size);

      
      // TODO: This code here _almost_ works, and would allow us to
      // read the pbuffer in without the extra copy through the
      // "buffer" object above.  Should be attempted again if this
      // ever becomes a performance problem.
      // google::protobuf::io::IstreamInputStream* is_istr = new google::protobuf::io::IstreamInputStream(&istr);
      // google::protobuf::io::CodedInputStream* c_istr = new google::protobuf::io::CodedInputStream(is_istr);
      // google::protobuf::io::CodedInputStream::Limit limit = c_istr->PushLimit(ib_size);
      // bool worked = m_record.ParseFromCodedStream(c_istr);
      // c_istr->PopLimit(limit);
      bool worked = m_record.ParseFromArray(buffer.get(),  ib_size);

      // Next we determine how many children we will have.
      uint8 num_children;
      istr.read( (char*)&num_children, sizeof(num_children) );

      for (uint8 i = 0; i < num_children; ++i) {
        uint8 child_id;
        istr.read( (char*)&child_id, sizeof(child_id) );
        
        // For each node that exists, we deserialize (in a depth first
        // manner), and then save the node as our child.
        boost::shared_ptr<TreeNode> node( new TreeNode(this, ElementT() ) );
        node->deserialize(istr);
        this->set_child(child_id, node);
      }
    }

  };

}} // namespace vw

#endif // __VW_PLATE_TREE_H__
