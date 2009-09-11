// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <cxxtest/TestSuite.h>
#include <boost/smart_ptr.hpp>

#include <vw/Plate/Tree.h>

using namespace std;
using namespace vw::platefile;

class TestTreeNode : public CxxTest::TestSuite {

public:
  typedef TreeNode<int> tree_type;

  void test_simple_tree_manipulation() {
    boost::shared_ptr<tree_type> root( new tree_type(0) );
    TS_ASSERT_EQUALS(root->value(), 0);

    // Insert some level 1 children
    root->insert(10, 0, 0, 1);
    root->insert(11, 1, 0, 1);
    root->insert(12, 0, 1, 1);
    root->insert(13, 1, 1, 1);

    TS_ASSERT_EQUALS(root->child(0)->value(), 10);
    TS_ASSERT_EQUALS(root->child(1)->value(), 11);
    TS_ASSERT_EQUALS(root->child(2)->value(), 12);
    TS_ASSERT_EQUALS(root->child(3)->value(), 13);

    // Insert one level 2 child
    root->insert(20, 3, 3, 2);
    root->insert(30, 3, 3, 3);
    TS_ASSERT_EQUALS(root->child(3)->child(3)->value(), 20);

    // Test an invalid index
    TS_ASSERT_THROWS(root->insert(10,2,1,1), IndexErr);
    TS_ASSERT_THROWS(root->insert(1,0,1,0), IndexErr);
  }

  void test_search() {

    boost::shared_ptr<tree_type> root( new tree_type(0) );
    TS_ASSERT_EQUALS(root->value(), 0);

    // Insert some level 1 children
    root->insert(10, 0, 0, 1);
    root->insert(11, 1, 0, 1);
    root->insert(12, 0, 1, 1);
    root->insert(13, 1, 1, 1);

    // Insert one level 2 child and one level 3 child
    root->insert(20, 3, 3, 2);
    root->insert(30, 3, 3, 3);

    // Now search for the nodes in level 3
    int result = root->search(3,3,2);
    TS_ASSERT_EQUALS(result, 20);

    // Check level 2
    result = root->search(3,3,3);
    TS_ASSERT_EQUALS(result, 30);

    // Check level 1
    result = root->search(0,0,1);
    TS_ASSERT_EQUALS(result, 10);
    result = root->search(1,0,1);
    TS_ASSERT_EQUALS(result, 11);
    result = root->search(0,1,1);
    TS_ASSERT_EQUALS(result, 12);
    result = root->search(1,1,1);
    TS_ASSERT_EQUALS(result, 13);
    
    // Check level 0
    result = root->search(0,0,1);

    // Test searches on nodes that should be invalid
    TS_ASSERT_THROWS(root->search(3,4,3), TileNotFoundErr);
    TS_ASSERT_THROWS(root->search(3,8,3), IndexErr);
    TS_ASSERT_THROWS(root->search(3,4,5), TileNotFoundErr);
    TS_ASSERT_THROWS(root->search(0,0,2), TileNotFoundErr);
  }

}; // class TestTree
