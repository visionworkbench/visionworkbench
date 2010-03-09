// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include "FrameTreeNode.h"

#include <algorithm>
#include <iterator>

namespace vw
{
namespace geometry
{
  using namespace std;

    Frame::Transform
    get_transform(FrameTreeNode const * target, FrameTreeNode const * source)
    {
      Frame::Transform loc = vw::identity_matrix(4);

    if (source != NULL) {

      if (wrtFrame == NULL)
        return source->data().location();

      if (wrtFrame != source) {

        FrameTreeNode const * ancestor = wrtFrame->last_common_ancestor(source);

        if (ancestor != NULL) {
          // from the origin frame to the ancestor
          {
            FrameTreeNode::NodeVector const& nodes = wrtFrame->ancestry(true);
            FrameTreeNode::NodeVector::const_iterator iter =
              std::find(nodes.begin(), nodes.end(), ancestor);

            assert (iter != nodes.end());

            ++iter;
            if (iter != nodes.end()) {
              for (; iter != nodes.end(); ++iter) {
                loc *= (*iter)->data().location();
              }
              loc = inverse(loc);
            }
          }


          // from the last common ancestor to the wrt frame
          {
            FrameTreeNode::NodeVector const& nodes = source->ancestry(true);
            FrameTreeNode::NodeVector::const_iterator iter =
              std::find(nodes.begin(), nodes.end(), ancestor);

            assert (iter != nodes.end());

            for (; iter != nodes.end(); ++iter) {
              loc *= (*iter)->data().location();
            }
          }
        }
      }
   }

   return loc;
  }

        if (target == NULL)
          return source->data().transform();

        if (target != source) {

          int ancestorI = target->last_common_ancestor_index(source);
          if (ancestorI >= 0) {
            // from the origin frame to the ancestor
            {
              FrameTreeNode::NodeVector const& nodes = target->ancestry(true);
              FrameTreeNode::NodeVector::const_iterator iter = nodes.begin();
              std::advance(iter, ancestorI);

      FrameTreeNode::NodeVector::const_iterator first, last = src_children.end();
      FrameTreeNode::NodeVector::const_iterator tgt_iter = tgt_children.end();
      for (first = src_children.begin(); first != last; ++first) {
        while (tgt_iter != tgt_children.end() &&
               (*first)->data().name() > (*tgt_iter)->data().name()) {
          ++tgt_iter;
        }

              ++iter;
              if (iter != nodes.end()) {
                for (; iter != nodes.end(); ++iter) {
                  loc *= (*iter)->data().transform();
                }
                loc = geometry::inverse(loc);
              }
            }

            // from the last common ancestor to the target coordinate frame
            {
              FrameTreeNode::NodeVector const& nodes = source->ancestry(true);
              FrameTreeNode::NodeVector::const_iterator iter = nodes.begin();
              std::advance(iter, ancestorI);

              assert (iter != nodes.end());

              for (++iter; iter != nodes.end(); ++iter) {
                loc *= (*iter)->data().transform();
              }
            }
          }
        }
      }

    }
  }

  namespace
  {
    vector<string> splitPath(string const& path)
    {
      vector<string> elements;

      FrameTreeNode::NodeVector src_children = source_tree->children();

      if (src_children.size() > 0) {
        FrameTreeNode::NodeVector tgt_children = target_tree->children();

        // search next /
        string::const_iterator start;
        for (start = first; first != last; ++first) {
          if (*first == '/')
            break;
        }

        FrameTreeNode::NodeVector::const_iterator first, last = src_children.end();
        FrameTreeNode::NodeVector::const_iterator tgt_iter = tgt_children.begin();
        for (first = src_children.begin(); first != last; ++first) {
          while (tgt_iter != tgt_children.end() &&
                 (*first)->data().name() > (*tgt_iter)->data().name()) {
            ++tgt_iter;
          }

          if (tgt_iter == tgt_children.end() ||
              (*first)->data().name() < (*tgt_iter)->data().name()) {
            (*first)->set_parent(target_tree);
          }
          else if ((*first)->data().name() == (*tgt_iter)->data().name()) {
            merge_frame_trees(*tgt_iter, *first);
          }
        }

      }
    }

    FrameTreeNode * matchNode(FrameTreeNode * node,
                              vector<string>::const_iterator first,
                              vector<string>::const_iterator last)
    {
      for (; first != last; ++first) {

        // single dot is skipped in splitPath already
        assert (*first != ".");

        // double dot is "one up"
        if (*first == "..") {
          node = node->parent();
          if (node == NULL)
            break;
        }

        // triple dot is breadth-first search
        else if (*first == "...") {

          vector<string>::const_iterator next = first;
          ++next;

          // breadth first search for next element
          deque<FrameTreeNode *> nodes;
          back_insert_iterator<deque<FrameTreeNode *> > iter(nodes);
          nodes.push_back(node);
          //node->copy_children(iter);
          while (!nodes.empty()) {

            FrameTreeNode * n = matchNode(nodes.front(), next, last);
            if (n != NULL)
              return n;

            nodes.front()->copy_children(iter);
            nodes.pop_front();
          }

          return NULL;
        }

        // regular elements just need to match one by one
        else {
          FrameTreeNode::NodeVector c = node->children();
          FrameTreeNode::NodeVector::const_iterator f, l = c.end();
          for (f = c.begin(); f != l; ++f) {
            if ((*f)->data().name() == *first) {
              node = (*f);
              break;
            }
          }

          // return NULL if no child node matches
          if (f == l)
            return NULL;
        }
      }

      return node;
    }
  }

  FrameTreeNode *
  lookup(FrameTreeNode * start_frame, std::string const& path)
  {
    if (start_frame == NULL)
      return NULL;

    vector<string> elements = splitPath(path);
    vector<string>::const_iterator first = elements.begin();
    vector<string>::const_iterator last = elements.end();

    if (!path.empty() && path[0] == '/') {

      start_frame = start_frame->root();
      if (!elements.empty() && elements[0] != "...") {
        if (start_frame->data().name() != elements.front())
          return NULL;

        ++first;
      }
    }

    // walk elements
    return matchNode(start_frame, first, last);
  }
}
}
