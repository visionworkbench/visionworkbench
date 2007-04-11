// __BEGIN_LICENSE__
// 
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
// 
// Copyright 2006 Carnegie Mellon University. All rights reserved.
// 
// This software is distributed under the NASA Open Source Agreement
// (NOSA), version 1.3.  The NOSA has been approved by the Open Source
// Initiative.  See the file COPYING at the top of the distribution
// directory tree for the complete NOSA document.
// 
// THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
// KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
// LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
// SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
// A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
// THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
// DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
// 
// __END_LICENSE__

/// \file PixelIterator.h
/// 
/// An STL-compliant iterator for general image views.
///
#ifndef __VW_IMAGE_PIXELITERATOR_H__
#define __VW_IMAGE_PIXELITERATOR_H__

#include <boost/iterator/iterator_facade.hpp>

#include <vw/Core/FundamentalTypes.h>

namespace vw {

  /// An STL-compliant iterator for accessing a general image view 
  /// as an ordered 1-dimensional container of pixels.
  ///
  /// The iterator first traverses rows, then columns, then planes.
  /// This iterator can be used in the generic standard template
  /// library algorithms.  Creating a fully STL-compliant iterator is
  /// actually quite a lot of work, so this class is a subclass of the
  /// boost iterator_facade<>, which takes care of defining the many
  /// types and member functions necessary for various STL iterator
  /// types.
  template <class ViewT>
  class PixelIterator : public boost::iterator_facade<PixelIterator<ViewT>,
                                                      typename ViewT::pixel_type,
                                                      boost::random_access_traversal_tag,
                                                      typename ViewT::result_type,
                                                      int64> {
      
    // This is required for boost::iterator_facade
    friend class boost::iterator_core_access;
    typedef typename ViewT::pixel_type pixel_type;

    // Private variables
    ViewT const* m_view_ptr;
    int32 m_width, m_height;
    int64 m_index, m_pixels_per_plane;

    // Testing equality and distance
    bool equal        (PixelIterator<ViewT> const& iter) const { return (m_index == iter.m_index); }
    int64 distance_to (PixelIterator<ViewT> const &iter) const { return iter.m_index - m_index;    }

    // Forward, backward, and random access movement
    void increment()      {  ++m_index;    }
    void decrement()      {  --m_index;    }
    void advance(int64 n) {  m_index += n; }

    // Dereferencing
    typename PixelIterator::reference dereference() const { 
      // Modulus arithmetic for random access iteratation
      int32 p = (int32)(m_index / m_pixels_per_plane);
      int32 r = (int32)(m_index % m_pixels_per_plane) / m_width;
      int32 c = (int32)(m_index % m_pixels_per_plane) % m_width;

      return (*m_view_ptr)(c,r,p);
    }

  public:
    // Constructors
    PixelIterator( ViewT const& view, int32 c, int32 r, int32 p=0 ) 
      : m_view_ptr(&view), m_width(view.cols()), m_height(view.rows()) {
      m_index = p*((int64) m_width*m_height) + r*m_width + c;
      m_pixels_per_plane = (int64) m_width * m_height;
    }
        
    explicit PixelIterator( ViewT const& view )
      :  m_view_ptr(&view), m_width(view.cols()), m_height(view.rows()), m_index(0) {
      m_pixels_per_plane = (int64) m_width * m_height;
    }
  };

} // namespace vw

#endif // __VW_IMAGE_PIXELITERATOR_H__
