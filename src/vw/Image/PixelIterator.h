// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
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
  /// The iterator first traverses from column to column, then from
  /// row to row, and finally from plane to plane.  This iterator can
  /// be used in the generic standard template library algorithms.
  /// Creating a fully STL-compliant iterator is actually quite a lot
  /// of work, so this class is a subclass of boost::iterator_facade,
  /// which takes care of defining the many types and member functions
  /// necessary for various STL iterator types.
  ///
  /// PixelIterator is a relatively slow mechanism for working with
  /// image pixels.  It stores a single offset into the image, and
  /// every time you dereference the iterator it uses modular
  /// arithmetic to compute the corresponding column, row, and plane
  /// number.  You might think it would be faster for the iterator to
  /// operate on those quantities directly, but then you end up having
  /// to do an equivalent amount of math every time you move the
  /// iterator---regardless of whether or not you end up dereferencing
  /// it---which turns out to be even worse in many applications.  If
  /// performance is a concern, you are generally better off working
  /// with pixel accessors rather than pixel iterators.
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
    ViewT const &m_view;
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
      return m_view(col(),row(),plane());
    }

  public:
    // Constructors
    PixelIterator( ViewT const& view, int32 c, int32 r, int32 p=0 )
      : m_view(view), m_width(view.cols()), m_height(view.rows()) {
      m_index = p*((int64) m_width*m_height) + r*m_width + c;
      m_pixels_per_plane = (int64) m_width * m_height;
    }

    PixelIterator():
      m_view(), m_width(-1), m_height(-1),
      m_index(-1), m_pixels_per_plane(-1)
    { }

    /* The compiler does not like the default copy constructor and
     * assignment operators when using the m_view variable, for some
     * reason, so we override them here.
    */
    PixelIterator( PixelIterator<ViewT> const& cpy ):
      m_view(cpy.m_view), m_width(cpy.m_width), m_height(cpy.m_height),
      m_index(cpy.m_index), m_pixels_per_plane(cpy.m_pixels_per_plane) { }

    PixelIterator& operator=( PixelIterator<ViewT> const &cpy ) {
      m_view = cpy.m_view;
      m_width = cpy.m_width;
      m_height = cpy.m_height;
      m_index = cpy.m_index;
      m_pixels_per_plane = cpy.m_pixels_per_plane;

      return *this;
    }

    explicit PixelIterator( ViewT const& view )
      :  m_view(view), m_width(view.cols()), m_height(view.rows()), m_index(0) {
      m_pixels_per_plane = (int64) m_width * m_height;
    }

    int32 col() const {
      return (int32)(m_index % m_pixels_per_plane) % m_width;
    }

    int32 row() const {
      return (int32)(m_index % m_pixels_per_plane) / m_width;
    }

    int32 plane() const {
      return (int32)(m_index / m_pixels_per_plane);
    }

  };

} // namespace vw

#endif // __VW_IMAGE_PIXELITERATOR_H__
