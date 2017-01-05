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


/// \file PixelIterator.h
///
/// An STL-compliant iterator for general image views.
///
#ifndef __VW_IMAGE_PIXELITERATOR_H__
#define __VW_IMAGE_PIXELITERATOR_H__

#include <boost/iterator/iterator_facade.hpp>

#include <vw/Core/FundamentalTypes.h>
#include <vw/Math/Vector.h>

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

  // TODO: Should probably consolidate with the BresenhamLine class.
  /// Iterates through positions in an image in an orthogonal or diagonal line.
  /// - This is a simplified version of the BresenhamLine class.
  /// - Tests for this class are located in TestPerPixelAccessorViews.cxx
  class PixelLineIterator  : public boost::iterator_facade<PixelLineIterator,
                                                      Vector2i,
                                                      boost::random_access_traversal_tag,
                                                      Vector2i,
                                                      int64> {
  public: // Definitions
  
    /// Specifies one of the eight cardinal directions
    enum Direction {TL, T, TR, L, R, BL, B, BR};
  
  public: // Functions
  
    PixelLineIterator(Vector2i start, Direction dir, Vector2i size)
      : m_col(start[0]), m_row(start[1]), m_width(size[0]), m_height(size[1]) {
      switch(dir) {
        case TL: m_dcol = -1; m_drow = -1; break;
        case T : m_dcol =  0; m_drow = -1; break;
        case TR: m_dcol =  1; m_drow = -1; break;
        case L : m_dcol = -1; m_drow =  0; break;
        case R : m_dcol =  1; m_drow =  0; break;
        case BL: m_dcol = -1; m_drow =  1; break;
        case B : m_dcol =  0; m_drow =  1; break;
        default: m_dcol =  1; m_drow =  1; break; // BR
      };
    }

    /// Return the current location along the line
    Vector2i operator*() const { return dereference();}

    int col() const {return m_col;}
    int row() const {return m_row;}

    /// Advance to the next pixel along the line
    void operator++() {increment();}

    /// Returns true if the iterator is still in the image bounds
    bool is_good() {
      return ((m_col >= 0) && (m_row >= 0) && (m_col < m_width) && (m_row < m_height));
    }

    /// Equality operator
    bool equal(PixelLineIterator const& iter) const { 
      return ((m_col == iter.m_col) && (m_row == iter.m_row)); 
    }
    
    /// Distance operator assumes the other iterator is on the same line.
    int64 distance_to(PixelLineIterator const &iter) const {
      if (m_dcol != 0)
        return (iter.m_col - m_col) / m_dcol;
      else
        return (iter.m_row - m_row) / m_drow;
    }

    // Forward, backward, and random access movement
    void increment()      { m_col +=   m_dcol;  m_row +=   m_drow; }
    void decrement()      { m_col -=   m_dcol;  m_row -=   m_drow; }
    void advance(int64 n) { m_col += n*m_dcol;  m_row += n*m_drow; }

    // Dereferencing
    Vector2i dereference() const { return Vector2f(m_col, m_row); }

    std::string to_string() const {
      std::stringstream s;
      s << "Position("<<m_col<<", "<<m_row<<"), Direction("<<m_dcol<<", "<<m_drow<<")";
      return s.str();
    }

  private: // Variables
    int m_col,   m_row;
    int m_width, m_height;
    int m_dcol,  m_drow;
    
  }; // End class PixelLineIterator

} // namespace vw

#endif // __VW_IMAGE_PIXELITERATOR_H__
